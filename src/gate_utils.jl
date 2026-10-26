export static_flux_coupling, get_qubit_subspace_couplings, build_H_S1, build_H_LC, sweep_couplings
export run_spectography_circuit, run_spectography_flux_sweep_2a
export flux_drive_operators
export cz_gate_matrix, get_νφ, unwrap_phase, subspace_leakage, static_ζ, get_νφ_point, get_νφ_sweep, run_cz_flux_sweep_2b
function correct_1q_phase_CZ(U, ψ_ge, ψ_eg, ψ_gg, ψ_ee)
    # set phase reference to phase of gg state
    ζ_gg = ψ_gg' * U * ψ_gg |> angle
    U_0 = exp(-im*ζ_gg) * U

    # get single-qubit Z rotation angles
    ζ_ge, ζ_eg = [ψ' * U_0 * ψ |> angle for ψ in (ψ_ge, ψ_eg)]

    # compose Z-correction matrix
    VZ = (ψ_gg * ψ_gg' +
        exp(-im*ζ_ge) * ψ_ge * ψ_ge' +
        exp(-im*ζ_eg) * ψ_eg * ψ_eg' +
        exp(-im * (ζ_ge + ζ_eg)) * ψ_ee * ψ_ee')
    return VZ * U_0
end

### Utils for Reproducing Native CZ Gate Results

function static_flux_coupling(
    sys::CompositeQuantumSystem,
    δϕA::Real,
    δϕB::Real,
    J::Real;
    fA_lbl::Symbol=:fA,
    fB_lbl::Symbol=:fB,
    )
    return J * (sys[fA_lbl, :φ] + δϕA) * (sys[fB_lbl, :φ] + δϕB)
end

"""
    flux_drive_operators(sys, J; fA_lbl=:fA, fB_lbl=:fB)

Operators `(ΦA_op, ΦB_op)` through which the flux detunings of the two fluxonia enter the Hamiltonian, i.e. ∂H/∂δϕA
and ∂H/∂δϕB, for the coupling `-J(φA + δϕA)(φB + δϕB)` of [`build_H_S1`](@ref):

    ΦA_op = EL_A φA - J φB,    ΦB_op = EL_B φB - J φA

The Hamiltonian is `sys.H + δϕA ΦA_op + δϕB ΦB_op`, up to constants. This is exact when the external flux enters
through the linear inductive terms.
"""
function flux_drive_operators(sys::CompositeQuantumSystem, J::Real; fA_lbl::Symbol=:fA, fB_lbl::Symbol=:fB)
    ELA, ELB = [sys[lbl].params.EL for lbl in (fA_lbl, fB_lbl)]
    ΦA_op = ELA * sys[fA_lbl, :φ] - J * sys[fB_lbl, :φ]
    ΦB_op = ELB * sys[fB_lbl, :φ] - J * sys[fA_lbl, :φ]
    return ΦA_op, ΦB_op
end

"""
    dynamical_flux_terms(sys, δϕtA, δϕtB, J; symmetric_drive=false, amplitudes=nothing, fA_lbl=:fA, fB_lbl=:fB)

Linear response of the Hamiltonian to time-dependent flux detunings `δϕtA(p, t)`, `δϕtB(p, t)` of the two fluxonia,
i.e. ∂H/∂δϕ contracted with the detunings, as a `QobjEvo`. The coupling enters as `-J(φA + δϕA)(φB + δϕB)`, see
[`build_H_S1`](@ref). This is exact when the external flux enters through the linear inductive terms.

By default the two qubits are driven by their own functions, as two time-dependent terms. If both are driven by the
same pulse shape, the terms can be merged into one, which is cheaper to integrate:
- `amplitudes=(ΩA, ΩB)`: both qubits follow one function (`δϕtA`, or `δϕtB` if `δϕtA` is `nothing`), scaled by `ΩA`
        for qubit A and `ΩB` for qubit B.
- `symmetric_drive=true`: shorthand for `amplitudes=(1, 1)`.
The two options are mutually exclusive, and the function of the unused qubit may be `nothing`.
"""
function dynamical_flux_terms(
    sys::CompositeQuantumSystem,
    δϕtA::Union{FNCA, Nothing}, # assumes takes parameters (p, t::Real)
    δϕtB::Union{FNCB, Nothing}, # assumes takes parameters (p, t::Real)
    J::Real;
    symmetric_drive::Bool=false,
    amplitudes::Union{Tuple{Real, Real}, Nothing}=nothing,
    fA_lbl::Symbol=:fA,
    fB_lbl::Symbol=:fB,
    ) where {FNCA<:Function, FNCB<:Function}
    symmetric_drive && !isnothing(amplitudes) &&
        throw(ArgumentError("`symmetric_drive` and `amplitudes` are mutually exclusive"))
    ΦAt_op, ΦBt_op = flux_drive_operators(sys, J; fA_lbl=fA_lbl, fB_lbl=fB_lbl)
    amps = symmetric_drive ? (1, 1) : amplitudes
    if !isnothing(amps)
        Ht = something(δϕtA, δϕtB) |> dr -> QobjEvo((amps[1] * ΦAt_op + amps[2] * ΦBt_op, dr))
    else
        Ht = mapreduce(+, (δϕtA, δϕtB), (ΦAt_op, ΦBt_op)) do δϕt, Φt_op
            isnothing(δϕt) ? 0 : QobjEvo((Φt_op, δϕt))
        end
    end
    return Ht
end

function get_qubit_subspace_couplings(
    sys::CompositeQuantumSystem,
    H_int::QuantumObject;
    bare_to_dressed_tracking::Symbol=:overlap
    )
    σA, σB = map(sys.subsys) do subsys
        g, e = [state(subsys, lbl; bare_to_dressed_tracking=bare_to_dressed_tracking)[2] for lbl in (:g, :e)]
        σx = g*e' + e*g'
        σz = g*g' - e*e'
        return (;σx=σx, σz=σz)
    end
    σxx = σA.σx ⊗ σB.σx
    σxz = σA.σx ⊗ σB.σz
    σzx = σA.σz ⊗ σB.σx
    σzz = σA.σz ⊗ σB.σz

    xx = tr(σxx * H_int) / 4
    xz = tr(σxz * H_int) / 4
    zx = tr(σzx * H_int) / 4
    zz = tr(σzz * H_int) / 4
    return (;xx=xx, xz=xz, zx=zx, zz=zz)
end

function build_H_S1(
    fA_params::Union{Dict, NamedTuple},
    fB_params::Union{Dict, NamedTuple},
    J::Real;
    qho_basis_nlev::Int=50,
    ebasis_nlev::Int=6,
    )
    # `fA_params`, `fB_params` must include the subsystem names, e.g. `name=:fA`
    fA, fB = map((fA_params, fB_params)) do params
        fluxonium(;params..., nlev=qho_basis_nlev) |>
            to_energy_basis |>
            sys -> truncate_sys(sys, ebasis_nlev)
    end
    sys_bare = fA ⊗ fB
    δϕA, δϕB = [subsys.params.δϕ for subsys in (fA, fB)]

    fA_cross_term = J * δϕB * fA.φ
    fB_cross_term = J * δϕA * fB.φ
    fA_tilted = fA - fA_cross_term
    fB_tilted = fB - fB_cross_term
    # coupling in the Hamiltonian is -J(φA + δϕA)(φB + δϕB) (opposed inductor orientations); the
    # single-qubit cross terms are moved into the tilted fluxonia, `H_int` is what remains
    H_int = -(
        static_flux_coupling(sys_bare, δϕA, δϕB, J; fA_lbl=fA.name, fB_lbl=fB.name)
        - (fA_cross_term ⊗ fB.I + fA.I ⊗ fB_cross_term)
    )
    return (fA_tilted ⊗ fB_tilted) + H_int, H_int
end

"""
    build_H_LC(fA_params, fB_params, J; J_C=0.0, kwargs...)

[`build_H_S1`](@ref) with an additional capacitive coupling `+J_C nA nB` between the charge operators of the two
fluxonia, i.e. an inductive and capacitive coupling `-J φφ + J_C nn`. This is the joint sign flip of the convention
`+g_L φφ - g_C nn` of 2Q.ipynb, which leaves ζzz unchanged, with `J = g_L` and `J_C = g_C`. `kwargs` are passed to
`build_H_S1`.

Like `build_H_S1`, it returns `(sys, H_int)`, with `H_int` including the capacitive term. Use
`H_builder_kwargs=(J_C=…,)` to pass `J_C` to the sweeps, with `first ∘ build_H_LC` as `H_builder` where only the
system is needed.
"""
function build_H_LC(
    fA_params::Union{Dict, NamedTuple},
    fB_params::Union{Dict, NamedTuple},
    J::Real;
    J_C::Real=0.0,
    kwargs...
    )
    sys, H_int = build_H_S1(fA_params, fB_params, J; kwargs...)
    H_C = J_C * sys[fA_params[:name], :n, fB_params[:name], :n]
    return sys + H_C, H_int + H_C
end

function sweep_couplings(
    δϕs::AbstractArray{<:Real},
    fA_params::Union{NamedTuple, Dict},
    fB_params::Union{NamedTuple, Dict},
    H_builder::Hfnc,
    J::Real;
    qho_basis_nlev::Int=50,
    ebasis_nlev::Int=6,
    bare_to_dressed_tracking::Symbol=:overlap,
    ) where {Hfnc<:Function}
    # shift the flux offset δϕ of both fluxonia by the same amount; `fluxonium`'s default offset is π
    base_params = map(p -> (; p..., δϕ=get(p, :δϕ, π)), (fA_params, fB_params))
    results = map(δϕs) do δϕ
        δfA_p, δfB_p = map(p -> (; p..., δϕ=p.δϕ + δϕ), base_params)
        δsys, H_int = H_builder(δfA_p, δfB_p, J;
                                qho_basis_nlev=qho_basis_nlev, ebasis_nlev=ebasis_nlev)
        return get_qubit_subspace_couplings(δsys, H_int;
                                            bare_to_dressed_tracking=bare_to_dressed_tracking)
    end
    return StructArray(results)
end

"""
    run_spectography_circuit(sys, H_dr, tlist; dr_params=nothing, integrator_pbar=false, state_tracking=:overlap, fA_lbl=:fA, sesolve_kwargs...)

Evolve the symmetric superposition (|gg⟩+|ge⟩)/√2 of the dressed states of `sys` under `sys.H + H_dr`
and return the population that has been transferred to qubit A's excited state, |⟨eg|ψ⟩|² + |⟨ee|ψ⟩|²,
at the end of `tlist`. The state is evolved directly with `sesolve`, which is much cheaper than computing
the full propagator, since only P(t_end) for a single initial state is needed.

# Arguments
- `sys`: system with the static Hamiltonian `sys.H`, whose dressed states are labelled `:gg, :ge, :eg, :ee`.
- `H_dr`: time-dependent drive Hamiltonian, a `QobjEvo` whose coefficient functions take `(p, t)`.
- `tlist`: times passed to `sesolve`; the result is evaluated at `tlist[end]`.
- `dr_params`: parameters `p` handed to the coefficient functions of `H_dr`.
- `integrator_pbar`: forwarded to `sesolve`. Leave it `false` when calling from a sweep that has its own
        progress bar, otherwise the per-call bars dominate the output.
- `state_tracking`: `:overlap` (default) or `:spectrum`, see `state`.
- `fA_lbl`: label of the driven subsystem (currently unused, kept for API symmetry with the other helpers).
- `sesolve_kwargs`: forwarded to `sesolve`, e.g. `abstol`, `reltol`.
"""
function run_spectography_circuit(
    sys::QuantumSystem,
    H_dr::QobjEvo,
    tlist::AbstractArray;
    dr_params::Union{PP, Nothing} = nothing,
    integrator_pbar::Bool=false,
    state_tracking::Symbol=:overlap,
    fA_lbl::Symbol=:fA,
    sesolve_kwargs...
    ) where {PP}
    gg, ge, eg, ee = (state(sys, lbl; bare_to_dressed_tracking=state_tracking)[2]
    for lbl in [:gg, :ge, :eg, :ee])
    ψ_0 = (1/sqrt(2)) * (gg + ge)
    # evolve the initial state directly instead of computing the full propagator
    sol = sesolve(sys.H + H_dr, ψ_0, tlist;
                  progress_bar=Val(integrator_pbar), params=dr_params, sesolve_kwargs...)
    ψ_final = sol.states[end]
    # total population in fA excited state
    P = abs2(eg' * ψ_final) + abs2(ee' * ψ_final)  # population in the |eg⟩ and |ee⟩ states
    # return probability value
    return P
end

"""
    run_spectography_flux_sweep_2a(fA_params, fB_params, tlist, env, δϕB_list, ω_list, J, Ω_ge;
        H_builder=first ∘ build_H_S1, H_builder_kwargs=nothing, dr_params=nothing, δϕA=-2π*0.042,
        integrator_pbar=false, showprogress=true, state_tracking=:overlap, sesolve_kwargs...)

Spectroscopy sweep of Fig. 2(a): drive qubit A's g↔e transition with `Ω_ge * n_A * env(p, t) * cos(ω t)`
and record the population transferred out of (|gg⟩+|ge⟩)/√2 (see [`run_spectography_circuit`](@ref)) on a
grid of drive frequencies and flux detunings of qubit B.

# Arguments
- `fA_params`, `fB_params`: `fluxonium` parameters of the two qubits. They must include the subsystem
        `name` and the flux offset `δϕ` (π at the sweet spot). `δϕ` is shifted by `δϕA` and the entries of
        `δϕB_list`.
- `tlist`: `[0, T]`, the drive acts on this window and the population is evaluated at `T`.
- `env`: pulse envelope `env(p, t)` without the carrier, e.g. wrapping `cosine_ramp_envelope` as
        `env(p, t) = env0(t)`.
- `δϕB_list`: detunings of qubit B from its sweet spot, as flux angles 2πδΦ/Φ0.
- `ω_list`: angular drive frequencies. `ω` is added to `p`, overriding any `ω` in `dr_params`, so `env`
        may read `p.ω`.
- `J`: inductive coupling, passed to `H_builder`. For the default builder it is the magnitude, the
        Hamiltonian contains `-J(φA + δϕA)(φB + δϕB)`, see [`build_H_S1`](@ref).
- `Ω_ge`: drive strength multiplying the charge operator of qubit A.
- `H_builder`: function `H_builder(pA, pB, J; H_builder_kwargs...)` returning the composite system `sys`,
        whose subsystems are named as in `fA_params` and `fB_params`. The default `first ∘ build_H_S1` keeps
        the system from the `(sys, H_int)` tuple returned by [`build_H_S1`](@ref). Unlike the builder of
        [`sweep_couplings`](@ref), it must not return a tuple.
- `H_builder_kwargs`: NamedTuple of keyword arguments for `H_builder`, e.g.
        `(qho_basis_nlev=50, ebasis_nlev=6)` for `build_H_S1`. Defaults to none.
- `dr_params`: NamedTuple of extra parameters for `env`.
- `δϕA`: detuning of qubit A from its sweet spot (default 0.458 Φ0, i.e. -0.042 Φ0 from Φ0/2).
- `showprogress`: show a progress bar over the `length(ω_list) * length(δϕB_list)` grid points.
- `integrator_pbar`, `state_tracking`, `sesolve_kwargs`: see [`run_spectography_circuit`](@ref).

# Returns
- Matrix of transferred populations, indexed `[i_ω, j_ϕB]`.
"""
function run_spectography_flux_sweep_2a(
    fA_params::Union{Dict, NamedTuple},
    fB_params::Union{Dict, NamedTuple},
    tlist::AbstractVector{<:Real},
    env::ENVFNC, # assumes function accepts (p,t::Real)
    δϕB_list::AbstractVector{<:Real},
    ω_list::AbstractVector{<:Real},
    J::Real,
    Ω_ge::Real;
    H_builder::HFUNC=(first ∘ build_H_S1), # assumes function accepts (fA_params, fB_params, J; H_builder_kwargs... )
    H_builder_kwargs::HPP=nothing,
    dr_params::PP=nothing,
    δϕA::Real=-2π * 0.042,
    integrator_pbar::Bool=false,
    showprogress::Bool=true,
    state_tracking::Symbol=:overlap,
    sesolve_kwargs...
        ) where {ENVFNC<:Function, HFUNC<:Function, PP, HPP}
    fA_lbl = fA_params[:name]
    δϕ0A, δϕ0B = fA_params[:δϕ], fB_params[:δϕ]
    pA = (; fA_params..., δϕ=δϕ0A + δϕA)
    P = Matrix{Float64}(undef, length(ω_list), length(δϕB_list))
    pbar = Progress(length(P); desc="Spectroscopy sweep: ", output=stdout, showspeed=true, enabled=showprogress)
    Threads.@threads for j in eachindex(δϕB_list)
        pB = (; fB_params..., δϕ=δϕ0B + δϕB_list[j])
        Hkwargs = something(H_builder_kwargs, (;))
        sys = H_builder(pA, pB, J; Hkwargs...)
        dr(p, t) = env(p, t) * cos(p.ω * t)
        H_dr = QobjEvo((Ω_ge * sys[fA_lbl, :n], dr)) # same for all ω, which enters through `p`
        for (i, ω) in enumerate(ω_list)
            p = something(dr_params, (;)) |> pp -> merge(pp, (;ω))
            P[i, j] = run_spectography_circuit(sys, H_dr, tlist; dr_params=p, integrator_pbar=integrator_pbar,
                                               state_tracking=state_tracking, fA_lbl=fA_lbl, sesolve_kwargs...)
            next!(pbar)
        end
    end
    finish!(pbar)
    return P
end

#function get_subspace_propagator(
    #Hevo::AbstractQuantumObject,
    #H0::QuantumObject,
    #tlist::AbstractVector{<:Real},
    #subspace_nlev::Int=4;
    #dr_params::DP=nothing,
    #progress_bar::Bool=false,
    #sesolve_kwargs...
        #) where {PP, DP}
    #_, ψs, _ = eigenstates(H0)
    #ψs_sub = ψs[1:subspace_nlev]
    #ψs_comp = ψs[subspace_nlev+1 : end]
    #Vt = [0 * qeye_like(Hevo, Val(true)) for _ in tlist]
    #V_comp = mapreduce(ψ -> ψ*ψ', +, ψs_comp)
    #for ψr in ψs_sub
        #Vψr = sesolve(Hevo, ψr, tlist; progress_bar=Val(progress_bar), params=dr_params, sesolve_kwargs...).states
        #for ψl in ψs_sub
            #@. Vt += (Ref(ψl') * Uψr) * Ref(ψl * ψr')
        #end
    #end
    #return Vt .+ Ref(V_comp)
#end
#


### CZ gate from flux pulses (Fig. 2b)

"""
    cz_gate_matrix(sys, τ, t_ramp, ΩA, ΩB, J; ω=0.0, ramp_profile=:tanh, t_rise=nothing, abstol=1e-10,
        reltol=1e-10, bare_to_dressed_tracking=:overlap, fA_lbl=:fA, fB_lbl=:fB)

Evolve the dressed computational states `gg, ge, eg, ee` of the idle system `sys` under one CZ flux pulse of
total length `τ` and return the matrix of the evolution in the computational basis.

The flux of qubit A (B) is pulsed by `ΩA * f(t)` (`ΩB * f(t)`), where `f` is
`full_pulse(ω, t_ramp, τ - 2 * t_ramp; ramp_profile, t_rise)`, i.e. a square pulse with ramps for `ω = 0`. The
pulsed Hamiltonian is `sys.H` plus the linear response to the flux, see [`dynamical_flux_terms`](@ref), which is
exact when the external flux enters through the linear inductive terms of the fluxonia. The two drives share the
pulse and are merged into one term (`amplitudes`), which is cheaper to integrate.

# Arguments
- `sys`: idle composite system, including the static coupling.
- `τ`, `t_ramp`: total pulse length and ramp length.
- `ΩA`, `ΩB`: flux amplitudes, as flux angles 2πδΦ/Φ0.
- `J`: inductive coupling, see [`build_H_S1`](@ref).
- `ω`: carrier frequency, `0` for the square pulses of Fig. 2(b).
- `ramp_profile`, `t_rise`: see [`full_pulse`](@ref).
- `abstol`, `reltol`: tolerances of `sesolve`. They are tight, since the lab-frame phase reaches ~10³ rad over a
        few hundred ns, and looser values give phase errors comparable to the conditional phases.
- `bare_to_dressed_tracking`: `:overlap` (default) or `:spectrum`, see `state`. The sweeps of Fig. 2(b) cross a
        level crossing of the qubits, where the energy ordering of `:spectrum` swaps the labels.
- `fA_lbl`, `fB_lbl`: labels of the two fluxonia in `sys`.

# Returns
- `(C, sts)`: `C[l, r] = ⟨ψ_l|U(τ)|ψ_r⟩` with `ψ ∈ (gg, ge, eg, ee)`, and `sts = (; gg, ge, eg, ee)` the idle
        dressed states.
"""
function cz_gate_matrix(
    sys::CompositeQuantumSystem,
    τ::Real,
    t_ramp::Real,
    ΩA::Real,
    ΩB::Real,
    J::Real;
    ω::Real=0.0,
    ramp_profile::Symbol=:tanh,
    t_rise::Maybe{Real}=nothing,
    abstol::Real=1e-10,
    reltol::Real=1e-10,
    bare_to_dressed_tracking::Symbol=:overlap,
    fA_lbl::Symbol=:fA,
    fB_lbl::Symbol=:fB,
    )
    pulse = full_pulse(ω, t_ramp, τ - 2 * t_ramp; ramp_profile=ramp_profile, t_rise=t_rise)
    # both qubits share the pulse, so the drive is a single time-dependent term, which is cheaper per time step
    H_dr = dynamical_flux_terms(sys, (p, t) -> pulse(t), nothing, J; amplitudes=(ΩA, ΩB),
                                fA_lbl=fA_lbl, fB_lbl=fB_lbl)
    H = sys.H + H_dr
    lbls = (:gg, :ge, :eg, :ee)
    ψs = [state(sys, lbl; bare_to_dressed_tracking=bare_to_dressed_tracking)[2] for lbl in lbls]
    # only the evolved computational states are needed, not the full propagator
    Uψs = [sesolve(H, ψ, [0.0, τ]; progress_bar=Val(false), abstol=abstol, reltol=reltol).states[end] for ψ in ψs]
    C = [ψl' * Uψr for ψl in ψs, Uψr in Uψs]
    return C, NamedTuple{lbls}(Tuple(ψs))
end

"""
    subspace_leakage(C)

Average population that leaves the 4-dimensional computational subspace, `1 - Tr(C†C)/4`, for the
matrix `C` of [`cz_gate_matrix`](@ref).
"""
subspace_leakage(C::AbstractMatrix{<:Number}) = 1 - real(sum(abs2, C)) / 4

"""
    get_νφ(C, τ)

Conditional phase ϕ = arg(c_ee c_gg / (c_eg c_ge)) with `c_k = C[k, k]` the diagonal of the matrix `C` of
[`cz_gate_matrix`](@ref), i.e. the difference of qubit A's Ramsey phases with B in |e⟩ vs |g⟩ (ideal π/2 pulses).
Single-qubit dynamic phases cancel; ϕ = -∫ζ dt with ζ = E_ee - E_eg - E_ge + E_gg. ϕ is only known mod 2π, see
[`get_νφ_sweep`](@ref) for the unwrapped version.

# Returns
- `(ϕ [rad], ϕ/τ [rad/µs], leakage 1-|c_k|² per state)`, the leakage as a NamedTuple `(; gg, ge, eg, ee)`.
"""
function get_νφ(C::AbstractMatrix{<:Number}, τ::Real)
    c = (; gg=C[1, 1], ge=C[2, 2], eg=C[3, 3], ee=C[4, 4])
    ϕ = angle(c.ee * c.gg / (c.eg * c.ge))
    leak = map(x -> 1 - abs2(x), c)
    return ϕ, 1e3 * ϕ / τ, leak
end

"""
    unwrap_phase(ϕs)

Unwrap the phases `ϕs` (mod 2π) into a continuous sequence, assuming successive entries differ by less than π.
"""
function unwrap_phase(ϕs)
    out = collect(float.(ϕs))
    for i in 2:length(out)
        out[i] = out[i-1] + rem2pi(ϕs[i] - ϕs[i-1], RoundNearest)
    end
    return out
end

"""
    static_ζ(fA_params, fB_params, δϕA, δϕB, J; H_builder=first ∘ build_H_S1, H_builder_kwargs=nothing,
        bare_to_dressed_tracking=:overlap)

Signed static ZZ rate ζ = E_ee - E_eg - E_ge + E_gg (rad/ns) with qubit A (B) held at the flux detuning `δϕA`
(`δϕB`) from the offset `δϕ` in its parameters, as flux angles 2πδΦ/Φ0. ζ = 4 g_zz.

`H_builder` and `H_builder_kwargs` build the system from the shifted parameters, see
[`run_spectography_flux_sweep_2a`](@ref).
"""
function static_ζ(
    fA_params::Union{Dict, NamedTuple},
    fB_params::Union{Dict, NamedTuple},
    δϕA::Real,
    δϕB::Real,
    J::Real;
    H_builder::HFUNC=(first ∘ build_H_S1),
    H_builder_kwargs::HPP=nothing,
    bare_to_dressed_tracking::Symbol=:overlap,
    ) where {HFUNC<:Function, HPP}
    pA = (; fA_params..., δϕ=fA_params[:δϕ] + δϕA)
    pB = (; fB_params..., δϕ=fB_params[:δϕ] + δϕB)
    sys = H_builder(pA, pB, J; something(H_builder_kwargs, (;))...)
    E = map(l -> real(energy(sys, l; bare_to_dressed_tracking=bare_to_dressed_tracking)),
            (; gg=:gg, ge=:ge, eg=:eg, ee=:ee))
    return E.ee - E.eg - E.ge + E.gg
end

"""
    get_νφ_sweep(sys, τs, t_ramp, ΩA, ΩB, J; ω=0.0, ramp_profile=:tanh, t_rise=nothing, abstol=1e-10,
        reltol=1e-10, bare_to_dressed_tracking=:overlap, fA_lbl=:fA, fB_lbl=:fB)

Run [`cz_gate_matrix`](@ref) for each total pulse length τ in `τs` and unwrap ϕ(τ). Requires |vϕ|·step(τs) < π and
|ϕ(τs[1])| < π (keep `τs[1]` short), and `τs[1] ≥ 2 t_ramp`.

# Returns
- `(vϕ [rad/µs] = ϕ(τ_max)/τ_max as in Eq. (3), plateau slope dϕ/dτ [rad/µs], unwrapped ϕs, max non-adiabatic
        error over all τ and computational states)`.
"""
function get_νφ_sweep(
    sys::CompositeQuantumSystem,
    τs::AbstractVector{<:Real},
    t_ramp::Real,
    ΩA::Real,
    ΩB::Real,
    J::Real;
    ω::Real=0.0,
    ramp_profile::Symbol=:tanh,
    t_rise::Maybe{Real}=nothing,
    abstol::Real=1e-10,
    reltol::Real=1e-10,
    bare_to_dressed_tracking::Symbol=:overlap,
    fA_lbl::Symbol=:fA,
    fB_lbl::Symbol=:fB,
    )
    res = map(τs) do τ
        C, _ = cz_gate_matrix(sys, τ, t_ramp, ΩA, ΩB, J; ω=ω, ramp_profile=ramp_profile, t_rise=t_rise,
                              abstol=abstol, reltol=reltol, bare_to_dressed_tracking=bare_to_dressed_tracking,
                              fA_lbl=fA_lbl, fB_lbl=fB_lbl)
        get_νφ(C, τ)
    end
    ϕs = unwrap_phase(first.(res))
    leak = maximum(maximum(r[3]) for r in res)
    # linear fit with intercept: the ramps only add a τ-independent offset
    _, slope = hcat(ones(length(τs)), collect(τs)) \ ϕs
    return 1e3 * ϕs[end] / τs[end], 1e3 * slope, ϕs, leak
end

"""
    get_νφ_point(fA_params, fB_params, sys_idle, τ, t_ramp, ΩA, ΩB, J; ramp_profile=:tanh, t_rise=nothing,
        abstol=1e-10, reltol=1e-10, H_builder=first ∘ build_H_S1, H_builder_kwargs=nothing,
        bare_to_dressed_tracking=:overlap)

One Ramsey-type CZ experiment (Fig. 2(b)): flux pulses of amplitude `ΩA`, `ΩB` (flux angles) and total length `τ`
(ns) on both qubits, square pulses (ω = 0). `sys_idle` is the idle system built from `fA_params`, `fB_params`.

[`get_νφ`](@ref) only gives ϕ mod 2π, and a 200 ns pulse accumulates up to ~60 rad. The 2π branch is chosen as the
one closest to the static estimate -ζ·(pulse area), with pulse area = τ - t_ramp for both ramp profiles and ζ from
[`static_ζ`](@ref). The static estimate agrees with the dynamics to <1% (checked at δΦ = -0.06 Φ0), far below the
π needed to pick the wrong branch; `branch_err` (the wrapped distance to the static estimate, rad) flags points
where that assumption fails -- treat |branch_err| ≳ π/2 as unreliable.

# Returns
- `(vϕ [rad/µs, paper convention vϕ = 4 g_zz = -ϕ/τ], max non-adiabatic error over the computational states,
        branch_err, average leakage out of the computational subspace)`.
"""
function get_νφ_point(
    fA_params::Union{Dict, NamedTuple},
    fB_params::Union{Dict, NamedTuple},
    sys_idle::CompositeQuantumSystem,
    τ::Real,
    t_ramp::Real,
    ΩA::Real,
    ΩB::Real,
    J::Real;
    ramp_profile::Symbol=:tanh,
    t_rise::Maybe{Real}=nothing,
    abstol::Real=1e-10,
    reltol::Real=1e-10,
    H_builder::HFUNC=(first ∘ build_H_S1),
    H_builder_kwargs::HPP=nothing,
    bare_to_dressed_tracking::Symbol=:overlap,
    ) where {HFUNC<:Function, HPP}
    C, _ = cz_gate_matrix(sys_idle, τ, t_ramp, ΩA, ΩB, J; ramp_profile=ramp_profile, t_rise=t_rise,
                          abstol=abstol, reltol=reltol, bare_to_dressed_tracking=bare_to_dressed_tracking,
                          fA_lbl=fA_params[:name], fB_lbl=fB_params[:name])
    ϕ, _, err = get_νφ(C, τ)
    ζ = static_ζ(fA_params, fB_params, ΩA, ΩB, J; H_builder=H_builder, H_builder_kwargs=H_builder_kwargs,
                 bare_to_dressed_tracking=bare_to_dressed_tracking)
    ϕ_ref = -ζ * (τ - t_ramp)
    branch_err = rem2pi(ϕ - ϕ_ref, RoundNearest)
    ϕ_unwrapped = ϕ_ref + branch_err
    return -1e3 * ϕ_unwrapped / τ, maximum(err), branch_err, subspace_leakage(C)
end

"""
    run_cz_flux_sweep_2b(fA_params, fB_params, δϕA_list, δϕB_list, J; τ=200.0, t_ramp=12.0, ramp_profile=:tanh,
        t_rise=nothing, abstol=1e-10, reltol=1e-10, H_builder=first ∘ build_H_S1, H_builder_kwargs=nothing,
        bare_to_dressed_tracking=:overlap, showprogress=true)

Fig. 2(b)/(c): accumulation speed of the conditional phase, vϕ = ϕ/τ (Eq. 3, paper sign convention vϕ = 4 g_zz),
from the dynamical CZ experiment on a grid of flux detunings `δϕA_list`, `δϕB_list` of the two qubits. Defaults
follow the paper: τ = 200 ns square pulses with tanh ramps of rise time t_r = t_ramp/6 = 2 ns. Each point uses
[`get_νφ_point`](@ref).

# Arguments
- `fA_params`, `fB_params`: `fluxonium` parameters of the two qubits, including their `name` and the flux offset
        `δϕ` (π at the sweet spot), see [`run_spectography_flux_sweep_2a`](@ref).
- `δϕA_list`, `δϕB_list`: detunings from the sweet spot, as flux angles 2πδΦ/Φ0.
- `J`: inductive coupling, passed to `H_builder`, see [`build_H_S1`](@ref).
- `H_builder`, `H_builder_kwargs`: build the composite system from the (shifted) parameters, see
        [`run_spectography_flux_sweep_2a`](@ref). The pulses act through [`dynamical_flux_terms`](@ref), so the
        external flux must enter through the linear inductive terms of the Hamiltonian.
- `showprogress`: show a progress bar over the grid.
- The remaining arguments are described in [`cz_gate_matrix`](@ref).

# Returns
- `(vϕ [rad/µs], max non-adiabatic error, branch_err [rad], avg. leakage)`, each indexed `[i_A, j_B]`.
"""
function run_cz_flux_sweep_2b(
    fA_params::Union{Dict, NamedTuple},
    fB_params::Union{Dict, NamedTuple},
    δϕA_list::AbstractVector{<:Real},
    δϕB_list::AbstractVector{<:Real},
    J::Real;
    τ::Real=200.0,
    t_ramp::Real=12.0,
    ramp_profile::Symbol=:tanh,
    t_rise::Maybe{Real}=nothing,
    abstol::Real=1e-10,
    reltol::Real=1e-10,
    H_builder::HFUNC=(first ∘ build_H_S1),
    H_builder_kwargs::HPP=nothing,
    bare_to_dressed_tracking::Symbol=:overlap,
    showprogress::Bool=true,
    ) where {HFUNC<:Function, HPP}
    sys_idle = H_builder(fA_params, fB_params, J; something(H_builder_kwargs, (;))...)
    # fill the dressed-state lookup cache before threading, so threads only read it
    foreach(lbl -> state(sys_idle, lbl; bare_to_dressed_tracking=bare_to_dressed_tracking), (:gg, :ge, :eg, :ee))
    vϕ_mtrx = Matrix{Float64}(undef, length(δϕA_list), length(δϕB_list))
    err_mtrx, branch_mtrx, leak_mtrx = similar(vϕ_mtrx), similar(vϕ_mtrx), similar(vϕ_mtrx)
    pbar = Progress(length(vϕ_mtrx); desc="CZ flux sweep: ", output=stdout, showspeed=true, enabled=showprogress)
    Threads.@threads for idx in CartesianIndices(vϕ_mtrx)
        i, j = Tuple(idx)
        vϕ_mtrx[i, j], err_mtrx[i, j], branch_mtrx[i, j], leak_mtrx[i, j] = get_νφ_point(
            fA_params, fB_params, sys_idle, τ, t_ramp, δϕA_list[i], δϕB_list[j], J;
            ramp_profile=ramp_profile, t_rise=t_rise, abstol=abstol, reltol=reltol, H_builder=H_builder,
            H_builder_kwargs=H_builder_kwargs, bare_to_dressed_tracking=bare_to_dressed_tracking,
        )
        next!(pbar)
    end
    finish!(pbar)
    return vϕ_mtrx, err_mtrx, branch_mtrx, leak_mtrx
end
