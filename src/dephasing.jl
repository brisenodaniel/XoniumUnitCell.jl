export thermal_energy, qubit_pink_noise, instantaneous_flux_expectations, flux_harmonic_differentials
export run_flux_modulation_sweep, floquet_instantaneous_check, dephasing_rates

"""
    thermal_energy(T)

Thermal energy `k_B T / ħ` as an angular frequency in rad/ns, for a temperature `T` in kelvin, with the rounded
constants `k_B = 1.38e-23 J/K`, `h = 6.63e-34 J s` used in the reference implementation. `thermal_energy(15e-3)`
equals the default `kT` of [`S_pink_noise`](@ref).
"""
thermal_energy(T::Real) = 2π * 1.38e-23 * T / 6.63e-34 * 1e-9

"""
    qubit_pink_noise(params; kT=thermal_energy(1.5e-3), qho_basis_nlev=50, ebasis_nlev=6, bare_to_dressed_tracking=:overlap, kwargs...)

[`PinkNoise`](@ref) of a single fluxonium with `fluxonium` parameters `params` (including its flux offset `δϕ`, π
at the sweet spot), built from its inductive and charging energies and the flux matrix element `φ_ge = ⟨g|φ|e⟩`.
`kT` is the thermal energy in rad/ns, by default that of 1.5 mK, see [`thermal_energy`](@ref). `kwargs` are passed to
[`noise_amplitudes`](@ref).
"""
function qubit_pink_noise(
    params::Union{Dict, NamedTuple};
    kT::Real=thermal_energy(1.5e-3),
    qho_basis_nlev::Int=50,
    ebasis_nlev::Int=6,
    bare_to_dressed_tracking::Symbol=:overlap,
    kwargs...
    )
    qubit = fluxonium(; params..., nlev=qho_basis_nlev) |>
        to_energy_basis |>
        sys -> truncate_sys(sys, ebasis_nlev)
    _, g = state(qubit, :g; bare_to_dressed_tracking=bare_to_dressed_tracking)
    _, e = state(qubit, :e; bare_to_dressed_tracking=bare_to_dressed_tracking)
    return PinkNoise(; El=params[:EL], Ec=params[:EC], φ_ge=g' * qubit.φ * e, kT=kT, kwargs...)
end

"""
    instantaneous_flux_expectations(Φ, H0, V, φ_op, ψ0, θs)

Follow the instantaneous eigenstates of `H0 + 2πΦ cos(θ) V` connected to the columns of `ψ0` (e.g. the idle dressed
states gg, ge, eg, ee) over one drive period, by maximal overlap between neighbouring samples. `H0`, `V`, `φ_op` are
dense matrices. Tracking starts at θ = π/2, where the detuning is 0. For drives far slower than the qubit
frequencies the Floquet states follow these eigenstates, and `⟨φ_op⟩` depends on `t` only through θ = ωt.

# Arguments
- `Φ`: modulation amplitude in units of Φ0.
- `H0`, `V`: static Hamiltonian and drive operator, see [`flux_drive_operators`](@ref).
- `φ_op`: operator whose expectation value is recorded.
- `ψ0`: matrix with the tracked initial states as columns.
- `θs`: equally spaced drive phases over one period, without the endpoint. The tracking needs them dense enough
        that the eigenvectors change little between samples.

# Returns
- `(φt, ψt)`: the `length(θs) × size(ψ0, 2)` matrix of `⟨φ_op⟩` and the tracked eigenvectors, both indexed by `θs`.
"""
function instantaneous_flux_expectations(
    Φ::Real,
    H0::AbstractMatrix,
    V::AbstractMatrix,
    φ_op::AbstractMatrix,
    ψ0::AbstractMatrix,
    θs::AbstractVector{<:Real},
    )
    nθ = length(θs)
    φt = Matrix{Float64}(undef, nθ, size(ψ0, 2))
    ψt = Vector{Matrix{ComplexF64}}(undef, nθ)
    ψprev = ψ0
    for k in circshift(1:nθ, -(nθ ÷ 4))
        Vecs = eigen(Hermitian(H0 + 2π * Φ * cos(θs[k]) * V)).vectors
        ovlp = abs2.(ψprev' * Vecs)
        idx = [argmax(ovlp[l, :]) for l in axes(ovlp, 1)]
        allunique(idx) || @warn "Ambiguous eigenstate tracking at Φ = $Φ, θ = $(θs[k]); increase nθ"
        ψprev = Vecs[:, idx]
        ψt[k] = ψprev
        φt[k, :] = real.(diag(ψprev' * φ_op * ψprev))
    end
    return φt, ψt
end

"""
    flux_harmonic_differentials(φt, θs; pairs=((1, 3), (2, 4)))

Period average and first-harmonic amplitude of `⟨φ⟩(θ)` for the columns of `φt` (ordered gg, ge, eg, ee), combined
as the mean over `pairs = ((a, b), …)` of `|φ[b] - φ[a]|`. The default pairs are the qubit-A transitions gg→eg and
ge→ee. The first harmonic is the complex amplitude φ1 in `⟨φ⟩(θ) ≈ φ0 + Re(φ1 e^{iθ}) + …`.

# Returns
- `(combined period average, combined first-harmonic amplitude)`.
"""
function flux_harmonic_differentials(φt::AbstractMatrix, θs::AbstractVector{<:Real}; pairs=((1, 3), (2, 4)))
    nθ = length(θs)
    φ0 = vec(sum(φt; dims=1)) ./ nθ
    φ1 = vec(sum(φt .* cis.(-θs); dims=1)) .* (2 / nθ)
    comb(φ) = sum(abs(φ[b] - φ[a]) for (a, b) in pairs) / length(pairs)
    return comb(φ0), comb(φ1)
end

"""
    run_flux_modulation_sweep(fA_params, fB_params, Φlist, J; δϕ_offsets=(0.005, -0.005), nθ=400, H_builder=first ∘ build_H_S1,
        H_builder_kwargs=nothing, bare_to_dressed_tracking=:overlap, showprogress=true)

Flux modulation of qubit A, `δϕA(t) = 2πΦ cos(ωt)` with qubit B at its offset, for each amplitude `Φ` in `Φlist`
(units of Φ0), in the adiabatic limit where the Floquet states follow the instantaneous eigenstates. For each `Φ`
this records `⟨φA⟩` along the eigenstates connected to gg, ge, eg, ee over one period, see
[`instantaneous_flux_expectations`](@ref), and its period average and first harmonic, see
[`flux_harmonic_differentials`](@ref). All of these are independent of the drive frequency.

# Arguments
- `fA_params`, `fB_params`: `fluxonium` parameters of the qubits, including `name` and the flux offset `δϕ`, see
        [`run_spectography_flux_sweep_2a`](@ref).
- `Φlist`: modulation amplitudes in units of Φ0.
- `J`: inductive coupling, passed to `H_builder`, see [`build_H_S1`](@ref).
- `δϕ_offsets`: shifts of `δϕ` of qubits A and B (rad) for the coupled system that is modulated, e.g. to lift
        the symmetry of the sweet spot.
- `nθ`: number of drive phases per period.
- `H_builder`, `H_builder_kwargs`: build the composite system from the shifted parameters, see
        [`run_spectography_flux_sweep_2a`](@ref).
- `bare_to_dressed_tracking`, `showprogress`: see [`state`](@ref) and [`run_spectography_flux_sweep_2a`](@ref).

# Returns
- NamedTuple `(sys, θs, φt_sweep, ∂0, ∂1)`: the coupled system, the drive phases, the vector over `Φlist` of the
        `nθ × 4` matrices of `⟨φA⟩`, and the vectors over `Φlist` of the period average and first harmonic.
"""
function run_flux_modulation_sweep(
    fA_params::Union{Dict, NamedTuple},
    fB_params::Union{Dict, NamedTuple},
    Φlist::AbstractVector{<:Real},
    J::Real;
    δϕ_offsets::Tuple{Real, Real}=(0.005, -0.005),
    nθ::Int=400,
    H_builder::HFUNC=(first ∘ build_H_S1),
    H_builder_kwargs::HPP=nothing,
    bare_to_dressed_tracking::Symbol=:overlap,
    showprogress::Bool=true,
    ) where {HFUNC<:Function, HPP}
    fA_lbl, fB_lbl = fA_params[:name], fB_params[:name]
    pA = (; fA_params..., δϕ=fA_params[:δϕ] + δϕ_offsets[1])
    pB = (; fB_params..., δϕ=fB_params[:δϕ] + δϕ_offsets[2])
    sys = H_builder(pA, pB, J; something(H_builder_kwargs, (;))...)
    ΦA_op, _ = flux_drive_operators(sys, J; fA_lbl=fA_lbl, fB_lbl=fB_lbl)
    H0, V, φ_op = Matrix(sys.H.data), Matrix(ΦA_op.data), Matrix(sys[fA_lbl, :φ].data)
    # fill the dressed-state lookup cache before threading, so threads only read it
    ψ0 = reduce(hcat, (state(sys, lbl; bare_to_dressed_tracking=bare_to_dressed_tracking)[2].data
                       for lbl in (:gg, :ge, :eg, :ee)))
    θs = 2π .* (0:nθ-1) ./ nθ
    φt_sweep = Vector{Matrix{Float64}}(undef, length(Φlist))
    ∂0, ∂1 = zeros(length(Φlist)), zeros(length(Φlist))
    pbar = Progress(length(Φlist); desc="Flux modulation sweep: ", output=stdout, showspeed=true, enabled=showprogress)
    Threads.@threads for j in eachindex(Φlist)
        φt, _ = instantaneous_flux_expectations(Φlist[j], H0, V, φ_op, ψ0, θs)
        φt_sweep[j] = φt
        ∂0[j], ∂1[j] = flux_harmonic_differentials(φt, θs)
        next!(pbar)
    end
    finish!(pbar)
    return (; sys, θs, φt_sweep, ∂0, ∂1)
end

"""
    floquet_instantaneous_check(sys, J, ω, Φ; nθ=400, nT=20, δλ=0.5, fA_lbl=:fA, fB_lbl=:fB,
        bare_to_dressed_tracking=:overlap, floquet_kwargs...)

Compare the instantaneous eigenstates used by [`run_flux_modulation_sweep`](@ref) with the tracked Floquet states of
`sys` under the modulation `δϕA(t) = 2πΦ cos(ωt)`, for drive frequency `ω` (rad/ns) and amplitude `Φ` (units of Φ0).
This is best done at the largest amplitude and fastest drive, where non-adiabatic corrections are largest.

The Floquet states are sampled at `nT` times over one period (`nθ` must be a multiple of `nT`), and are tracked
with the arguments `δλ` and `floquet_kwargs` of `FloquetBasis`, e.g. `n_basis_states=10`.

# Returns
- NamedTuple with, for each of the states gg, ge, eg, ee: `min_overlap`, the smallest `|⟨inst|Floquet⟩|²`, and
        `max_Δφ`, the largest `|Δ⟨φA⟩|`; and `Floquet` and `instantaneous` NamedTuples `(∂0, ∂1)` with the
        differentials of [`flux_harmonic_differentials`](@ref) from the `nT` samples.
"""
function floquet_instantaneous_check(
    sys::CompositeQuantumSystem,
    J::Real,
    ω::Real,
    Φ::Real;
    nθ::Int=400,
    nT::Int=20,
    δλ::Real=0.5,
    fA_lbl::Symbol=:fA,
    fB_lbl::Symbol=:fB,
    bare_to_dressed_tracking::Symbol=:overlap,
    floquet_kwargs...
    )
    nθ % nT == 0 || throw(ArgumentError("`nθ` must be a multiple of `nT`"))
    lbls = (:gg, :ge, :eg, :ee)
    ΦA_op, ΦB_op = flux_drive_operators(sys, J; fA_lbl=fA_lbl, fB_lbl=fB_lbl)
    φA = sys[fA_lbl, :φ]
    T = 2π / ω
    ts = T .* (0:nT-1) ./ nT
    drA = (p, t) -> 2π * Φ * cos(ω * t)
    drB = (p, t) -> 0.0
    fb = FloquetBasis(sys.H, [(drA, ΦA_op), (drB, ΦB_op)], T, ts; δλ=δλ, nT=nT, progress_bar=Val(false), floquet_kwargs...)
    # instantaneous eigenstates at the same drive phases θ = ωt
    θs = 2π .* (0:nθ-1) ./ nθ
    ψ0 = reduce(hcat, (state(sys, lbl; bare_to_dressed_tracking=bare_to_dressed_tracking)[2].data for lbl in lbls))
    φt_inst, ψt_inst = instantaneous_flux_expectations(Φ, Matrix(sys.H.data), Matrix(ΦA_op.data),
                                                       Matrix(φA.data), ψ0, θs)
    ks = 1 .+ (0:nT-1) .* (nθ ÷ nT)
    φt_fl = Matrix{Float64}(undef, nT, length(lbls))
    ovlp = similar(φt_fl)
    for (n, t) in enumerate(ts), (l, lbl) in enumerate(lbls)
        ψ = state(sys, fb, lbl, t; bare_to_dressed_tracking=bare_to_dressed_tracking)[2]
        φt_fl[n, l] = real(expect(φA, ψ))
        ovlp[n, l] = abs2(dot(ψt_inst[ks[n]][:, l], ψ.data))
    end
    fl = flux_harmonic_differentials(φt_fl, θs[ks])
    inst = flux_harmonic_differentials(φt_inst[ks, :], θs[ks])
    return (;
        min_overlap=NamedTuple{lbls}(Tuple(minimum(ovlp[:, l]) for l in eachindex(lbls))),
        max_Δφ=NamedTuple{lbls}(Tuple(maximum(abs, φt_fl[:, l] .- φt_inst[ks, l]) for l in eachindex(lbls))),
        Floquet=(; ∂0=fl[1], ∂1=fl[2]),
        instantaneous=(; ∂0=inst[1], ∂1=inst[2]),
    )
end

"""
    dephasing_rates(φt_sweep, ω_list, S; pairs=((1, 3), (2, 4)))

Estimated dephasing rate Γφ of qubit A due to the noise PSD `S` (e.g. a [`PinkNoise`](@ref)) for each drive
frequency `ω` in `ω_list` (rad/ns; a vector or tuple) and each modulation amplitude of `φt_sweep` from
[`run_flux_modulation_sweep`](@ref). For each `(a, b)` in `pairs` (columns of `φt` ordered gg, ge, eg, ee) the
weights are those of [`calc_gs`](@ref) with `b` the upper and `a` the lower state, and the rates of
[`calc_Γφ`](@ref) are averaged. The default pairs are the qubit-A transitions gg→eg and ge→ee. The weights
do not depend on `ω`, which enters only through `S(kω)`.

# Returns
- Matrix of rates in 1/ns, indexed `[i_ω, j_Φ]`.
"""
function dephasing_rates(φt_sweep::AbstractVector{<:AbstractMatrix}, ω_list, S; pairs=((1, 3), (2, 4)))
    Γφ = Matrix{Float64}(undef, length(ω_list), length(φt_sweep))
    for (j, φt) in enumerate(φt_sweep), (i, ω) in enumerate(ω_list)
        Γφ[i, j] = sum(calc_Γφ(φt[:, b], φt[:, a], S, ω) for (a, b) in pairs) / length(pairs)
    end
    return Γφ
end
