# XoniumUnitCell.jl API

Reference for the functions and types defined in `src/`. For setup and the notebooks, see [README.md](README.md).

```julia
using XoniumUnitCell    # exports everything listed below, and re-exports QuantumSystems
```

## Conventions

- **Units.** Energies and frequencies are in GHz, times in ns, and ħ = 1. Many functions take **angular** quantities (`2π·…`): `EC`, `EL`, `EJ`, `J`, drive frequencies `ω`, and `kT` are all in rad/ns where a function says so. Convert with [`tounit_ħ`](#tounit_ħ--tounit_h). Check each docstring.
- **Flux.** Flux offsets and detunings (`δϕ`, `δϕA`, `Φ`, …) are **flux angles** `2πδΦ/Φ0`, except where a function says "units of Φ0". `δϕ = π` is the sweet spot.
- **Subsystem names.** Composite systems are addressed by subsystem `name` (`:fA`, `:fB` by default), e.g. `sys[:fA, :φ]`. Parameter sets passed to the sweep functions must contain `name`.
- **Dressed-state tracking.** Functions with `bare_to_dressed_tracking` / `state_tracking` choose how bare labels (`:g, :e`, `:gg, :ge, :eg, :ee`) are mapped to dressed states: `:overlap` (default) or `:spectrum`.
- **Coupling sign.** TODO. See `notebooks/parametric_CZ_notes.md`.
- **Parameter containers.** `fA_params`, `fB_params` are `NamedTuple`s (or `Dict`s) of `fluxonium` keywords, e.g. `(EC=…, EL=…, EJ=…, δϕ=π, name=:fA)`.

## Contents

1. [QuantumSystems (re-exported)](#quantumsystems-re-exported)
2. [Circuit Hamiltonians](#circuit-hamiltonians) (`circuits.jl`)
3. [Circuit graphs](#circuit-graphs) (`circuit_graphs.jl`)
4. [Unit conversions](#unit-conversions) (`unit_conversions.jl`)
5. [Coupled-system construction and sweeps](#coupled-system-construction-and-sweeps) (`gate_utils.jl`)
6. [Spectroscopy (Fig. 2a)](#spectroscopy-fig-2a)
7. [CZ gate from flux pulses (Fig. 2b)](#cz-gate-from-flux-pulses-fig-2b)
8. [Pulse envelopes](#pulse-envelopes) (`pulses.jl`)
9. [Noise models and dephasing rates](#noise-models-and-dephasing-rates) (`env_psd.jl`)
10. [Flux modulation and dephasing](#flux-modulation-and-dephasing) (`dephasing.jl`)
11. [Plotting](#plotting) (`plotting.jl`)
12. [Not exported](#not-exported)

---

## QuantumSystems (re-exported)

`XoniumUnitCell` re-exports the names of [`QuantumSystems.jl`](https://github.com/brisenodaniel/QuantumSystems.jl) (`SimpleQuantumSystem`, `CompositeQuantumSystem`, `⊗`, `state`, `energy`, `to_energy_basis`, `truncate_sys`, `Maybe`, …). TODO: document, or link to that package's docs.

---

## Circuit Hamiltonians

Each constructor returns a `SimpleQuantumSystem` with named operators and a `params` record.

### `fluxonium`

```julia
fluxonium(; EC, EL, EJ, δϕ=π, n_g=0.0, nlev=40, Φext_incos=false, name=:Fluxonium)
```

Fluxonium in the harmonic-oscillator basis, with

    H = 4 EC (n - n_g)² + ½ EL (φ + δϕ)² - EJ cos φ

If `Φext_incos=true` the external flux is moved into the cosine instead: `½ EL φ² - EJ cos(φ - δϕ)`.

- `EC`, `EL`, `EJ`: charging, inductive and Josephson energies. Use angular units (`2π·GHz`) for dynamics in ns.
- `δϕ`: external flux angle (default `π`, the sweet spot).
- `n_g`: offset charge. `nlev`: dimension of the oscillator basis.
- **Returns** a system with operators `φ`, `n`, `H` and params `(EC, EL, EJ, δϕ, Φext_incos)`.

```julia
fl = fluxonium(; EC=2π*1.0, EL=2π*1.0, EJ=2π*4.5, nlev=30)
```

### `transmon`

```julia
transmon(; EC, EJ, ng=0.0, ncut=20, name=:Transmon)
```

Transmon in the charge basis, `H = 4 EC (n - ng)² - EJ cos φ`, with charge states `-ncut:ncut` (dimension `2·ncut + 1`).

- **Returns** a system with operators `n`, `cosφ`, `H` and params `(EC, EJ, ng)`.

### `split_transmon`

```julia
split_transmon(; EC, EJ1, EJ2, ϕext=0.0, ng=0.0, ncut=20, name=:Transmon)
```

TODO: Hamiltonian.

- **Returns** a system with operators `n`, `cosφ`, `H` and params `(EC, EJ1, EJ2, ϕext, ng, ncut)`.

### `φn_operators`

```julia
φn_operators(; EC, EJ=nothing, EL=nothing, nlev=10) -> (φ, n)
```

Flux and charge operators in the oscillator basis, `φ = φ_zpf·x`, `n = n_zpf·p`, with `φ_zpf = (8Z)^¼`, `n_zpf = (8Z)^-¼` and `Z = EC/EL` (or `EC/EJ` if `EL` is `nothing`). At least one of `EJ`, `EL` is required, otherwise an `ArgumentError` is thrown.

### `tounit_ħ` / `tounit_h`

```julia
tounit_ħ(params, convert_keys=(:EC, :EJ, :EL))
tounit_h(params, convert_keys=(:EC, :EJ, :EL))
```

Multiply (`tounit_ħ`, GHz → rad/ns) or divide (`tounit_h`) the entries of `params` whose key is in `convert_keys` by `2π`. Other entries are unchanged.

**Returns a tuple of `key => value` pairs**, not a `NamedTuple`. Wrap it to get one:

```julia
p = NamedTuple(tounit_ħ((EC=1.0, EJ=4.5, EL=1.0)))   # (EC = 2π, EJ = 9π, EL = 2π)
```

---

## Circuit graphs

A small toolkit for circuit quantization: describe a circuit as nodes and edges, get its capacitance and inverse-inductance matrices, change basis, and convert matrix elements to circuit parameters with [Unit conversions](#unit-conversions). Capacitances are in fF and inductances in nH. Values may be numeric or symbolic (`Symbolics.Num`).

### `CNode`

```julia
CNode(label::Symbol, Cg; Lg=nothing)
```

A node with capacitance `Cg` (fF) to ground and an optional inductor `Lg` (nH) to ground. `Lg` must be nonzero.

### `CEdge`

```julia
CEdge(a::Symbol, b::Symbol; C=nothing, L=nothing)
CEdge((a, b); C=nothing, L=nothing)
```

A branch between nodes `a` and `b` with an optional capacitor `C` (fF) and an optional inductor `L` (nH) in parallel. At least one is required, and `L` must be nonzero.

### `CGraph`

```julia
CGraph(nodes, edges=[]; tol=1e-9, simplify=identity)
CGraph(labels, C, L_inv; tol=1e-9, simplify=identity)
```

A circuit with fields `nodes`, `edges`, `C` (capacitance matrix), `L_inv` (inverse-inductance matrix), `tol` and `simplify`. The matrices are stored in full, in node order.

- **From nodes and edges:** each node adds `Cg` and `1/Lg` to its diagonal entry, and each edge adds `C` and `1/L` to both diagonals and subtracts it from the two off-diagonals. Parallel edges are merged, with capacitances adding and inductances combining as `1/L = 1/L₁ + 1/L₂`. Node labels must be unique, and self-loops, unknown nodes and empty edges are rejected.
- **From matrices:** nodes and edges are read back from `C` and `L_inv`. Edge values are the negated off-diagonals and ground values are the row sums. After a change of basis the result can contain negative capacitances and inductances, so it need not be a physical circuit.
- **Symmetry:** `C` and `L_inv` must be symmetric to within `tol`, otherwise an `ArgumentError` is thrown.
- **`simplify`:** applied to symbolic quantities only, never to numbers. A good choice for symbolic work is `x -> Symbolics.simplify(x; expand=true)`. It is kept by the graph and reused by the operations below.

```julia
g = CGraph([CNode(:a, 10.0; Lg=5.0), CNode(:b, 10.0)], [CEdge(:a, :b; C=2.0)])
g.C        # [12.0 -2.0; -2.0 12.0]
```

**Operations**

| Expression | Result |
|---|---|
| `M * g * Q` | New `CGraph` with `C → M*C*Q`, `L_inv → M*L_inv*Q`, for square matrices `M`, `Q`. A change of basis `φ = Vψ` is `V' * g * V`. |
| `g * v` | A `CVector` holding `(g.C * v, g.L_inv * v)`. |
| `v' * g` | A `CAdjoint` holding `(v' * g.C, v' * g.L_inv)`. |
| `a' * g * b` | `(C = a' * g.C * b, L = a' * g.L_inv * b)`. |
| `g1 == g2` | Same labels and matrices (structural equality for symbolic entries). |

### `CVector`, `CAdjoint`

```julia
CVector(C, L_inv)
CAdjoint(C, L_inv)
```

Wrappers for the pair of products `g * v` and `v' * g`. They multiply with plain vectors and adjoints like the pair they hold and return a `NamedTuple` `(C = …, L = …)`: `u' * x` and `x * w` give pairs of scalars, `x * u'` and `w * x` give pairs of matrices. The adjoint of a `CVector` is a `CAdjoint` and vice versa.

### `basis_vectors`

```julia
basis_vectors(graph) -> NamedTuple
```

Unit vectors keyed by node label, in node order: `basis_vectors(g).A1` has a `1` at the position of `:A1`. Unpack with `(; A1, B2) = basis_vectors(g)`.

### `relabel`

```julia
relabel(graph, labels::AbstractVector{Symbol}) -> CGraph
```

Same matrices with the `k`-th node renamed to `labels[k]`. Use it to name the new variables after `V' * g * V`.

---

## Unit conversions

Units: capacitances in fF, inductances in nH, critical currents in nA, energies as E/h in GHz, couplings in MHz.

| Function | Converts | To |
|---|---|---|
| `CPW_freq2cap(freq; Z0=50)` | frequency `freq` (GHz) and impedance `Z0` (Ω) of a λ/4 CPW resonator | capacitance to ground (fF) |
| `Cinv2E_C(Cinv)` | inverse-capacitance diagonal element (1/fF) | charging energy `E_C/h` (GHz), `E_C = e²/2·Cinv` |
| `Cinv2J(Cinv)` | inverse-capacitance off-diagonal element (1/fF) | qubit-qubit coupling `J/h` (MHz), `J = 4e²·Cinv` |
| `Cinv2g(Cinv)` | inverse-capacitance off-diagonal element (1/fF) | qubit-resonator coupling `g/2π` (MHz) |
| `Linv2E_L(Linv)` | inverse-inductance element (1/nH) | inductive energy `E_L/h` (GHz), `E_L = (ħ/2e)²·Linv` |
| `L2E_L(L)` | inductance (nH) | inductive energy `E_L/h` (GHz) |
| `Ic2E_J(Ic)` | critical current (nA) | Josephson energy `E_J/h` (GHz) |
| `E_J2Ic(EJ)` | Josephson energy `E_J/h` (GHz) | critical current (nA) |

```julia
Cinv2E_C(1/100)    # charging energy for 100 fF, ≈ 0.194 GHz
Ic2E_J(30.0)       # ≈ 14.9 GHz
```

---

## Coupled-system construction and sweeps

### `build_H_S1`

```julia
build_H_S1(fA_params, fB_params, J; qho_basis_nlev=50, ebasis_nlev=6) -> (sys, H_int)
```

Two fluxonia coupled inductively, `-J(φA + δϕA)(φB + δϕB)`. Each fluxonium is built with `fluxonium(; params..., nlev=qho_basis_nlev)`, converted to its energy basis and truncated to `ebasis_nlev` levels. The single-qubit cross terms `J δϕB φA` and `J δϕA φB` are moved into the (tilted) fluxonia, and `H_int` is what remains.

- `fA_params`, `fB_params`: `fluxonium` keywords. They must include `name`, and the flux offset `δϕ` of each qubit is read from them.
- **Returns** `(sys, H_int)`: the composite system including the coupling, and the interaction Hamiltonian.

### `build_H_LC`

```julia
build_H_LC(fA_params, fB_params, J; J_C=0.0, kwargs...) -> (sys, H_int)
```

[`build_H_S1`](#build_h_s1) plus a capacitive coupling `+J_C nA nB`, i.e. an inductive and capacitive coupling `-J φφ + J_C nn`. This is the joint sign flip of the `+g_L φφ - g_C nn` convention of `2Q.ipynb`, which leaves ζzz unchanged, with `J = g_L` and `J_C = g_C`. `kwargs` go to `build_H_S1`. `H_int` includes the capacitive term. Where a sweep needs only the system, use `first ∘ build_H_LC` and pass `J_C` via `H_builder_kwargs=(J_C=…,)`.

### `static_flux_coupling`

```julia
static_flux_coupling(sys, δϕA, δϕB, J; fA_lbl=:fA, fB_lbl=:fB)
```

The operator `J (φA + δϕA)(φB + δϕB)` on the composite system `sys`. (Note the sign: `build_H_S1` negates it.)

### `get_qubit_subspace_couplings`

```julia
get_qubit_subspace_couplings(sys, H_int; bare_to_dressed_tracking=:overlap) -> (xx, xz, zx, zz)
```

Projects `H_int` onto the Pauli products of the two qubits' dressed `g`/`e` states, `σσ' = tr(σ⊗σ' · H_int)/4`, for `σ, σ' ∈ {x, z}`. `sys` must be a composite of two subsystems.

### `sweep_couplings`

```julia
sweep_couplings(δϕs, fA_params, fB_params, H_builder, J;
                qho_basis_nlev=50, ebasis_nlev=6, bare_to_dressed_tracking=:overlap) -> StructArray
```

Shift the flux offset of **both** fluxonia by each `δϕ` in `δϕs` (relative to their `δϕ`, default `π`), build the system with `H_builder`, and compute the subspace couplings.

- `H_builder(pA, pB, J; qho_basis_nlev, ebasis_nlev)` must return a `(sys, H_int)` tuple, e.g. [`build_H_S1`](#build_h_s1). (The builders of the sweeps below must return the system only.)
- **Returns** a `StructArray` with fields `xx`, `xz`, `zx`, `zz`, one entry per `δϕ`.

### `flux_drive_operators`

```julia
flux_drive_operators(sys, J; fA_lbl=:fA, fB_lbl=:fB) -> (ΦA_op, ΦB_op)
```

The operators `∂H/∂δϕA` and `∂H/∂δϕB` through which flux detunings enter the Hamiltonian:

    ΦA_op = EL_A φA - J φB,    ΦB_op = EL_B φB - J φA

so the Hamiltonian is `sys.H + δϕA ΦA_op + δϕB ΦB_op` up to constants. Exact when the external flux enters through the linear inductive terms.

---

## Spectroscopy (Fig. 2a)

### `run_spectography_circuit`

```julia
run_spectography_circuit(sys, H_dr, tlist; dr_params=nothing, integrator_pbar=false,
                         state_tracking=:overlap, fA_lbl=:fA, sesolve_kwargs...) -> Float64
```

Evolve `(|gg⟩ + |ge⟩)/√2` under `sys.H + H_dr` with `sesolve` and return the population transferred to qubit A's excited state, `|⟨eg|ψ⟩|² + |⟨ee|ψ⟩|²`, at `tlist[end]`.

- `sys`: system whose dressed states are labelled `:gg, :ge, :eg, :ee`.
- `H_dr`: drive Hamiltonian, a `QobjEvo` whose coefficient functions take `(p, t)`.
- `dr_params`: the `p` passed to those functions.
- `integrator_pbar`: per-call progress bar for `sesolve`. Leave it `false` inside sweeps.
- `fA_lbl` is currently unused. Extra keywords (`abstol`, `reltol`, …) go to `sesolve`.

### `run_spectography_flux_sweep_2a`

```julia
run_spectography_flux_sweep_2a(fA_params, fB_params, tlist, env, δϕB_list, ω_list, J, Ω_ge;
    H_builder=first ∘ build_H_S1, H_builder_kwargs=nothing, dr_params=nothing, δϕA=-2π*0.042,
    integrator_pbar=false, showprogress=true, state_tracking=:overlap, sesolve_kwargs...) -> Matrix
```

Spectroscopy of Fig. 2(a): drive qubit A's g↔e transition with `Ω_ge · nA · env(p, t) · cos(ω t)` and record the transferred population on a grid of drive frequencies and flux detunings of qubit B. **Multithreaded** over `δϕB_list`.

| Argument | Meaning |
|---|---|
| `fA_params`, `fB_params` | `fluxonium` parameters including `name` and `δϕ`. `δϕ` is shifted by `δϕA` (qubit A) and the entries of `δϕB_list` (qubit B). |
| `tlist` | `[0, T]`. The population is evaluated at `T`. |
| `env` | Envelope `env(p, t)` without the carrier, e.g. wrapping [`cosine_ramp_envelope`](#cosine_ramp_envelope). |
| `δϕB_list` | Detunings of qubit B from its sweet spot (flux angles). |
| `ω_list` | Angular drive frequencies. `ω` is added to `p` (overriding any `ω` in `dr_params`), so `env` may read `p.ω`. |
| `J` | Inductive coupling passed to `H_builder`. |
| `Ω_ge` | Drive strength multiplying `nA`. |
| `H_builder` | `H_builder(pA, pB, J; H_builder_kwargs...)` returning the system **only** (not a tuple). |
| `H_builder_kwargs` | `NamedTuple` of keywords for `H_builder`, e.g. `(qho_basis_nlev=50, ebasis_nlev=6)`. |
| `dr_params` | `NamedTuple` of extra parameters for `env`. |
| `δϕA` | Detuning of qubit A from its sweet spot. |

**Returns** a matrix of transferred populations indexed `[i_ω, j_δϕB]`.

---

## CZ gate from flux pulses (Fig. 2b)

### `cz_gate_matrix`

```julia
cz_gate_matrix(sys, τ, t_ramp, ΩA, ΩB, J; ω=0.0, ramp_profile=:tanh, t_rise=nothing,
               abstol=1e-10, reltol=1e-10, bare_to_dressed_tracking=:overlap,
               fA_lbl=:fA, fB_lbl=:fB) -> (C, sts)
```

Evolve the dressed computational states `gg, ge, eg, ee` of the idle system `sys` under one flux pulse of total length `τ`. Qubit A (B) is pulsed by `ΩA·f(t)` (`ΩB·f(t)`), with `f = full_pulse(ω, t_ramp, τ - 2·t_ramp; ramp_profile, t_rise)`. The two drives share the pulse and are merged into one time-dependent term.

- `ΩA`, `ΩB`: flux amplitudes (flux angles). `ω = 0` gives the square pulses of Fig. 2(b).
- `abstol`, `reltol`: tight by default, since the lab-frame phase reaches ~10³ rad and looser values give phase errors comparable to the conditional phases.
- **Returns** `(C, sts)`: `C[l, r] = ⟨ψ_l|U(τ)|ψ_r⟩` for `ψ ∈ (gg, ge, eg, ee)`, and `sts = (; gg, ge, eg, ee)`, the idle dressed states.

### `subspace_leakage`

```julia
subspace_leakage(C) -> Float64
```

Average population leaving the computational subspace, `1 - Tr(C†C)/4`, for the matrix `C` of `cz_gate_matrix`.

### `get_νφ`

```julia
get_νφ(C, τ) -> (ϕ, νϕ, leak)
```

Conditional phase `ϕ = arg(c_ee c_gg / (c_eg c_ge))` with `c_k = C[k, k]`. Single-qubit phases cancel. `ϕ` is only known mod 2π.

- **Returns** `ϕ` (rad), `ϕ/τ` in rad/µs, and the per-state leakage `1 - |c_k|²` as a `NamedTuple (; gg, ge, eg, ee)`.

### `unwrap_phase`

```julia
unwrap_phase(ϕs) -> Vector{Float64}
```

Unwrap phases mod 2π into a continuous sequence, assuming successive entries differ by less than π.

### `static_ζ`

```julia
static_ζ(fA_params, fB_params, δϕA, δϕB, J; H_builder=first ∘ build_H_S1,
         H_builder_kwargs=nothing, bare_to_dressed_tracking=:overlap) -> Float64
```

Signed static ZZ rate `ζ = E_ee - E_eg - E_ge + E_gg` (rad/ns) with the qubits at flux detunings `δϕA`, `δϕB` from their `δϕ` offsets. `ζ = 4 g_zz`.

### `get_νφ_sweep`

```julia
get_νφ_sweep(sys, τs, t_ramp, ΩA, ΩB, J; ω=0.0, ramp_profile=:tanh, t_rise=nothing,
             abstol=1e-10, reltol=1e-10, bare_to_dressed_tracking=:overlap,
             fA_lbl=:fA, fB_lbl=:fB) -> (νϕ, slope, ϕs, leak)
```

Run `cz_gate_matrix` for each pulse length in `τs` and unwrap `ϕ(τ)`. Requires `|vϕ|·step(τs) < π`, `|ϕ(τs[1])| < π` (keep `τs[1]` short) and `τs[1] ≥ 2·t_ramp`.

- **Returns** `νϕ` (`ϕ(τ_max)/τ_max`, rad/µs), the plateau slope `dϕ/dτ` (rad/µs, from a linear fit with intercept), the unwrapped `ϕs`, and the largest leakage over all `τ` and states.

### `get_νφ_point`

```julia
get_νφ_point(fA_params, fB_params, sys_idle, τ, t_ramp, ΩA, ΩB, J; ramp_profile=:tanh,
             t_rise=nothing, abstol=1e-10, reltol=1e-10, H_builder=first ∘ build_H_S1,
             H_builder_kwargs=nothing, bare_to_dressed_tracking=:overlap)
    -> (vϕ, max_err, branch_err, leakage)
```

One Ramsey-type CZ experiment with square pulses (ω = 0). The 2π branch of `ϕ` is chosen closest to the static estimate `-ζ·(τ - t_ramp)` using [`static_ζ`](#static_ζ). `|branch_err|` (rad) flags points where that assumption fails, and `|branch_err| ≳ π/2` is unreliable.

- **Returns** `vϕ` in rad/µs (paper convention `vϕ = 4 g_zz = -ϕ/τ`), the largest non-adiabatic error over the computational states, `branch_err`, and the average leakage.

### `run_cz_flux_sweep_2b`

```julia
run_cz_flux_sweep_2b(fA_params, fB_params, δϕA_list, δϕB_list, J; τ=200.0, t_ramp=12.0,
    ramp_profile=:tanh, t_rise=nothing, abstol=1e-10, reltol=1e-10,
    H_builder=first ∘ build_H_S1, H_builder_kwargs=nothing,
    bare_to_dressed_tracking=:overlap, showprogress=true) -> (vϕ, err, branch_err, leakage)
```

Fig. 2(b)/(c): the conditional-phase accumulation speed `vϕ = ϕ/τ` on a grid of flux detunings of the two qubits, using [`get_νφ_point`](#get_νφ_point) for each point. Defaults follow the paper: `τ = 200` ns with tanh ramps of rise time `t_ramp/6 = 2` ns. **Multithreaded** over the grid. The idle system is built once, and its dressed-state cache is filled before threading.

- `H_builder` must return the system only. The pulses act through the linear response of `flux_drive_operators`, so the external flux must enter through the linear inductive terms.
- **Returns** four matrices indexed `[i_A, j_B]`: `vϕ` (rad/µs), the maximum non-adiabatic error, `branch_err` (rad) and the average leakage.

---

## Pulse envelopes

All envelope functions return a function of time `t`.

### `cosine_ramp_envelope`

```julia
cosine_ramp_envelope(t_ramp, t_plateau) -> f(t)
```

Flat-top envelope with raised-cosine ramps: 0 → 1 on `[0, t_ramp]`, 1 on the plateau, 1 → 0 on the last `t_ramp`, and 0 outside `[0, 2·t_ramp + t_plateau]`. The envelope is C¹ at every boundary.

### `tanh_ramp_envelope`

```julia
tanh_ramp_envelope(t_ramp, t_plateau; t_rise=t_ramp/6) -> f(t)
```

Flat-top envelope with tanh edges (Fig. 2(c) of the paper, `t_r = 2 ns`), with the same footprint `T = 2·t_ramp + t_plateau` as the cosine envelope and edges centered at `t_ramp/2` and `T - t_ramp/2`. It is shifted and rescaled so that it is exactly 0 at `t = 0, T` and exactly 1 at `T/2`.

### `full_pulse`

```julia
full_pulse(ω, t_ramp, t_plateau; phase=0.0, ramp_profile=:cosine, t_rise=nothing) -> f(t)
```

`f(t) = envelope(t)·cos(ω t + phase)`, where the envelope is `cosine_ramp_envelope` (`ramp_profile=:cosine`) or `tanh_ramp_envelope` (`ramp_profile=:tanh`, with `t_rise` defaulting to `t_ramp/6`). `ω` is an angular frequency, `phase` is in radians, and any other `ramp_profile` throws an `ArgumentError`. To match the paper's `t_r = 2 ns`, use `ramp_profile=:tanh, t_ramp=12.0`.

```julia
f = full_pulse(2π*0.1, 12.0, 176.0; ramp_profile=:tanh)
f(100.0)
```

---

## Noise models and dephasing rates

Pink-noise PSDs and the dephasing rate of a driven qubit (`env_psd.jl`). Frequencies are angular (rad/ns). The default thermal energy is `kT = KT_15mK` (T = 15 mK).

### `noise_amplitudes`

```julia
noise_amplitudes(El, Ec, φ_ge; tndC=1.1e-6, df=1.8e-6) -> (Ad, Af)
```

Amplitudes of the dielectric-loss and 1/f flux noise terms for a qubit with inductive energy `El`, charging energy `Ec` and flux matrix element `φ_ge = ⟨g|φ|e⟩`: `Af = df·El·|φ_ge|` and `Ad = π²·tndC·|φ_ge|²/Ec`.

### `S_pink_noise`

```julia
S_pink_noise(w, Ad, Af; kT=KT_15mK)
```

Noise power spectral density at angular frequency `w`:

    S(w) = Af²·T + α·Ad/T²,    T = |2π/w|,    α = |coth(w/2kT) + 1| / 2

1/f flux noise plus a thermally weighted dielectric-loss term. `w = 0` is not regularized (the result is `Inf` or `NaN`), so handle the DC component separately. Broadcast over arrays: `S_pink_noise.(ws, Ad, Af)`.

### `PinkNoise`

```julia
PinkNoise(Ad, Af; kT=KT_15mK)
PinkNoise(; El, Ec, φ_ge, kT=KT_15mK, tndC, df)
```

A callable noise model. `S = PinkNoise(...)`, then `S(w)` equals `S_pink_noise(w, S.Ad, S.Af; kT=S.kT)`, and `S.(ws)` works on arrays. The keyword constructor computes `(Ad, Af)` with `noise_amplitudes`.

### `default_pink_noise`

```julia
default_pink_noise(; EJ=4.5, EL=1.0, EC=1.0, kT=KT_15mK, kwargs...) -> PinkNoise
```

`PinkNoise` for a fluxonium at its sweet spot with energies in GHz (converted to angular units), with `φ_ge` taken between its two lowest eigenstates. `kwargs` go to `noise_amplitudes`.

### `calc_gs`

```julia
calc_gs(φ_upper, φ_lower) -> Vector
```

Harmonic weights `gₖ` of the flux coupling of a Floquet doublet. The inputs are `⟨φ⟩(tₙ)` for the upper and lower state at `N` equally spaced times `tₙ = nT/N`, `n = 0…N-1` over one period (**without** the endpoint `T`). Returns the `N` weights of `g(t) = (⟨φ⟩_upper - ⟨φ⟩_lower)/2` with the convention of `numpy.fft.ifft`, ordered by harmonic index `0, 1, …, -1`. `g₀` is the period average.

### `calc_Γφ`

```julia
calc_Γφ(gs, S, Af, w_base; g0_centered=false)
calc_Γφ(gs, S::PinkNoise, w_base; g0_centered=false)
calc_Γφ(φ_upper, φ_lower, S::PinkNoise, w_base)
```

Dephasing rate of a driven qubit:

    Γφ = 8·Af·|g₀| + Σ_{k≠0} |gₖ|²·S(k·w_base)

- `gs`: harmonic weights from `calc_gs`, in `0, 1, …, -1` order, or `-N÷2 … N÷2-1` if `g0_centered`.
- `S`: noise PSD, a function of angular frequency that accepts scalars and arrays (e.g. a `PinkNoise`).
- `Af`: amplitude of the zero-frequency noise (taken from `S.Af` for a `PinkNoise`).
- `w_base = 2π/T`: angular frequency of the drive period.
- The third method computes `gs` from the `φ` samples with `calc_gs`.

---

## Flux modulation and dephasing

Flux modulation of qubit A, `δϕA(t) = 2πΦ cos(ωt)`, in the adiabatic limit where Floquet states follow the instantaneous eigenstates (`dephasing.jl`). The workflow is `run_flux_modulation_sweep` → `dephasing_rates`, with `floquet_instantaneous_check` to validate the adiabatic assumption.

### `thermal_energy`

```julia
thermal_energy(T) -> Float64
```

`k_B T/ħ` in rad/ns for a temperature `T` in kelvin, using the rounded constants `k_B = 1.38e-23`, `h = 6.63e-34`. `thermal_energy(15e-3)` equals `KT_15mK`.

### `qubit_pink_noise`

```julia
qubit_pink_noise(params; kT=thermal_energy(1.5e-3), qho_basis_nlev=50, ebasis_nlev=6,
                 bare_to_dressed_tracking=:overlap, kwargs...) -> PinkNoise
```

`PinkNoise` of a single fluxonium with `fluxonium` parameters `params` (including its flux offset `δϕ`), from its `EL`, `EC` and the matrix element `⟨g|φ|e⟩`. `kT` is in rad/ns (default 1.5 mK). `kwargs` go to `noise_amplitudes`.

### `instantaneous_flux_expectations`

```julia
instantaneous_flux_expectations(Φ, H0, V, φ_op, ψ0, θs) -> (φt, ψt)
```

Follow the instantaneous eigenstates of `H0 + 2πΦ cos(θ) V` that are connected to the columns of `ψ0` over one drive period, by maximal overlap between neighbouring samples, starting at `θ = π/2`. `H0`, `V` and `φ_op` are **dense matrices**.

- `Φ`: modulation amplitude in units of Φ0. `θs`: equally spaced drive phases over one period without the endpoint, dense enough that the eigenvectors change little between samples (a warning is emitted if tracking is ambiguous).
- **Returns** the `length(θs) × size(ψ0, 2)` matrix of `⟨φ_op⟩` and the vector of tracked eigenvector matrices.

### `flux_harmonic_differentials`

```julia
flux_harmonic_differentials(φt, θs; pairs=((1, 3), (2, 4))) -> (∂0, ∂1)
```

Period average and first-harmonic amplitude of `⟨φ⟩(θ)` for the columns of `φt` (ordered gg, ge, eg, ee), combined as the mean over `pairs = ((a, b), …)` of `|φ[b] - φ[a]|`. The default pairs are the qubit-A transitions gg→eg and ge→ee.

### `run_flux_modulation_sweep`

```julia
run_flux_modulation_sweep(fA_params, fB_params, Φlist, J; δϕ_offsets=(0.005, -0.005), nθ=400,
    H_builder=first ∘ build_H_S1, H_builder_kwargs=nothing,
    bare_to_dressed_tracking=:overlap, showprogress=true) -> NamedTuple
```

For each amplitude `Φ` in `Φlist` (units of Φ0), record `⟨φA⟩` along the eigenstates connected to gg, ge, eg, ee over one period, and its period average and first harmonic. All of these are independent of the drive frequency. **Multithreaded** over `Φlist`.

- `δϕ_offsets`: shifts of `δϕ` for qubits A and B (rad), e.g. to lift the symmetry of the sweet spot.
- `nθ`: number of drive phases per period.
- **Returns** `(sys, θs, φt_sweep, ∂0, ∂1)`: the coupled system, the drive phases, the vector over `Φlist` of `nθ × 4` matrices of `⟨φA⟩`, and the vectors over `Φlist` of the period average and first harmonic.

### `floquet_instantaneous_check`

```julia
floquet_instantaneous_check(sys, J, ω, Φ; nθ=400, nT=20, δλ=0.5, fA_lbl=:fA, fB_lbl=:fB,
    bare_to_dressed_tracking=:overlap, floquet_kwargs...) -> NamedTuple
```

Compare the instantaneous eigenstates with the tracked Floquet states of `sys` under `δϕA(t) = 2πΦ cos(ωt)` (`ω` in rad/ns). Best done at the largest amplitude and fastest drive, where non-adiabatic corrections are largest. The Floquet states are sampled at `nT` times per period, and `nθ` must be a multiple of `nT`. `floquet_kwargs` (e.g. `n_basis_states=10`) go to `FloquetBasis`.

- **Returns** a `NamedTuple` with `min_overlap` (smallest `|⟨inst|Floquet⟩|²`) and `max_Δφ` (largest `|Δ⟨φA⟩|`), each a `NamedTuple` over `gg, ge, eg, ee`, plus `Floquet` and `instantaneous`, each `(; ∂0, ∂1)`.

### `dephasing_rates`

```julia
dephasing_rates(φt_sweep, ω_list, S; pairs=((1, 3), (2, 4))) -> Matrix
```

Estimated dephasing rate of qubit A from the noise PSD `S` (e.g. a `PinkNoise`) for each drive frequency in `ω_list` (rad/ns) and each amplitude of `φt_sweep` (from `run_flux_modulation_sweep`). For each `(a, b)` in `pairs` the weights are those of `calc_gs` (with `b` upper and `a` lower), and the rates of `calc_Γφ` are averaged. The weights do not depend on `ω`, which enters only through `S(kω)`.

- **Returns** a matrix of rates in 1/ns indexed `[i_ω, j_Φ]`.

---

## Plotting

Requires `CairoMakie`, which `XoniumUnitCell` loads.

### `zz_heatmap_figure`

```julia
zz_heatmap_figure(Z; x, y, xlabel, ylabel, title="", colorbar_label, colormap=Reverse(:Blues),
                  clims=(1e-2, 1e2), dashed_curve=nothing, device_point=nothing,
                  device_label="Our device") -> (fig, ax)
```

Log-scaled heatmap of a ZZ map `Z` over the grid `x × y`, with a colorbar. The labels default to `J₁₂ (MHz)`, `J_c (MHz)` and `ζ (GHz)`. The color range is taken from `extrema(Z)` (a constant `Z` is widened so the log scale is drawable), so `Z` must be positive, and `clims` is currently unused.

- `dashed_curve = (x_curve, y_curve)`: draws a dashed black "Minimum ZZ" curve.
- `device_point = (x0, y0)`: marks a point labelled `device_label`.
- **Returns** the `Figure` and the `Axis`.

---

## Not exported

- `dynamical_flux_terms(sys, δϕtA, δϕtB, J; symmetric_drive=false, amplitudes=nothing, fA_lbl, fB_lbl)` (`gate_utils.jl`): linear response `∂H/∂δϕ` to time-dependent flux detunings, as a `QobjEvo`. Used by `cz_gate_matrix`. Call as `XoniumUnitCell.dynamical_flux_terms`.
