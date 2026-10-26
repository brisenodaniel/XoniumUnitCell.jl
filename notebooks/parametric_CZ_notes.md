# Parametric CZ: working notes

Summary of the work on `parametric_CZ.ipynb`, reproducing Ma et al., arXiv:2308.16040
("Native approach to controlled-Z gates in inductively coupled fluxonium qubits").
Units: GHz and ns with ħ = 1, unless stated otherwise.

## 1. Hamiltonian conventions vs. the paper (Fig. 2 mismatch)

- The code uses the fluxonium Hamiltonian
  `H = 4EC n² + ½EL(φ + φext)² − EJ cos φ` (`src/circuits.jl`). Of the candidate conventions,
  this is the only one consistent with the paper: it reproduces the qubit-crossing points at
  δΦ ≈ ±0.0614 Φ0 and J·φge·φge ≈ 19.3 MHz. Eq. (S1), taken literally, is internally inconsistent.
- φ_zpf = (8EC/EL)^(1/4) in `φn_operators` was checked and is correct.
- Coupling: inductive, `J(φA + φextA)(φB + φextB)`. The sign of J must be **negative**
  (`J = -3.5·2π·1e-3`) to match the sign structure of Figs. 2(a) and 2(c).
- Remaining discrepancy: a rigid offset of ≈ 35 MHz in qubit-A frequencies (fA(0.458 Φ0) ≈ 1.026 GHz
  in simulation vs ≈ 1.06 GHz in the paper). It is not explained by any convention choice and is
  probably a parameter difference.
- Flux detuning from the sweet spot enters exactly (up to a c-number) as
  `H(δ) = H_idle + δA(EL_A φA + J φB) + δB(EL_B φB + J φA)` (`get_perturbed_system`, `cz_drive_op`).

## 2. Conditional phase (Fig. 2(c))

- Conditional phase: `ϕ = arg(c_ee c_gg / (c_eg c_ge))` with `c_k = ⟨k|U|k⟩`.
  Statically, ϕ = −∫ζ dt with ζ = E_ee − E_eg − E_ge + E_gg = 4 g_zz.
  Paper convention: v_ϕ = −ϕ/τ.
- ϕ is unwrapped onto the 2π branch closest to the static estimate −ζ(τ − t_ramp);
  `branch_err` records the distance to that estimate, as a diagnostic.
- Pulses: `full_pulse(ω, t_ramp, t_plateau; ramp_profile=:cosine|:tanh, t_rise)`. The tanh ramp
  (rise time t_r = t_ramp/6 = 2 ns, as in the paper) gives a non-adiabatic error of ~2e-4.
  A 2 ns cosine ramp gives ~1e-2.
- `get_propagator(sys, tlist, H_int, dr; abstol, reltol)` returns the computational-subspace
  propagator P·U(t)·P from four `sesolve` runs. Verified against the full propagator to 6e-15.
- Leakage: `L = 1 − Tr(U†U)/4`, the average population leaving the computational subspace.
  Per-state leakage is `1 − |c_k|²`; the worst case comes from the singular values of PUP.
- Notebook functions: `cz_propagator`, `get_νφ`, `unwrap_phase`, `static_ζ`, `get_νφ_point`,
  `figure_2b_sweep` (threaded over the δΦA × δΦB grid). The superseded functions
  `get_driven_flux_ops`, `get_CZ_hamiltonian` and `run_CZ` are commented out, not deleted.
- Check at (δΦA, δΦB) = (−0.06, −0.06): plateau slope 228.44 rad/µs (dynamic) vs 228.48 rad/µs
  (static −ζ). The new path is ~18× faster than the old one.

## 3. QuantumToolbox fork and Floquet tools

- `Project.toml` takes QuantumToolbox from the fork
  `https://github.com/brisenodaniel/QuantumToolbox.jl.git#main` (currently commit `586220f`).
- Fork additions used here:
  - `track_quasienergies(; H0, H_dr, T, nT, δλ, min_overlap, ...)`: ramps the drive λ = 0 → 1 and
    matches Floquet states to the static eigenstates of `H0` by overlap.
  - `FloquetBasis(H0, H_dr, T, precompute; δλ, nT, min_overlap, ...)`: entry `i` is connected to the
    `i`-th eigenstate of `H0`, and the quasienergies are unfolded onto the static branches.
  - Pass `progress_bar=Val(false)` to silence the `sesolve` progress bar. This keyword is also stored
    for later propagator calls.

## 4. `state(sys, fb::FloquetBasis, lbl, t=0.0; ...)` (`src/QuantumSystems.jl`)

- Resolves `lbl` with the usual `state` rules (`bare_to_dressed_tracking = :spectrum | :overlap`).
  It returns `(ε, ψ(t))`: the tracked quasienergy, and the Floquet state ψ(t) = U(t,0) ψ(0).
  This is the state, not the periodic mode u(t).
- Throws on a dimension mismatch. Warns if |⟨static|ψ_fl(0)⟩|² < `min_overlap`, which flags an
  untracked basis or a basis built from a different `H0`.
- `cache=true` stores the micromotion propagator in `fb`. Also adds `energy(sys, fb, lbl)`.
- Tests: the "state with FloquetBasis" testset in `test/test_CompositeQuantumSystem.jl`.
  Full suite: 43 / 69 / 6 / 2 passing.

## 5. Flux-dispersion differentials under sinusoidal modulation (Fig. 3)

Goal: the sensitivity of qubit A to flux under a modulation Φ_A(t) = Φ cos(ωt) about Φ0/2, with
qubit B parked at Φ0/2, for ω/2π = 50, 5 and 0.5 MHz.

- Quantities, from ⟨φA⟩ in the states gg, ge, eg, ee:
  - `∂ωs`: the period-averaged differential (|⟨φA⟩₀,eg − ⟨φA⟩₀,gg| + |⟨φA⟩₀,ee − ⟨φA⟩₀,ge|)/2.
  - `∂λ_1`: the same combination of the complex first-harmonic amplitudes,
    ⟨φA⟩(t) ≈ φ₀ + Re(φ₁ e^{iωt}).
  - Multiply by EL_A to convert to ∂ω/∂(2πΦ_A); the J term (~1%) is neglected.
- `∂ωs` vanishes, up to ~1e-12, by symmetry: modulating about the sweet spot averages ∂ω/∂Φ to
  zero (paper, Fig. 3(a)). T2 instead depends on the spectral content of ∂ω/∂Φ(t), which is what
  `∂λ_1` (and, more completely, the filter function in Eqs. S9–S13) captures.
- Full Floquet approach: too slow. The 36-level space reaches 167 rad/ns, so one period solve takes
  ~6 s at 50 MHz and the cost grows ∝ T. That is ~220 s per point at 5 MHz and roughly 1 h per point
  at 0.5 MHz.
- Current approach (cell 20): **instantaneous eigenstates**.
  - Diagonalize H(δ) at `nθ = 400` phases θ = ωt, with δ = 2πΦ cos θ.
  - Track gg/ge/eg/ee by overlap with the previous sample, starting from δ = 0.
  - In this adiabatic limit the results do not depend on ω, so all frequency rows are identical.
- Consistency check (end of cell 20):
  - Builds a tracked `FloquetBasis` with δλ = 0.25 at the largest Φ and the fastest ω (50 MHz).
  - Prints the minimum overlaps |⟨inst|Floquet⟩|², the max |Δ⟨φA⟩|, and both `∂ωs` and `∂λ_1`
    from each method.
- Cell 21 plots `∂ωs` and `∂λ_1` side by side vs Φ, one curve per drive frequency.
- Open issues:
  - Earlier Floquet runs with δλ = 1 gave a matched overlap of 0.469 at Φ = 0.1, where qubit A
    sweeps through its crossing with qubit B (~±0.06 Φ0). The labels near this crossing need care.
  - With few time samples, higher harmonics alias onto the first harmonic.
  - The new cell 20 has not been run yet.

## 6. Pending / ideas

- Run the full Fig. 2(c) grid (`figure_2b_sweep`; ~1.5 min estimated).
- Fig. 3(b) T2: use a filter function with 1/f flux noise, S_Φ = A²/(2|f|) with A ≈ 10.6 µΦ0/√Hz,
  applied to the time series ∂ω/∂Φ(t). Alternatively, Monte Carlo over noise realizations.
- Notebook backups: in the session scratchpad (`parametric_CZ.backup*.ipynb`).
