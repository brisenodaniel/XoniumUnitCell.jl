export CPW_freq2cap, Cinv2E_C, Cinv2J, Cinv2g, Linv2E_L, L2E_L, Ic2E_J, E_J2Ic

# Units: capacitances in fF, inductances in nH, critical currents in nA, energies as E/h in GHz
# (couplings in MHz). The matrices of a `CGraph` use fF and 1/nH.

# Physical constants (CODATA 2018 exact SI values)
const _e_charge = 1.602176634e-19  # C
const _h_planck = 6.62607015e-34   # J s
const _ħ_planck = _h_planck / 2π   # J s

"""
    CPW_freq2cap(freq; Z0=50)

Capacitance (fF) to ground of a λ/4 CPW resonator with frequency `freq` (GHz, ωr/2π) and
impedance `Z0` (Ω).
"""
CPW_freq2cap(freq; Z0=50) = 1e15 / (8 * freq * 1e9 * Z0)

"""
    Cinv2E_C(Cinv)

Convert an inverse-capacitance matrix element `Cinv` (1/fF) to a charging energy `E_C/h` (GHz),
`E_C = e²/2 · Cinv`. For a mode, `Cinv` is the diagonal element of the inverse of the capacitance
matrix in the basis of that mode.
"""
Cinv2E_C(Cinv::Real) = 1e-9 * 0.5 * _e_charge^2 * Cinv * 1e15 / _h_planck

"""
    Cinv2J(Cinv)

Convert an inverse-capacitance matrix element `Cinv` (1/fF) to a qubit-qubit capacitive coupling
`J/h` (MHz), `J = 4e² · Cinv`. `Cinv` is the off-diagonal element of the inverse capacitance
matrix between the two modes.
"""
Cinv2J(Cinv::Real) = 1e-6 * 4 * _e_charge^2 * Cinv * 1e15 / _h_planck

"""
    Cinv2g(Cinv)

Convert an inverse-capacitance matrix element `Cinv` (1/fF) to a qubit-resonator coupling
`g/2π` (MHz).
"""
Cinv2g(Cinv::Real) = 1e-6 * _e_charge * sqrt(π * _ħ_planck / 100) * Cinv * 1e15 / _h_planck

"""
    Linv2E_L(Linv)

Convert an inverse-inductance matrix element `Linv` (1/nH) to an inductive energy `E_L/h` (GHz),
`E_L = (ħ/2e)² · Linv`, so that the potential is `½ E_L φ²` for the phase `φ`. For a mode, `Linv`
is the diagonal element of the inverse-inductance matrix in the basis of that mode, e.g.
`graph.L_inv[k, k]`. One nanohenry corresponds to about 163.46 GHz.
"""
Linv2E_L(Linv::Real) = (_ħ_planck / (2 * _e_charge))^2 * Linv / _h_planck

"""
    L2E_L(L)

Convert an inductance `L` (nH) to an inductive energy `E_L/h` (GHz), see [`Linv2E_L`](@ref).
"""
L2E_L(L::Real) = Linv2E_L(inv(L))

"""
    Ic2E_J(Ic)

Convert a junction critical current `Ic` (nA) to a Josephson energy `E_J/h` (GHz),
`E_J = ħ Ic / 2e`. One nanoampere corresponds to about 0.4967 GHz.
"""
Ic2E_J(Ic::Real) = 1e-18 * Ic / (4π * _e_charge)

"""
    E_J2Ic(EJ)

Convert a Josephson energy `E_J/h` (GHz) to a critical current `Ic` (nA), the inverse of
[`Ic2E_J`](@ref).
"""
E_J2Ic(EJ::Real) = EJ * 4π * _e_charge / 1e-18
