

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
