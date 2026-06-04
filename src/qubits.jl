export φn_operators, fluxonium, transmon

function φn_operators(;EC::Real, EJ::Maybe{Real}=nothing, EL::Maybe{Real}=nothing, nlev::Int=10)
    isnothing(EJ) && isnothing(EL) && throw(
        ArgumentError("No potential term. Both `EJ` and `EL` cannot be `nothing`.")
    )
    Z = isnothing(EL) ? EC/EJ : EC/EL
    φ_zpf = (8.0 * Z)^(1/4)
    n_zpf = (8.0 * Z)^(-1/4)
    φ = φ_zpf * QuantumToolbox.position(nlev)
    n = n_zpf * QuantumToolbox.momentum(nlev)
    return φ, n
end


function fluxonium(;
                   EC::Real,
                   EL::Real,
                   EJ::Real,
                   δϕ::Real = 0.0,
                   n_g::Real=0.0,
                   nlev::Int=40,
                   name::Symbol=:Fluxonium)
    φ, n = φn_operators(;EC=EC, EL=EL, nlev=nlev)
    Φ = π - δϕ
    H = 4.0 * EC * (n-n_g)^2 + # Capacitive Energy
        (1.0/2.0) * EL * (φ + Φ)^2 - # Linear inductor energy
        EJ * cos(φ) # JJ energy
    ops = (;φ=φ, n=n, H=H,)
    params = (;EC=EC, EL=EL, EJ=EJ, δϕ=δϕ,)
    return SimpleQuantumSystem(ops; params=params, name=name)
end

function transmon(;
                  EC::Real,
                  EJ::Real,
                  n_g::Real=0.0,
                  nlev::Int=40,
                  name::Symbol=:Transmon)
    φ, n = φn_operators(;EC=EC, EJ=EJ, nlev=nlev)
    H = 4.0 * EC * (n-n_g)^2 - # Capacitive energy
        EJ * cos(φ) # JJ energy
    ops = (;φ=φ, n=n, H=H)
    params = (;EC=EC, EJ=EJ)
    return SimpleQuantumSystem(ops; params=params, name=name)
end
