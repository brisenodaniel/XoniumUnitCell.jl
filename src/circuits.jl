export φn_operators, fluxonium, transmon, split_transmon, transmon_nbasis
export tounit_ħ, tounit_h

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

function cosφ_operator(n::AbstractQuantumObject, ϕext::Float64=0.0)
    dim = size(n)[1]
    cosφ = 0.5 * (
    exp(-im*ϕext)*diagm(1=>ones(dim-1)) +
        exp(im*ϕext)*diagm(-1=>ones(dim-1))
    )
    return Qobj(cosφ, type=n.type, dims=n.dims)
end


function fluxonium(;
                   EC::Real,
                   EL::Real,
                   EJ::Real,
                   δϕ::Real = 0.0,
                   n_g::Real=0.0,
                   nlev::Int=40,
                   Φext_incos::Bool=false,
                   name::Symbol=:Fluxonium)
    φ, n = φn_operators(;EC=EC, EL=EL, nlev=nlev)
    Φ = π - δϕ
    Φ_lin, Φ_cos = Φext_incos ? (0.0 , Φ) : (Φ, 0.0)
    H = 4.0 * EC * (n-n_g)^2 + # Capacitive Energy
        (1.0/2.0) * EL * (φ - Φ_lin)^2 - # Linear inductor energy
        EJ * cos(φ - Φ_cos) # JJ energy
    ops = (;φ=φ, n=n, H=H,)
    params = (;EC=EC, EL=EL, EJ=EJ, δϕ=δϕ, Φext_incos=Φext_incos)
    return SimpleQuantumSystem(ops; params=params, name=name)
end


function transmon_hbasis(;
                  EC::Real,
                  EJ::Real,
                  ng::Real=0.0,
                  nlev::Int=40,
                  name::Symbol=:Transmon,
    )
    φ, n = φn_operators(;EC=EC, EJ=EJ, nlev=nlev)
    H = 4.0 * EC * (n-ng)^2 - # Capacitive energy
        EJ * cos(φ) # JJ energy
    ops = (;φ=φ, n=n, H=H,)
    params = (;EC=EC, EJ=EJ, ng=ng)
    return SimpleQuantumSystem(ops; params=params, name=name)
end

function transmon(;
                 EC::Real,
                 EJ::Real,
                 ng::Real=0.0,
                 ncut::Int=20,
                 name::Symbol=:Transmon,
                         )
    nlev = 2*ncut + 1
    n = collect(0:nlev-1) .- ncut |> Diagonal |> Qobj
    cosφ = cosφ_operator(n)
    H = 4.0 * EC * (n - ng)^2 - # capacitive energy
        EJ * cosφ # JJ energy
    ops = (;n=n, cosφ=cosφ, H=H)
    params = (; EC=EC, EJ=EJ, ng=ng)
    return SimpleQuantumSystem(ops; params=params, name=name)
end

function split_transmon(;
                        EC::Real,
                        EJ1::Real,
                        EJ2::Real,
                        ϕext::Real=0.0,
                        ng::Real=0.0,
                        ncut::Int=20,
                        name::Symbol=:Transmon
                        )
    trans = transmon(;EC=EC, EJ=EJ1, ng=ng, ncut=ncut, name=name)
    H = trans.H - EJ2*cosφ_operator(trans[:n], ϕext)
    ops = (;trans.ops..., H=H)
    params = (;EC=EC, EJ1=EJ1, EJ2=EJ2, ϕext=ϕext, ng=ng, ncut=ncut)
    return SimpleQuantumSystem(ops; params=params, name=name)
end

function split_transmon_hbasis(;
                        EC::Real,
                        EJ1::Real,
                        EJ2::Real,
                        n_g::Real=0.0,
                        Φ_ext::Real=0.0,
                        nlev::Int=40,
                        name::Symbol=:Transmon,
                        )
    trans = transmon(; EC=EC, EJ=EJ1, n_g=n_g, nlev=nlev, name=name)
    return trans + (-1* EJ2*cos(trans[:φ] - Φ_ext))
end



function tounit_ħ(
    params::Union{NamedTuple, Dict},
    convert_keys::Union{Tuple, AbstractVector}=(:EC, :EJ, :EL)
    )
    return map(keys(params), values(params)) do k, v
        k ∈ convert_keys ? k => 2π*v : k => v
    end
end

function tounit_h(
    params::Union{NamedTuple, Dict},
    convert_keys::Union{Tuple, AbstractVector}=(:EC, :EJ, :EL)
    )
    return map(keys(params), values(params)) do k, v
        k ∈ convert_keys ? k => v/2π : k => v
    end
end
