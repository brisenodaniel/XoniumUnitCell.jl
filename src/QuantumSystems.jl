abstract type AbstractQuantumSystem end
const Qsys = AbstractQuantumSystem
const QuantumSystem = AbstractQuantumSystem
const StateIndex = Union{Int, AbstractVector{<:Int}, Symbol}
const Maybe{T} = Union{T, Nothing} where T
const _STATE_SYMBOLS = (; g=1, e=2, f=3, h=4) # Julia indexing

export AbstractQuantumSystem, Qsys, QuantumSystem, StateIndex, Maybe
export SimpleQuantumSystem, CompositeQuantumSystem, Qsys, tensor_dims
export state, energy, to_energy_basis, truncate_sys, I, state_idx

struct SimpleQuantumSystem{
    Ops<:NamedTuple,
    PT<:NamedTuple} <: AbstractQuantumSystem
    ops::Ops
    params::Maybe{PT}
    name::Maybe{Symbol}

    function SimpleQuantumSystem(
        ops::Ops;
        params::Maybe{PT} = nothing,
        name::Maybe{Symbol} = nothing
        ) where {
        Ops<:NamedTuple,
        PT<:NamedTuple
        }
        err = _check_ops(ops)
        isnothing(err) ? nothing : throw(err)
        return new{Ops, PT}(ops, params, name)
    end
end

struct CompositeQuantumSystem{
    Ops<:NamedTuple,
    Systms<:NamedTuple,
    PT<:NamedTuple,} <: AbstractQuantumSystem

    subsys::Systms
    ops::Ops
    params::Maybe{PT}
    U_basis :: Maybe{AbstractQuantumObject}
    dressed_dim :: Maybe{Int}
    name :: Maybe{Symbol}
    state_lookup::AbstractArray{Union{Int64, Missing}}

    function CompositeQuantumSystem(
        subsys::Systms,
        ops::Ops;
        params::Maybe{PT}=nothing,
        U_basis::Maybe{AbstractQuantumObject}=nothing,
        dressed_dim::Maybe{Int}=nothing,
        name::Maybe{Symbol}=nothing,
        ) where{
            Ops<:NamedTuple,
            Systms<:NamedTuple,
            PT<:NamedTuple,
        }
        err = _check_subsys(subsys)
        err = isnothing(err) ? _check_ops(ops) : err
        isnothing(err) ? nothing : throw(err)
        n_sys = length(subsys)
        tdims = (size(sys.H)[1] for (_, sys) in pairs(subsys)) |> collect
        state_lookup = Array{Union{Int64, Missing}, n_sys}(missing, tdims...)
        return new{Ops, Systms, PT}(subsys, ops, params, U_basis, dressed_dim, name, state_lookup)
    end
end

"""
    _require_name(sys)

Throw if `sys.name === nothing`, otherwise return `sys.name`.
"""
_require_name(sys::AbstractQuantumSystem) =
    isnothing(sys.name) ?
        throw(ArgumentError("All subsystems must have a non-nothing `name`")) :
        sys.name


"""
    _tensor_identity(H)

Return identity operator on the Hilbert space of operator `H`.
"""
_tensor_identity(H::AbstractQuantumObject) = qeye_like(H)


"""
    _lift_left(op, Hright)

Return op ⊗ I
"""
_lift_left(op::AbstractQuantumObject, Hright::AbstractQuantumObject) =
    tensor(op, _tensor_identity(Hright))


"""
    _lift_right(Hleft, op)

Return I ⊗ op
"""
_lift_right(Hleft::AbstractQuantumObject, op::AbstractQuantumObject) =
    tensor(_tensor_identity(Hleft), op)

function _lift_local_operator(
    sys::CompositeQuantumSystem,
    subkey::Symbol,
    opkey::Symbol,
)
    # --- sanity checks ---
    hasproperty(sys.subsys, subkey) ||
        throw(KeyError("No subsystem named $subkey"))

    subsys = sys.subsys[subkey]

    hasproperty(subsys.ops, opkey) ||
        throw(KeyError("Subsystem $subkey has no operator $opkey"))

    # --- build tensor factors in correct order ---
    factors = map(keys(sys.subsys)) do name
        if name === subkey
            subsys.ops[opkey]
        else
            qeye(sys.subsys[name])
        end
    end

    # tensor in order
    O = reduce(tensor, factors)

    # --- apply dressed basis if present ---
    if sys.U_basis !== nothing
        U = sys.U_basis
        O = U' * O * U
    end

    return O
end

function _lift_multi_operator(
    sys::CompositeQuantumSystem,
    names::NTuple{N,Symbol},
) where {N}
    # --- must be pairs ---
    iseven(N) ||
        throw(ArgumentError("Expected (subsystem, operator) pairs"))

    # --- parse pairs ---
    pairs = Pair{Symbol,Symbol}[
        names[i] => names[i+1] for i in 1:2:N
    ]

    # --- forbid duplicate subsystems ---
    subs = first.(pairs)
    length(unique(subs)) == length(subs) ||
        throw(ArgumentError(
            "Each subsystem may appear at most once in operator lifting"
        ))

    # --- validate subsystems and operators ---
    for (subkey, opkey) in pairs
        hasproperty(sys.subsys, subkey) ||
            throw(ArgumentError("No subsystem named $subkey"))

        hasproperty(sys.subsys[subkey].ops, opkey) ||
            throw(ArgumentError(
                "Subsystem $subkey has no operator $opkey"
            ))
    end

    # --- map requested local operators ---
    local_ops = Dict{Symbol,AbstractQuantumObject}(
        subkey => sys.subsys[subkey].ops[opkey]
        for (subkey, opkey) in pairs
    )

    # --- build tensor factors in subsystem order ---
    factors = map(keys(sys.subsys)) do name
        if haskey(local_ops, name)
            local_ops[name]
        else
           qeye(sys.subsys[name])
        end
    end

    # --- tensor product ---
    O = reduce(tensor, factors)

    # --- apply dressed basis if present ---
    if sys.U_basis !== nothing
        U = sys.U_basis
        O = U' * O * U
    end

    return O
end

function _check_ops_addition(
    sys::AbstractQuantumSystem,
    ops::NamedTuple,
)
    for (k, op) in pairs(ops)
        op isa AbstractQuantumObject ||
            throw(ArgumentError(
                "Value for key $k must be an AbstractQuantumObject"
            ))

        sys.H.dimensions == op.dimensions ||
            throw(DimensionMismatch(
                "Operator $k has incompatible dimensions with system Hamiltonian"
            ))
    end
end

function _check_ops(ops::NamedTuple)
    msg = nothing
    if !_all_qobj(ops)
        msg = "All elements of `ops` must be QuantumObjects"
    elseif !(values(ops) |> same_dims)
        msg = "All QuantumObjects in `ops` must have the same dimension"
    elseif !hasproperty(ops, :H)
        msg = "`ops` must contain Hamiltonian operator with key `:H`"
    end
    return isnothing(msg) ? msg : ArgumentError(msg)
end

function _check_subsys(subsys::NamedTuple)
    msg = nothing
    if isempty(subsys) || !all(
        sys -> sys isa AbstractQuantumSystem, subsys
        )
        msg = "Parameter `subsys` must be a non-empty NamedTuple with " *
            "`AbstractQuantumSystem` values"
    end
    return isnothing(msg) ? msg : ArgumentError(msg)
end

function _state_index(lbl::StateIndex, sys::AbstractQuantumSystem)
    if lbl isa Int
        idx = [lbl + 1,]
    elseif lbl isa AbstractVector{<:Int}
        idx = lbl .+ 1
    else
        _lbl = (collect ∘ String)(lbl)
        idx =  map(_lbl) do l_n
            l_n = Symbol(l_n)
            return _STATE_SYMBOLS[l_n]
        end
    end
    _check_idx(idx, sys)
    return idx
end

function _check_idx(lbl::AbstractVector{<:Int}, sys::AbstractQuantumSystem)
    tdims = length(lbl) == 1 ? (size(sys.H)[1],) : tensor_dims(sys)
    length(tdims) == length(lbl) || throw(ArgumentError(
        "Number of characters or integers in `lbl` must match number of subsystems in `sys`"
    ))
    idx_checks = map(lbl, tdims) do l_n, t_n
        l_n ∉ 1:t_n
    end
    bad_idx = findall(idx_checks)
    return isempty(bad_idx) || throw(ArgumentError(
        "State indices at `lbl` indices $bad_idx lie outside of the range allowed  " *
            "by the Hilbert space dimension of their subsystem "
    ))
end

function same_dims(ops_vals)
    if isempty(ops_vals)
        return true
    end
    dim0 = first(ops_vals).dimensions
    return all(q -> q.dimensions == dim0, ops_vals)
end

function _all_qobj(ops::NamedTuple)
    if isempty(ops)
        return false
    end
    return all(q -> q isa AbstractQuantumObject, ops)
end



function Qsys(
    ops::Ops;
    params::Maybe{PT}=nothing,
    name::Maybe{Symbol}=nothing,
    ) where {
        Ops <: NamedTuple,
        PT <: NamedTuple
    }
    return SimpleQuantumSystem(ops; params=params, name=name)
end

function Qsys(
    sys::SimpleQuantumSystem,
    name::Maybe{Symbol}=nothing
    )
    new_ops = (;sys.ops...)
    new_params = isnothing(sys.params) ? nothing : (;sys.params...)
    name = isnothing(name) ? sys.name : name
    return SimpleQuantumSystem(new_ops; new_params=params, name=name)
end

function tensor_dims(sys::AbstractQuantumSystem)
    if sys isa SimpleQuantumSystem
        return collect(sys.H.dims)
    else
        return ((tensor_dims(subsys)[1] for subsys in sys.subsys)...,)
    end
end

function QuantumToolbox.qeye(sys::QuantumSystem)
    return qeye_like(sys.H)
end

function I(sys::QuantumSystem)
    return qeye(sys)
end

function Base.getproperty(sys::AbstractQuantumSystem, key::Symbol)
    if hasfield(typeof(sys), key)
        return Core.getfield(sys, key)
    elseif hasproperty(Core.getfield(sys, :ops), key)
        return Core.getfield(sys, :ops)[key]
    elseif sys isa CompositeQuantumSystem &&
        hasproperty(Core.getfield(sys, :subsys), key)
        return Core.getfield(sys, :subsys)[key]
    elseif key==:I
        return qeye(sys)
    else
        throw(FieldError(typeof(sys), key))
    end
end

function Base.getindex(sys::SimpleQuantumSystem, key::Symbol)
    return Core.getfield(sys, :ops)[key]
end

function Base.getindex(
    sys::CompositeQuantumSystem,
    key::Symbol,
)
    if hasproperty(sys.ops, key)
        return sys.ops[key]
    elseif hasproperty(sys.subsys, key)
        return sys.subsys[key]
    else
        throw(KeyError("No operator or subsystem named $key"))
    end
end


function Base.getindex(
    sys::CompositeQuantumSystem,
    subkey::Symbol,
    opkey::Symbol,
)
    return _lift_local_operator(sys, subkey, opkey)
end


function Base.getindex(
    sys::CompositeQuantumSystem,
    keys::Vararg{Symbol},
)
    return _lift_multi_operator(sys, keys)
end


function to_energy_basis(
    sys::SimpleQuantumSystem;
    set_energy_zero::Union{Bool, Float64}=false,
    etol::Union{Bool, Float64}=false)
    λ, _, U = eigenstates(sys.H)
    U = Qobj(U, dims=sys.H.dims)
    e0 = set_energy_zero === true ? λ[1] :
        set_energy_zero isa Float64 ? set_energy_zero : 0.0
    new_ops = (;(k=> U'*op*U for (k, op) in pairs(sys.ops))...)
    if e0 != 0.0
        H = new_ops[:H] - e0
        new_ops = (; new_ops..., H=H)
    end
    sys = SimpleQuantumSystem(new_ops; params=sys.params, name=sys.name)
    if etol isa Float64
        return tidyup(sys, etol)
    elseif etol
        return tidyup(sys)
    else
        return sys
    end
end

function to_energy_basis(
    sys::CompositeQuantumSystem;
    set_energy_zero::Union{Bool, Float64} = true,
    etol::Union{Bool, Float64} = false,
)
    # --- 1. Diagonalize composite Hamiltonian ---
    λ, _, U = eigenstates(sys.H)

    # Wrap U as a QuantumObject so algebra stays in QuantumToolbox
    Uq = QuantumObject(U, sys.H.type, sys.H.dims)

    # --- 2. Similarity transform of Hamiltonian ---
    H_diag = Uq' * sys.H * Uq

    # --- 3. Optional energy shift ---
    if set_energy_zero !== false
        e0 = set_energy_zero === true ? (minimum ∘ real)(λ) : set_energy_zero
        H_diag = H_diag - e0
    end

    # --- 4. Build new ops ---
    new_ops = (; H = H_diag)

    # --- 5. Construct dressed system ---
    dressed = CompositeQuantumSystem(
        sys.subsys,
        new_ops;
        params       = sys.params,
        U_basis      = Uq,
        dressed_dim  = length(λ),
        name         = sys.name,
    )

    # --- 6. Optional tidyup ---
    if etol isa Float64
        return tidyup(dressed, etol)
    elseif etol === true
        return tidyup(dressed)
    else
        return dressed
    end
end

function truncate_sys(sys::SimpleQuantumSystem, N::Int)
    new_ops = (;(k => truncate_Qobj(op, N) for (k, op) in pairs(sys.ops))...)
    return SimpleQuantumSystem(new_ops; params=sys.params, name=sys.name)
end

function QuantumToolbox.tidyup(sys::SimpleQuantumSystem, etol::Maybe{Float64}=nothing)
    etol = isnothing(etol) ? QuantumToolbox.settings.tidyup_tol : etol
    new_ops = (;(k=> tidyup(op, etol) for (k, op) in pairs(sys.ops))...)
    return SimpleQuantumSystem(new_ops; params=sys.params, name=sys.name)
end

function QuantumToolbox.tidyup(
    sys::CompositeQuantumSystem,
    etol::Maybe{Float64} = nothing,
)
    # Resolve tolerance
    tol = isnothing(etol) ? QuantumToolbox.settings.tidyup_tol : etol

    # --- 1. Tidy composite operators ---
    new_ops = (; (
        k => tidyup(op, tol)
        for (k, op) in pairs(sys.ops)
    )...)

    # --- 2. Tidy subsystems recursively ---
    new_subsys = (; (
        name => tidyup(subsys, tol)
        for (name, subsys) in pairs(sys.subsys)
    )...)

    # --- 3. Tidy basis transform if present ---
    new_U_basis = isnothing(sys.U_basis) ? nothing : tidyup(sys.U_basis, tol)

    # --- 4. Rebuild CompositeQuantumSystem ---
    return CompositeQuantumSystem(
        new_subsys,
        new_ops;
        params      = sys.params,
        U_basis     = new_U_basis,
        dressed_dim = sys.dressed_dim,
        name        = sys.name,
    )
end

function truncate_Qobj(q::AbstractQuantumObject, N::Int)
    return QuantumObject(q.data[1:N, 1:N], q.type, (N,))
end

@memoize LRU{Any, Any}(maxsize=100) function QuantumToolbox.eigenstates(sys::AbstractQuantumSystem)
    eigenstates(sys.H)
end

@doc raw"""
    `state(sys::SimpleQuantumSystem, lbl::StateIndex)`
Given a fock-space eigenstate label, returns the corresponding
eigenstate of the system Hamiltonian.

# Arguments
- `sys`: Simple QuantumSystem for which to return the Hamiltonian eigenstate
- `lbl`: The state label of the desired system eigenstate. Must be a positive integer, a symbol,
         or a vector of Ints of length 1.
- `bare_to_dressed_tracking`: Either `:overlap` or `:spectrum` (default). Accepted for API symmetry with
         `CompositeQuantumSystem` and validated, but otherwise ignored since the eigenstates of a
         simple system are already energy-ordered.

# Notes:
- If `lbl` is a positive integer `n`, then the function will return the `nth` eigenstate of the
  system Hamiltonian where `n=0` corresponds to the ground state.
- If `lbl` is a Symbol, it must be a single character in the set ``\{g, e, f, h\}``, which will
  correspond to the state indices ``(0,1,2,3)``.
- If `lbl` is a symbol outside of the allowed set, more than a single character, or an integer not in
.  the interval `[0, size(sys.H[1]))`, the function will throw an `ArgumentError`

# Returns:
- The eigenstate of `sys.H` corresponding to the state label `lbl`.
 """
function state(sys::SimpleQuantumSystem, lbl::StateIndex; bare_to_dressed_tracking::Symbol=:spectrum)
    _check_tracking(bare_to_dressed_tracking)
    λ, ψ, _ = eigenstates(sys.H)
    idx = _state_index(lbl, sys)[1]
    return λ[idx], ψ[idx]
end

@doc raw"""
    ``state(sys::CompositeQuantumSystem, lbl::StateIndex)

Given a fock-space eigenstate label, returns the corresponding eigenstate
of the system Hamiltonian.

# Arguments
- `sys`: CompositeQuantumSystem for which to return the Hamiltonian eigenstate.
- `lbl`: The state label of the desired system eigenstate.
- `bare`: Boolean flag which sets whether to return the eigenstate in the bare product
          basis or the dressed eigenbasis of `sys.H`.
- `bare_to_dressed_tracking`: Keyword, either `:overlap` or `:spectrum` (default). Only used when
          `bare==false` and `lbl` is a multi-subsystem label.
        - `:overlap`: the dressed eigenstate with maximal overlap with the bare product state.
        - `:spectrum`: the dressed eigenstate whose energy rank equals the energy rank of the
          bare product state in the bare spectrum (e.g. if `:ge` is the 2nd excited bare state,
          the 2nd excited eigenstate of `sys.H` is returned). Degenerate bare levels are ordered
          by label order (first subsystem index varies fastest), so the result is somewhat
          arbitrary for them; prefer `:overlap` there.

# Notes
- For a `sys` with `N` subsystems, `lbl` must be a positive integer, a Symbol of `N` characters,
  or a `Vector{Int, N}`.
- If `lbl` is a single positive integer `k`, function will return `eigenstates(sys.H)[2][k]`
- If `lbl` is a `Vector{Int, N}` or a `Symbol` of `N` characters in the set ``\{g,e,f,h\}``:
      - When `bare==true`, will return the kronecker product
             ``math
             \hat{U}\left(\bigotimes_{n=0}^N |l[n]\rangle\right) \hat{U}^\dagger
             ``
             where ``l_n`` is the `nth` element or character in `lbl`, `|l_n\rangle` is the corresponding eigenstate of
             `nth` subsystem, and `U` is the basis transform to the basis defined by `sys.U_basis`. If `isnothing(sys.U_basis)`,
             then ``U=1``.
        - When `bare==false`, will return the eigenstate of `sys.H` with the maximal value of
             ``math
             |\langle l | k\rangle|
             ``
             where ``| l \langle`` is the output of `state(sys, lbl, true)` and `|k\rangle` is any eigenstate of
             sys.H`.

# Returns
- The dressed or bare eigenenergy and eigenstate pair ``(\lambda, \psi)`` of `sys` specified by `lbl`.
"""
function state(
    sys::CompositeQuantumSystem,
    lbl::StateIndex,
    bare::Bool=false;
    bare_to_dressed_tracking::Symbol=:spectrum,
)
    _check_tracking(bare_to_dressed_tracking)
    if bare
        idx = _state_index(lbl, sys)
        return _bare_state(sys, idx)
    else
        idx = state_idx(sys, lbl; bare_to_dressed_tracking)
        λ, ψ, _ = eigenstates(sys.H)
        return λ[idx], ψ[idx]
    end
end

function _check_tracking(mode::Symbol)
    mode in (:overlap, :spectrum) || throw(ArgumentError(
        "`bare_to_dressed_tracking` must be `:overlap` or `:spectrum`, got `:$mode`"
    ))
    return mode
end

function _bare_state(sys::CompositeQuantumSystem, idx::AbstractVector{<:Int})

    if length(idx) == 1
        i = idx[1]
        H = tensor((sub.H for sub in sys.subsys)...)
        λ, ψ, _ = eigenstates(H)
        return λ[i], ψ[i]
    end
    ens_states = map(sys.subsys, idx) do sub, i
        λ_sub, ψ_sub, _ = eigenstates(sub.H)
        return λ_sub[i], ψ_sub[i]
    end
    λ, ψ = first.(ens_states) |> prod, last.(ens_states) |> x -> tensor(x...)
    U = sys.U_basis
    if isnothing(U)
        return λ, ψ
    else
        return λ, U'*ψ
    end
end

function state_idx(
    sys::CompositeQuantumSystem,
    lbl::StateIndex,
    bare::Bool=false;
    bare_to_dressed_tracking::Symbol=:spectrum,
)
    _check_tracking(bare_to_dressed_tracking)
    idx = _state_index(lbl, sys)
    if length(idx) == 1
        return idx[1]
    elseif !bare && bare_to_dressed_tracking == :spectrum
        return _spectrum_idx(sys, idx)
    elseif !bare && !ismissing(sys.state_lookup[idx...])
        return sys.state_lookup[idx...]
    else
        H = bare ? tensor((sub.H for sub in sys.subsys)...) : sys.H
        _, ψ, _ = eigenstates(H)
        _, ψbare = _bare_state(sys, idx)
        best_idx = abs.(ψbare' * ψn for ψn in ψ) |> argmax
        if !bare
            sys.state_lookup[idx...] = best_idx
        end
        return best_idx
    end
end

# Energy rank of the bare product state `idx` among all bare product states.
function _spectrum_idx(sys::CompositeQuantumSystem, idx::AbstractVector{<:Int})
    sub_λ = [real.(first(eigenstates(sub.H))) for sub in sys.subsys]
    # first subsystem index varies fastest, matching the kron ordering of `tensor`
    inds = CartesianIndices(Tuple(length.(sub_λ)))
    energies = [sum(sub_λ[n][I[n]] for n in eachindex(sub_λ)) for I in inds]
    order = sortperm(vec(energies))  # stable: ties broken by label order
    target = LinearIndices(inds)[idx...]
    return findfirst(==(target), order)
end

energy(sys::CompositeQuantumSystem, lbl::StateIndex, bare::Bool=false; kwargs...) =
    state(sys, lbl, bare; kwargs...)[1]
energy(sys::SimpleQuantumSystem, lbl::StateIndex; kwargs...) = state(sys, lbl; kwargs...)[1]

## Infix Operators #######################

function QuantumToolbox.:⊗(A::SimpleQuantumSystem, B::SimpleQuantumSystem)
    # --- 1. Require names ---
    isnothing(A.name) && throw(ArgumentError(
        "Left-hand SimpleQuantumSystem must have a non-nothing `name`"
    ))
    isnothing(B.name) && throw(ArgumentError(
        "Right-hand SimpleQuantumSystem must have a non-nothing `name`"
    ))

    nameA = A.name
    nameB = B.name

    # --- 2. Build subsystem NamedTuple (order matters!) ---
    subsys = (; nameA => A, nameB => B)

    # --- 3. Build composite Hamiltonian ---
    H = tensor(A.H, qeye(B)) + tensor(qeye(A), B.H)

    ops = (; H = H)

    # --- 4. Build params NamedTuple ---
    params = (; nameA => A.params, nameB => B.params)

    # --- 5. Construct CompositeQuantumSystem ---
    return CompositeQuantumSystem(
        subsys,
        ops;
        params = params,
        name = Symbol(nameA, :_, nameB),
    )
end

function QuantumToolbox.:⊗(A::SimpleQuantumSystem, B::CompositeQuantumSystem)
    nameA = _require_name(A)

    # Subsystems (A comes first)
    subsys = (; nameA => A, B.subsys...)

    # Hamiltonian
    H = tensor(A.H, qeye(B)) + tensor(qeye(A), B.H)

    ops = (; H = H)

    # Params
    params = (; nameA => A.params, B.params...)

    return CompositeQuantumSystem(
        subsys,
        ops;
        params = params,
        name = Symbol(nameA, :_, B.name),
    )
end

function QuantumToolbox.:⊗(A::CompositeQuantumSystem, B::SimpleQuantumSystem)
    nameB = _require_name(B)

    subsys = (; A.subsys..., nameB => B)

    H = tensor(A.H, qeye(B)) + tensor(qeye(A), B.H)

    ops = (; H = H)

    params = (; A.params..., nameB => B.params)

    return CompositeQuantumSystem(
        subsys,
        ops;
        params = params,
        name = Symbol(A.name, :_, nameB),
    )
end

function QuantumToolbox.:⊗(A::CompositeQuantumSystem, B::CompositeQuantumSystem)
    # Ensure no name collisions
    overlap = intersect(keys(A.subsys), keys(B.subsys))
    isempty(overlap) || throw(ArgumentError(
        "Subsystem name collision in composite product: $(collect(overlap))"
    ))

    subsys = (; A.subsys..., B.subsys...)

    H = tensor(A.H, qeye(B)) + tensor(qeye(A), B.H)

    ops = (; H = H)

    params = (; A.params..., B.params...)

    return CompositeQuantumSystem(
        subsys,
        ops;
        params = params,
        name = Symbol(A.name, :_, B.name),
    )
end

function Base.:+(sys::SimpleQuantumSystem, Op::AbstractQuantumObject)
    # Dimension check (QuantumToolbox will throw if incompatible,
    # but we fail early with a clearer message)
    sys.H.dimensions == Op.dimensions ||
        throw(DimensionMismatch(
            "Operator dimensions do not match system Hamiltonian"
        ))

    new_H = sys.H + Op
    new_ops = (; sys.ops..., H = new_H)

    return SimpleQuantumSystem(
        new_ops;
        params = sys.params,
        name   = sys.name,
    )
end

function Base.:+(sys::CompositeQuantumSystem, Op::AbstractQuantumObject)
    sys.H.dimensions == Op.dimensions ||
        throw(DimensionMismatch(
            "Operator dimensions do not match composite Hamiltonian"
        ))

    new_H = sys.H + Op
    new_ops = (; sys.ops..., H = new_H)

    return CompositeQuantumSystem(
        sys.subsys,
        new_ops;
        params      = sys.params,
        U_basis     = sys.U_basis,
        dressed_dim = sys.dressed_dim,
        name        = sys.name,
    )
end

function Base.:+(sys::SimpleQuantumSystem, ops::NamedTuple)
    _check_ops_addition(sys, ops)

    new_ops = merge(sys.ops, ops)

    return SimpleQuantumSystem(
        new_ops;
        params = sys.params,
        name   = sys.name,
    )
end

function Base.:+(sys::CompositeQuantumSystem, ops::NamedTuple)
    _check_ops_addition(sys, ops)

    new_ops = merge(sys.ops, ops)

    return CompositeQuantumSystem(
        sys.subsys,
        new_ops;
        params      = sys.params,
        U_basis     = sys.U_basis,
        dressed_dim = sys.dressed_dim,
        name        = sys.name,
    )
end

Base.:+(ops::NamedTuple, sys::AbstractQuantumSystem) = sys + ops
Base.:+(Op::AbstractQuantumObject, sys::AbstractQuantumSystem) = sys + Op

Base.:-(sys::AbstractQuantumSystem, Op::AbstractQuantumObject) = sys +(-1*Op)
