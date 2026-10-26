export CNode, CEdge, CGraph, CVector, CAdjoint, basis_vectors, relabel

"""
    CNode(label, Cg; Lg=nothing)

A circuit node named `label` with capacitance `Cg` (fF) to ground and, optionally, an inductor
`Lg` (nH) to ground (`nothing` means none). `Lg` must be nonzero, since zero would be a short.
"""
struct CNode{CT<:Real, LT<:Real}
    label::Symbol
    Cg::CT
    Lg::Union{LT, Nothing}
end
function CNode(label::Symbol, Cg::Real; Lg::Union{Real, Nothing}=nothing)
    return CNode{typeof(Cg), _elt_type(Lg, Cg)}(label, Cg, Lg)
end

"""
    CEdge(nodes; C=nothing, L=nothing)
    CEdge(a, b; C=nothing, L=nothing)

A branch between the nodes labeled `a` and `b`, with an optional capacitor `C` (fF) and an
optional inductor `L` (nH) in parallel. `nothing` means the element is absent, so
`CEdge(a, b; C=10.0)` is a pure capacitor and `CEdge(a, b; L=5.0)` a pure inductor. An edge needs
at least one of the two. `L` must be nonzero (`L = 0` would be a short). It is usually positive,
but negative values are allowed.
"""
struct CEdge{CT<:Real, LT<:Real}
    nodes::Tuple{Symbol, Symbol}
    C::Union{CT, Nothing}
    L::Union{LT, Nothing}
end
# the type parameter of an absent element follows the present one
_elt_type(x, y) = isnothing(x) ? (isnothing(y) ? Float64 : typeof(y)) : typeof(x)
function CEdge(nodes::Tuple{Symbol, Symbol}; C::Union{Real, Nothing}=nothing,
               L::Union{Real, Nothing}=nothing)
    return CEdge{_elt_type(C, L), _elt_type(L, C)}(nodes, C, L)
end
CEdge(a::Symbol, b::Symbol; C=nothing, L=nothing) = CEdge((a, b); C, L)

"""
    CGraph(nodes, edges=[]; tol=1e-9, simplify=identity)
    CGraph(labels, C, L_inv; tol=1e-9, simplify=identity)

A circuit, described by nodes and edges and by the two matrices that follow from them: the
capacitance matrix `C` and the inverse-inductance matrix `L_inv`. The fields are

- `nodes::Vector{CNode}` and `edges::Vector{CEdge}`, with at most one edge between two nodes;
- `C` and `L_inv`, stored in full, with rows and columns in the order of `nodes`, so the `k`-th
  node is the `k`-th basis vector;
- `tol`, the tolerance used to decide that a matrix is symmetric and that an entry is zero;
- `simplify`, a function applied to symbolic quantities (see below).

`C` and `L_inv` must be symmetric to within `tol`, otherwise an `ArgumentError` is thrown. A matrix
that passes the check is symmetrized exactly. Symbolic entries cannot be evaluated, so they are only
checked if a `simplify` function is given: `C[i, j] - C[j, i]` must then simplify to zero.

**Symbolic values.** `simplify` (default `identity`) is applied to every symbolic matrix entry that
is computed, to the ground values and edge values that are read back from the matrices, and to the
differences in the symmetry check. It is never applied to numbers. Without it, products such as
`V' * graph * V` leave expressions like `2/Lj - 2/Lj` that are not recognised as zero, so the nodes
and edges read back from them contain spurious entries. A good choice is
`x -> Symbolics.simplify(x; expand=true)`. The function is stored in the graph and reused by
`M * graph * Q`, `relabel` and the other operations that build a new graph.

**From nodes and edges.** Each `CNode` contributes its ground capacitance `Cg` and `1/Lg` to the
diagonal, and each `CEdge` adds `C` (and `1/L`) to both diagonal entries and subtracts it from the
two off-diagonal entries. Parallel edges between the same two nodes are merged: capacitances add and
inductances combine as `1/L = 1/L₁ + 1/L₂`. Node labels must be unique, edges may not be
self-loops (a ground capacitor belongs in the node's `Cg`), must connect nodes in the list, must
have a capacitor or an inductor, and no inductance may be zero. Numeric types are promoted to one
type, so values may be symbolic (e.g. `Symbolics.Num`).

**From matrices.** The nodes and edges are read back from `C` and `L_inv`: `C_ij = -C[i, j]` and
`1/L_ij = -L_inv[i, j]` on the edge between nodes `i` and `j`, and the ground values are the row
sums, `Cg_i = Σⱼ C[i, j]` and `1/Lg_i = Σⱼ L_inv[i, j]`. Entries within `tol` of zero (exactly zero
for symbolic values) mean no edge, no inductor, or no ground capacitance. The result can contain
negative capacitances and inductances, e.g. after a change of basis, so it need not be a physical
circuit.

    M * graph * Q

For square matrices `M` and `Q` this returns a new graph with `C → M*C*Q` and `L_inv → M*L_inv*Q`,
the same labels, `tol` and `simplify`, and with the nodes and edges read back from the new matrices.
The result must be symmetric to within `tol`. A change of basis `φ = Vψ` is `V' * graph * V`.

    graph * v
    v' * graph

For a vector `v` these return a [`CVector`](@ref) and a [`CAdjoint`](@ref), which hold the two
products side by side.
"""
struct CGraph{T<:Real}
    nodes::Vector{CNode{T, T}}
    edges::Vector{CEdge{T, T}}
    C::Matrix{T}
    L_inv::Matrix{T}
    tol::Float64
    simplify::Function

    function CGraph{T}(nodes, edges, C, L_inv, tol, simplify) where {T<:Real}
        labels = [nd.label for nd in nodes]
        allunique(labels) || throw(ArgumentError("duplicate node labels in $labels"))
        n = length(labels)
        for (name, M) in (("C", C), ("L_inv", L_inv))
            size(M) == (n, n) || throw(DimensionMismatch(
                "`$name` has size $(size(M)) but the graph has $n nodes, so it must be ($n, $n)"))
        end
        tol >= 0 || throw(ArgumentError("`tol` must not be negative, got $tol"))
        return new{T}(
            collect(CNode{T, T}, nodes), collect(CEdge{T, T}, edges),
            Matrix{T}(_checked_symmetric(C, tol, "C", simplify)),
            Matrix{T}(_checked_symmetric(L_inv, tol, "L_inv", simplify)),
            Float64(tol), simplify,
        )
    end
end

# only concrete numbers can be checked; symbolic values are left alone
_is_concrete(x) = x isa Union{AbstractFloat, Integer, Rational}
_zero_L(x) = _is_concrete(x) && iszero(x)

# `simplify` is only for symbolic values: numbers are returned as they are
_simp(simplify, x) = _is_concrete(x) ? x : simplify(x)
_simp_all(simplify, A) = simplify === identity ? A : map(x -> _simp(simplify, x), A)

# `A` itself if it is symmetric; its symmetrization if it is symmetric to within `tol`. Symbolic
# matrices are only checked if there is a `simplify` function to show that differences are zero.
function _checked_symmetric(A, tol, name, simplify=identity)
    if !_is_concrete(zero(eltype(A)))
        simplify === identity && return A
        for j in axes(A, 2), i in 1:j-1
            d = _simp(simplify, A[i, j] - A[j, i])
            isequal(d, 0) || throw(ArgumentError(
                "`$name` is not symmetric: $name[$i, $j] - $name[$j, $i] = $d does not simplify " *
                "to zero"))
        end
        return A
    end
    asym = maximum(abs, A .- A'; init=zero(real(eltype(A))))
    asym <= tol || throw(ArgumentError(
        "`$name` is not symmetric: the largest |$name[i, j] - $name[j, i]| is $asym, above " *
        "tol = $tol"))
    iszero(asym) && return A
    return (A .+ A') ./ 2
end

# an entry counts as zero if it is within `tol` of zero, or literally zero for symbolic values
_is_zero(x, tol) = _is_concrete(x) ? abs(x) <= tol : isequal(x, 0)

_edge_key(e::CEdge) = e.nodes[1] < e.nodes[2] ? e.nodes : (e.nodes[2], e.nodes[1])

_parallel_C(::Nothing, ::Nothing) = nothing
_parallel_C(::Nothing, C) = C
_parallel_C(C, ::Nothing) = C
_parallel_C(C1, C2) = C1 + C2

_parallel_L(::Nothing, ::Nothing) = nothing
_parallel_L(::Nothing, L) = L
_parallel_L(L, ::Nothing) = L
function _parallel_L(L1, L2)
    s = L1 + L2
    # equal and opposite inductors cancel to an open circuit, which is "no inductor"
    _is_concrete(s) && iszero(s) && return nothing
    return L1 * L2 / s
end

function _check_edge(labels, e::CEdge)
    a, b = e.nodes
    a == b && throw(ArgumentError(
        "edge $(a)–$(b) is a self-loop; use the node's `Cg` for a capacitor to ground"))
    for l in (a, b)
        l in labels || throw(ArgumentError("edge $(a)–$(b) refers to unknown node :$(l)"))
    end
    isnothing(e.C) && isnothing(e.L) && throw(ArgumentError(
        "edge $(a)–$(b) has neither a capacitor `C` nor an inductor `L`"))
    _zero_L(e.L) && throw(ArgumentError("edge $(a)–$(b) has zero inductance (a short)"))
    return nothing
end

_convert_opt(::Type{T}, x) where {T} = isnothing(x) ? nothing : convert(T, x)

function CGraph(nodes::AbstractVector{<:CNode}, edges::AbstractVector{<:CEdge}=CEdge{Bool,Bool}[];
                tol::Real=1e-9, simplify::Function=identity)
    isempty(nodes) && throw(ArgumentError("a CGraph needs at least one node"))
    labels = [nd.label for nd in nodes]
    allunique(labels) || throw(ArgumentError("duplicate node labels in $labels"))
    for nd in nodes
        _zero_L(nd.Lg) && throw(ArgumentError(
            "node :$(nd.label) has zero ground inductance (a short)"))
    end
    foreach(e -> _check_edge(labels, e), edges)

    # one numeric type for everything, inverse inductances included
    types = DataType[]
    for nd in nodes
        push!(types, typeof(nd.Cg))
        isnothing(nd.Lg) || push!(types, typeof(inv(nd.Lg)))
    end
    for e in edges
        isnothing(e.C) || push!(types, typeof(e.C))
        isnothing(e.L) || push!(types, typeof(inv(e.L)))
    end
    T = promote_type(types...)

    new_nodes = [CNode{T, T}(nd.label, nd.Cg, _convert_opt(T, nd.Lg)) for nd in nodes]
    merged = CEdge{T, T}[]
    index = Dict{Tuple{Symbol, Symbol}, Int}()
    for e in edges
        key = _edge_key(e)
        c, l = _convert_opt(T, e.C), _convert_opt(T, e.L)
        if haskey(index, key)
            old = merged[index[key]]
            merged[index[key]] = CEdge{T, T}(old.nodes, _parallel_C(old.C, c), _parallel_L(old.L, l))
        else
            push!(merged, CEdge{T, T}(e.nodes, c, l))
            index[key] = length(merged)
        end
    end

    idx = Dict(label => i for (i, label) in enumerate(labels))
    n = length(labels)
    C, L_inv = zeros(T, n, n), zeros(T, n, n)
    for (i, nd) in enumerate(new_nodes)
        C[i, i] += nd.Cg
        isnothing(nd.Lg) || (L_inv[i, i] += inv(nd.Lg))
    end
    for e in merged
        i, j = idx[e.nodes[1]], idx[e.nodes[2]]
        for (M, w) in ((C, e.C), (L_inv, isnothing(e.L) ? nothing : inv(e.L)))
            isnothing(w) && continue
            M[i, j] -= w
            M[j, i] -= w
            M[i, i] += w
            M[j, j] += w
        end
    end
    return CGraph{T}(new_nodes, merged, _simp_all(simplify, C), _simp_all(simplify, L_inv), tol,
                     simplify)
end

function CGraph(labels::AbstractVector{Symbol}, C::AbstractMatrix{<:Real},
                L_inv::AbstractMatrix{<:Real}; tol::Real=1e-9, simplify::Function=identity)
    T = promote_type(eltype(C), eltype(L_inv))
    T <: Integer && (T = Float64)             # inverse inductances are generally not integers
    n = length(labels)
    allunique(labels) || throw(ArgumentError("duplicate node labels in $(collect(labels))"))
    for (name, M) in (("C", C), ("L_inv", L_inv))
        size(M) == (n, n) || throw(DimensionMismatch(
            "`$name` has size $(size(M)) but the graph has $n nodes, so it must be ($n, $n)"))
    end
    # the nodes and edges are read from the simplified, symmetrized matrices
    C = _checked_symmetric(_simp_all(simplify, C), tol, "C", simplify)
    L_inv = _checked_symmetric(_simp_all(simplify, L_inv), tol, "L_inv", simplify)

    nodes = map(1:n) do i
        Cg = _simp(simplify, sum(@view C[i, :]))
        invLg = _simp(simplify, sum(@view L_inv[i, :]))
        CNode{T, T}(labels[i], _is_zero(Cg, tol) ? zero(T) : Cg,
                    _is_zero(invLg, tol) ? nothing : _simp(simplify, inv(invLg)))
    end
    edges = CEdge{T, T}[]
    for i in 1:n, j in i+1:n
        c, invl = -C[i, j], -L_inv[i, j]
        czero, lzero = _is_zero(c, tol), _is_zero(invl, tol)
        (czero && lzero) && continue
        push!(edges, CEdge{T, T}((labels[i], labels[j]), czero ? nothing : c,
                                 lzero ? nothing : _simp(simplify, inv(invl))))
    end
    return CGraph{T}(nodes, edges, C, L_inv, tol, simplify)
end

# `==` on numbers, and structural equality (`isequal`) where that gives a symbolic expression
function _entries_equal(A, B)
    size(A) == size(B) || return false
    return all(zip(A, B)) do (x, y)
        r = x == y
        return r isa Bool ? r : isequal(x, y)
    end
end
function Base.:(==)(a::CGraph, b::CGraph)
    return [nd.label for nd in a.nodes] == [nd.label for nd in b.nodes] &&
           _entries_equal(a.C, b.C) && _entries_equal(a.L_inv, b.L_inv)
end

function Base.:*(M::AbstractMatrix{<:Real}, graph::CGraph, Q::AbstractMatrix{<:Real})
    labels = [nd.label for nd in graph.nodes]
    return CGraph(labels, M * graph.C * Q, M * graph.L_inv * Q; tol=graph.tol,
                  simplify=graph.simplify)
end

"""
    CVector(C, L_inv)

The result of `graph * v` for a vector `v`: a wrapper around the two column vectors
`C = graph.C * v` and `L_inv = graph.L_inv * v`. It multiplies with plain vectors and adjoints
like the pair of vectors it holds, and the two results come back as a `NamedTuple` `(C = …, L = …)`:

- `u' * x` is `(C = u' * x.C, L = u' * x.L_inv)`, a pair of scalars.
- `x * u'` is `(C = x.C * u', L = x.L_inv * u')`, a pair of outer-product matrices.

`x'` is the [`CAdjoint`](@ref) of the two adjoint vectors.
"""
struct CVector{A<:AbstractVector, B<:AbstractVector}
    C::A
    L_inv::B
end

"""
    CAdjoint(C, L_inv)

The result of `v' * graph` for a vector `v`: a wrapper around the two row vectors
`C = v' * graph.C` and `L_inv = v' * graph.L_inv` (each an `Adjoint` of a vector). It multiplies
with plain vectors and adjoints like the pair it holds, returning a `NamedTuple` `(C = …, L = …)`:

- `x * w` is `(C = x.C * w, L = x.L_inv * w)`, a pair of scalars.
- `w * x` is `(C = w * x.C, L = w * x.L_inv)`, a pair of outer-product matrices.

So `a' * graph * b` is `(C = a' * graph.C * b, L = a' * graph.L_inv * b)`. `x'` is the
[`CVector`](@ref) of the two adjoint row vectors.
"""
struct CAdjoint{A<:Adjoint{<:Any, <:AbstractVector}, B<:Adjoint{<:Any, <:AbstractVector}}
    C::A
    L_inv::B
end

Base.adjoint(x::CVector) = CAdjoint(x.C', x.L_inv')
Base.adjoint(x::CAdjoint) = CVector(x.C', x.L_inv')
function Base.:(==)(a::T, b::T) where {T<:Union{CVector, CAdjoint}}
    return _entries_equal(a.C, b.C) && _entries_equal(a.L_inv, b.L_inv)
end

# `v'` of a vector is also an `AbstractMatrix`, so the element type is constrained to keep these
# methods strictly more specific than the matrix products above
const _RealAdjoint = Adjoint{<:Real, <:AbstractVector{<:Real}}

Base.:*(graph::CGraph, v::AbstractVector{<:Real}) = CVector(graph.C * v, graph.L_inv * v)
Base.:*(v::_RealAdjoint, graph::CGraph) = CAdjoint(v * graph.C, v * graph.L_inv)

Base.:*(u::_RealAdjoint, x::CVector) = (C=u * x.C, L=u * x.L_inv)
Base.:*(x::CAdjoint, w::AbstractVector{<:Real}) = (C=x.C * w, L=x.L_inv * w)
Base.:*(x::CVector, u::_RealAdjoint) = (C=x.C * u, L=x.L_inv * u)
Base.:*(w::AbstractVector{<:Real}, x::CAdjoint) = (C=w * x.C, L=w * x.L_inv)

"""
    basis_vectors(graph::CGraph) -> NamedTuple

The basis vectors of `graph`, as a `NamedTuple` of `Vector{Int}`s keyed by the node labels. The
`k`-th node of `graph.nodes` corresponds to the `k`-th unit vector, so `basis_vectors(g).A1` has
a `1` in the position of node `:A1` in `g.nodes` and zeros elsewhere. The keys follow the order of
`g.nodes`, and `(; A1, B2) = basis_vectors(g)` unpacks them into variables.
"""
function basis_vectors(graph::CGraph)
    n = length(graph.nodes)
    vectors = ntuple(k -> (v = zeros(Int, n); v[k] = 1; v), n)
    return NamedTuple{Tuple(nd.label for nd in graph.nodes)}(vectors)
end

"""
    relabel(graph::CGraph, labels::AbstractVector{Symbol}) -> CGraph

A new graph with the same matrices as `graph`, where the `k`-th node is renamed to `labels[k]`, in
the nodes and in the edges. `labels` must contain one unique label per node. The order of the
nodes, and so the basis, does not change. This is how to name the new variables after a change of
basis such as `V' * graph * V`.
"""
function relabel(graph::CGraph{T}, labels::AbstractVector{Symbol}) where {T}
    length(labels) == length(graph.nodes) || throw(DimensionMismatch(
        "got $(length(labels)) labels for a graph with $(length(graph.nodes)) nodes"))
    allunique(labels) || throw(ArgumentError("`labels` has duplicate names: $(collect(labels))"))
    new_label = Dict(nd.label => label for (nd, label) in zip(graph.nodes, labels))
    nodes = [CNode{T, T}(new_label[nd.label], nd.Cg, nd.Lg) for nd in graph.nodes]
    edges = [CEdge{T, T}((new_label[e.nodes[1]], new_label[e.nodes[2]]), e.C, e.L)
             for e in graph.edges]
    return CGraph{T}(nodes, edges, graph.C, graph.L_inv, graph.tol, graph.simplify)
end
