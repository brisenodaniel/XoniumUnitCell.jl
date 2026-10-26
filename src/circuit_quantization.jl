export CPW_freq2cap, Cinv2E_C, Cinv2J, Cinv2g
export CNode, CEdge, CGraph, node_index, add_node!, add_edge!
export cap_matrix, inv_inductance_matrix, to_matrix, ModeBasis, sum_diff_transform, transform_cap_matrix
export GraphVector, to_vector, node_vectors, reorder_nodes, relabel, transform, mode_index, diff_mode_indices, truncate_to_diff_modes, Cinv_params

# Physical constants (CODATA 2018 exact SI values)
const _e_charge = 1.602176634e-19  # C
const _h_planck = 6.62607015e-34   # J s
const _ħ_planck = _h_planck / 2π   # J s


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
but graphs produced by [`transform`](@ref) can contain negative values.
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

_edge_key(a::Symbol, b::Symbol) = a < b ? (a, b) : (b, a)
_edge_key(e::CEdge) = _edge_key(e.nodes...)

# convert an optional value, keeping `nothing`
_convert_opt(::Type{T}, x) where {T} = isnothing(x) ? nothing : convert(T, x)
function _edge_as(::Type{E}, e::CEdge) where {E}
    return CEdge{E,E}(e.nodes, _convert_opt(E, e.C), _convert_opt(E, e.L))
end

_node_as(::Type{N}, nd::CNode) where {N} = CNode{N,N}(nd.label, nd.Cg, _convert_opt(N, nd.Lg))

# capacitors in parallel add; inductors in parallel combine as 1/L = 1/L₁ + 1/L₂
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

function _merge_edges(e1::CEdge{E,E}, e2::CEdge) where {E}
    return CEdge{E,E}(
        e1.nodes,
        _convert_opt(E, _parallel_C(e1.C, e2.C)),
        _convert_opt(E, _parallel_L(e1.L, e2.L)),
    )
end

# only concrete numbers can be checked; symbolic values are left alone
_is_concrete(x) = x isa Union{AbstractFloat, Integer, Rational}
_zero_L(x) = _is_concrete(x) && iszero(x)

"""
    CGraph(nodes, edges=[])

Graph of a circuit: the nodes carry the capacitance (and optionally an inductance) to ground, the
edges the coupling capacitors and inductors. The order of `nodes` fixes the basis: node `k` is
the `k`-th basis vector, so it is row and column `k` of every matrix built from the graph. Use
[`reorder_nodes`](@ref) to change it.

Node labels must be unique. Edges may not be self-loops (a ground capacitor belongs in the node's
`Cg`), must connect nodes in the graph, must have a capacitor or an inductor, and no inductance
may be zero. Parallel branches between the same two nodes are merged:
capacitances add and inductances combine in parallel. Numeric types are promoted to one type for the nodes and one for the
edges (shared by `C` and `L`), so a graph may hold symbolic values (e.g. `Symbolics.Num`).
`g[label]` returns the node `label`. `g.C` and `g.L_inv` return the capacitance and
inverse-inductance matrices in the order of `g.nodes`, see [`to_matrix`](@ref). They are
recomputed on every access, so store the result if you need it repeatedly.
"""
struct CGraph{N<:Real,E<:Real}
    nodes::Vector{CNode{N,N}}
    edges::Vector{CEdge{E,E}}

    function CGraph{N,E}(nodes, edges) where {N<:Real,E<:Real}
        labels = [nd.label for nd in nodes]
        allunique(labels) || throw(ArgumentError("duplicate node labels in $labels"))
        for nd in nodes
            _zero_L(nd.Lg) && throw(ArgumentError(
                "node :$(nd.label) has zero ground inductance (a short)"))
        end
        merged = CEdge{E,E}[]
        index = Dict{Tuple{Symbol,Symbol},Int}()
        for e in edges
            _check_edge(labels, e)
            key = _edge_key(e)
            if haskey(index, key)
                merged[index[key]] = _merge_edges(merged[index[key]], e)
            else
                push!(merged, _edge_as(E, e))
                index[key] = length(merged)
            end
        end
        return new{N,E}(collect(CNode{N,N}, nodes), merged)
    end
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

function CGraph(nodes::AbstractVector{<:CNode}, edges::AbstractVector{<:CEdge}=CEdge{Bool,Bool}[])
    isempty(nodes) && throw(ArgumentError("a CGraph needs at least one node"))
    node_types = DataType[]
    for nd in nodes
        push!(node_types, typeof(nd.Cg))
        isnothing(nd.Lg) || push!(node_types, typeof(nd.Lg))
    end
    N = promote_type(node_types...)
    edge_types = DataType[]
    for e in edges
        isnothing(e.C) || push!(edge_types, typeof(e.C))
        isnothing(e.L) || push!(edge_types, typeof(e.L))
    end
    E = isempty(edge_types) ? N : promote_type(edge_types...)
    return CGraph{N,E}(
        [_node_as(N, nd) for nd in nodes],
        [_edge_as(E, e) for e in edges],
    )
end

Base.getindex(g::CGraph, label::Symbol) = g.nodes[node_index(g, label)]

"""
    node_index(g, label)

Row/column of node `label` in `cap_matrix(g)`.
"""
function node_index(g::CGraph, label::Symbol)
    idx = findfirst(nd -> nd.label == label, g.nodes)
    isnothing(idx) && throw(KeyError(label))
    return idx
end

"""
    add_node!(g, node)

Add `node` to `g`, replacing the existing node with the same label if there is one. The
values are converted to the node type of `g`. Returns `g`.
"""
function add_node!(g::CGraph{N}, node::CNode) where {N}
    _zero_L(node.Lg) && throw(ArgumentError(
        "node :$(node.label) has zero ground inductance (a short)"))
    new = _node_as(N, node)
    idx = findfirst(nd -> nd.label == node.label, g.nodes)
    isnothing(idx) ? push!(g.nodes, new) : (g.nodes[idx] = new)
    return g
end

"""
    add_edge!(g, edge)

Add a branch between two nodes already in `g`. If the two nodes are already connected, the
branches are in parallel: capacitances add and inductances combine as `1/L = 1/L₁ + 1/L₂`. The
values are converted to the edge type of `g`. Returns `g`.
"""
function add_edge!(g::CGraph{N,E}, edge::CEdge) where {N,E}
    _check_edge([nd.label for nd in g.nodes], edge)
    key = _edge_key(edge)
    idx = findfirst(e -> _edge_key(e) == key, g.edges)
    if isnothing(idx)
        push!(g.edges, _edge_as(E, edge))
    else
        g.edges[idx] = _merge_edges(g.edges[idx], edge)
    end
    return g
end

"""
    CPW_freq2cap(freq; Z0=50)

Capacitance (fF) to ground of a λ/4 CPW resonator with frequency `freq` (GHz, ωr/2π) and
impedance `Z0` (Ω).
"""
CPW_freq2cap(freq; Z0=50) = 1e15 / (8 * freq * 1e9 * Z0)

"""
    Cinv2E_C(Cinv)

Convert an inverse-capacitance matrix element `Cinv` (1/fF) to a charging energy `E_C/h` (GHz).
"""
Cinv2E_C(Cinv) = 1e-9 * 0.5 * _e_charge^2 * Cinv * 1e15 / _h_planck

"""
    Cinv2J(Cinv)

Convert an inverse-capacitance matrix element `Cinv` (1/fF) to a qubit-qubit capacitive
coupling `J/h` (MHz).
"""
Cinv2J(Cinv::Real) = 1e-6 * 4 * _e_charge^2 * Cinv * 1e15 / _h_planck

"""
    Cinv2g(Cinv)

Convert an inverse-capacitance matrix element `Cinv` (1/fF) to a qubit-resonator coupling
`g/2π` (MHz).
"""
Cinv2g(Cinv::Real) = 1e-6 * _e_charge * sqrt(π * _ħ_planck / 100) * Cinv * 1e15 / _h_planck

"""
    reorder_nodes(graph::CGraph, order::AbstractVector{Symbol}) -> CGraph

A new graph with the same nodes and edges as `graph`, with the nodes in the order `order`. `order`
must list every node label exactly once. Since the order of the nodes is the basis of the
matrices of the graph, `reorder_nodes(g, order).C == g.C[perm, perm]` for the corresponding
permutation `perm`.
"""
function reorder_nodes(graph::CGraph, order::AbstractVector{Symbol})
    labels = [nd.label for nd in graph.nodes]
    allunique(order) || throw(ArgumentError("`order` has duplicate node labels: $(collect(order))"))
    absent, unknown = setdiff(labels, order), setdiff(order, labels)
    (isempty(absent) && isempty(unknown)) || throw(ArgumentError(
        "`order` must list every node of the graph exactly once; missing $(absent), unknown $(unknown)"))
    nodes = [graph[label] for label in order]
    return typeof(graph)(nodes, graph.edges)
end

"""
    relabel(graph::CGraph, labels::AbstractVector{Symbol}) -> CGraph

A new graph with the same capacitances, inductances, connections and node order as `graph`, where
the `k`-th node is renamed to `labels[k]`. `labels` must contain one unique label per node. Use
[`reorder_nodes`](@ref) to change the order instead. `GraphVector`s made for the old labels no
longer match the new graph, so rebuild them with [`to_vector`](@ref).
"""
function relabel(graph::CGraph{N,E}, labels::AbstractVector{Symbol}) where {N,E}
    length(labels) == length(graph.nodes) || throw(DimensionMismatch(
        "got $(length(labels)) labels for a graph with $(length(graph.nodes)) nodes"))
    allunique(labels) || throw(ArgumentError("`labels` has duplicate names: $(collect(labels))"))
    new_label = Dict(nd.label => label for (nd, label) in zip(graph.nodes, labels))
    nodes = [CNode{N,N}(new_label[nd.label], nd.Cg, nd.Lg) for nd in graph.nodes]
    edges = [CEdge{E,E}((new_label[e.nodes[1]], new_label[e.nodes[2]]), e.C, e.L)
             for e in graph.edges]
    return CGraph{N,E}(nodes, edges)
end

# Add each edge's weight `w(e)` as a graph Laplacian: -w off the diagonal, +w on both diagonals.
function _add_laplacian!(M, g::CGraph, idx, w)
    for e in g.edges
        we = w(e)
        isnothing(we) && continue
        i, j = idx[e.nodes[1]], idx[e.nodes[2]]
        M[i, j] -= we
        M[j, i] -= we
        M[i, i] += we
        M[j, j] += we
    end
    return M
end

"""
    cap_matrix(g::CGraph) -> C

Maxwell capacitance matrix of the circuit `g`. Off-diagonals are `-C_ij`; each diagonal entry is
the node's ground capacitance plus the capacitors attached to it. Edges without a capacitor do
not contribute. Rows and columns follow the order of `g.nodes`.
"""
function cap_matrix(g::CGraph{N,E}) where {N,E}
    idx = Dict(nd.label => i for (i, nd) in enumerate(g.nodes))
    C = zeros(promote_type(N, E), length(g.nodes), length(g.nodes))
    for nd in g.nodes
        C[idx[nd.label], idx[nd.label]] = nd.Cg
    end
    return _add_laplacian!(C, g, idx, e -> e.C)
end

"""
    inv_inductance_matrix(g::CGraph) -> Linv

Inverse-inductance matrix `L⁻¹` (1/nH) of the circuit `g`, built like the capacitance matrix:
off-diagonals are `-1/L_ij` and each diagonal entry is the sum of `1/L` of the inductors
attached to the node, including its inductor to ground `Lg`. Edges without an inductor do not
contribute. Nodes are ordered as in [`cap_matrix`](@ref). This is not the inverse of a matrix: it is singular whenever no inductor
connects to ground, so it is assembled directly.
"""
function inv_inductance_matrix(g::CGraph{N,E}) where {N,E}
    idx = Dict(nd.label => i for (i, nd) in enumerate(g.nodes))
    Linv = zeros(typeof(inv(one(promote_type(N, E)))), length(g.nodes), length(g.nodes))
    for nd in g.nodes
        isnothing(nd.Lg) || (Linv[idx[nd.label], idx[nd.label]] += inv(nd.Lg))
    end
    return _add_laplacian!(Linv, g, idx, e -> isnothing(e.L) ? nothing : inv(e.L))
end

"""
    to_matrix(circuit_graph::CGraph) -> C, Linv

The capacitance matrix and the inverse-inductance matrix of `circuit_graph`, see
[`cap_matrix`](@ref) and [`inv_inductance_matrix`](@ref). Rows and columns follow the order of
`circuit_graph.nodes`.
"""
function to_matrix(circuit_graph::CGraph)
    return cap_matrix(circuit_graph), inv_inductance_matrix(circuit_graph)
end

# `g.C` and `g.L_inv` are computed from the nodes and edges on each access
function Base.getproperty(g::CGraph, name::Symbol)
    name === :C && return cap_matrix(g)
    name === :L_inv && return inv_inductance_matrix(g)
    return getfield(g, name)
end
Base.propertynames(::CGraph, private::Bool=false) = (:nodes, :edges, :C, :L_inv)

"""
    graph * v
    v' * graph

Apply the matrices of the graph to a vector, returning `(C = …, L = …)` with the capacitance
matrix in `C` and the inverse-inductance matrix in `L`: `graph * v` is
`(C = graph.C * v, L = graph.L_inv * v)` and `v' * graph` is
`(C = v' * graph.C, L = v' * graph.L_inv)`. The length of `v` must be the number of nodes, in the
order of `graph.nodes`, see [`to_vector`](@ref). The matrices are only built once per call.
"""
function Base.:*(graph::CGraph, v::AbstractVector)
    C, Linv = to_matrix(graph)
    return (C=C * v, L=Linv * v)
end

function Base.:*(v::Adjoint{<:Any, <:AbstractVector}, graph::CGraph)
    C, Linv = to_matrix(graph)
    return (C=v * C, L=v * Linv)
end

"""
    GraphVector(data, basis)

A vector together with the basis it is expressed in: `data` holds the numbers, a `Vector` or its
`Adjoint` (a row vector), and `basis` the node labels of the graph, in the order of `graph.nodes`,
that the entries of `data` refer to. Build one with [`to_vector`](@ref). `v'` is the `GraphVector`
with `data'`, so `v.data` is a plain `Vector` again after two adjoints.

Products check that the bases agree and throw an `ArgumentError` otherwise (for example after
[`reorder_nodes`](@ref)):
- `a * b` for two `GraphVector`s is `a.data * b.data`, so `a' * b` is the scalar product.
- `a * v` and `v' * a` with a plain vector `v` are `a.data * v` and `v' * a.data`.
- `graph * a` and `a' * graph` apply the matrices of the graph, as for plain vectors, and give
  `(C = …, L = …)` of `GraphVector`s. Multiplying that by another `GraphVector` on the right (or,
  for `graph * a`, on the left) completes a bilinear form, so `a' * graph * b` is
  `(C = a' * graph.C * b, L = a' * graph.L_inv * b)` as scalars.
"""
struct GraphVector{D<:AbstractVecOrMat}
    data::D
    basis::Vector{Symbol}

    function GraphVector(data::D, basis::AbstractVector{Symbol}) where {D<:AbstractVecOrMat}
        length(data) == length(basis) || throw(DimensionMismatch(
            "`data` has $(length(data)) entries but the basis has $(length(basis)) nodes"))
        return new{D}(data, collect(Symbol, basis))
    end
end
Base.adjoint(v::GraphVector) = GraphVector(v.data', v.basis)
Base.:(==)(a::GraphVector, b::GraphVector) = a.data == b.data && a.basis == b.basis

"""
    to_vector(node::Union{Symbol, CNode}, graph::CGraph) -> GraphVector

The basis vector of `node` in `graph`: a `GraphVector` whose `data` is the unit vector with a `1` in
the position of the node in `graph.nodes`, so that it picks out row and column `k` of the matrices
of the graph. A `CNode` is identified by its label alone. Throws a `KeyError` if the node is not
in the graph. Use `.data` for the plain `Vector{Int}`. There is no method without a graph, since
the position and length of the vector depend on it.
"""
function to_vector(node::Union{Symbol, CNode}, graph::CGraph)
    label = node isa CNode ? node.label : node
    data = zeros(Int, length(graph.nodes))
    data[node_index(graph, label)] = 1
    return GraphVector(data, [nd.label for nd in graph.nodes])
end

"""
    node_vectors(graph::CGraph) -> NamedTuple

The basis vector of every node of `graph`, as a `NamedTuple` of [`GraphVector`](@ref)s keyed by
the node labels, in the order of `graph.nodes`. `node_vectors(g).A1` is `to_vector(:A1, g)`, and
`(; A1, B2) = node_vectors(g)` unpacks them.
"""
function node_vectors(graph::CGraph)
    labels = Tuple(nd.label for nd in graph.nodes)
    return NamedTuple{labels}(Tuple(to_vector(label, graph) for label in labels))
end

function _check_basis(v::GraphVector, graph::CGraph)
    v.basis == [nd.label for nd in graph.nodes] || throw(ArgumentError(
        "vector is expressed in the basis $(v.basis), but the graph's nodes are " *
        "$([nd.label for nd in graph.nodes]); rebuild it with `to_vector`"))
    return nothing
end

function Base.:*(a::GraphVector, b::GraphVector)
    a.basis == b.basis || throw(ArgumentError(
        "vectors are expressed in different bases: $(a.basis) and $(b.basis)"))
    return a.data * b.data
end
Base.:*(a::GraphVector, v::AbstractVector) = a.data * v
Base.:*(v::Adjoint{<:Any, <:AbstractVector}, a::GraphVector) = v * a.data

function Base.:*(graph::CGraph, v::GraphVector)
    _check_basis(v, graph)
    C, Linv = to_matrix(graph)
    return (C=GraphVector(C * v.data, v.basis), L=GraphVector(Linv * v.data, v.basis))
end

function Base.:*(v::GraphVector, graph::CGraph)
    _check_basis(v, graph)
    C, Linv = to_matrix(graph)
    return (C=GraphVector(v.data * C, v.basis), L=GraphVector(v.data * Linv, v.basis))
end

const _GraphVectorPair = NamedTuple{(:C, :L), <:Tuple{GraphVector, GraphVector}}
Base.:*(x::_GraphVectorPair, v::GraphVector) = (C=x.C * v, L=x.L * v)
Base.:*(v::GraphVector, x::_GraphVectorPair) = (C=v * x.C, L=v * x.L)

# zero test for matrix entries: relative tolerance for floats, exact otherwise (e.g. symbolic)
function _zero_test(M; rtol=1e-10)
    if eltype(M) <: AbstractFloat
        tol = rtol * maximum(abs, M; init=zero(eltype(M)))
        return x -> abs(x) <= tol
    end
    return x -> isequal(x, 0)
end

"""
    transform(graph::CGraph, T, relabel=nothing; simplify=identity) -> CGraph

Change of variables `φ = T ψ` (the paper's `V`) applied to `graph`: the capacitance matrix `C`
and inverse-inductance matrix `L⁻¹` become `K = T' C T` and `M = T' L⁻¹ T`, and the result is
returned as a new `CGraph` whose nodes are the new variables, named by `relabel` in the order of
the columns of `T`.

- The rows of `T` follow the order of `graph.nodes` (see [`reorder_nodes`](@ref)).
- `relabel` names the columns of `T`. By default the nodes keep their labels, in the same order,
  so the new node `k` has the label of the old node `k`. `T` need not be square, which restricts
  the circuit to a subspace of the new variables, but then `relabel` is required because the
  old labels no longer fit.
- The graph is read back from `K` and `M` as in [`cap_matrix`](@ref) and
  [`inv_inductance_matrix`](@ref): `C_ij = -K_ij` and `1/L_ij = -M_ij` on the edges, and the
  ground values are the row sums, `Cg_i = Σⱼ K_ij` and `1/Lg_i = Σⱼ M_ij`. Zero entries mean no
  edge or no inductor. Floats are compared to zero relative to the largest entry of the matrix.
  For symbolic values, entries count as zero only if they are literally zero, so pass a function
  such as `x -> Symbolics.simplify(x; expand=true)` as `simplify`; it is applied to every matrix entry and to the row sums.

The result is a container for the transformed matrices and not generally a physical circuit:
capacitances and inductances can be negative. `to_matrix(transform(g, T, relabel))` returns
`(T' C T, T' L⁻¹ T)`.

Note that [`ModeBasis`](@ref) stores the opposite convention, `ψ = T φ`, so its `T` is the inverse
of the one taken here.
"""
function transform(
    graph::CGraph,
    T::AbstractMatrix{<:Real},
    relabel::Union{AbstractVector{Symbol}, Nothing}=nothing;
    simplify=identity,
)
    n, m = size(T)
    n == length(graph.nodes) || throw(DimensionMismatch(
        "`T` has $n rows but the graph has $(length(graph.nodes)) nodes"))
    if isnothing(relabel)
        m == n || throw(ArgumentError(
            "`T` maps to $m new variables but the graph has $n nodes; pass `relabel` to name them"))
        names = [nd.label for nd in graph.nodes]
    else
        names = collect(relabel)
    end
    length(names) == m || throw(DimensionMismatch(
        "`relabel` has $(length(names)) names but `T` has $m columns"))
    allunique(names) || throw(ArgumentError("`relabel` has duplicate names: $names"))

    C, Linv = to_matrix(graph)
    K = simplify.(T' * C * T)
    M = simplify.(T' * Linv * T)
    isC0, isL0 = _zero_test(K), _zero_test(M)

    nodes = map(1:m) do k
        Cg = simplify(sum(@view K[k, :]))
        invLg = simplify(sum(@view M[k, :]))
        CNode(names[k], Cg; Lg=isL0(invLg) ? nothing : inv(invLg))
    end
    edges = CEdge[]
    for i in 1:m, j in i+1:m
        Cij = isC0(K[i, j]) ? nothing : -K[i, j]
        Lij = isL0(M[i, j]) ? nothing : simplify(-inv(M[i, j]))
        (isnothing(Cij) && isnothing(Lij)) || push!(edges, CEdge(names[i], names[j]; C=Cij, L=Lij))
    end
    return CGraph(nodes, edges)
end

"""
    ModeBasis

Change of basis `T` from node fluxes to (sum, difference) modes, with the `names` of the node
pairs. Modes are ordered `[Σ₁, δ₁, Σ₂, δ₂, …]` following `names`. Built by
[`sum_diff_transform`](@ref).
"""
struct ModeBasis
    T::Matrix{Int}
    names::Vector{Symbol}
end

function _pair_position(mb::ModeBasis, name::Symbol)
    k = findfirst(==(name), mb.names)
    isnothing(k) && throw(KeyError(name))
    return k
end

"""
    mode_index(mb, name; kind=:diff)

Index of the `:diff` (default) or `:sum` mode of the pair `name` in the mode basis `mb`.
"""
function mode_index(mb::ModeBasis, name::Symbol; kind::Symbol=:diff)
    kind in (:diff, :sum) || throw(ArgumentError("`kind` must be :diff or :sum, got :$kind"))
    k = _pair_position(mb, name)
    return kind == :diff ? 2k : 2k - 1
end

"""
    diff_mode_indices(mb)

Indices of all difference modes in `mb`.
"""
diff_mode_indices(mb::ModeBasis) = collect(2:2:2length(mb.names))

"""
    sum_diff_transform(g, pairs) -> ModeBasis

Sum/difference modes of the node pairs in `pairs`, a vector of `name => (a, b)`, e.g.
`[:A => (:A1, :A2), :B => (:B1, :B2)]`. The sum mode is `φ_a + φ_b` and the difference mode
`φ_a - φ_b`. Every node of `g` must appear in exactly one pair. The columns of the transform follow
the order of `g.nodes`.
"""
function sum_diff_transform(
    g::CGraph,
    pairs::AbstractVector{<:Pair{Symbol,<:Tuple{Symbol,Symbol}}},
)
    order = [nd.label for nd in g.nodes]
    used = collect(Iterators.flatten(last.(pairs)))
    allunique(used) || throw(ArgumentError("a node appears in more than one pair: $used"))
    Set(used) == Set(order) || throw(ArgumentError(
        "pairs must cover every node of the graph exactly once"))
    allunique(first.(pairs)) || throw(ArgumentError("pair names must be unique"))
    T = zeros(Int, length(order), length(order))
    for (k, (_, (a, b))) in enumerate(pairs)
        ia, ib = findfirst(==(a), order), findfirst(==(b), order)
        T[2k-1, ia] = T[2k-1, ib] = 1
        T[2k, ia] = 1
        T[2k, ib] = -1
    end
    return ModeBasis(T, collect(first.(pairs)))
end

"""
    transform_cap_matrix(C, T)
    transform_cap_matrix(C, mb::ModeBasis)

Capacitance matrix `K = T⁻ᵀ C T⁻¹` in the mode basis `φ_mode = T φ_node`, from
`½ φ̇ₙᵀ C φ̇ₙ` with `φₙ = T⁻¹ φ_mode`.
"""
function transform_cap_matrix(C, T::AbstractMatrix)
    Tinv = inv(Rational{Int}.(T))  # exact, so symbolic entries stay rational
    return Tinv' * C * Tinv
end
transform_cap_matrix(C, mb::ModeBasis) = transform_cap_matrix(C, mb.T)

"""
    truncate_to_diff_modes(K, mb)

Drop the sum modes from `K` *before* inverting. This is not a good approximation: the sum
modes still load the difference modes, so the truncated `inv(K)` underestimates `E_C` and
overestimates `J`. Invert the full `K` and read off the difference-mode block instead.
"""
function truncate_to_diff_modes(K, mb::ModeBasis)
    idx = diff_mode_indices(mb)
    return K[idx, idx]
end

"""
    Cinv_params(Kinv, a, b)
    Cinv_params(Kinv, mb, a, b; truncated=false)

Charging energies `E_C/h` (GHz) of modes `a` and `b` and their capacitive coupling `J/h` (MHz),
read off the inverse capacitance matrix `Kinv` (1/fF), as `(; EC_a, EC_b, J_ab)`. The modes are
given by index, or by the pair names `a`, `b` of `mb` (difference modes). Pass `truncated=true`
if `Kinv` is the inverse of [`truncate_to_diff_modes`](@ref).
"""
function Cinv_params(Kinv, a::Integer, b::Integer)
    return (; EC_a=Cinv2E_C(Kinv[a, a]), EC_b=Cinv2E_C(Kinv[b, b]), J_ab=Cinv2J(Kinv[a, b]))
end

function Cinv_params(Kinv, mb::ModeBasis, a::Symbol, b::Symbol; truncated::Bool=false)
    idx(name) = truncated ? _pair_position(mb, name) : mode_index(mb, name)
    return Cinv_params(Kinv, idx(a), idx(b))
end
