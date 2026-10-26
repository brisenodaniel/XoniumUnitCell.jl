@testset "circuit_quantization" begin
    nodes = [CNode(:A1, 20.0), CNode(:A2, 80.0), CNode(:B1, 20.0), CNode(:B2, 80.0)]
    edges = [CEdge(:A1, :A2; C=10.0), CEdge(:A2, :B2; C=5.0), CEdge(:B1, :B2; C=10.0)]
    g = CGraph(nodes, edges)

    @testset "graph" begin
        @test g[:A2].Cg == 80.0
        @test node_index(g, :B1) == 3
        @test_throws KeyError node_index(g, :Z)
        @test_throws ArgumentError CGraph([CNode(:a, 1.0), CNode(:a, 2.0)])
        @test_throws ArgumentError CGraph(nodes, [CEdge(:A1, :A1; C=1.0)])
        @test_throws ArgumentError CGraph(nodes, [CEdge(:A1, :Z; C=1.0)])
        # parallel capacitors add, in either orientation
        gp = CGraph(nodes, [CEdge(:A1, :A2; C=4.0), CEdge(:A2, :A1; C=6.0)])
        @test only(gp.edges).C == 10.0
        # no edges gives a graph of singletons
        @test cap_matrix(CGraph(nodes)) == Diagonal([20.0, 80.0, 20.0, 80.0])
    end

    @testset "mutation" begin
        gm = CGraph(nodes, edges)
        add_node!(gm, CNode(:A2, 70.0))          # replaces
        @test length(gm.nodes) == 4 && gm[:A2].Cg == 70.0
        add_node!(gm, CNode(:C, 5.0))             # appends
        @test length(gm.nodes) == 5
        add_edge!(gm, CEdge(:A1, :A2; C=1.0))       # merges with the existing 10.0
        @test length(gm.edges) == 3
        @test only(filter(e -> XoniumUnitCell._edge_key(e) == (:A1, :A2), gm.edges)).C == 11.0
        @test_throws ArgumentError add_edge!(gm, CEdge(:A1, :Z; C=1.0))
        # the original graph is untouched
        @test g[:A2].Cg == 80.0 && length(g.nodes) == 4
    end

    @testset "inductive edges" begin
        @test CEdge(:A1, :A2; C=10.0).L === nothing
        @test CEdge(:A1, :A2; L=5.0).C === nothing
        @test CEdge(:A1, :A2; C=10.0, L=5).L === 5      # own type parameter, not converted to C's
        @test CEdge(:A1, :A2; L=5.0) isa CEdge{Float64,Float64}
        @test CEdge(:A1, :A2) isa CEdge                  # constructible, but rejected by a graph
        @test_throws ArgumentError CGraph(nodes, [CEdge(:A1, :A2)])
        # one numeric type for C and L in the graph; `nothing` is preserved
        gl = CGraph(nodes, [CEdge(:A1, :A2; C=10, L=5.0), CEdge(:A2, :B2; C=5.0), CEdge(:B1, :B2; L=2)])
        @test gl isa CGraph{Float64,Float64}
        @test gl.edges[1].C === 10.0 && gl.edges[1].L === 5.0
        @test gl.edges[2].L === nothing
        @test gl.edges[3].C === nothing && gl.edges[3].L === 2.0
        # parallel inductors combine as 1/L = 1/L1 + 1/L2, and capacitors add
        gp = CGraph(nodes, [CEdge(:A1, :A2; C=4.0, L=6.0), CEdge(:A2, :A1; C=6.0, L=3.0)])
        @test only(gp.edges).C == 10.0
        @test only(gp.edges).L ≈ 2.0
        # an absent element in one branch leaves the other branch's value unchanged
        gq = CGraph(nodes, [CEdge(:A1, :A2; C=4.0, L=6.0), CEdge(:A2, :A1; C=6.0)])
        @test only(gq.edges).C == 10.0 && only(gq.edges).L == 6.0
        gadd = CGraph(nodes, [CEdge(:A1, :A2; C=4.0)])
        add_edge!(gadd, CEdge(:A1, :A2; L=8.0))
        @test only(gadd.edges).L == 8.0 && only(gadd.edges).C == 4.0
        add_edge!(gadd, CEdge(:A1, :A2; L=8.0))
        @test only(gadd.edges).L == 4.0
        # L = 0 would be a short; negative L is allowed (it appears in transformed graphs)
        @test_throws ArgumentError CGraph(nodes, [CEdge(:A1, :A2; L=0.0)])
        @test_throws ArgumentError add_edge!(gadd, CEdge(:A1, :A2; L=0))
        @test only(CGraph(nodes, [CEdge(:A1, :A2; L=-1.0)]).edges).L == -1.0
        # equal and opposite inductors in parallel cancel to an open circuit
        gc = CGraph(nodes, [CEdge(:A1, :A2; C=1.0, L=3.0), CEdge(:A1, :A2; L=-3.0)])
        @test only(gc.edges).L === nothing && only(gc.edges).C == 1.0
        # inductors do not enter the capacitance matrix; pure inductors are skipped
        Cl = cap_matrix(gl)
        @test Cl[1, 2] == -10.0
        @test Cl[3, 4] == 0.0
    end

    @testset "to_matrix" begin
        gi = CGraph(nodes, [CEdge(:A1, :A2; C=10.0, L=5.0), CEdge(:A2, :B2; C=5.0),
                            CEdge(:B1, :B2; L=2.0), CEdge(:A2, :B1; L=4.0)])
        Cm, Linv = to_matrix(gi)
        @test Cm == cap_matrix(gi)
        @test Linv == inv_inductance_matrix(gi)
        @test Linv ≈ Linv'
        @test Linv[1, 2] ≈ -1 / 5
        @test Linv[2, 2] ≈ 1 / 5 + 1 / 4          # A2 touches the 5 nH and 4 nH inductors
        @test Linv[3, 4] ≈ -1 / 2
        @test Linv[1, 4] == 0.0                   # A1 and B2 are not connected by an inductor
        @test all(abs.(sum(Linv; dims=2)) .< 1e-12)   # no inductor to ground: singular
        # no inductors at all gives a zero matrix
        @test iszero(inv_inductance_matrix(CGraph(nodes, [CEdge(:A1, :A2; C=1.0)])))
    end

    @testset "reorder_nodes" begin
        gi = CGraph(nodes, [CEdge(:A1, :A2; C=10.0, L=5.0), CEdge(:A2, :B2; C=5.0),
                            CEdge(:B1, :B2; L=2.0), CEdge(:A2, :B1; L=4.0)])
        Cm, Linv = to_matrix(gi)
        ord = [:B2, :A1, :B1, :A2]
        perm = [node_index(gi, l) for l in ord]
        gr = reorder_nodes(gi, ord)
        @test [nd.label for nd in gr.nodes] == ord
        @test gr.edges == gi.edges
        # reordering is a symmetric permutation of both matrices
        @test gr.C == Cm[perm, perm]
        @test gr.L_inv == Linv[perm, perm]
        @test node_index(gr, :B2) == 1
        # the original graph is untouched and the same order is a no-op
        @test [nd.label for nd in gi.nodes] == [:A1, :A2, :B1, :B2]
        @test to_matrix(reorder_nodes(gi, [:A1, :A2, :B1, :B2])) == (Cm, Linv)
        @test typeof(gr) == typeof(gi)

        @test_throws ArgumentError reorder_nodes(gi, [:A1, :A2, :B1])
        @test_throws ArgumentError reorder_nodes(gi, [:A1, :A2, :B1, :B1])
        @test_throws ArgumentError reorder_nodes(gi, [:A1, :A2, :B1, :Z])
        @test_throws ArgumentError reorder_nodes(gi, [:A1, :A2, :B1, :B2, :Z])

        # the mode transform follows the node order, so physics is unchanged by reordering
        pairs = [:A => (:A1, :A2), :B => (:B1, :B2)]
        mb0 = sum_diff_transform(gi, pairs)
        K0 = transform_cap_matrix(Cm, mb0)
        mb = sum_diff_transform(gr, pairs)
        @test mb.T != mb0.T
        Kp = transform_cap_matrix(gr.C, mb)
        pp, p0 = Cinv_params(inv(Kp), mb, :A, :B), Cinv_params(inv(K0), mb0, :A, :B)
        @test all(pp[k] ≈ p0[k] for k in keys(p0))
    end

    @testset "relabel" begin
        gi = CGraph([CNode(:A1, 20.0; Lg=3.0), CNode(:A2, 80.0), CNode(:B1, 20.0), CNode(:B2, 80.0)],
                    [CEdge(:A1, :A2; C=10.0, L=5.0), CEdge(:A2, :B2; C=5.0), CEdge(:B1, :B2; L=2.0)])
        new = [:q1, :q2, :q3, :q4]
        gn = relabel(gi, new)
        @test [nd.label for nd in gn.nodes] == new
        # same numbers in the same positions, so identical matrices
        @test gn.C == gi.C && gn.L_inv == gi.L_inv
        @test gn[:q1].Cg == 20.0 && gn[:q1].Lg == 3.0 && gn[:q4].Lg === nothing
        @test Set(XoniumUnitCell._edge_key(e) for e in gn.edges) ==
              Set([(:q1, :q2), (:q2, :q4), (:q3, :q4)])
        e24 = only(filter(e -> XoniumUnitCell._edge_key(e) == (:q2, :q4), gn.edges))
        @test e24.C == 5.0 && e24.L === nothing
        e12 = only(filter(e -> XoniumUnitCell._edge_key(e) == (:q1, :q2), gn.edges))
        @test e12.C == 10.0 && e12.L == 5.0
        @test typeof(gn) == typeof(gi)
        # the original graph is untouched
        @test [nd.label for nd in gi.nodes] == [:A1, :A2, :B1, :B2]
        # labels may be permuted among the existing ones, which is not a reordering
        gs = relabel(gi, [:A2, :A1, :B1, :B2])
        @test gs.C == gi.C
        @test gs[:A2].Cg == 20.0
        # vectors of the old graph no longer match
        @test_throws ArgumentError to_vector(:A1, gi)' * gn
        @test (to_vector(:q2, gn)' * gn * to_vector(:q4, gn)).C == gi.C[2, 4]

        @test_throws DimensionMismatch relabel(gi, [:a, :b, :c])
        @test_throws ArgumentError relabel(gi, [:a, :b, :c, :c])
    end

    @testset "to_vector" begin
        gv = CGraph(nodes, [CEdge(:A1, :A2; C=10.0)])
        v = to_vector(:A1, gv)
        @test v isa GraphVector
        @test v.data == [1, 0, 0, 0] && v.basis == [:A1, :A2, :B1, :B2]
        @test to_vector(:B1, gv).data == [0, 0, 1, 0]
        @test to_vector(gv[:B2], gv).data == [0, 0, 0, 1]
        @test to_vector(gv[:B2], gv) == to_vector(:B2, gv)
        @test_throws MethodError to_vector(gv[:B2])          # no method without a graph
        # follows the node order, so it changes under reorder_nodes
        gr = reorder_nodes(gv, [:B2, :A1, :B1, :A2])
        @test to_vector(:A1, gr).data == [0, 1, 0, 0]
        @test to_vector(:B2, gr).data == [1, 0, 0, 0]
        @test to_vector(:A1, gr).basis == [:B2, :A1, :B1, :A2]
        # the k-th basis vector picks out row/column k of the matrices
        @test to_vector(:A2, gr).data' * gr.C * to_vector(:A2, gr).data == gr.C[4, 4]
        @test hcat((to_vector(nd, gr).data for nd in gr.nodes)...) == LinearAlgebra.I(4)
        @test_throws KeyError to_vector(:Z, gv)
        @test_throws KeyError to_vector(CNode(:Z, 1.0), gv)
    end

    @testset "node_vectors" begin
        gv = CGraph(nodes, [CEdge(:A1, :A2; C=10.0, L=5.0), CEdge(:A2, :B2; C=5.0)])
        nv = node_vectors(gv)
        @test nv isa NamedTuple && keys(nv) == (:A1, :A2, :B1, :B2)
        @test all(v isa GraphVector for v in nv)
        @test nv.A1 == to_vector(:A1, gv) && nv.B2.data == [0, 0, 0, 1]
        @test hcat((v.data for v in nv)...) == LinearAlgebra.I(4)
        # follows the node order, and unpacks into variables
        gr = reorder_nodes(gv, [:B2, :A1, :B1, :A2])
        @test keys(node_vectors(gr)) == (:B2, :A1, :B1, :A2)
        (; A2, B2) = node_vectors(gr)
        @test (A2' * gr * B2).C == gr.C[4, 1]
        @test_throws ArgumentError nv.A1' * gr * nv.B2
    end

    @testset "GraphVector algebra" begin
        gm = CGraph(nodes, [CEdge(:A1, :A2; C=10.0, L=5.0), CEdge(:A2, :B2; C=5.0),
                            CEdge(:B1, :B2; L=2.0)])
        a, b = to_vector(:A2, gm), to_vector(:B2, gm)
        # the adjoint adjoints the data
        @test a'.data == a.data' && a'.basis == a.basis
        @test a''.data == a.data && a'' == a
        @test a'.data isa Adjoint
        @test_throws DimensionMismatch GraphVector([1, 2], [:a, :b, :c])
        # scalar products
        @test a' * a == 1 && a' * b == 0
        @test a * b' == a.data * b.data'                      # outer product is a Matrix
        w = [0.0, 1.0, 1.0, -1.0]
        @test a' * w == w[2] && w' * a == w[2]
        @test_throws MethodError a * a                       # column times column
        @test_throws DimensionMismatch a' * w[1:3]
        # graph products carry the basis along
        r = a' * gm
        @test r.C isa GraphVector && r.L isa GraphVector
        @test r.C.data == gm.C[2, :]' && r.C.basis == a.basis
        c = gm * a
        @test c.C.data == gm.C[:, 2] && c.L.data == gm.L_inv[:, 2]
        # bilinear forms pick out matrix elements, with either association
        elem = a' * gm * b
        @test elem == (C=gm.C[2, 4], L=gm.L_inv[2, 4])
        @test (a' * gm) * b == elem == a' * (gm * b)
        @test (b' * gm * a).C == gm.C[4, 2]
        @test (a' * gm * a).C == gm.C[2, 2]
        # and generic vectors through the plain-vector products
        @test (a' * gm).C * w == gm.C[2, :]' * w
        # mismatched bases are rejected, e.g. after reordering
        gr = reorder_nodes(gm, [:B2, :A1, :B1, :A2])
        @test_throws ArgumentError a' * gr
        @test_throws ArgumentError gr * a
        ar, br = to_vector(:A2, gr), to_vector(:B2, gr)
        @test (ar' * gr * br).C == gr.C[4, 1]
        @test (ar' * gr * br).C == (a' * gm * b).C      # same physics in either basis
        @test_throws ArgumentError a' * br
    end

    @testset "ground inductance" begin
        @test CNode(:a, 1.0).Lg === nothing
        @test CNode(:a, 1.0; Lg=2).Lg === 2
        @test_throws ArgumentError CGraph([CNode(:a, 1.0; Lg=0.0)])
        gn = CGraph([CNode(:a, 1.0; Lg=4.0), CNode(:b, 2.0)], [CEdge(:a, :b; C=1.0, L=2.0)])
        @test gn isa CGraph{Float64,Float64}
        @test gn[:a].Lg === 4.0 && gn[:b].Lg === nothing
        Cn, Ln = to_matrix(gn)
        @test Ln ≈ [1/4+1/2 -1/2; -1/2 1/2]
        @test Cn ≈ [2.0 -1.0; -1.0 3.0]
        add_node!(gn, CNode(:b, 2.0; Lg=5))
        @test gn[:b].Lg === 5.0
        @test_throws ArgumentError add_node!(gn, CNode(:b, 2.0; Lg=0))
    end

    @testset "transform" begin
        gi = CGraph(nodes, [CEdge(:A1, :A2; C=10.0, L=5.0), CEdge(:A2, :B2; C=5.0),
                            CEdge(:B1, :B2; L=2.0), CEdge(:A2, :B1; L=4.0)])
        Cm, Linv = to_matrix(gi)

        # a generic invertible T: the matrices of the new graph are T'CT and T'L⁻¹T
        T = [1.0 0.5 0.0 0.2; 0.0 1.0 0.3 0.0; 0.4 0.0 1.0 0.1; 0.0 0.2 0.0 1.0]
        gt = transform(gi, T)
        # without `relabel` the nodes keep their labels
        @test [nd.label for nd in gt.nodes] == [:A1, :A2, :B1, :B2]
        Ct, Lt = to_matrix(gt)
        @test Ct ≈ T' * Cm * T
        @test Lt ≈ T' * Linv * T

        # custom names, and a reordered graph (the rows of T follow the node order)
        names = [:w, :x, :y, :z]
        gro = reorder_nodes(gi, [:B2, :A1, :B1, :A2])
        gp = transform(gro, T, names)
        Cp, Lp = to_matrix(gro)
        Ctp, Ltp = to_matrix(gp)
        @test [nd.label for nd in gp.nodes] == names
        @test Ctp ≈ T' * Cp * T
        @test Ltp ≈ T' * Lp * T
        # with the default the reordered labels are kept, in the reordered order
        @test [nd.label for nd in transform(gro, T).nodes] == [:B2, :A1, :B1, :A2]

        # the identity reproduces the graph
        gid = transform(gi, Matrix(1.0LinearAlgebra.I, 4, 4))
        @test cap_matrix(gid) ≈ Cm && inv_inductance_matrix(gid) ≈ Linv
        by_pair(g) = Dict(XoniumUnitCell._edge_key(e) => (e.C, e.L) for e in g.edges)
        orig, back = by_pair(gi), by_pair(gid)
        @test keys(orig) == keys(back)
        @test all(isnothing.(orig[k]) == isnothing.(back[k]) for k in keys(orig))

        # a restriction to a subspace is allowed, but needs names
        gr = transform(gi, T[:, 1:2], [:p, :q])
        @test length(gr.nodes) == 2
        @test cap_matrix(gr) ≈ T[:, 1:2]' * Cm * T[:, 1:2]
        @test_throws ArgumentError transform(gi, T[:, 1:2])

        @test_throws DimensionMismatch transform(gi, T[1:3, :])
        @test_throws DimensionMismatch transform(gi, T, [:a, :b])
        @test_throws ArgumentError transform(gi, T, [:a, :a, :b, :c])

        # sum/difference modes of a fluxonium: all of the inductance sits on the difference mode
        gf = CGraph([CNode(:a, 20.0), CNode(:b, 80.0)], [CEdge(:a, :b; C=10.0, L=5.0)])
        mb = sum_diff_transform(gf, [:A => (:a, :b)])
        V = inv(mb.T)                          # φ = Vψ, the inverse of ModeBasis.T
        gs = transform(gf, V, [:sumA, :diffA])
        @test cap_matrix(gs) ≈ transform_cap_matrix(cap_matrix(gf), mb)
        @test gs[:sumA].Lg === nothing         # sum mode has no potential energy
        @test gs[:diffA].Lg ≈ 5.0              # difference mode: L to ground
        @test only(gs.edges).L === nothing     # and no inductive coupling between them
    end

    @testset "matrix accessors" begin
        gi = CGraph(nodes, [CEdge(:A1, :A2; C=10.0, L=5.0), CEdge(:A2, :B2; C=5.0)])
        @test gi.C == cap_matrix(gi)
        @test gi.L_inv == inv_inductance_matrix(gi)
        @test (gi.C, gi.L_inv) == to_matrix(gi)
        @test gi.nodes === getfield(gi, :nodes) && gi.edges === getfield(gi, :edges)
        @test :C in propertynames(gi) && :L_inv in propertynames(gi)
        # reflects later changes to the graph
        add_edge!(gi, CEdge(:B1, :B2; C=2.0, L=1.0))
        @test gi.C[3, 4] == -2.0 && gi.L_inv[3, 4] == -1.0
        @test_throws ErrorException gi.C = zeros(4, 4)
    end

    @testset "graph-vector products" begin
        gm = CGraph(nodes, [CEdge(:A1, :A2; C=10.0, L=5.0), CEdge(:A2, :B2; C=5.0)])
        v = [1.0, -2.0, 0.5, 3.0]
        w = [0.0, 1.0, 1.0, -1.0]
        @test gm * v == (C=gm.C * v, L=gm.L_inv * v)
        @test v' * gm == (C=v' * gm.C, L=v' * gm.L_inv)
        @test keys(gm * v) == (:C, :L)
        # a basis vector picks out a row (or column) of each matrix
        e = to_vector(:A2, gm).data
        @test (gm * e).C == gm.C[:, 2] && (gm * e).L == gm.L_inv[:, 2]
        @test (e' * gm).C == gm.C[2, :]'
        # the bilinear form v' C w
        @test (v' * gm).C * w ≈ v' * gm.C * w
        @test_throws DimensionMismatch gm * v[1:3]
        @test_throws DimensionMismatch v[1:3]' * gm
    end

    C = cap_matrix(g)
    @test C ≈ C'
    @test C[1, 2] == -10.0
    @test C[2, 2] == 80.0 + 10.0 + 5.0

    pairs = [:A => (:A1, :A2), :B => (:B1, :B2)]
    @test_throws ArgumentError sum_diff_transform(g, [:A => (:A1, :A2)])
    @test_throws ArgumentError sum_diff_transform(g, [:A => (:A1, :A2), :B => (:A1, :B2)])
    mb = sum_diff_transform(g, pairs)
    @test mb.T == [1 1 0 0; 1 -1 0 0; 0 0 1 1; 0 0 1 -1]
    @test mode_index(mb, :B) == 4
    @test mode_index(mb, :B; kind=:sum) == 3
    K = transform_cap_matrix(C, mb)
    @test K ≈ K'

    # reference values from the original Python notebook (xonium-taichi-couplers.py)
    p = Cinv_params(inv(K), mb, :A, :B)
    @test p.EC_a ≈ 0.7398708283158923
    @test p.EC_b ≈ 0.7398708283158923
    @test p.J_ab ≈ 41.10393490643847

    # truncating the sum modes first underestimates E_C and overestimates J
    pt = Cinv_params(inv(truncate_to_diff_modes(K, mb)), mb, :A, :B; truncated=true)
    @test pt.EC_a < p.EC_a
    @test pt.J_ab > p.J_ab

    # node order does not matter: a permuted graph gives the same physics
    gperm = CGraph(nodes[[3, 1, 4, 2]], edges)
    mbp = sum_diff_transform(gperm, pairs)
    pp = Cinv_params(inv(transform_cap_matrix(cap_matrix(gperm), mbp)), mbp, :A, :B)
    @test pp.J_ab ≈ p.J_ab
end
