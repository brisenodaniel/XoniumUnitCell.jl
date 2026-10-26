# Tests for src/circuit_graphs.jl, loaded into its own module so they are independent of the
# copy of the file that the package includes.
module CircuitGraphsTests
using Test, LinearAlgebra, Symbolics
include(joinpath(@__DIR__, "..", "src", "circuit_graphs.jl"))

labels_of(g) = [nd.label for nd in g.nodes]
edge_dict(g) = Dict(CircuitGraphsTests._edge_key(e) => e for e in g.edges)

@testset "circuit_graphs" begin
    nodes = [CNode(:A1, 20.0), CNode(:A2, 80.0), CNode(:B1, 20.0), CNode(:B2, 80.0)]
    edges = [CEdge(:A1, :A2; C=10.0, L=5.0), CEdge(:A2, :B2; C=5.0),
             CEdge(:B1, :B2; L=2.0), CEdge(:A2, :B1; L=4.0)]
    g = CGraph(nodes, edges)
    V4 = [1.0 0.5 0.0 0.2; 0.0 1.0 0.3 0.0; 0.4 0.0 1.0 0.1; 0.0 0.2 0.0 1.0]

    @testset "nodes and edges" begin
        @test CNode(:a, 1.0).Lg === nothing
        @test CNode(:a, 1.0; Lg=2).Lg === 2
        @test CEdge(:a, :b; C=1.0).L === nothing
        @test CEdge(:a, :b; L=5.0).C === nothing
        @test CEdge(:a, :b; L=5.0) isa CEdge{Float64,Float64}
        @test CEdge(:a, :b; C=10.0, L=5).L === 5
    end

    @testset "building from nodes and edges" begin
        @test g isa CGraph{Float64}
        @test labels_of(g) == [:A1, :A2, :B1, :B2]
        @test g.tol == 1e-9
        # the nodes and edges are stored as given
        @test g.nodes[2].Cg == 80.0 && all(nd.Lg === nothing for nd in g.nodes)
        @test length(g.edges) == 4
        @test edge_dict(g)[(:A1, :A2)].C == 10.0 && edge_dict(g)[(:A1, :A2)].L == 5.0
        @test edge_dict(g)[(:B1, :B2)].C === nothing && edge_dict(g)[(:B1, :B2)].L == 2.0
        # and give the matrices
        @test g.C[1, 2] == g.C[2, 1] == -10.0
        @test g.C[1, 1] == 20.0 + 10.0
        @test g.C[2, 2] == 80.0 + 10.0 + 5.0
        @test g.C[1, 3] == 0.0
        @test g.L_inv[1, 2] ≈ -1 / 5
        @test g.L_inv[2, 2] ≈ 1 / 5 + 1 / 4          # A2 touches the 5 nH and 4 nH inductors
        @test g.L_inv[3, 4] ≈ -1 / 2
        @test g.L_inv[3, 3] ≈ 1 / 2 + 1 / 4
        @test g.L_inv[1, 4] == 0.0
        @test all(abs.(sum(g.L_inv; dims=2)) .< 1e-12)   # no inductor to ground: singular
        @test g.C == g.C' && g.L_inv == g.L_inv'

        # no edges gives diagonal matrices, no inductors gives a zero L_inv
        g0 = CGraph(nodes)
        @test g0.C == Diagonal([20.0, 80.0, 20.0, 80.0]) && iszero(g0.L_inv) && isempty(g0.edges)

        # parallel edges are merged: capacitances add, inductances combine in parallel
        gp = CGraph(nodes, [CEdge(:A1, :A2; C=4.0, L=6.0), CEdge(:A2, :A1; C=6.0, L=3.0)])
        @test gp.C[1, 2] == -10.0 && gp.L_inv[1, 2] ≈ -(1 / 6 + 1 / 3)
        @test length(gp.edges) == 1
        @test only(gp.edges).C == 10.0 && only(gp.edges).L ≈ 2.0
        gq = CGraph(nodes, [CEdge(:A1, :A2; C=4.0, L=6.0), CEdge(:A2, :A1; C=6.0)])
        @test only(gq.edges).C == 10.0 && only(gq.edges).L == 6.0
        # inductor to ground on a node
        gn = CGraph([CNode(:a, 1.0; Lg=4.0), CNode(:b, 2.0)], [CEdge(:a, :b; C=1.0, L=2.0)])
        @test gn.L_inv ≈ [1/4+1/2 -1/2; -1/2 1/2] && gn.C ≈ [2.0 -1.0; -1.0 3.0]
        @test gn.nodes[1].Lg === 4.0 && gn.nodes[2].Lg === nothing
        # negative inductance is allowed, and equal and opposite inductors cancel
        @test CGraph(nodes, [CEdge(:A1, :A2; L=-1.0)]).L_inv[1, 2] == 1.0
        gc = CGraph(nodes, [CEdge(:A1, :A2; C=1.0, L=3.0), CEdge(:A1, :A2; L=-3.0)])
        @test iszero(gc.L_inv) && only(gc.edges).L === nothing && only(gc.edges).C == 1.0
        # numeric types are promoted, with inverse inductances taken into account
        @test CGraph([CNode(:a, 1), CNode(:b, 2)], [CEdge(:a, :b; C=1, L=2)]) isa CGraph{Float64}
        @test CGraph([CNode(:a, 1), CNode(:b, 2)], [CEdge(:a, :b; C=1)]) isa CGraph{Int}

        # tolerance
        @test CGraph(nodes, edges; tol=1e-6).tol == 1e-6
        @test_throws ArgumentError CGraph(nodes, edges; tol=-1.0)

        @test_throws ArgumentError CGraph([CNode(:a, 1.0), CNode(:a, 2.0)])
        @test_throws ArgumentError CGraph(CNode[])
        @test_throws ArgumentError CGraph(nodes, [CEdge(:A1, :A1; C=1.0)])
        @test_throws ArgumentError CGraph(nodes, [CEdge(:A1, :Z; C=1.0)])
        @test_throws ArgumentError CGraph(nodes, [CEdge(:A1, :A2)])
        @test_throws ArgumentError CGraph(nodes, [CEdge(:A1, :A2; L=0.0)])
        @test_throws ArgumentError CGraph([CNode(:a, 1.0; Lg=0.0)])
    end

    @testset "building from matrices" begin
        # a graph read back from its own matrices has the same nodes and edges
        gb = CGraph(labels_of(g), g.C, g.L_inv)
        @test gb isa CGraph{Float64} && labels_of(gb) == labels_of(g)
        @test [nd.Cg for nd in gb.nodes] ≈ [nd.Cg for nd in g.nodes]
        @test all(nd.Lg === nothing for nd in gb.nodes)
        eg, eb = edge_dict(g), edge_dict(gb)
        @test keys(eg) == keys(eb)
        for k in keys(eg)
            @test isnothing(eb[k].C) == isnothing(eg[k].C) && isnothing(eb[k].L) == isnothing(eg[k].L)
            isnothing(eg[k].C) || @test eb[k].C ≈ eg[k].C
            isnothing(eg[k].L) || @test eb[k].L ≈ eg[k].L
        end
        # ground capacitance and inductance are the row sums
        C = [3.0 -1.0; -1.0 2.0]
        L = [0.5 -0.25; -0.25 0.25]
        gm = CGraph([:p, :q], C, L)
        @test gm.nodes[1].Cg == 2.0 && gm.nodes[2].Cg == 1.0
        @test gm.nodes[1].Lg ≈ 4.0 && gm.nodes[2].Lg === nothing
        @test only(gm.edges).nodes == (:p, :q) && only(gm.edges).C == 1.0 && only(gm.edges).L ≈ 4.0
        C[1, 1] = 99.0                              # the graph owns its matrices
        @test gm.C[1, 1] == 3.0
        @test CGraph([:p, :q], [1 -1; -1 2], [1.5 0; 0 1]) isa CGraph{Float64}
        # entries within tol are zero
        gz = CGraph([:p, :q], [1.0 1e-12; 1e-12 1.0], zeros(2, 2))
        @test isempty(gz.edges) && all(nd.Lg === nothing for nd in gz.nodes)
        @test_throws DimensionMismatch CGraph([:p, :q], ones(2, 3), ones(2, 2))
        @test_throws DimensionMismatch CGraph([:p, :q], ones(2, 2), ones(3, 3))
        @test_throws ArgumentError CGraph([:p, :p], ones(2, 2), ones(2, 2))
        @test gm == CGraph([:p, :q], [3.0 -1.0; -1.0 2.0], L)
        @test gm != CGraph([:q, :p], [3.0 -1.0; -1.0 2.0], L)
    end

    @testset "symmetry" begin
        sym = [1.0 2.0; 2.0 4.0]
        # exactly asymmetric matrices are rejected, whichever of the two is asymmetric
        @test_throws ArgumentError CGraph([:p, :q], [1.0 2.0; 3.0 4.0], sym)
        @test_throws ArgumentError CGraph([:p, :q], sym, [0.0 1.0; -1.0 5.0])
        err = try CGraph([:p, :q], [1.0 2.0; 3.0 4.0], sym) catch e e end
        @test occursin("C", err.msg) && occursin("1.0", err.msg) && occursin("1.0e-9", err.msg)
        # small asymmetry within tol is accepted and symmetrized exactly
        near = [1.0 2.0; 2.0 + 1e-10 4.0]
        gs = CGraph([:p, :q], near, sym)
        @test gs.C == gs.C' && gs.C[1, 2] ≈ 2.0 + 5e-11
        # above the default tolerance it is rejected, unless tol is raised
        far = [1.0 2.0; 2.0 + 1e-8 4.0]
        @test_throws ArgumentError CGraph([:p, :q], far, sym)
        @test CGraph([:p, :q], far, sym; tol=1e-6).C == CGraph([:p, :q], far, sym; tol=1e-6).C'
        @test CGraph([:p, :q], far, sym; tol=1e-6).tol == 1e-6
        # an exactly symmetric matrix is stored unchanged
        @test CGraph([:p, :q], sym, sym).C == sym
    end

    @testset "M * graph * Q" begin
        # the congruence V' * graph * V
        V = [1.0 0.5 0.0 0.2; 0.0 1.0 0.3 0.0; 0.4 0.0 1.0 0.1; 0.0 0.2 0.0 1.0]
        gt = V' * g * V
        @test gt isa CGraph{Float64} && labels_of(gt) == labels_of(g)
        @test gt.C ≈ V' * g.C * V && gt.L_inv ≈ V' * g.L_inv * V
        @test gt.C == gt.C' && gt.L_inv == gt.L_inv'
        # the inverse capacitance transforms with the inverse matrices
        @test inv(gt.C) ≈ inv(V) * inv(g.C) * inv(V)'
        # the original graph is untouched
        @test g.C == CGraph(nodes, edges).C && labels_of(g) == [:A1, :A2, :B1, :B2]
        # different matrices on the two sides are fine if the result is symmetric
        S = Diagonal([2.0, 1.0, 1.0, 3.0])
        @test (S * g * S).C ≈ S * g.C * S
        @test ((2 * Matrix(1.0LinearAlgebra.I, 4, 4)) * g * Matrix(1.0LinearAlgebra.I, 4, 4)).C ≈ 2 * g.C

        # the nodes and edges are updated to match the new matrices
        @test CGraph(gt.nodes, gt.edges).C ≈ gt.C
        @test CGraph(gt.nodes, gt.edges).L_inv ≈ gt.L_inv
        @test gt.nodes[1].Cg ≈ sum(gt.C[1, :]) && gt.edges != g.edges
        egt = edge_dict(gt)
        @test egt[(:A1, :A2)].C ≈ -gt.C[1, 2]
        @test egt[(:A1, :A2)].L ≈ -1 / gt.L_inv[1, 2]
        # and a transformed graph can be read back from its matrices again
        @test CGraph(labels_of(gt), gt.C, gt.L_inv) == gt

        # sum and difference modes of a fluxonium: the inductance sits on the difference mode
        gf = CGraph([CNode(:a, 20.0), CNode(:b, 80.0)], [CEdge(:a, :b; C=10.0, L=5.0)])
        Vf = inv(Rational{Int}.([1 1; 1 -1]))           # φ = Vψ for ψ = (φ_a + φ_b, φ_a - φ_b)
        gs = relabel(Vf' * gf * Vf, [:sum, :diff])
        @test gs.nodes[1].Lg === nothing                 # the sum mode has no potential energy
        @test gs.nodes[2].Lg ≈ 5.0                       # the difference mode: all of the inductor
        @test gs.L_inv ≈ [0.0 0.0; 0.0 1 / 5] atol = 1e-12
        @test all(isnothing(e.L) for e in gs.edges)      # and no inductive coupling
        @test gs.C ≈ Vf' * gf.C * Vf
        # the tolerance is inherited
        @test (V' * CGraph(nodes, edges; tol=1e-6) * V).tol == 1e-6

        # the result must be symmetric to within tol
        M = [1.0 1.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 2.0 0.0; 1.0 0.0 0.0 1.0]
        @test_throws ArgumentError M * g * Matrix(1.0LinearAlgebra.I, 4, 4)
        @test_throws ArgumentError Matrix(1.0LinearAlgebra.I, 4, 4) * g * M
        @test (M * g * M').C ≈ M * g.C * M'
        # numerical noise in an otherwise symmetric product is accepted
        noisy = V' * g * (V .+ 1e-14)
        @test noisy.C == noisy.C'
        # the result must still have one row per node
        @test_throws DimensionMismatch ones(2, 4) * g * ones(4, 2)
        @test_throws DimensionMismatch ones(3, 3) * g * ones(3, 3)
        # the two-sided products with a single matrix are not defined
        @test_throws MethodError ones(4, 4) * g
        @test_throws MethodError g * ones(4, 4)
    end

    @testset "basis_vectors" begin
        bv = basis_vectors(g)
        @test bv isa NamedTuple && keys(bv) == (:A1, :A2, :B1, :B2)
        @test bv.A1 == [1, 0, 0, 0] && bv.B1 == [0, 0, 1, 0] && bv.B2 == [0, 0, 0, 1]
        @test all(v isa Vector{Int} for v in bv)
        @test hcat(values(bv)...) == I(4)
        # fresh vectors on each call
        bv.A1[1] = 7
        @test basis_vectors(g).A1[1] == 1 && bv.A2 == [0, 1, 0, 0]
        # the order of the nodes decides the correspondence
        gr = CGraph(nodes[[4, 1, 3, 2]], [CEdge(:A1, :A2; C=10.0)])
        br = basis_vectors(gr)
        @test keys(br) == (:B2, :A1, :B1, :A2)
        @test br.B2 == [1, 0, 0, 0] && br.A1 == [0, 1, 0, 0] && br.A2 == [0, 0, 0, 1]
        (; A2, B2) = br
        @test A2' * B2 == 0 && A2' * A2 == 1
        # a basis vector picks out an entry or a row of the matrices
        (; A1, B2) = basis_vectors(g)
        @test A1' * g.C * B2 == g.C[1, 4]
        @test g.C * A1 == g.C[:, 1]
        @test basis_vectors(CGraph([CNode(:x, 1.0)])) == (x=[1],)
        @test size(br.A1) == (4,) && size(br.A1') == (1, 4)
    end

    @testset "relabel" begin
        new = [:q1, :q2, :q3, :q4]
        gn = relabel(g, new)
        @test gn isa CGraph{Float64} && labels_of(gn) == new
        @test gn.C == g.C && gn.L_inv == g.L_inv && gn.tol == g.tol
        @test keys(basis_vectors(gn)) == (:q1, :q2, :q3, :q4)
        @test basis_vectors(gn).q3 == basis_vectors(g).B1
        # the edges are relabeled too, and still match the matrices
        @test Set(keys(edge_dict(gn))) == Set([(:q1, :q2), (:q2, :q4), (:q3, :q4), (:q2, :q3)])
        @test edge_dict(gn)[(:q1, :q2)].C == 10.0 && edge_dict(gn)[(:q1, :q2)].L == 5.0
        @test CGraph(gn.nodes, gn.edges).C ≈ gn.C
        # a new graph that owns its matrices, and the original is untouched
        @test gn.C !== g.C
        gn.C[1, 1] = -1.0
        @test g.C[1, 1] == 30.0 && labels_of(g) == [:A1, :A2, :B1, :B2]
        # existing labels may be permuted among the nodes, which is not a reordering
        gs = relabel(g, [:A2, :A1, :B1, :B2])
        @test gs.C == g.C && labels_of(gs) == [:A2, :A1, :B1, :B2]
        # naming the variables after a change of basis
        V = [1.0 1.0 0.0 0.0; 1.0 -1.0 0.0 0.0; 0.0 0.0 1.0 1.0; 0.0 0.0 1.0 -1.0]
        gm = relabel(V' * g * V, [:sumA, :diffA, :sumB, :diffB])
        @test labels_of(gm) == [:sumA, :diffA, :sumB, :diffB] && gm.C ≈ V' * g.C * V
        @test all(e.nodes[1] in labels_of(gm) && e.nodes[2] in labels_of(gm) for e in gm.edges)
        @test_throws DimensionMismatch relabel(g, [:a, :b, :c])
        @test_throws DimensionMismatch relabel(g, [:a, :b, :c, :d, :e])
        @test_throws ArgumentError relabel(g, [:a, :b, :c, :c])
    end

    @testset "CVector and CAdjoint" begin
        v = [1.0, -2.0, 0.5, 3.0]
        w = [0.0, 1.0, 1.0, -1.0]
        cv = g * v
        @test cv isa CVector && cv.C == g.C * v && cv.L_inv == g.L_inv * v
        ca = v' * g
        @test ca isa CAdjoint && ca.C == v' * g.C && ca.L_inv == v' * g.L_inv
        @test ca.C isa Adjoint && ca.L_inv isa Adjoint
        # integer vectors
        @test (g * [1, 0, 0, 0]).C == g.C[:, 1]
        @test ([1, 0, 0, 0]' * g).C == g.C[1, :]'

        # vector/adjoint interactions: u' * x and x * w
        u = [2.0, 0.0, 1.0, -1.0]
        @test u' * cv == (C=u' * g.C * v, L=u' * g.L_inv * v)
        @test ca * w == (C=v' * g.C * w, L=v' * g.L_inv * w)
        @test keys(u' * cv) == (:C, :L) && u' * cv isa NamedTuple
        @test (u' * cv).C isa Float64
        # outer products
        @test cv * u' == (C=(g.C * v) * u', L=(g.L_inv * v) * u')
        @test w * ca == (C=w * (v' * g.C), L=w * (v' * g.L_inv))
        @test (cv * u').C isa Matrix
        # associativity: either way round gives the same bilinear form
        @test (u' * g) * v == u' * (g * v) == (C=u' * g.C * v, L=u' * g.L_inv * v)
        @test u' * g * v == (C=u' * g.C * v, L=u' * g.L_inv * v)
        # matrix elements through the basis vectors
        (; A1, B2) = basis_vectors(g)
        @test A1' * g * B2 == (C=g.C[1, 4], L=g.L_inv[1, 4])
        @test B2' * g * A1 == (C=g.C[4, 1], L=g.L_inv[4, 1])
        @test A1' * (g * B2) == A1' * g * B2

        # adjoints of the wrappers
        @test cv' isa CAdjoint && cv'.C == (g.C * v)'
        @test ca' isa CVector && ca'.L_inv == (v' * g.L_inv)'
        @test cv'' == cv && ca'' == ca
        @test cv' * w == (C=(g.C * v)' * w, L=(g.L_inv * v)' * w)
        @test cv == g * v && cv != g * w && ca == v' * g && ca != w' * g

        # not defined: wrapper times wrapper, and wrong shapes
        @test_throws MethodError ca * cv
        @test_throws MethodError cv * cv
        @test_throws MethodError ca * ca
        @test_throws DimensionMismatch g * v[1:3]
        @test_throws DimensionMismatch v[1:3]' * g
        @test_throws DimensionMismatch u[1:3]' * cv
        @test_throws DimensionMismatch ca * w[1:3]

        # with a pair of vectors and the three-way matrix product written next to each other
        gx = CGraph([:p, :q], [1.0 2.0; 2.0 4.0], [0.0 1.0; 1.0 5.0])
        e1, e2 = [1, 0], [0, 1]
        @test e1' * gx * e2 == (C=2.0, L=1.0) && e2' * gx * e1 == (C=2.0, L=1.0)
        @test (gx * e2).C == [2.0, 4.0] && (e2' * gx).C == [2.0 4.0]

        # symbolic entries
        @variables Cg1 Cg2 Ca Lj
        gy = CGraph([CNode(:a, Cg1), CNode(:b, Cg2)], [CEdge(:a, :b; C=Ca, L=Lj)])
        r = [1, 0]' * gy * [0, 1]
        @test r isa NamedTuple && isequal(r.C, -Ca) && isequal(r.L, -1 / Lj)
        @test (gy * [1, 0]) isa CVector
    end

    @testset "simplify hooks" begin
        @variables Cg1 Cg2 Ca Lj
        simp = x -> Symbolics.simplify(x; expand=true)
        build(; kw...) = CGraph([CNode(:a, Cg1), CNode(:b, Cg2)], [CEdge(:a, :b; C=Ca, L=Lj)]; kw...)
        V = [1//2 1//2; 1//2 -1//2]            # φ = Vψ for ψ = (φ_a + φ_b, φ_a - φ_b)

        gs = build(; simplify=simp)
        @test gs.simplify === simp && build().simplify === identity
        gt = V' * gs * V
        @test gt.simplify === simp                              # inherited by the product
        @test relabel(gt, [:s, :d]).simplify === simp           # and by relabel
        # with simplify, the nodes and edges read back from the transformed matrices are clean
        @test gt.nodes[1].Lg === nothing                         # the sum mode has no inductor
        @test isequal(gt.nodes[2].Lg, Lj)                        # the difference mode has all of it
        @test all(isnothing(e.L) for e in gt.edges)
        @test isequal(simp(gt.nodes[1].Cg - Cg1 / 2), 0)
        @test isequal(simp(gt.nodes[2].Cg - (Ca + Cg1 / 2)), 0)
        @test isequal(simp(only(gt.edges).C - (Cg2 - Cg1) / 4), 0)
        @test isequal(simp(gt.L_inv[2, 2] - 1 / Lj), 0) && isequal(gt.L_inv[1, 2], 0)
        # the stored matrices are simplified too, and the graph can be rebuilt from its parts
        @test isequal(gt.L_inv[1, 1], 0)
        @test isequal(simp(CGraph(gt.nodes, gt.edges; simplify=simp).C[2, 2] - gt.C[2, 2]), 0)

        # without it, mathematically zero entries are not recognised and show up as spurious
        # inductors and edges
        g0 = V' * build() * V
        @test g0.simplify === identity
        @test g0.nodes[1].Lg !== nothing

        # symbolic matrices are checked for symmetry when there is a simplify function
        sym = [Cg1 Ca; Ca Cg2]
        asym = [Cg1 Ca; Cg2 Cg2]
        zero2 = Num.(zeros(2, 2))
        @test CGraph([:a, :b], sym, zero2; simplify=simp) isa CGraph{Num}
        @test CGraph([:a, :b], [Cg1 Ca + Lj - Lj; Ca Cg2], zero2; simplify=simp) isa CGraph{Num}
        err = try CGraph([:a, :b], asym, zero2; simplify=simp) catch e e end
        @test err isa ArgumentError && occursin("not symmetric", err.msg) && occursin("C[1, 2]", err.msg)
        @test_throws ArgumentError CGraph([:a, :b], sym, [Lj Ca; Cg1 Lj]; simplify=simp)
        # ... and not without one, because unsimplified differences cannot be shown to vanish
        @test CGraph([:a, :b], asym, zero2) isa CGraph{Num}

        # numbers are never passed to simplify, and are not affected by it
        calls = Ref(0)
        counting = x -> (calls[] += 1; x)
        gn = CGraph(nodes, edges; simplify=counting)
        gnt = V4' * gn * V4
        @test calls[] == 0 && gnt.simplify === counting
        @test gnt.C ≈ V4' * g.C * V4
        @test CGraph(labels_of(gn), gn.C, gn.L_inv; simplify=counting).C == gn.C && calls[] == 0
        # a numeric graph is the same with or without it
        @test gn == g
    end

    @testset "symbolic" begin
        @variables Cg1 Cg2 Ca Lj
        gs = CGraph([CNode(:a, Cg1), CNode(:b, Cg2)], [CEdge(:a, :b; C=Ca, L=Lj)])
        @test gs isa CGraph{Num}
        @test isequal(gs.C[1, 2], -Ca) && isequal(gs.C[1, 1], Cg1 + Ca)
        @test isequal(gs.L_inv[1, 2], -1 / Lj)
        @test isequal(gs.nodes[1].Cg, Cg1) && isequal(only(gs.edges).C, Ca)
        # symbolic matrices cannot be checked for symmetry, and are stored as computed
        V = [1 1; 1 -1]
        gt = V' * gs * V
        @test gt isa CGraph{Num} && labels_of(gt) == [:a, :b]
        @test all(isequal.(gt.C, V' * gs.C * V))
        @test all(isequal.(gt.L_inv, V' * gs.L_inv * V))
        @test length(gt.nodes) == 2
        @test gs == CGraph([:a, :b], gs.C, gs.L_inv)
        @test gs != gt
        @test isequal(Symbolics.simplify(gt.L_inv[1, 1]; expand=true), 0)
    end
end
end # module
