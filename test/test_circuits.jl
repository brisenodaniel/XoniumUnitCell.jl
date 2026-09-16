@testset "Single Qubit Circuits" begin

    @testset "Aashish Clerk 2023 PRA Fluxonium" begin
        # define parameters for Method I fluxonia A, B
        pA = (EJ=5.5, EC=1.5, EL=1.0) |> tounit_ħ
        pB = (EJ=5.7, EC=1.2, EL=1.0) |> tounit_ħ
        fA, fB = (fluxonium(;p...) for p in (pA, pB))
        ΔA, ΔB = map((fA, fB)) do fl
            g, e = energy.(Ref(fl), [:g, :e])
            e-g
        end
        @test ΔA/2π ≈ 0.606 atol=0.0005
        @test ΔB/2π ≈ 0.354 atol=0.0005

        # define parameters for method 2
        pA = (EJ=4.5, EC=1.8, EL=1.5) |> tounit_ħ
        pB = (EJ=3.5, EC=1.1, EL=1.0) |> tounit_ħ
        fA, fB = (fluxonium(;p...) for p in (pA, pB))
        ΔA, ΔB = map((fA, fB)) do fl
            g, e = energy.(Ref(fl), [:g, :e])
            e-g
        end
        @test ΔA/2π ≈ 1.79 atol=0.005
        @test ΔB/2π ≈ 0.86 atol=0.005
    end

    # SCqubits agreement
    @testset "SCQubits agreement" begin
        # transmon
        ptrans = (; EJ= 30.02, EC=1.2, ng=0.3, ncut=31) |>
            p -> tounit_ħ(p, (:EC, :EJ, :EL))
        trans = transmon(;ptrans...)
        Δtrans = energy(trans, :e) - energy(trans, :g) |> abs

        #tunable transmon
        ptune = (; EC=0.32, EJ1=3.4, EJ2=13, ϕext=0.0) |>
            p -> tounit_ħ(p, (:EC, :EJ1, :EJ2))
        split_trans = split_transmon(;ptune...)
        Δsplit_trans = energy(split_trans, :e) - energy(split_trans, :g) |> abs

        @test  Δtrans/2π ≈ 15.6686330484788 atol=1e-9
        @test Δsplit_trans/2π ≈ 6.14118624 atol=1e-9
    end
end

@testset "Multi Qubit Circuits" begin
    @testset "MIT FTF 2023 Device A" begin
        pf1 = (EC=1.41, EL=0.80, EJ=6.27) |> tounit_ħ
        pf2 = (EC=1.30, EL=0.59, EJ=5.71) |> tounit_ħ
        ptrans = (EC=0.32, EJ1=3.4, EJ2=13, ϕext=0.0) |>
            p -> tounit_ħ(p, (:EC, :EJ1, :EJ2))
        f1, f2 = fluxonium(;pf1..., name=:f1), fluxonium(;pf2..., name=:f2)
        cplr = split_transmon(;ptrans..., name=:cplr )
        # diagonalize subsystems, set to same energy reference
        e_ref = energy.([f1, f2, cplr], :g) .|> real |> minimum
        f1, f2, cplr = [to_energy_basis(sys; etol=1e-9, set_energy_zero=true) |>
            x -> truncate_sys(x, 6)
            for sys in (f1, f2, cplr)]
        ftf = f1 ⊗ cplr ⊗ f2
        # define couplings
        n1nc = ftf[:f1, :n, :cplr, :n]
        n2nc = ftf[:cplr, :n, :f2, :n]
        n1n2 = ftf[:f1, :n, :f2, :n]
        Jc = (570, 560, 125) .* (1e-3 * 2π)
        H_int = Jc[1]*n1nc + Jc[2]*n2nc + Jc[3]*n1n2
        ftf = ftf + H_int
        # get ground state energy
        g_en = energy(ftf, :ggg)
        # first fluxonium gap
        e1 = energy(ftf, :egg) - g_en |> real
        # second fluxonium gap
        e2 = energy(ftf, :gge) - g_en |> real
        # check with scqubit values
        @test e1/2π ≈ 0.29301234 atol=1e-6
        @test e2/2π ≈ 0.21770244 atol=1e-6
    end

end
