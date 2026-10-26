using Symbolics
@testset "unit conversions" begin
    # independent SI values
    e, h = 1.602176634e-19, 6.62607015e-34
    ħ = h / 2π

    @testset "capacitance" begin
        # E_C = e²/2C: 1 fF is 19.37 GHz
        @test Cinv2E_C(1.0) ≈ 1e-9 * e^2 / (2 * 1e-15) / h
        @test Cinv2E_C(1.0) ≈ 19.37 rtol = 1e-3
        @test Cinv2E_C(1 / 20.0) ≈ Cinv2E_C(1.0) / 20                  # linear in 1/C
        # J = 4e² Cinv, so J (MHz) is 8000 times E_C (GHz)
        @test Cinv2J(0.3) ≈ 8000 * Cinv2E_C(0.3)
        @test Cinv2J(1.0) ≈ 1e-6 * 4e15 * e^2 / h
        @test Cinv2g(1.0) ≈ 1e-6 * e * sqrt(π * ħ / 100) * 1e15 / h
        @test Cinv2E_C(0.0) == 0 && Cinv2J(-0.5) == -Cinv2J(0.5)       # sign is kept
        # reference values of the original notebook: ground 20/80/20/80 fF, couplers 10/5/10 fF,
        # in the sum/difference coordinates ψ = Mφ, so K = M⁻ᵀ C M⁻¹. E_C = 0.73987 GHz, J = 41.104 MHz
        M = [1 1 0 0; 1 -1 0 0; 0 0 1 1; 0 0 1 -1]
        C = [30.0 -10.0 0.0 0.0; -10.0 95.0 0.0 -5.0; 0.0 0.0 30.0 -10.0; 0.0 -5.0 -10.0 95.0]
        Mi = inv(Rational{Int}.(M))
        Kinv = inv(Mi' * C * Mi)
        @test Cinv2E_C(Kinv[2, 2]) ≈ 0.7398708283158923
        @test Cinv2E_C(Kinv[4, 4]) ≈ 0.7398708283158923
        @test Cinv2J(Kinv[2, 4]) ≈ 41.10393490643847
        # resonator
        @test CPW_freq2cap(7.0) ≈ 1e15 / (8 * 7e9 * 50)
        @test CPW_freq2cap(7.0; Z0=75) ≈ CPW_freq2cap(7.0) * 50 / 75
        @test CPW_freq2cap(2.0) ≈ 2 * CPW_freq2cap(4.0)
    end

    @testset "inductance" begin
        # E_L = (Φ0/2π)²/L: 1 nH is 163.46 GHz
        @test Linv2E_L(1.0) ≈ 1e-9 * (ħ / (2e))^2 / 1e-9 / h
        @test Linv2E_L(1.0) ≈ 163.46 rtol = 1e-4
        @test Linv2E_L(0.2) ≈ Linv2E_L(1.0) / 5
        @test L2E_L(5.0) ≈ Linv2E_L(1 / 5.0)
        @test L2E_L(5.0) ≈ 163.46 / 5 rtol = 1e-4
        @test L2E_L(2.0) ≈ L2E_L(1.0) / 2
        # the L_inv element of a mode gives its E_L directly: a 5 nH fluxonium inductor
        g = CGraph([CNode(:a, 20.0), CNode(:b, 80.0)], [CEdge(:a, :b; C=10.0, L=5.0)])
        Vi = [1//2 1//2; 1//2 -1//2]                                     # φ = Vψ, ψ_Δ = φ_a - φ_b
        gd = Vi' * g * Vi
        @test Linv2E_L(gd.L_inv[2, 2]) ≈ L2E_L(5.0)
        @test Linv2E_L(gd.L_inv[1, 1]) ≈ 0 atol = 1e-12                  # the sum mode has none
    end

    @testset "Josephson energy" begin
        # E_J = ħ Ic / 2e: 1 nA is 0.4967 GHz
        @test Ic2E_J(1.0) ≈ ħ * 1e-9 / (2e) / h * 1e-9
        @test Ic2E_J(1.0) ≈ 0.4967 rtol = 1e-3
        @test Ic2E_J(30.0) ≈ 30 * Ic2E_J(1.0)
        @test E_J2Ic(Ic2E_J(12.5)) ≈ 12.5
        @test Ic2E_J(E_J2Ic(4.0)) ≈ 4.0
        # E_J and E_L share the flux quantum: E_J · E_L-normalized product = Ic · (Φ0/2π)
        φ0 = ħ / (2e)
        @test Ic2E_J(1.0) * L2E_L(1.0) ≈ (φ0 * 1e-9 / h * 1e-9) * (φ0^2 / 1e-9 / h * 1e-9) rtol = 1e-12
    end

    @testset "types" begin
        @test Cinv2E_C(1) isa Float64 && Linv2E_L(2) isa Float64 && Ic2E_J(3) isa Float64
        @test Cinv2E_C(1.0f0) isa Real
        @test_throws MethodError Cinv2E_C("1")
        # symbolic values pass through
        @variables Cinv Linv Ic
        @test Cinv2E_C(Cinv) isa Num && Linv2E_L(Linv) isa Num && Ic2E_J(Ic) isa Num
        @test isequal(Symbolics.simplify(Cinv2E_C(2Cinv) - 2 * Cinv2E_C(Cinv); expand=true), 0)
    end
end
