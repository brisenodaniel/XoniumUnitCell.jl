fftshift_like(x) = circshift(x, length(x) ÷ 2)

@testset "Environmental PSD" begin
    Ad, Af = 2.0e-7, 3.0e-6
    kT = XoniumUnitCell.KT_15mK

    @testset "S_pink_noise" begin
        w = 2π * 0.5
        T = abs(2π / w)
        α = abs(1 / tanh(w / (2kT)) + 1) / 2
        @test S_pink_noise(w, Ad, Af) ≈ Af^2 * T + α * Ad / T^2
        # 1/f term is even in w, the thermal weight α is not
        @test S_pink_noise(w, Ad, 0.0) > S_pink_noise(-w, Ad, 0.0)
        @test S_pink_noise(w, 0.0, Af) ≈ S_pink_noise(-w, 0.0, Af)
        @test S_pink_noise(w, 0.0, Af) ≈ Af^2 * 2π / w
        @test !isfinite(S_pink_noise(0.0, Ad, Af))
    end

    @testset "PinkNoise" begin
        S = PinkNoise(Ad, Af)
        ws = 2π .* [-0.3, 0.1, 0.7]
        @test S.(ws) ≈ S_pink_noise.(ws, Ad, Af)
        @test S(ws[2]) == S_pink_noise(ws[2], Ad, Af)
        Ad2, Af2 = noise_amplitudes(1.0, 2.0, 0.5 + 0im)
        @test Af2 ≈ 1.8e-6 * 1.0 * 0.5
        @test Ad2 ≈ π^2 * 1.1e-6 * 0.25 / 2.0
        S2 = PinkNoise(; El=1.0, Ec=2.0, φ_ge=0.5)
        @test (S2.Ad, S2.Af) == (Ad2, Af2)
    end

    @testset "default_pink_noise" begin
        S = default_pink_noise()
        @test S.Ad > 0 && S.Af > 0
        @test isfinite(S(2π * 0.1))
    end

    @testset "calc_gs / calc_Γφ" begin
        S = PinkNoise(Ad, Af)
        w = 2π * 0.05
        N = 40
        θ = 2π .* (0:N-1) ./ N
        # ⟨φ⟩_upper - ⟨φ⟩_lower = 2a + 2b cos θ  =>  g₀ = a, g_{±1} = b/2
        a, b = 0.7, 0.3
        φu, φl = a .+ b .* cos.(θ), -a .- b .* cos.(θ)
        gs = calc_gs(φu, φl)
        @test gs[1] ≈ a
        @test gs[2] ≈ b / 2 && gs[end] ≈ b / 2
        @test maximum(abs, gs[3:end-1]) < 1e-12
        expected = 8 * Af * a + (b / 2)^2 * (S(w) + S(-w))
        @test calc_Γφ(gs, S, w) ≈ expected
        @test calc_Γφ(gs, S, Af, w) ≈ expected
        @test calc_Γφ(φu, φl, S, w) ≈ expected
        @test calc_Γφ(fftshift_like(gs), S, w; g0_centered=true) ≈ expected
        # a harmonic-1 signal is weighted with S(+w) for e^{-iθ} dependence: sign convention of numpy's ifft
        gs_sin = calc_gs(cos.(θ .+ 0.4), -cos.(θ .+ 0.4))
        @test calc_Γφ(gs_sin, S, w) ≈ (1 / 2)^2 * (S(w) + S(-w))
        @test_throws DimensionMismatch calc_gs(φu, φl[1:end-1])
    end
end
