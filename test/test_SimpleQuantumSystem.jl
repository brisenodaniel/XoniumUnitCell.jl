using Test
using LinearAlgebra
using QuantumToolbox
using XoniumUnitCell

@testset "SimpleQuantumSystem" begin

    # -----------------------------
    # Helper: small fluxonium system
    # -----------------------------
    function small_fluxonium(; name=:f)
        fluxonium(
            EC = 1.0,
            EL = 1.0,
            EJ = 5.0,
            nlev = 6,
            name = name,
        )
    end

    # -----------------------------
    # Constructor behavior
    # -----------------------------
    @testset "constructor" begin
        sys = small_fluxonium()

        @test sys isa SimpleQuantumSystem
        @test sys.name == :f
        @test hasproperty(sys.ops, :H)
        @test ishermitian(sys.H)
        @test size(sys.H) == (6, 6)

        # Missing Hamiltonian should error
        φ, n = φn_operators(EC=1.0, EL=1.0, nlev=4)
        @test_throws ArgumentError SimpleQuantumSystem((φ=φ, n=n))
    end

    # -----------------------------
    # tensor_dims
    # -----------------------------
    @testset "tensor_dims" begin
        sys = small_fluxonium()
        @test tensor_dims(sys) == [6]
    end

    # -----------------------------
    # state and energy
    # -----------------------------
    @testset "state / energy" begin
        sys = small_fluxonium()

        E0, ψ0 = state(sys, :g)
        E1, ψ1 = state(sys, :e)

        @test E1 > E0
        @test norm(ψ0) ≈ 1.0 atol=1e-10
        @test norm(ψ1) ≈ 1.0 atol=1e-10

        @test energy(sys, :g) ≈ E0
        @test energy(sys, 1) ≈ E1
    end

    # -----------------------------
    # to_energy_basis
    # -----------------------------
    @testset "to_energy_basis" begin
        sys = small_fluxonium()
        λ = eigenenergies(sys.H)

        dressed = to_energy_basis(sys; set_energy_zero=true)

        # Hamiltonian should be diagonal
        @test isapprox(norm(dressed.H.data - Diagonal(diag(dressed.H))), 0.0; atol=1e-10)

        # Ground state energy set to zero
        @test abs(real(dressed.H[1,1])) ≤ 1e-12

        # Spectrum preserved up to shift
        λ_ = diag(dressed.H)
        @test isapprox(λ_, λ .- λ[1], atol=1e-12)
    end

    # -----------------------------
    # truncate_sys
    # -----------------------------
    @testset "truncate_sys" begin
        sys = small_fluxonium() |> to_energy_basis
        sys2 = truncate_sys(sys, 4)

        @test size(sys2.H) == (4, 4)
        @test tensor_dims(sys2) == [4]

        # Truncation preserves ordering
        λ = eigenenergies(sys.H)[1:4]
        λ .-= λ[1]
        λ2 = eigenenergies(sys2.H)
        λ2 .-= λ2[1]
        @test λ2 ≈ λ
    end

    # -----------------------------
    # Addition: H + operator
    # -----------------------------
    @testset "operator addition (+)" begin
        sys = small_fluxonium()

        # Add a small perturbation
        δH = 0.01 * sys[:n]^2
        sys′ = sys + δH

        @test sys′ isa SimpleQuantumSystem
        @test ishermitian(sys′.H)
        @test norm(sys′.H - (sys.H + δH)) ≤ 1e-12

        # Dimension mismatch should error
        bad_op = qeye_like(sys.H)[1:4, 1:4] |> Qobj
        @test_throws DimensionMismatch sys + bad_op
    end


    @testset "state / energy vs direct eigenstates" begin
        sys = small_fluxonium()

        # Direct diagonalization of the Hamiltonian operator
        λ_raw, ψ_raw, _ = QuantumToolbox.eigenstates(sys.H)

        # Ground and first excited states
        for (lbl, idx) in zip((:g, :e, :f), (1, 2, 3))
            E, ψ = state(sys, lbl)

            @test E ≈ λ_raw[idx]
            @test ψ ≈ ψ_raw[idx]
            @test norm(ψ) ≈ 1.0 atol=1e-12

            @test energy(sys, lbl) ≈ λ_raw[idx]
        end
    end

    @testset "operator addition vs direct Hamiltonian addition" begin
        sys = small_fluxonium()

        # Define perturbation purely at operator level
        δH = 0.05 * sys[:n]^2

        # Direct operator-level Hamiltonian
        H_direct = sys.H + δH

        # System-level addition
        sys′ = sys + δH


        @test sys′.H ≈ H_direct
        @test ishermitian(sys′.H)

        # Spectra must match
        λ1 = eigenenergies(H_direct)
        λ2 = eigenenergies(sys′.H)

        @test λ1 ≈ λ2
    end

    @testset "to_energy_basis vs manual operator transform" begin
        sys = small_fluxonium()

        # Raw diagonalization
        λ, _, U = QuantumToolbox.eigenstates(sys.H)
        Uq = QuantumObject(U, sys.H.type, sys.H.dims)

        # Manual energy-zeroed Hamiltonian
        H_manual = Uq' * sys.H * Uq
        H_manual -= λ[1]

        dressed = to_energy_basis(sys; set_energy_zero=true)

        # --- Hamiltonian checks ---
        @test dressed.H ≈ H_manual
        @test isapprox(norm(dressed.H.data - Diagonal(diag(dressed.H))), 0.0; atol=1e-12)
        @test abs(real(dressed.H[1,1])) ≤ 1e-12

        # --- Operator checks ---
        for opkey in (:φ, :n)
            O_manual = Uq' * sys[opkey] * Uq
            @test dressed[opkey] ≈ O_manual
        end

        # --- Spectrum consistency ---
        @test diag(dressed.H) ≈ (λ .- λ[1])
    end
end
