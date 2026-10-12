@testset "CompositeQuantumSystem" begin

    # -------------------------------------------------
    # Helpers: two small subsystems
    # -------------------------------------------------
    function small_fluxonium(name; ΔEL=0.0)
        fluxonium(
            EC = 1.0,
            EL = 1.0 + ΔEL,
            EJ = 5.0,
            nlev = 5,
            name = name,
        ) |> sys -> to_energy_basis(sys; etol=1e-9)
    end

    A = small_fluxonium(:A)
    B = small_fluxonium(:B, ΔEL=0.05)

    # -------------------------------------------------
    # Tensor product operator (⊗)
    # -------------------------------------------------
    @testset "tensor product ⊗" begin
        AB = A ⊗ B

        @test AB isa CompositeQuantumSystem
        @test keys(AB.subsys) == (:A, :B)

        # Manual Hamiltonian tensor product
        H_manual = tensor(A.H, qeye(B)) + tensor(qeye(A), B.H)

        @test AB.H ≈ H_manual
        @test ishermitian(AB.H)
    end

    AB = A ⊗ B

    # -------------------------------------------------
    # Operator lifting: single-subsystem
    # -------------------------------------------------
    @testset "operator lifting (single)" begin
        # Manual lifted operators
        φA_manual = tensor(A[:φ], qeye(B))
        φB_manual = tensor(qeye(A), B[:φ])

        @test AB[:A, :φ] ≈ φA_manual
        @test AB[:B, :φ] ≈ φB_manual
    end

    # -------------------------------------------------
    # Operator lifting: multi-operator
    # -------------------------------------------------
    @testset "operator lifting (multi)" begin
        # Manual construction
        O_manual = tensor(A[:φ], B[:n])

        @test AB[:A, :φ, :B, :n] ≈ O_manual
    end

    # -------------------------------------------------
    # Bare state consistency
    # -------------------------------------------------
    @testset "bare product states" begin
        # Bare subsystem states
        λA, ψA = state(A, :g)
        λB, ψB = state(B, :e)
        λB0, _ = state(B, :g)


        # Manual bare product state
        ψ_manual = tensor(ψA, ψB)
        λ_manual = λA * λB

        # Composite bare state
        λ, ψ = state(AB, :ge, true)
        λ0, _ = state(AB, :gg, true)

        @test λ  ≈ λ_manual
        @test ψ ≈ ψ_manual
        @test norm(ψ) ≈ 1.0 atol=1e-12
    end

    # -------------------------------------------------
    # Dressed vs bare state consistency
    # -------------------------------------------------
    @testset "dressed vs bare state mapping" begin
        # Bare product state
        _, ψ_bare = state(AB, :ge, true)

        # Dressed eigenstates
        λ, ψs, _ = QuantumToolbox.eigenstates(AB.H)

        # Find dressed eigenstate with maximal overlap
        overlaps = [abs(ψ_bare' * ψ) for ψ in ψs]
        idx = argmax(overlaps)

        λ_dressed, ψ_dressed = state(AB, :ge, false)

        @test λ_dressed ≈ λ[idx]
        @test ψ_dressed ≈ ψs[idx]
    end

    # -------------------------------------------------
    # Basis transformation correctness
    # -------------------------------------------------
    @testset "to_energy_basis (composite)" begin
        # Manual diagonalization
        λ, _, U = QuantumToolbox.eigenstates(AB.H)
        Uq = QuantumObject(U, AB.H.type, AB.H.dims)

        # Manual transformed operators
        φA_manual = Uq' * tensor(A[:φ], qeye(B)) * Uq
        φB_manual = Uq' * tensor(qeye(A), B[:φ]) * Uq

        dressed = to_energy_basis(AB; set_energy_zero=true)

        # Hamiltonian checks
        @test dressed.H ≈ ((Uq' * AB.H * Uq) - λ[1])
        @test abs(real(dressed.H[1,1])) ≤ 1e-12
        @test isapprox(norm(dressed.H.data - Diagonal(diag(dressed.H))), 0.0; atol=1e-12)

        # Operator checks
        @test dressed[:A, :φ] ≈ φA_manual
        @test dressed[:B, :φ] ≈ φB_manual
    end


    @testset "bare_to_dressed_tracking" begin
        AB = A ⊗ B
        @test state(AB, :gg; bare_to_dressed_tracking=:spectrum)[1] ≈ first(eigenstates(AB.H))[1]
        @test state(AB, :ge) == state(AB, :ge; bare_to_dressed_tracking=:spectrum)
        @test state_idx(AB, :ge; bare_to_dressed_tracking=:overlap) isa Int
        # weakly coupled: both modes agree on low-lying states
        for lbl in (:gg, :ge, :eg, :ee)
            @test state_idx(AB, lbl; bare_to_dressed_tracking=:overlap) == state_idx(AB, lbl)
        end
        # spectrum rank matches bare energy ordering
        λd = first(eigenstates(AB.H))
        @test energy(AB, :ee; bare_to_dressed_tracking=:spectrum) ≈
              λd[state_idx(AB, :ee; bare_to_dressed_tracking=:spectrum)]
        @test_throws ArgumentError state(AB, :gg; bare_to_dressed_tracking=:bogus)
        @test state(A, :e; bare_to_dressed_tracking=:spectrum) == state(A, :e)
        @test_throws ArgumentError state(A, :e; bare_to_dressed_tracking=:bogus)
    end
end
