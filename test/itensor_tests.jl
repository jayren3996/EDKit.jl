@testset "ITensor Integration" begin
    @testset "MPS/Vector Conversion" begin
        for L in 2:4, d in 2:3
            v = randn(ComplexF64, d^L) |> normalize
            s = siteinds(d, L)
            ψ = vec2mps(v, s)
            @test mps2vec(ψ) ≈ v
        end
    end

    @testset "Product MPS And Operator Conversion" begin
        s = siteinds(2, 4)
        states = [[1.0, 0.0], [0.0, 1.0], [1 / sqrt(2), 1 / sqrt(2)], [1.0, 0.0]]
        ψp = EDKit.productstate(s, states)
        @test norm(mps2vec(ψp) - kron(states...)) < 1e-12

        mat = randn(ComplexF64, 4, 4)
        op = mat2op(mat, s[1], s[2])
        @test op2mat(op, s[1], s[2]) ≈ mat
    end

    @testset "Symmetry-Resolved mps2vec" begin
        L = 6
        B = basis(L = L, N = L ÷ 2, k = 0, p = 1)
        P = sector_embedding(B)
        v_sector = randn(ComplexF64, size(B, 1)) |> normalize
        ψ_full = P * v_sector
        ψ_mps = vec2mps(ψ_full, siteinds(B.B, length(B)))
        @test mps2vec(ψ_mps, B) ≈ v_sector atol = 1e-10

        B2 = basis(L = L, N = L ÷ 2, k = 0, p = 1, a = 2)
        P2 = sector_embedding(B2)
        v2 = randn(ComplexF64, size(B2, 1)) |> normalize
        ψ2 = vec2mps(P2 * v2, siteinds(B2.B, length(B2)))
        @test mps2vec(ψ2, B2) ≈ v2 atol = 1e-10
    end

    @testset "Pauli Basis Utilities" begin
        for L in 1:3
            coeffs = randn(4^L)
            @test pauli_list(pauli(coeffs)) ≈ coeffs
        end

        H = randn(ComplexF64, 4, 4) |> Hermitian |> Array
        ρ = randn(ComplexF64, 4, 4)
        ρ = ρ * ρ'
        ρ ./= tr(ρ)
        @test pauli(commutation_mat(H) * pauli_list(ρ)) ≈ -1im * (H * ρ - ρ * H)

        Lm = randn(ComplexF64, 4, 4)
        @test pauli(dissipation_mat(Lm) * pauli_list(ρ)) ≈ Lm * ρ * Lm' - (Lm' * Lm * ρ + ρ * Lm' * Lm) / 2
    end

    @testset "PMPS And MPO Conversion" begin
        s = siteinds("S=1/2", 3)
        ps = siteinds("Pauli", 3)
        ψ = vec2mps(randn(ComplexF64, 2^3) |> normalize, s)
        pmps = mps2pmps(ψ, ps)
        mpo = pmps2mpo(pmps, s)
        pmpo = mpo2pmpo(mpo, ps)
        @test length(pmps) == 3
        @test length(mpo) == 3
        @test length(pmpo) == 3
    end

    @testset "TEBD And DAOE Smoke Tests" begin
        s = siteinds("S=1/2", 4)
        ps = siteinds("Pauli", 4)
        h2 = [Array(spin((0.1, "xx"), (0.2, "zz"))) for _ in 1:3]
        gates = tebd4(h2, s, 0.05)
        @test !isempty(gates)

        ψ = vec2mps(randn(ComplexF64, 2^4) |> normalize, s)
        idgates = [mat2op(Matrix(I, 4, 4), s[i], s[i+1]) for i in 1:3]
        before = mps2vec(ψ)
        tebd_n!(ψ, idgates, reverse(idgates), cutoff = 1e-12, maxdim = 16)
        @test norm(mps2vec(ψ)) ≈ 1.0 atol = 1e-10
        @test abs(abs(dot(before, mps2vec(ψ))) - 1) < 1e-8

        D = daoe(ps, 2, 0.3)
        @test length(D) == 4

        ψI = productMPS(ps, fill("I", 4))
        ψX = productMPS(ps, ["X", "I", "I", "I"])
        ψXX = productMPS(ps, ["X", "X", "I", "I"])
        @test inner(ψI, apply(D, ψI)) ≈ 1.0
        @test inner(ψX, apply(D, ψX)) ≈ 1.0
        @test inner(ψXX, apply(D, ψXX)) ≈ exp(-0.3)

        FD = fdaoe(ps, 2, 0.3)
        ψZ = productMPS(ps, ["Z", "I", "I", "I"])
        ψZZ = productMPS(ps, ["Z", "Z", "I", "I"])
        ψXZ = productMPS(ps, ["X", "Z", "I", "I"])
        @test inner(ψZ, apply(FD, ψZ)) ≈ 1.0
        @test inner(ψZZ, apply(FD, ψZZ)) ≈ 1.0
        @test inner(ψXX, apply(FD, ψXX)) ≈ 1.0
        @test inner(ψXZ, apply(FD, ψXZ)) ≈ exp(-0.6)
    end
end
