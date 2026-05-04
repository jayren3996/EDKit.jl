@testset "Core Operators And Helpers" begin
    @test Array(spin("X")) ≈ [0 1; 1 0]
    @test Array(spin("Y")) ≈ [0 -1im; 1im 0]
    @test Array(spin("xx")) ≈ kron(Array(spin("x")), Array(spin("x")))
    @test Array(spin((1.0, "xx"), (1.0, "yy"), (0.5, "zz"))) ≈ Array(spin("xx") + spin("yy") + 0.5 * spin("zz"))

    L = 4
    mat = spin((1.0, "xx"), (0.3, "zz"))
    inds = [2, 3]
    H = operator(mat, inds, L)
    @test Array(H) ≈ kron_embed(mat, inds, L)
    @test sparse(H) ≈ sparse(Array(H))

    bond = spin((1.0, "xx"), (1.0, "yy"))
    H1 = trans_inv_operator(bond, 2, L)
    H2 = operator(fill(bond, L), [[i, mod(i, L) + 1] for i in 1:L], L)
    @test Array(H1) ≈ Array(H2)

    v = randn(ComplexF64, 2^L)
    @test H1 * v ≈ Array(H1) * v
    @test EDKit.mul(H1, v) ≈ Array(H1) * v

    accumulated = ones(ComplexF64, size(H1, 1))
    mul!(accumulated, H1, v)
    @test accumulated ≈ ones(ComplexF64, size(H1, 1)) + Array(H1) * v

    scaled = ones(ComplexF64, size(H1, 1))
    mul!(scaled, H1, v, 2, 3)
    @test scaled ≈ 2 * (Array(H1) * v) .+ 3

    M = randn(ComplexF64, size(H1, 2), 3)
    @test H1 * M ≈ Array(H1) * M
    @test EDKit.mul(H1, M) ≈ Array(H1) * M

    accumulated_matrix = ones(ComplexF64, size(H1, 1), size(M, 2))
    mul!(accumulated_matrix, H1, M)
    @test accumulated_matrix ≈ ones(ComplexF64, size(H1, 1), size(M, 2)) + Array(H1) * M

    scaled_matrix = ones(ComplexF64, size(H1, 1), size(M, 2))
    mul!(scaled_matrix, H1, M, 2, 3)
    @test scaled_matrix ≈ 2 * (Array(H1) * M) .+ 3

    Br = TranslationalBasis(L = L, k = 0, N = 2, threaded = false)
    Hr = trans_inv_operator(bond, 2, Br)
    Hrd = Array(Hr)
    vr = randn(ComplexF64, size(Hr, 2))
    @test Hr * vr ≈ Hrd * vr
    @test EDKit.mul(Hr, vr) ≈ Hrd * vr

    scaled_reduced = ones(ComplexF64, size(Hr, 1))
    mul!(scaled_reduced, Hr, vr, 2, 3)
    @test scaled_reduced ≈ 2 * (Hrd * vr) .+ 3

    Mr = randn(ComplexF64, size(Hr, 2), 3)
    @test Hr * Mr ≈ Hrd * Mr
    @test EDKit.mul(Hr, Mr) ≈ Hrd * Mr

    scaled_reduced_matrix = ones(ComplexF64, size(Hr, 1), size(Mr, 2))
    mul!(scaled_reduced_matrix, Hr, Mr, 2, 3)
    @test scaled_reduced_matrix ≈ 2 * (Hrd * Mr) .+ 3

    @testset "Operator application agrees with dense matrices across basis types" begin
        Lsmall = 4
        Nsmall = 2
        conserving_bond = spin((1.0, "+-"), (1.0, "-+"), (0.2, "zz"))
        basis_cases = [
            TensorBasis(L = Lsmall, base = 2),
            ProjectedBasis(L = Lsmall, N = Nsmall, threaded = false),
            ParityBasis(L = Lsmall, N = Nsmall, p = 1, threaded = false),
            FlipBasis(L = Lsmall, N = Nsmall, p = 1, threaded = false),
            ParityFlipBasis(L = Lsmall, N = Nsmall, p = 1, z = 1, threaded = false),
            TranslationalBasis(L = Lsmall, N = Nsmall, k = 0, threaded = false),
            TranslationalBasis(L = Lsmall, N = Nsmall, k = 1, threaded = false),
            TranslationParityBasis(L = Lsmall, N = Nsmall, k = 0, p = 1, threaded = false),
            TranslationFlipBasis(L = Lsmall, N = Nsmall, k = 0, p = 1, threaded = false),
            basis(L = Lsmall, N = Nsmall, k = 0, p = 1, z = 1, threaded = false),
        ]

        for Bcase in basis_cases
            iszero(size(Bcase, 1)) && continue
            Hcase = trans_inv_operator(conserving_bond, 2, Bcase)
            Hdense = Array(Hcase)
            vcase = randn(ComplexF64, size(Hcase, 2))
            Mcase = randn(ComplexF64, size(Hcase, 2), 3)

            @test Hcase * vcase ≈ Hdense * vcase
            @test EDKit.mul(Hcase, vcase) ≈ Hdense * vcase

            y = ones(ComplexF64, size(Hcase, 1))
            mul!(y, Hcase, vcase)
            @test y ≈ ones(ComplexF64, size(Hcase, 1)) + Hdense * vcase

            y_scaled = ones(ComplexF64, size(Hcase, 1))
            mul!(y_scaled, Hcase, vcase, 2, 3)
            @test y_scaled ≈ 2 * (Hdense * vcase) .+ 3

            @test Hcase * Mcase ≈ Hdense * Mcase
            @test EDKit.mul(Hcase, Mcase) ≈ Hdense * Mcase

            Y = ones(ComplexF64, size(Hcase, 1), size(Mcase, 2))
            mul!(Y, Hcase, Mcase)
            @test Y ≈ ones(ComplexF64, size(Hcase, 1), size(Mcase, 2)) + Hdense * Mcase

            Y_scaled = ones(ComplexF64, size(Hcase, 1), size(Mcase, 2))
            mul!(Y_scaled, Hcase, Mcase, 2, 3)
            @test Y_scaled ≈ 2 * (Hdense * Mcase) .+ 3

            clear_sparse_cache!()
            @test sparse!(Hcase) ≈ sparse(Hdense)
            @test Hcase * Mcase ≈ Hdense * Mcase
            @test EDKit.mul(Hcase, Mcase) ≈ Hdense * Mcase
            clear_sparse_cache!()
        end
    end

    B = TensorBasis(L = L, base = 2)
    state = productstate([0, 1, 0, 1], B)
    @test count(!iszero, state) == 1
    @test state[index([0, 1, 0, 1], base = 2)] == 1

    bell = normalize([1.0, 0.0, 0.0, 1.0])
    @test ent_S(bell, [1], TensorBasis(L = 2, base = 2)) ≈ log(2)
    @test ent_S(bell, [1], 2) ≈ log(2)

    E = [0.0, 1.0, 3.0, 6.0]
    @test gapratio(E) ≈ [0.5, 2 / 3]
    @test meangapratio(E) ≈ 7 / 12
end
