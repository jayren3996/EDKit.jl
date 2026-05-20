@testset "Basis Constructors And Symmetry Decomposition" begin
    L = 6
    Nhalf = L ÷ 2

    pxp(v) = all(v[i] == 0 || v[i + 1] == 0 for i in 1:length(v)-1)
    BP = ProjectedBasis(L = L, f = pxp)
    @test size(BP, 1) == 21

    BN = ProjectedBasis(L = L, N = Nhalf)
    @test size(BN, 1) == binomial(L, Nhalf)

    @testset "Fixed-N uses EDKit digit convention" begin
        dgt = zeros(Int, 4)

        Bproj = ProjectedBasis(L = 4, N = 1, threaded = false)
        @test size(Bproj, 1) == 4
        for idx in Bproj.I
            change!(dgt, idx; base = 2)
            @test count(==(0), dgt) == 1
            @test sum(dgt) == 3
        end

        shift = [4, 1, 2, 3]
        Bsym = basis(; L = 4, N = 1, symmetries = [(shift, 0)], threaded = false)
        for idx in Bsym.I
            change!(dgt, idx; base = 2)
            @test count(==(0), dgt) == 1
            @test sum(dgt) == 3
        end
    end

    XXZ = spin((1.0, "xx"), (1.0, "yy"), (0.4, "zz"))
    full_vals = eigvals(Hermitian(trans_inv_operator(XXZ, 2, L))) |> sort

    vals_N = Float64[]
    for N in 0:L
        B = basis(L = L, N = N)
        append!(vals_N, eigvals(Hermitian(trans_inv_operator(XXZ, 2, B))))
    end
    @test sort(vals_N) ≈ full_vals

    vals_k = Float64[]
    for k in 0:L-1
        B = basis(L = L, k = k)
        append!(vals_k, eigvals(Hermitian(trans_inv_operator(XXZ, 2, B))))
    end
    @test sort(vals_k) ≈ full_vals

    Bk = basis(L = L, N = Nhalf, k = 0)
    vals_kp = Float64[]
    for p in (1, -1)
        B = basis(L = L, N = Nhalf, k = 0, p = p)
        append!(vals_kp, eigvals(Hermitian(trans_inv_operator(XXZ, 2, B))))
    end
    @test sort(vals_kp) ≈ eigvals(Hermitian(trans_inv_operator(XXZ, 2, Bk))) |> sort

    vals_kz = Float64[]
    for z in (1, -1)
        B = basis(L = L, N = Nhalf, k = 1, z = z)
        append!(vals_kz, eigvals(Hermitian(trans_inv_operator(XXZ, 2, B))))
    end
    @test sort(vals_kz) ≈ eigvals(Hermitian(trans_inv_operator(XXZ, 2, basis(L = L, N = Nhalf, k = 1)))) |> sort

    vals_a2 = Float64[]
    for k in 0:(L ÷ 2 - 1)
        B = TranslationalBasis(L = L, k = k, a = 2)
        append!(vals_a2, eigvals(Hermitian(trans_inv_operator(XXZ, 2, B))))
    end
    @test sort(vals_a2) ≈ full_vals

    @testset "Fixed-charge small_N paths match full scans" begin
        BP_small = ProjectedBasis(L = 8, N = 3, small_N = true)
        BP_scan = ProjectedBasis(L = 8, N = 3, threaded = false)
        @test BP_small.I == BP_scan.I

        BT_small = TranslationalBasis(L = 8, N = 3, k = 2, small_N = true)
        BT_scan = TranslationalBasis(L = 8, N = 3, k = 2, threaded = false)
        @test BT_small.I == BT_scan.I
        @test BT_small.R ≈ BT_scan.R

        BP3_small = ProjectedBasis(L = 4, N = 4, base = 3, small_N = true)
        BP3_scan = ProjectedBasis(L = 4, N = 4, base = 3, threaded = false)
        @test BP3_small.I == BP3_scan.I

        BT3_small = TranslationalBasis(L = 4, N = 4, k = 0, base = 3, small_N = true)
        BT3_scan = TranslationalBasis(L = 4, N = 4, k = 0, base = 3, threaded = false)
        @test BT3_small.I == BT3_scan.I
        @test BT3_small.R ≈ BT3_scan.R
    end

    @testset "Threaded basis constructors match single-thread scans" begin
        basis_makers = [
            threaded -> ProjectedBasis(L = 5, N = 2, threaded = threaded),
            threaded -> ProjectedBasis(L = 2, f = _ -> true, threaded = threaded),
            threaded -> ParityBasis(L = 5, N = 2, p = 1, threaded = threaded),
            threaded -> FlipBasis(L = 6, N = 3, p = 1, threaded = threaded),
            threaded -> ParityFlipBasis(L = 6, N = 3, p = 1, z = 1, threaded = threaded),
            threaded -> TranslationalBasis(L = 6, N = 3, k = 1, threaded = threaded),
            threaded -> TranslationParityBasis(L = 6, N = 3, k = 0, p = 1, threaded = threaded),
            threaded -> TranslationFlipBasis(L = 6, N = 3, k = 1, p = 1, threaded = threaded),
            threaded -> basis(L = 6, N = 3, k = 0, p = 1, threaded = threaded),
        ]

        for make_basis in basis_makers
            assert_threaded_basis_matches_serial(make_basis)
        end
    end

    @testset "Fixed-charge edge sectors are explicit" begin
        Bempty = ProjectedBasis(L = 4, N = 5, threaded = false)
        @test size(Bempty, 1) == 0
        @test isempty(Bempty.I)
        @test index(Bempty, [0, 0, 0, 0]) == (0, 1)
        @test_throws ErrorException index(Bempty, [0, 0, 0, 0]; check = true)

        Bvac = ProjectedBasis(L = 4, N = 0, threaded = false)
        Bfull = ProjectedBasis(L = 4, N = 4, threaded = false)
        @test Bvac.I == [16]
        @test Bfull.I == [1]

        BT_empty = TranslationalBasis(L = 4, N = 5, k = 0, threaded = false)
        @test size(BT_empty, 1) == 0
        @test isempty(BT_empty.I)

        BT_vac = TranslationalBasis(L = 4, N = 0, k = 0, threaded = false)
        @test BT_vac.I == [16]
        @test BT_vac.R ≈ [4.0]
    end

    @testset "Sector embeddings are isometries" begin
        for B in (
            TensorBasis(L = 4, base = 2),
            ProjectedBasis(L = 4, N = 2, threaded = false),
            TranslationalBasis(L = 4, N = 2, k = 0, threaded = false),
            basis(L = 4, N = 2, k = 0, p = 1, threaded = false),
        )
            assert_sector_embedding_isometry(B)
        end
    end

    @testset "TranslationFlipBasis copy preserves fields" begin
        B = TranslationFlipBasis(L = 4, k = 0, p = 1)
        Bcopy = copy(B)
        @test Bcopy isa TranslationFlipBasis
        @test Bcopy.M == B.M
        @test Bcopy.B == B.B
        @test Bcopy.I == B.I
    end

    @testset "TranslationFlipBasis default predicate is no-op" begin
        B_default = TranslationFlipBasis(L = 6, k = 1, p = 1)
        B_true = TranslationFlipBasis(L = 6, k = 1, p = 1, f = _ -> true)
        B_nothing = TranslationFlipBasis(L = 6, k = 1, p = 1, f = nothing)
        @test B_default.I == B_true.I == B_nothing.I
        @test B_default.R ≈ B_true.R
        @test B_default.R ≈ B_nothing.R
    end

    @testset "TranslationFlipBasis requires even L" begin
        # Even L is supported.
        @test_nowarn TranslationFlipBasis(L = 4, k = 0, p = 1)
        @test_nowarn TranslationFlipBasis(L = 6, k = 0, p = 1)

        # Odd L must error — flip-parity check uses L÷2 integer division
        # and is not correct for odd L.
        @test_throws AssertionError TranslationFlipBasis(L = 5, k = 0, p = 1)
        @test_throws AssertionError TranslationFlipBasis(L = 7, k = 0, p = 1)
    end
end
