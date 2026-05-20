@testset "Core Operators And Helpers" begin
    @test Array(spin("X")) ≈ [0 1; 1 0]
    @test Array(spin("Y")) ≈ [0 -1im; 1im 0]
    @test Array(spin("X")) ≈ 2 .* Array(spin("x"))
    @test Array(spin("Y")) ≈ 2 .* Array(spin("y"))
    @test Array(spin("Z")) ≈ 2 .* Array(spin("z"))
    @test Array(spin("xx")) ≈ kron(Array(spin("x")), Array(spin("x")))
    @test Array(spin("XX")) ≈ 4 .* Array(spin("xx"))
    @test Array(spin("YY")) ≈ 4 .* Array(spin("yy"))
    @test Array(spin("ZZ")) ≈ 4 .* Array(spin("zz"))
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

    assert_operator_matches_dense(H1; rng=stable_rng(101))

    Br = TranslationalBasis(L = L, k = 0, N = 2, threaded = false)
    Hr = trans_inv_operator(bond, 2, Br)
    assert_operator_matches_dense(Hr; rng=stable_rng(102))

    @testset "Operator application agrees with dense matrices across basis types" begin
        Lsmall = 4
        Nsmall = 2
        conserving_bond = spin((1.0, "+-"), (1.0, "-+"), (0.2, "zz"))
        custom_shift = [2, 3, 4, 1]
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
            basis(L = Lsmall, N = Nsmall, symmetries = [(custom_shift, 0)], threaded = false),
        ]

        for (case_id, Bcase) in enumerate(basis_cases)
            iszero(size(Bcase, 1)) && continue
            Hcase = trans_inv_operator(conserving_bond, 2, Bcase)
            assert_operator_matches_dense(Hcase; rng=stable_rng(200 + case_id))
        end
    end

    @testset "Abelian operator application does not mutate group state" begin
        Lsmall = 4
        custom_shift = [2, 3, 4, 1]
        Bsym = basis(L = Lsmall, N = 2, symmetries = [(custom_shift, 0)], threaded = false)
        Hsym = trans_inv_operator(spin((1.0, "+-"), (1.0, "-+"), (0.2, "zz")), 2, Bsym)
        v = random_vector(size(Hsym, 2); rng=stable_rng(301))
        M = random_matrix(size(Hsym, 2), 2; rng=stable_rng(302))

        Bsym.G.s .= 2
        initial_state = copy(Bsym.G.s)
        Hsym * v
        @test Bsym.G.s == initial_state
        EDKit.mul(Hsym, v)
        @test Bsym.G.s == initial_state
        Hsym * M
        @test Bsym.G.s == initial_state
        EDKit.mul(Hsym, M)
        @test Bsym.G.s == initial_state
    end

    @testset "Explicit sparse operator construction matches dense reference" begin
        Lsmall = 4
        Nsmall = 2
        conserving_bond = spin((1.0, "+-"), (1.0, "-+"), (0.2, "zz"))
        custom_shift = [2, 3, 4, 1]
        sparse_basis_cases = [
            TensorBasis(L = Lsmall, base = 2),
            ProjectedBasis(L = Lsmall, N = Nsmall, threaded = false),
            TranslationalBasis(L = Lsmall, N = Nsmall, k = 0, threaded = false),
            basis(L = Lsmall, N = Nsmall, k = 0, p = 1, z = 1, threaded = false),
            basis(L = Lsmall, N = Nsmall, symmetries = [(custom_shift, 0)], threaded = false),
        ]

        for (case_id, Bcase) in enumerate(sparse_basis_cases)
            iszero(size(Bcase, 1)) && continue
            Hcase = trans_inv_operator(conserving_bond, 2, Bcase)
            dense = Array(Hcase)
            Mcase = random_matrix(size(Hcase, 2), 2; rng=stable_rng(400 + case_id))
            @test sparse(Hcase) ≈ sparse(dense)
            clear_sparse_cache!()
            @test sparse!(Hcase) ≈ sparse(dense)
            @test Hcase * Mcase ≈ dense * Mcase
            clear_sparse_cache!()
        end

        Hduplicate = operator([spin("z"), spin("z")], [[1], [2]], TensorBasis(L = 2, base = 2))
        @test sparse(Hduplicate) ≈ sparse(Array(Hduplicate))

        Hlarge_duplicate = operator(
            [spin("z"), 2 * spin("z"), -0.5 * spin("z")],
            [[1], [1], [1]],
            TensorBasis(L = 4, base = 2),
        )
        @test eltype(sparse(Hlarge_duplicate)) == eltype(Hlarge_duplicate)
        @test sparse(Hlarge_duplicate) ≈ sparse(Array(Hlarge_duplicate))
    end

    @testset "TensorBasis base-2 operator application matches independent dense embedding" begin
        Lfast = 5
        Bfast = TensorBasis(L = Lfast, base = 2)
        mats = [
            spin((1.0, "xx"), (0.7, "yy"), (0.3, "zz")),
            spin((0.4, "x"), (0.2im, "y"), (-0.1, "z")),
            spin((0.6, "+-"), (-0.5, "-+"), (0.25, "zz")),
            spin((0.2, "xx"), (-0.4, "yy"), (0.15, "zz")),
        ]
        inds = [[1, 2], [4], [5, 1], [3, 1]]
        Hfast = operator(mats, inds, Bfast)
        dense_ref = sum(dense_local_embedding(mats[i], inds[i], Lfast; base = 2) for i in eachindex(mats))
        vfast = random_vector(size(Hfast, 2); rng=stable_rng(501))
        Mfast = random_matrix(size(Hfast, 2), 4; rng=stable_rng(502))

        @test Hfast * vfast ≈ dense_ref * vfast
        @test EDKit.mul(Hfast, vfast) ≈ dense_ref * vfast

        clear_sparse_cache!()
        @test Hfast * Mfast ≈ dense_ref * Mfast
        @test EDKit.mul(Hfast, Mfast) ≈ dense_ref * Mfast

        target = similar(vfast, size(Hfast, 1))
        mul!(target, Hfast, vfast, 1, 0)
        @test target ≈ dense_ref * vfast

        mtarget = similar(Mfast, size(Hfast, 1), size(Mfast, 2))
        mul!(mtarget, Hfast, Mfast, 1, 0)
        @test mtarget ≈ dense_ref * Mfast

        @test sparse(Hfast) ≈ sparse(dense_ref)
        @test Array(Hfast) ≈ dense_ref

        @test sparse!(Hfast) ≈ sparse(dense_ref)
        @test Hfast * Mfast ≈ dense_ref * Mfast
        @test EDKit.mul(Hfast, Mfast) ≈ dense_ref * Mfast
        clear_sparse_cache!()

        mul!(target, Hfast, vfast, 1, 0)
        @test (@allocated mul!(target, Hfast, vfast, 1, 0)) < sizeof(Int) * length(Bfast)
        @test (@allocated mul!(mtarget, Hfast, Mfast, 1, 0)) < sizeof(Int) * length(Bfast)
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

    @test isnan(meangapratio(Float64[]))
    @test isnan(meangapratio([0.0]))
    @test isnan(meangapratio([0.0, 1.0]))
    @test meangapratio([0.0, 1.0, 3.0]) ≈ 0.5

    # SPIN_CACHE is concretely typed so cache hits are type-stable.
    @test valtype(EDKit.SPIN_CACHE) === SparseMatrixCSC{ComplexF64, Int}
    @test spin("X") isa SparseMatrixCSC{ComplexF64, Int}
    @test spin("xx") isa SparseMatrixCSC{ComplexF64, Int}
    @test Matrix(spin("X")) ≈ ComplexF64[0 1; 1 0]
    @test Matrix(spin("Y")) ≈ ComplexF64[0 -im; im 0]
    @test Matrix(spin("Z")) ≈ ComplexF64[1 0; 0 -1]

    @testset "_SPARSE_CACHE evicts oldest entries past the 8-operator cap" begin
        clear_sparse_cache!()
        # 9 distinct operators (different bases means distinct objectid).
        ops = [trans_inv_operator(spin("zz"), 2, TensorBasis(L = L, base = 2)) for L in 2:10]
        for opt in ops
            sparse!(opt)
            @test EDKit._cached_sparse(opt) !== nothing
        end
        # First operator should now have been evicted by the LRU.
        @test EDKit._cached_sparse(ops[1]) === nothing
        @test EDKit._cached_sparse(ops[end]) !== nothing
        clear_sparse_cache!()
        @test EDKit._cached_sparse(ops[end]) === nothing
    end
end
