@testset "Linear Maps And Algorithms" begin
    L = 4
    N = 2

    Bfull = TensorBasis(L = L, base = 2)
    ψfull = randn(ComplexF64, size(Bfull, 1))
    T = DoubleBasis(Bfull, Bfull)
    Sfull = symmetrizer(T)
    @test Sfull * ψfull ≈ ψfull
    @test Sfull ≈ Matrix{ComplexF64}(I, size(Bfull, 1), size(Bfull, 1))

    Bpar = ParityBasis(L = L, p = 1, base = 2)
    Tpar = DoubleBasis(Bpar, Bfull)
    Ppar = sector_embedding(Bpar)
    Spar = symmetrizer(Tpar)
    @test Spar ≈ Ppar'
    @test Tpar(ψfull) ≈ Spar * ψfull

    ψproj = Ppar * (Spar * ψfull)
    @test ψproj ≈ Ppar * (Ppar' * ψfull)

    bases = [
        ("tensor", TensorBasis(L = L, base = 2), randn(ComplexF64, size(Bfull, 1))),
        ("projected", ProjectedBasis(L = L, N = N, base = 2, threaded = false), randn(ComplexF64, binomial(L, N))),
        ("abelian", basis(L = L, N = N, k = 0, p = 1, base = 2, threaded = false), nothing),
        ("parity", ParityBasis(L = L, p = 1, base = 2, threaded = false), nothing),
        ("flip", FlipBasis(L = L, p = 1, N = N, base = 2, threaded = false), nothing),
        ("parity-flip", ParityFlipBasis(L = L, p = 1, z = 1, N = N, base = 2, threaded = false), nothing),
        ("translation", TranslationalBasis(L = L, k = 0, N = N, base = 2, threaded = false), nothing),
        ("translation-parity", TranslationParityBasis(L = L, k = 0, p = 1, N = N, base = 2, threaded = false), nothing),
        ("translation-flip", TranslationFlipBasis(L = L, k = 0, p = 1, N = N, base = 2, threaded = false), nothing),
    ]
    bases = [(name, B, isnothing(v) ? randn(ComplexF64, size(B, 1)) : v) for (name, B, v) in bases]

    @testset "DoubleBasis fast action matches symmetrizer for all basis types" begin
        for (target_name, Btarget, _) in bases
            for (source_name, Bsource, v) in bases
                T = DoubleBasis(Btarget, Bsource)
                @test T(v) ≈ symmetrizer(T) * v atol = 1e-10 rtol = 1e-10
            end
        end
    end

    @test_throws ArgumentError DoubleBasis(TensorBasis(L = 3, base = 2), TensorBasis(L = 4, base = 2))
    @test_throws ArgumentError DoubleBasis(TensorBasis(L = 4, base = 2), TensorBasis(L = 4, base = 3))

    mat = rand(2, 2) |> Hermitian
    Hrect = trans_inv_operator(mat, 1, Tpar) |> Array
    @test size(Hrect) == (size(Bpar, 1), size(Bfull, 1))

    Bn = ProjectedBasis(L = L, N = N, base = 2, threaded = false)
    Bnp = ProjectedBasis(L = L, N = N + 1, base = 2, threaded = false)
    Tup = DoubleBasis(Bnp, Bn)
    Splus_rect = operator(fill(spin("+"), L), collect(1:L), Tup) |> Array
    Splus_full = operator(fill(spin("+"), L), collect(1:L), Bfull) |> Array
    Psrc = symmetrizer(DoubleBasis(Bn, Bfull))
    Ptgt = symmetrizer(DoubleBasis(Bnp, Bfull))
    @test size(Splus_rect) == (size(Bnp, 1), size(Bn, 1))
    @test Splus_rect ≈ Ptgt * Splus_full * Psrc'

    @testset "DoubleBasis between momentum sectors respects conservation" begin
        Bk0 = TranslationalBasis(L = L, k = 0, base = 2, threaded = false)
        Bk1 = TranslationalBasis(L = L, k = 1, base = 2, threaded = false)
        Tk01 = DoubleBasis(Bk1, Bk0)

        # Sizes for the rectangular inter-sector operator.
        @test size(Tk01) == (size(Bk1, 1), size(Bk0, 1))

        # A translation-invariant operator built on the DoubleBasis must
        # produce a zero matrix between distinct momentum sectors — this
        # is the defining physical property of momentum conservation.
        H_xx = trans_inv_operator(spin((1.0, "xx"), (1.0, "yy")), 2, Tk01) |> Array
        @test size(H_xx) == (size(Bk1, 1), size(Bk0, 1))
        @test norm(H_xx) < 1e-12

        # The diagonal-sector check guards against a degenerate "always zero"
        # DoubleBasis implementation: same operator on Tk00 = DoubleBasis(Bk0, Bk0)
        # must reproduce the within-sector operator and be non-zero.
        Tk00 = DoubleBasis(Bk0, Bk0)
        H_xx_00 = trans_inv_operator(spin((1.0, "xx"), (1.0, "yy")), 2, Tk00) |> Array
        H_xx_within = trans_inv_operator(spin((1.0, "xx"), (1.0, "yy")), 2, Bk0) |> Array
        @test H_xx_00 ≈ H_xx_within atol = 1e-12
        @test norm(H_xx_00) > 1e-6
    end

    @testset "Momentum-sector inter-block reconstruction matches ⟨k₁|H|k₀⟩" begin
        # Build the canonical momentum eigenstates |k, j⟩ in the full Hilbert
        # space, then verify that the symmetrizer-based reconstruction
        # Pk₁ * H * Pk₀' reproduces the exact ⟨k₁,p|H|k₀,j⟩ matrix element
        # for a single-site σ_x operator that violates translation symmetry.
        L_mom = 4
        Bfull_mom = TensorBasis(L = L_mom, base = 2)
        Bk0 = TranslationalBasis(L = L_mom, k = 0, base = 2, threaded = false)
        Bk1 = TranslationalBasis(L = L_mom, k = 1, base = 2, threaded = false)

        function momentum_state(Bk, j, Lloc)
            v = zeros(ComplexF64, 2^Lloc)
            rep = Bk.I[j] - 1
            mask = (1 << Lloc) - 1
            cur, n = rep, 1
            while true
                v[cur + 1] += Bk.C[n]
                cur = ((cur >> 1) | (cur << (Lloc - 1))) & mask
                cur == rep && break
                n += 1
            end
            v ./ norm(v)
        end

        σx_full = Array(operator(ComplexF64[0 1; 1 0], [1], Bfull_mom))
        ground = zeros(ComplexF64, size(Bk1, 1), size(Bk0, 1))
        for j in 1:size(Bk0, 1), p in 1:size(Bk1, 1)
            ψj = momentum_state(Bk0, j, L_mom)
            ψp = momentum_state(Bk1, p, L_mom)
            ground[p, j] = ψp' * σx_full * ψj
        end

        Pk0 = symmetrizer(DoubleBasis(Bk0, Bfull_mom))
        Pk1 = symmetrizer(DoubleBasis(Bk1, Bfull_mom))
        recon = Pk1 * σx_full * Pk0'
        @test recon ≈ ground atol = 1e-12
    end

    @testset "DoubleBasis symmetrizer is physically correct on momentum sectors" begin
        # The earlier test `T(v) ≈ symmetrizer(T) * v` only validates that the
        # two derivations agree — it does not pin down the convention. The
        # three checks below pin physical correctness directly, so a
        # regression that conjugates basis_embedding (and T(v) in tandem)
        # would be caught here even though it passes the consistency test.
        L_orth = 6
        for k in 0:L_orth÷2
            Bk = TranslationalBasis(L = L_orth, k = k, base = 2, threaded = false)
            iszero(size(Bk, 1)) && continue
            S = symmetrizer(DoubleBasis(Bk, Bk))
            @test S ≈ I atol = 1e-10  # orthonormality within a sector
        end
        for k1 in 0:L_orth-1, k2 in (k1+1):L_orth-1
            B1 = TranslationalBasis(L = L_orth, k = k1, base = 2, threaded = false)
            B2 = TranslationalBasis(L = L_orth, k = k2, base = 2, threaded = false)
            (iszero(size(B1, 1)) || iszero(size(B2, 1))) && continue
            S = symmetrizer(DoubleBasis(B1, B2))
            @test norm(S) < 1e-10  # distinct momentum sectors are orthogonal
        end
        # Projection identity: a momentum eigenstate built in Bfull projects
        # exactly to the corresponding unit vector in Bk-coordinates.
        Bfull_orth = TensorBasis(L = L_orth, base = 2)
        for k in 0:L_orth-1
            Bk = TranslationalBasis(L = L_orth, k = k, base = 2, threaded = false)
            iszero(size(Bk, 1)) && continue
            for j in 1:min(3, size(Bk, 1))
                v = zeros(ComplexF64, 2^L_orth)
                rep = Bk.I[j] - 1
                mask = (1 << L_orth) - 1
                cur, n = rep, 1
                while true
                    v[cur + 1] += Bk.C[n]
                    cur = ((cur >> 1) | (cur << (L_orth - 1))) & mask
                    cur == rep && break
                    n += 1
                end
                v ./= norm(v)
                c = DoubleBasis(Bk, Bfull_orth)(v)
                expected = zeros(ComplexF64, size(Bk, 1))
                expected[j] = 1
                @test c ≈ expected atol = 1e-10
            end
        end
    end

    σz = Array(spin("Z"))
    dm = densitymatrix([1.0, 0.0])
    @test expectation(σz, dm) ≈ 1.0

    H = zeros(2, 2)
    jumps = [sqrt(0.2) * [0.0 1.0; 0.0 0.0]]
    lb = lindblad(H, jumps)
    dm2 = lb(dm, 0.1, order = 4)
    @test tr(dm2.ρ) ≈ 1.0 atol = 1e-8
    @test all(isreal, eigvals(Hermitian(dm2.ρ)))

    ql = quadraticlindblad(zeros(4, 4), zeros(4, 2), [zeros(4, 4)])
    cm = covariancematrix([1, 0])
    cm2 = ql(cm, 0.1, order = 3)
    @test size(cm2.Γ) == (4, 4)

    A = [1.0 0.0; 0.0 -1.0]
    B = zeros(2, 2)
    M = majoranaform(A, B)
    @test size(M) == (4, 4)

    ops = [spin("X"), spin("Z")]
    vec = [1.0, 0.0]
    cmat = covmat(ops, vec)
    @test size(cmat) == (2, 2)
    sol = qimsolve(ops, vec)
    @test size(sol, 1) == 2

    @testset "qimsolve returns all null directions when all variances are small" begin
        null_ops = [zeros(2, 2), zeros(2, 2)]
        null_sol = qimsolve(null_ops, vec; tol = 1e-7)
        @test size(null_sol) == (2, 2)
    end

    @testset "simplify handles columns that eliminate to zero" begin
        # Two linearly dependent columns: elimination zeroes the second one.
        h = Float64[1.0 1.0;
                    0.0 0.0]
        result = EDKit.simplify(copy(h))
        @test result[:, 1] ≈ [1.0, 0.0]
        @test all(iszero, result[:, 2])

        # qimsolve with linearly dependent ops should not crash.
        dep_ops = [Float64[1 0; 0 0], Float64[1 0; 0 0]]
        dep_sol = qimsolve(dep_ops, [1.0, 0.0])
        @test size(dep_sol, 1) == 2
    end
end
