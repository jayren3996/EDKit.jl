@testset "Adaptive Krylov Time Evolution" begin
    Random.seed!(11)

    # Build a small Heisenberg-like chain used in several tests below.
    function _xxz_operator(L, B)
        bond = spin((1.0, "xx"), (1.0, "yy"), (0.7, "zz"))
        trans_inv_operator(bond, 1:2, B)
    end

    @testset "Dense reference — Operator path" begin
        L = 6
        B = TensorBasis(L = L, base = 2)
        H = _xxz_operator(L, B)
        Hd = Hermitian(Array(H))
        ψ0 = randn(ComplexF64, size(B, 1)); normalize!(ψ0)

        # Single-time agreement with dense exponential.
        for t in (0.1, 0.5, 1.7, -0.0)
            ψref = exp(-1im * t * Array(Hd)) * ψ0
            ψkry = timeevolve(H, ψ0, float(t); tol = 1e-12, m_init = 25, m_max = 50)
            @test norm(ψkry - ψref) < 1e-10
            @test isapprox(norm(ψkry), 1.0; atol = 1e-10)
        end

        # Multi-time: reuse diagnostics checked separately below.
        ts = collect(range(0.0, 3.0; length = 31))
        ψs, diag = timeevolve(H, ψ0, ts;
                              tol = 1e-12, m_init = 25, m_max = 50,
                              return_diagnostics = true)
        for (i, t) in enumerate(ts)
            ψref = exp(-1im * t * Array(Hd)) * ψ0
            @test norm(ψs[:, i] - ψref) < 1e-10
        end
        @test all(isapprox(norm(ψs[:, i]), 1.0; atol = 1e-10) for i in 1:length(ts))
        @test diag.total_times_served == length(ts)
    end

    @testset "AbstractMatrix Hamiltonian path" begin
        N = 40
        A = randn(ComplexF64, N, N); A = (A + A') / 2  # Hermitian dense
        ψ0 = randn(ComplexF64, N); normalize!(ψ0)
        ψref = exp(-1im * 0.8 * A) * ψ0
        ψkry = timeevolve(A, ψ0, 0.8; tol = 1e-12, m_init = 20, m_max = 35)
        @test norm(ψkry - ψref) < 1e-10

        # Also works with Hermitian wrapper and sparse form.
        Ah = Hermitian(A)
        @test norm(timeevolve(Ah, ψ0, 0.8; tol = 1e-12, m_init = 20, m_max = 35) - ψref) < 1e-10
        As = sparse(A)
        @test norm(timeevolve(As, ψ0, 0.8; tol = 1e-12, m_init = 20, m_max = 35) - ψref) < 1e-10
    end

    @testset "Symmetry sector compatibility" begin
        L = 8
        # Half-filling k=0 sector of an XXZ chain.
        Bk = TranslationalBasis(L = L, k = 0, N = L ÷ 2, base = 2, threaded = false)
        bond = spin((1.0, "+-"), (1.0, "-+"), (0.5, "zz"))
        H = trans_inv_operator(bond, 1:2, Bk)
        Hd = Hermitian(Array(H))

        ψ0 = randn(ComplexF64, size(Bk, 1)); normalize!(ψ0)
        ts = [0.3, 0.9, 1.4]
        ψs = timeevolve(H, ψ0, ts; tol = 1e-11, m_init = 20, m_max = 40)
        for (i, t) in enumerate(ts)
            ψref = exp(-1im * t * Array(Hd)) * ψ0
            @test norm(ψs[:, i] - ψref) < 1e-9
            @test isapprox(norm(ψs[:, i]), 1.0; atol = 1e-9)
        end
    end

    @testset "Reuse actually happens" begin
        # One Lanczos basis should serve many closely-spaced output times
        # without any restarts or extensions.
        L = 6
        B = TensorBasis(L = L, base = 2)
        H = _xxz_operator(L, B)
        ψ0 = randn(ComplexF64, size(B, 1)); normalize!(ψ0)

        ts = collect(range(0.0, 0.4; length = 25))
        _, diag = timeevolve(H, ψ0, ts;
                             tol = 1e-10, m_init = 25, m_max = 50,
                             return_diagnostics = true)
        @test diag.basis_builds == 1
        @test diag.restarts == 0
        @test diag.total_times_served == length(ts)
        # The matvec budget must stay tied to the single basis build.
        @test diag.matvecs ≤ diag.max_dim_used
    end

    @testset "Restart happens on a long interval" begin
        # A small basis forced on a long interval must restart before reaching
        # the final requested time. We keep the tolerance loose enough that the
        # test mostly measures the adaptive controller rather than the
        # underlying numerics.
        L = 8
        B = TensorBasis(L = L, base = 2)
        H = _xxz_operator(L, B)
        Hd = Hermitian(Array(H))
        ψ0 = randn(ComplexF64, size(B, 1)); normalize!(ψ0)

        ts = collect(range(0.0, 20.0; length = 11))
        ψs, diag = timeevolve(H, ψ0, ts;
                              tol = 1e-10,
                              m_init = 6, m_max = 6,        # tiny basis
                              extend_basis = false,
                              return_diagnostics = true)
        for (i, t) in enumerate(ts)
            ψref = exp(-1im * t * Array(Hd)) * ψ0
            @test norm(ψs[:, i] - ψref) < 1e-6
        end
        @test diag.restarts ≥ 1
        @test diag.basis_builds ≥ 2
    end

    @testset "Extension path" begin
        # With extension allowed, a small initial basis should grow in place
        # rather than (or before) restarting.
        L = 6
        B = TensorBasis(L = L, base = 2)
        H = _xxz_operator(L, B)
        Hd = Hermitian(Array(H))
        ψ0 = randn(ComplexF64, size(B, 1)); normalize!(ψ0)

        ts = collect(range(0.0, 3.0; length = 15))
        ψs, diag = timeevolve(H, ψ0, ts;
                              tol = 1e-12,
                              m_init = 5, m_max = 40,
                              extend_step = 5,
                              extend_basis = true,
                              return_diagnostics = true)
        for (i, t) in enumerate(ts)
            ψref = exp(-1im * t * Array(Hd)) * ψ0
            @test norm(ψs[:, i] - ψref) < 1e-9
        end
        @test diag.basis_extensions ≥ 1
    end

    @testset "Explicit cache workflow and in-place variants" begin
        L = 6
        B = TensorBasis(L = L, base = 2)
        H = _xxz_operator(L, B)
        Hd = Hermitian(Array(H))
        ψ0 = randn(ComplexF64, size(B, 1)); normalize!(ψ0)

        cache = KrylovEvolutionCache(H, ψ0; tol = 1e-12, m_init = 20, m_max = 35)

        # Sequential forward calls share the cache and advance the anchor.
        ψa = timeevolve!(cache, 0.5)
        @test norm(ψa - exp(-1im * 0.5 * Array(Hd)) * ψ0) < 1e-10

        ψb = timeevolve!(cache, 1.5)
        @test norm(ψb - exp(-1im * 1.5 * Array(Hd)) * ψ0) < 1e-10

        # In-place single-time form.
        out = similar(ψ0)
        timeevolve!(out, H, ψ0, 0.7; tol = 1e-12, m_init = 20, m_max = 35)
        @test norm(out - exp(-1im * 0.7 * Array(Hd)) * ψ0) < 1e-10

        # Backwards time on the same cache errors (documented limitation).
        @test_throws ErrorException timeevolve!(cache, 0.1)

        # Zero norm input errors.
        @test_throws ErrorException KrylovEvolutionCache(H, zeros(ComplexF64, size(B, 1)))
    end

    @testset "Forward-only time rules" begin
        L = 6
        B = TensorBasis(L = L, base = 2)
        H = _xxz_operator(L, B)
        Hd = Hermitian(Array(H))
        ψ0 = randn(ComplexF64, size(B, 1)); normalize!(ψ0)

        # Negative single-time request is rejected.
        @test_throws ErrorException timeevolve(H, ψ0, -0.3; tol = 1e-12)

        # Any negative entry in ts is rejected by the stateless form.
        @test_throws ErrorException timeevolve(H, ψ0, [-0.1, 0.3]; tol = 1e-12)

        # Stateless multi-time form sorts internally and restores the caller's
        # order on the output columns.
        ts_sorted = collect(range(0.0, 2.0; length = 9))
        ψs_sorted = timeevolve(H, ψ0, ts_sorted; tol = 1e-12, m_init = 20, m_max = 40)
        perm = [5, 1, 9, 3, 7, 2, 8, 4, 6]
        ts_shuffled = ts_sorted[perm]
        ψs_shuffled = timeevolve(H, ψ0, ts_shuffled; tol = 1e-12, m_init = 20, m_max = 40)
        for (j, p) in enumerate(perm)
            @test isapprox(ψs_shuffled[:, j], ψs_sorted[:, p]; atol = 1e-10)
        end
        # And both still agree with the dense reference.
        for (j, p) in enumerate(perm)
            ψref = exp(-1im * ts_shuffled[j] * Array(Hd)) * ψ0
            @test norm(ψs_shuffled[:, j] - ψref) < 1e-9
        end

        # The stateful cache-based in-place form deliberately does *not* sort:
        # it rejects unsorted ts to avoid surprising the caller with anchor
        # moves between requests.
        cache = KrylovEvolutionCache(H, ψ0; tol = 1e-12, m_init = 20, m_max = 40)
        out = Matrix{ComplexF64}(undef, size(B, 1), length(ts_shuffled))
        @test_throws ErrorException timeevolve!(out, cache, ts_shuffled)
    end

    @testset "hermitian kwarg has been removed" begin
        # Setting `hermitian` is no longer accepted — the solver is
        # Hermitian-only by construction and the vestigial kwarg is gone.
        L = 4
        B = TensorBasis(L = L, base = 2)
        H = _xxz_operator(L, B)
        ψ0 = randn(ComplexF64, size(B, 1)); normalize!(ψ0)
        @test_throws MethodError KrylovEvolutionCache(H, ψ0; hermitian = true)
        @test_throws MethodError timeevolve(H, ψ0, 0.1; hermitian = false)
    end

    @testset "normalize_output default does not mask error" begin
        # With a deliberately under-sized Krylov basis and a very loose tol,
        # the reconstructed state should be *inaccurate*. Default
        # normalize_output=false exposes that inaccuracy; opting in to
        # normalize_output=true silently hides norm drift and is therefore
        # NOT equivalent to a correct solve.
        L = 6
        B = TensorBasis(L = L, base = 2)
        H = _xxz_operator(L, B)
        Hd = Hermitian(Array(H))
        ψ0 = randn(ComplexF64, size(B, 1)); normalize!(ψ0)

        # Default: norm should still be ~1 (Lanczos is unitary) and the
        # answer should match the dense reference for a reasonable setup.
        ψt = timeevolve(H, ψ0, 1.0; tol = 1e-12, m_init = 20, m_max = 40)
        @test abs(norm(ψt) - 1.0) < 1e-10
        @test norm(ψt - exp(-1im * 1.0 * Array(Hd)) * ψ0) < 1e-10

        # Opting in to renormalization still produces a unit-norm output but
        # must not be trusted for accuracy control.
        ψt_n = timeevolve(H, ψ0, 1.0; tol = 1e-12, m_init = 20, m_max = 40,
                          normalize_output = true)
        @test isapprox(norm(ψt_n), 1.0; atol = 1e-12)
    end

    @testset "Defect monitor stress: wide spectral spread" begin
        # A Hermitian operator with a deliberately wide spectrum and a long
        # propagation horizon: the defect monitor is a fast-oscillating
        # trigonometric sum whose frequency content scales with the spectral
        # spread. A naive fixed 9-sample check could miss peaks between
        # samples on long intervals, so the adaptive sample-count selector
        # has to kick in.
        Random.seed!(42)
        N = 60
        U = qr(randn(ComplexF64, N, N)).Q |> Matrix
        # Wide spectrum from -8 to +8 with randomly placed eigenvalues.
        λ = sort(16 .* rand(N) .- 8)
        H = U * Diagonal(λ) * U'
        H = (H + H') / 2
        ψ0 = randn(ComplexF64, N); normalize!(ψ0)

        # Large τ means many beat cycles of the monitor on each candidate
        # interval — this is the regime where adaptive densification matters.
        ts = collect(range(0.0, 12.0; length = 25))
        ψs = timeevolve(H, ψ0, ts; tol = 1e-10, m_init = 25, m_max = 50)
        for (i, t) in enumerate(ts)
            ψref = exp(-1im * t * H) * ψ0
            @test norm(ψs[:, i] - ψref) < 1e-8
            @test abs(norm(ψs[:, i]) - 1.0) < 1e-8
        end

        # Same test with a single shot at t = 12: one cache, possibly many
        # restarts, but the defect controller must still hit every requested
        # output time accurately without over-confident acceptance.
        ψfinal = timeevolve(H, ψ0, 12.0; tol = 1e-10, m_init = 20, m_max = 40)
        @test norm(ψfinal - exp(-1im * 12.0 * H) * ψ0) < 1e-8
    end
end
