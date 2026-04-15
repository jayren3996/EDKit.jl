@testset "Adaptive Lindblad Arnoldi Evolution" begin
    Random.seed!(17)

    function _explicit_liouvillian(H, jumps)
        T = promote_type(ComplexF64, eltype(H), map(eltype, jumps)...)
        d = size(H, 1)
        I_d = Matrix{T}(I, d, d)
        D = zeros(T, d, d)
        for L in jumps
            D .+= adjoint(L) * L
        end

        Lmat = -1im * (kron(I_d, H) - kron(transpose(H), I_d))
        for L in jumps
            Lmat .+= kron(conj(L), L)
        end
        Lmat .+= -0.5 * kron(I_d, D)
        Lmat .+= -0.5 * kron(transpose(D), I_d)
        Lmat
    end

    function _explicit_evolve(H, jumps, dm, t)
        d = size(dm.ρ, 1)
        ρvec = exp(t * _explicit_liouvillian(H, jumps)) * vec(dm.ρ)
        reshape(ρvec, d, d)
    end

    function _unitary_density_evolve(H, dm, t)
        U = exp(-1im * t * H)
        U * dm.ρ * U'
    end

    function _random_hermitian(rng, d)
        A = randn(rng, ComplexF64, d, d)
        (A + A') / 2
    end

    @testset "Single-time dense Lindblad path matches explicit Liouvillian" begin
        H = ComplexF64[0.2 0.1-0.05im; 0.1+0.05im -0.3]
        jumps = [
            sqrt(0.4) * ComplexF64[0 1; 0 0],
            sqrt(0.15) * ComplexF64[0 0; 1 0],
        ]
        lb = lindblad(H, jumps)
        ρ0 = densitymatrix(ComplexF64[0.3 + 0.2im, 0.9])
        t = 0.7

        ρt, diag = lindblad_timeevolve(lb, ρ0, t;
                                       tol = 1e-12,
                                       m_init = 4,
                                       m_max = 12,
                                       return_diagnostics = true,
                                       record_min_eig = true)

        ρref = _explicit_evolve(H, jumps, ρ0, t)
        @test norm(ρt.ρ - ρref) < 1e-10
        @test diag.total_times_served == 1
        @test length(diag.trace_drifts) == 1
        @test length(diag.hermiticity_drifts) == 1
        @test length(diag.min_hermitian_eigs) == 1
        @test diag.max_trace_drift ≈ diag.trace_drifts[1]
        @test diag.max_hermiticity_drift ≈ diag.hermiticity_drifts[1]
        @test diag.min_min_hermitian_eig ≈ diag.min_hermitian_eigs[1]
    end

    @testset "Multi-time reuse and diagnostics" begin
        H = ComplexF64[
            0.2 0.1 0.0 0.0;
            0.1 -0.3 0.15 0.0;
            0.0 0.15 0.4 0.12;
            0.0 0.0 0.12 -0.1
        ]
        jumps = [
            sqrt(0.2) * ComplexF64[
                0 1 0 0;
                0 0 0 0;
                0 0 0 1;
                0 0 0 0
            ],
        ]
        lb = lindblad(H, jumps)
        ψ0 = randn(ComplexF64, size(H, 1))
        normalize!(ψ0)
        ρ0 = densitymatrix(ψ0)
        ts = collect(range(0.0, 0.35; length = 8))

        ρs, diag = lindblad_timeevolve(lb, ρ0, ts;
                                       tol = 1e-11,
                                       m_init = 5,
                                       m_max = 14,
                                       return_diagnostics = true,
                                       record_min_eig = true)

        @test length(ρs) == length(ts)
        @test diag.total_times_served == length(ts)
        @test diag.basis_builds == 1
        @test diag.restarts == 0
        @test length(diag.trace_drifts) == length(ts)
        @test length(diag.hermiticity_drifts) == length(ts)
        @test length(diag.min_hermitian_eigs) == length(ts)
        @test diag.max_trace_drift == maximum(diag.trace_drifts)
        @test diag.max_hermiticity_drift == maximum(diag.hermiticity_drifts)
        @test diag.min_min_hermitian_eig == minimum(diag.min_hermitian_eigs)

        for (ρt, t) in zip(ρs, ts)
            ρref = _explicit_evolve(H, jumps, ρ0, t)
            @test norm(ρt.ρ - ρref) < 1e-9
        end
    end

    @testset "Sparse LiouvillianMap cache workflow" begin
        H = sparse(ComplexF64[0.0 0.2; 0.2 0.1])
        jumps = sparse.([sqrt(0.3) * ComplexF64[0 1; 0 0]])
        A = LiouvillianMap(H, jumps)
        ρ0 = densitymatrix(ComplexF64[0.0, 1.0])
        cache = LindbladArnoldiCache(A, ρ0;
                                     tol = 1e-12,
                                     m_init = 4,
                                     m_max = 10,
                                     record_min_eig = true)
        t = 0.5
        ρt = lindblad_timeevolve!(cache, t)
        ρref = _explicit_evolve(Matrix(H), Matrix.(jumps), ρ0, t)

        @test norm(ρt.ρ - ρref) < 1e-10
        @test cache.diagnostics.total_times_served == 1
        @test cache.diagnostics.liouvillian_applies ≥ cache.m
        @test length(cache.diagnostics.trace_drifts) == 1
    end

    @testset "Forward-only time rules" begin
        H = ComplexF64[0 1; 1 0]
        jumps = [sqrt(0.2) * ComplexF64[0 1; 0 0]]
        lb = lindblad(H, jumps)
        ρ0 = densitymatrix(ComplexF64[0.0, 1.0])

        @test_throws ErrorException lindblad_timeevolve(lb, ρ0, -0.1)
        @test_throws ErrorException lindblad_timeevolve(lb, ρ0, [-0.1, 0.2])

        cache = LindbladArnoldiCache(lb, ρ0; tol = 1e-10, m_init = 4, m_max = 10)
        lindblad_timeevolve!(cache, 0.2)
        @test_throws ErrorException lindblad_timeevolve!(cache, 0.1)
    end

    @testset "Pure Hamiltonian limit matches unitary density evolution" begin
        rng = MersenneTwister(23)
        H = _random_hermitian(rng, 3)
        lb = lindblad(H, Matrix{ComplexF64}[])
        ψ0 = randn(rng, ComplexF64, 3)
        normalize!(ψ0)
        ρ0 = densitymatrix(ψ0)

        t = 1.4
        ρt = lindblad_timeevolve(lb, ρ0, t; tol = 1e-12, m_init = 4, m_max = 9)
        ρref = _unitary_density_evolve(H, ρ0, t)
        ψt = timeevolve(H, ψ0, t; tol = 1e-12, m_init = 6, m_max = 12)
        ρψ = densitymatrix(ψt).ρ
        @test norm(ρt.ρ - ρref) < 1e-10
        @test norm(ρt.ρ - ρψ) < 1e-10

        ts = [0.0, 0.2, 0.7, 1.5]
        ρs, diag = lindblad_timeevolve(lb, ρ0, ts;
                                       tol = 1e-12,
                                       m_init = 4,
                                       m_max = 9,
                                       return_diagnostics = true)
        @test diag.total_times_served == length(ts)
        for (ρτ, τ) in zip(ρs, ts)
            @test norm(ρτ.ρ - _unitary_density_evolve(H, ρ0, τ)) < 1e-10
        end

        cache = LindbladArnoldiCache(lb, ρ0;
                                     tol = 1e-6,
                                     m_init = 3,
                                     m_max = 3,
                                     extend_basis = false)
        ρrestart = lindblad_timeevolve!(cache, 5.0)
        @test norm(ρrestart.ρ - _unitary_density_evolve(H, ρ0, 5.0)) < 5e-6
        @test cache.diagnostics.restarts ≥ 1
    end

    @testset "Generator preserves trace and Hermiticity" begin
        rng = MersenneTwister(91)
        H = _random_hermitian(rng, 3)
        jumps = [
            randn(rng, ComplexF64, 3, 3) / 5,
            randn(rng, ComplexF64, 3, 3) / 7,
        ]
        ρh = _random_hermitian(rng, 3)

        for A in (
            LiouvillianMap(H, jumps),
            LiouvillianMap(sparse(H), sparse.(jumps)),
        )
            Lρ = A * ρh
            @test abs(tr(Lρ)) < 1e-12
            @test norm(Lρ - Lρ') < 1e-12

            for _ in 1:3
                ψ = randn(rng, ComplexF64, 3)
                normalize!(ψ)
                ρ = densitymatrix(ψ).ρ
                @test abs(tr(A * ρ)) < 1e-12
            end
        end
    end

    @testset "Semigroup consistency for stateless and stateful evolution" begin
        rng = MersenneTwister(202)
        H = _random_hermitian(rng, 3)
        jumps = [
            ComplexF64[0 1 0; 0 0 0; 0 0 0] * sqrt(0.2),
            ComplexF64[0 0 0; 0 0 1; 0 0 0] * sqrt(0.15),
        ]
        lb = lindblad(H, jumps)
        ρ0 = densitymatrix(randn(rng, ComplexF64, 3))

        t1, t2 = 0.25, 0.6
        ρ_one = lindblad_timeevolve(lb, ρ0, t1 + t2; tol = 1e-11, m_init = 4, m_max = 9)
        ρ_seq = lindblad_timeevolve(lb,
                                    lindblad_timeevolve(lb, ρ0, t1; tol = 1e-11, m_init = 4, m_max = 9),
                                    t2;
                                    tol = 1e-11, m_init = 4, m_max = 9)
        @test norm(ρ_one.ρ - ρ_seq.ρ) < 1e-9

        cache = LindbladArnoldiCache(lb, ρ0; tol = 1e-11, m_init = 4, m_max = 9)
        lindblad_timeevolve!(cache, t1)
        ρ_cache = lindblad_timeevolve!(cache, t1 + t2)
        @test norm(ρ_one.ρ - ρ_cache.ρ) < 1e-9

        total = 2.0
        cache_restart = LindbladArnoldiCache(lb, ρ0;
                                             tol = 1e-6,
                                             m_init = 5,
                                             m_max = 5,
                                             extend_basis = false)
        ρ_stateful = lindblad_timeevolve!(cache_restart, 0.6)
        ρ_stateful = lindblad_timeevolve!(cache_restart, total)
        ρ_split = lindblad_timeevolve(lb,
                                      lindblad_timeevolve(lb, ρ0, 0.6;
                                                          tol = 1e-6,
                                                          m_init = 5,
                                                          m_max = 5,
                                                          extend_basis = false),
                                      total - 0.6;
                                      tol = 1e-6,
                                      m_init = 5,
                                      m_max = 5,
                                      extend_basis = false)
        ρ_total = lindblad_timeevolve(lb, ρ0, total;
                                      tol = 1e-6,
                                      m_init = 5,
                                      m_max = 5,
                                      extend_basis = false)
        @test norm(ρ_stateful.ρ - ρ_split.ρ) < 1e-6
        @test norm(ρ_stateful.ρ - ρ_total.ρ) < 1e-6
        @test cache_restart.diagnostics.restarts ≥ 1
    end

    @testset "Known steady state remains fixed" begin
        γ = 0.4
        jump = sqrt(γ) * ComplexF64[0 1; 0 0]
        lb = lindblad(zeros(ComplexF64, 2, 2), [jump])
        ρss = densitymatrix(ComplexF64[1.0, 0.0])
        ts = collect(range(0.0, 3.0; length = 9))

        ρs, diag = lindblad_timeevolve(lb, ρss, ts;
                                       tol = 1e-12,
                                       m_init = 2,
                                       m_max = 4,
                                       return_diagnostics = true,
                                       record_min_eig = true)

        for ρ in ρs
            @test norm(ρ.ρ - ρss.ρ) < 1e-12
        end
        @test maximum(diag.trace_drifts) < 1e-12
        @test maximum(diag.hermiticity_drifts) < 1e-12
        @test diag.min_min_hermitian_eig ≥ -1e-12
    end

    @testset "Dense and sparse backend parity" begin
        rng = MersenneTwister(44)
        H = _random_hermitian(rng, 4)
        jumps = [
            randn(rng, ComplexF64, 4, 4) / 6,
            randn(rng, ComplexF64, 4, 4) / 8,
        ]
        lb = lindblad(H, jumps)
        Ad = LiouvillianMap(H, jumps)
        As = LiouvillianMap(sparse(H), sparse.(jumps))
        ρ0 = densitymatrix(randn(rng, ComplexF64, 4))
        ts = [0.0, 0.15, 0.45, 0.8]

        ρ_lb, d_lb = lindblad_timeevolve(lb, ρ0, ts;
                                         tol = 1e-11,
                                         m_init = 5,
                                         m_max = 12,
                                         return_diagnostics = true,
                                         record_min_eig = true)
        ρ_d, d_d = lindblad_timeevolve(Ad, ρ0, ts;
                                       tol = 1e-11,
                                       m_init = 5,
                                       m_max = 12,
                                       return_diagnostics = true,
                                       record_min_eig = true)
        ρ_s, d_s = lindblad_timeevolve(As, ρ0, ts;
                                       tol = 1e-11,
                                       m_init = 5,
                                       m_max = 12,
                                       return_diagnostics = true,
                                       record_min_eig = true)

        for i in eachindex(ts)
            @test norm(ρ_lb[i].ρ - ρ_d[i].ρ) < 1e-10
            @test norm(ρ_lb[i].ρ - ρ_s[i].ρ) < 1e-10
        end
        @test d_lb.trace_drifts ≈ d_d.trace_drifts atol = 1e-12
        @test d_lb.trace_drifts ≈ d_s.trace_drifts atol = 1e-12
        @test d_lb.hermiticity_drifts ≈ d_d.hermiticity_drifts atol = 1e-12
        @test d_lb.hermiticity_drifts ≈ d_s.hermiticity_drifts atol = 1e-12
        @test d_lb.min_hermitian_eigs ≈ d_d.min_hermitian_eigs atol = 1e-12
        @test d_lb.min_hermitian_eigs ≈ d_s.min_hermitian_eigs atol = 1e-12
    end

    @testset "Tightening controls improves accuracy on a nontrivial model" begin
        rng = MersenneTwister(66)
        H = _random_hermitian(rng, 4)
        jumps = [
            ComplexF64[
                0 1 0 0;
                0 0 0 0;
                0 0 0 1;
                0 0 0 0
            ] * sqrt(0.2),
            randn(rng, ComplexF64, 4, 4) / 10,
        ]
        lb = lindblad(H, jumps)
        ρ0 = densitymatrix(randn(rng, ComplexF64, 4))
        ts = [0.4, 1.0, 1.7]

        ρ_lo = lindblad_timeevolve(lb, ρ0, ts;
                                   tol = 1e-2,
                                   m_init = 2,
                                   m_max = 2,
                                   extend_basis = false)
        ρ_hi = lindblad_timeevolve(lb, ρ0, ts;
                                   tol = 1e-12,
                                   m_init = 6,
                                   m_max = 12,
                                   extend_basis = true)

        err_lo = maximum(norm(ρ.ρ - _explicit_evolve(H, jumps, ρ0, t)) for (ρ, t) in zip(ρ_lo, ts))
        err_hi = maximum(norm(ρ.ρ - _explicit_evolve(H, jumps, ρ0, t)) for (ρ, t) in zip(ρ_hi, ts))
        @test err_hi < err_lo
        @test err_hi < 1e-9
    end
end
