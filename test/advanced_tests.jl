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
end
