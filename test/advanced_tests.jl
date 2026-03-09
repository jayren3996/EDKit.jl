@testset "Linear Maps And Algorithms" begin
    L = 4
    N = 2

    Bfull = TensorBasis(L = L, base = 2)
    T = DoubleBasis(Bfull, Bfull)
    ψfull = randn(ComplexF64, size(Bfull, 1))
    @test T(ψfull) ≈ ψfull

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
