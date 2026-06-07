# Phase 7 — Extensions regression/feature tests.
#
# These exercise the additive analysis/dynamics features (RDM + mutual
# information, steady-state Lindblad, backward/two-sided evolution) and the
# basis extensions (spinful fermions, mixed local dimensions). Each testset is
# self-contained and references a brute-force or analytic value, never the
# implementation under test.

using EDKit, Test, LinearAlgebra, SparseArrays, Random

# Independent reduced-density-matrix reference: build the amplitude tensor
# ψ[a, b] by direct assignment over full-space states (uses only change!/index,
# not schmidt), then contract ρ_A = ψ ψ†. Ordering matches rdm's subsystem
# bases so the matrices are directly comparable.
function _rdm_bruteforce(v, Ainds, L, base)
    Ainds = sort(collect(Ainds))
    Binds = [i for i in 1:L if !(i in Ainds)]
    B1 = TensorBasis(L=length(Ainds), base=base)
    B2 = TensorBasis(L=length(Binds), base=base)
    bfull = TensorBasis(L=L, base=base)
    dgt = zeros(Int, L)
    ψ = zeros(eltype(v), size(B1, 1), size(B2, 1))
    for i in 1:base^L
        EDKit.change!(bfull, i, dgt)
        ia = EDKit.index(B1, collect(dgt[Ainds]))[2]
        ib = EDKit.index(B2, collect(dgt[Binds]))[2]
        ψ[ia, ib] = v[i]
    end
    ψ * ψ'
end

# von Neumann entropy directly from a density matrix (independent of ent_S).
function _vn_entropy(ρ; cutoff=1e-12)
    λ = eigvals(Hermitian(Matrix(ρ)))
    -sum(x * log(x) for x in λ if x > cutoff)
end

@testset "schmidt-4: reduced density matrix" begin
    Random.seed!(2718)
    L, base = 6, 2
    b = TensorBasis(L=L, base=base)
    v = randn(ComplexF64, base^L); v ./= norm(v)

    for Ainds in ([1, 2, 3], [1, 2], [2, 4], [1, 3, 5], [4, 5, 6])
        ρ = rdm(v, Ainds, b)
        @test size(ρ) == (base^length(Ainds), base^length(Ainds))
        @test ρ ≈ ρ'                                  # Hermitian
        @test tr(ρ) ≈ 1                               # normalized
        @test minimum(real(eigvals(Hermitian(ρ)))) > -1e-10   # PSD
        @test ρ ≈ _rdm_bruteforce(v, Ainds, L, base)  # matches brute force
        # eigenvalues are squared Schmidt singular values
        @test sort(real(eigvals(Hermitian(ρ)))) ≈ sort(EDKit.ent_spec(v, Ainds, b) .^ 2)
        # von Neumann entropy from ρ matches ent_S
        @test _vn_entropy(ρ) ≈ ent_S(v, Ainds, b)
    end

    # Product state: ρ_A is pure (rank 1).
    up = [1.0, 0.0]; dn = [0.0, 1.0]
    vp = kron(up, dn, up, dn, up, dn)             # site 1 = up, ...
    ρp = rdm(vp, [1, 2], b)
    @test rank(ρp; rtol=1e-10) == 1
    @test _vn_entropy(ρp) < 1e-10

    # Bell pair on L=2: ρ_A = I/2.
    b2 = TensorBasis(L=2, base=2)
    bell = (kron(up, up) + kron(dn, dn)) / sqrt(2)
    @test rdm(bell, [1], b2) ≈ [0.5 0.0; 0.0 0.5]

    # Real input stays real.
    vr = randn(base^L); vr ./= norm(vr)
    @test eltype(rdm(vr, [1, 2, 3], b)) <: Real

    # Convenience overload inferring the basis from length.
    @test rdm(v, [1, 2, 3], L) ≈ rdm(v, [1, 2, 3], b)
end

@testset "schmidt-4: mutual information" begin
    Random.seed!(31415)
    L, base = 6, 2
    b = TensorBasis(L=L, base=base)
    v = randn(ComplexF64, base^L); v ./= norm(v)

    A = [1, 2]; C = [4, 5]
    I_AC = mutual_information(v, A, C, b)
    @test I_AC ≥ -1e-10                              # non-negative
    @test I_AC ≈ mutual_information(v, C, A, b)       # symmetric
    # I(A:C) = S(A) + S(C) - S(A∪C)
    @test I_AC ≈ ent_S(v, A, b) + ent_S(v, C, b) - ent_S(v, sort(vcat(A, C)), b)

    # Complementary bipartition of a pure state: I(A:Ā) = 2 S(A).
    Ā = [3, 4, 5, 6]
    @test mutual_information(v, [1, 2], Ā, b) ≈ 2 * ent_S(v, [1, 2], b)

    # Product state: zero mutual information between any disjoint regions.
    up = [1.0, 0.0]; dn = [0.0, 1.0]
    vp = kron(up, dn, up, dn, up, dn)
    @test abs(mutual_information(vp, [1, 2], [4, 5], b)) < 1e-10

    # Overlapping subsystems error.
    @test_throws Exception mutual_information(v, [1, 2, 3], [3, 4], b)

    # Convenience overload.
    @test mutual_information(v, A, C, L) ≈ I_AC
end

@testset "lindblad-6: steady state" begin
    # Single-qubit amplitude damping: H = 0, L = σ⁻ ⇒ ρ_ss = |0⟩⟨0|.
    σm = ComplexF64[0 1; 0 0]
    A1 = LiouvillianMap(zeros(ComplexF64, 2, 2), [σm])
    ρ1 = steadystate(A1)
    @test Array(ρ1) ≈ [1 0; 0 0]
    @test tr(ρ1.ρ) ≈ 1
    @test norm(ρ1.ρ - ρ1.ρ') < 1e-9                # Hermitian
    @test norm(A1 * ρ1.ρ) < 1e-9                    # 𝓛[ρ_ss] = 0

    # Two-qubit driven-dissipative system with a unique steady state.
    σx = ComplexF64[0 1; 1 0]; I2 = ComplexF64[1 0; 0 1]
    H = kron(σx, I2) + kron(I2, σx)
    jumps = [kron(σm, I2), kron(I2, σm)]
    A = LiouvillianMap(H, jumps)

    ρdense  = steadystate(A; method=:dense)
    ρkrylov = steadystate(A; method=:krylov)
    for ρ in (ρdense, ρkrylov)
        M = ρ.ρ
        @test tr(M) ≈ 1
        @test norm(M - M') < 1e-7                            # Hermitian
        @test minimum(real(eigvals(Hermitian(M)))) > -1e-7   # PSD
        @test norm(A * M) < 1e-6                             # 𝓛[ρ_ss] = 0
    end
    @test Array(ρdense) ≈ Array(ρkrylov)                      # methods agree

    # Cross-check against long-time evolution with the explicit propagator.
    lb = lindblad(H, jumps)
    dm = densitymatrix(Matrix{ComplexF64}(I, 4, 4) / 4)       # maximally mixed start
    for _ in 1:4000
        dm = lb(dm, 0.02)
    end
    normalize!(dm)
    @test Matrix(Array(ρdense)) ≈ dm.ρ atol=1e-4

    # Overloads: Lindblad and (H, jumps).
    @test Array(steadystate(lb)) ≈ Array(ρdense)
    @test Array(steadystate(H, jumps)) ≈ Array(ρdense)

    # Unknown method errors.
    @test_throws ArgumentError steadystate(A; method=:bogus)
end
