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

@testset "te-4: backward / two-sided evolution" begin
    Random.seed!(170)
    L = 8
    H = trans_inv_operator(spin((1.0, "xx"), (1.0, "yy"), (0.7, "zz")), 1:2, TensorBasis(L=L, base=2))
    N = 2^L
    Hm = Matrix(H)
    ψ0 = randn(ComplexF64, N); ψ0 ./= norm(ψ0)

    # Backward evolution to -t reproduces exp(+i t H) ψ0.
    for t in (0.5, 1.3, 2.0)
        ψ_back = timeevolve(H, ψ0, -t; tol=1e-11)
        @test ψ_back ≈ exp(+im * t * Hm) * ψ0 rtol=1e-7
    end

    # Forward then backward returns the original state (unitary round trip).
    ψf = timeevolve(H, ψ0, 2.5; tol=1e-11)
    ψr = timeevolve(H, ψf, -2.5; tol=1e-11)
    @test ψr ≈ ψ0 rtol=1e-7

    # Within one cache: evolve forward, then step backward to an earlier time.
    cache = KrylovEvolutionCache(H, ψ0; tol=1e-11)
    ψa = timeevolve!(cache, 1.0)
    @test ψa ≈ exp(-im * 1.0 * Hm) * ψ0 rtol=1e-7
    cache2 = KrylovEvolutionCache(H, ψ0; tol=1e-11)
    ψb = timeevolve!(cache2, -1.0)
    @test ψb ≈ exp(+im * 1.0 * Hm) * ψ0 rtol=1e-7

    # Multi-time backward: a vector of negative times, returned in input order.
    ts = [-0.4, -1.2, -0.8]
    ψs = timeevolve(H, ψ0, ts; tol=1e-11)
    for (k, t) in enumerate(ts)
        @test ψs[:, k] ≈ exp(-im * t * Hm) * ψ0 rtol=1e-7
    end

    # Mixing positive and negative target times in one call is rejected.
    @test_throws Exception timeevolve(H, ψ0, [-0.5, 0.5]; tol=1e-11)

    # OTOC-style two-sided round trip: e^{+iHt} W e^{-iHt} ψ with W diagonal.
    W = Diagonal(cis.(randn(N)))           # unitary "operator insertion"
    ψ_fwd = timeevolve(H, ψ0, 1.5; tol=1e-11)
    ψ_ins = W * ψ_fwd
    ψ_otoc = timeevolve(H, ψ_ins, -1.5; tol=1e-11)
    @test ψ_otoc ≈ exp(+im * 1.5 * Hm) * (W * (exp(-im * 1.5 * Hm) * ψ0)) rtol=1e-7
end

@testset "e1: spinful / multi-species fermions" begin
    # Sector dimensions.
    @test size(SpinfulFermionBasis(L=3, S=2, N=(2, 1)), 1) == binomial(3, 2) * binomial(3, 1)  # 9
    @test size(SpinfulFermionBasis(L=4, S=2, N=(2, 2)), 1) == binomial(4, 2)^2                  # 36
    @test size(SpinfulFermionBasis(L=3, S=2, N=3), 1) == binomial(6, 3)                         # total N=3 over 6 modes
    @test size(SpinfulFermionBasis(L=2, S=2), 1) == 2^4                                         # all sectors, 4 modes
    @test size(SpinfulFermionBasis(L=3, S=3, N=(1, 1, 1)), 1) == binomial(3, 1)^3               # 27

    # Mode map (blocked): mode(i,σ) = (σ-1)L + i.
    B = SpinfulFermionBasis(L=4, S=2)
    @test fermionmode(B, 1, 1) == 1
    @test fermionmode(B, 4, 1) == 4
    @test fermionmode(B, 1, 2) == 5
    @test fermionmode(B, 3, 2) == 7
    @test fermionmode(B, 2, :↑) == 2
    @test fermionmode(B, 2, :↓) == 6
    @test fermionmode(B, 2, :up) == fermionmode(B, 2, :↑)
    @test fermionmode(B, 2, :down) == fermionmode(B, 2, :↓)

    # Canonical anticommutation {c_p, c†_q} = δ_pq I on the full Fock space.
    Bfull = SpinfulFermionBasis(L=2, S=2)               # 4 modes, dim 16
    d = size(Bfull, 1)
    cs = [Array(fermion_operator("-", [p], Bfull)) for p in 1:4]   # annihilation
    cd = [Array(fermion_operator("+", [q], Bfull)) for q in 1:4]   # creation
    for p in 1:4, q in 1:4
        anti = cs[p] * cd[q] + cd[q] * cs[p]
        @test anti ≈ (p == q ? Matrix(I, d, d) : zeros(d, d))
    end

    # Spin-aware construction agrees with raw-mode construction.
    Bs = SpinfulFermionBasis(L=3, S=2, N=(2, 1))
    Braw = SpinlessFermionBasis(L=6, N=3, f = dgt -> sum(dgt[1:3]) == 2 && sum(dgt[4:6]) == 1)
    @test Bs.I == Braw.I
    op_spin = fermion_operator("+-", [(1, :↑), (2, :↑)], Bs)
    op_raw  = fermion_operator("+-", [1, 2], Bs)        # modes for (1,↑),(2,↑)
    @test Array(op_spin) ≈ Array(op_raw)

    # Number conservation: [H, N̂] = 0 on the full Fock space.
    Bf = SpinfulFermionBasis(L=2, S=2)
    Hf = hubbard(Bf; t=1.0, U=3.0, μ=0.7)
    Nop = reduce(+, [fermion_operator("n", [m], Bf) for m in 1:4])
    Hfm, Nfm = Array(Hf), Array(Nop)
    @test Hfm ≈ Hfm'
    @test Hfm * Nfm ≈ Nfm * Hfm

    # 2-site half-filled Hubbard: exact ground state E0 = (U - √(U²+16t²))/2.
    for (t, U) in [(1.0, 0.0), (1.0, 4.0), (0.7, 2.5)]
        B2 = SpinfulFermionBasis(L=2, S=2, N=(1, 1))
        H2 = hubbard(B2; t=t, U=U, boundary=:open)
        E0 = minimum(eigvals(Hermitian(Array(H2))))
        @test E0 ≈ (U - sqrt(U^2 + 16t^2)) / 2
    end

    # Errors.
    @test_throws Exception SpinfulFermionBasis(L=3, S=2, N=(4, 1))           # N↑ > L
    @test_throws Exception fermionmode(SpinfulFermionBasis(L=2, S=1), 1, :↓)  # :↓ needs S ≥ 2
    @test_throws Exception hubbard(SpinfulFermionBasis(L=3, S=3))            # hubbard needs S=2
end

@testset "e2: momentum-resolved spinless fermions" begin
    # Translation-invariant fermion Hamiltonian: hopping + nn interaction.
    fermH(B, V) = begin
        hop = trans_inv_fermion_operator("+-", [1, 2], B)
        -(Array(hop) + adjoint(Array(hop))) + V * Array(trans_inv_fermion_operator("nn", [1, 2], B))
    end

    for (L, N, V) in [(6, 3, 0.0), (6, 2, 0.8), (7, 3, 1.1), (8, 4, 0.5), (6, 4, 1.3)]
        # Full N-sector reference.
        Bfull = SpinlessFermionBasis(L=L, N=N)
        ref = sort(real(eigvals(Hermitian(fermH(Bfull, V)))))

        # Union of all momentum sectors.
        kvals = Float64[]
        dimsum = 0
        for k in 0:L-1
            Bk = TranslationalFermionBasis(L=L, N=N, k=k)
            dimsum += size(Bk, 1)
            size(Bk, 1) == 0 && continue
            append!(kvals, real(eigvals(Hermitian(Matrix(fermH(Bk, V))))))
        end
        @test dimsum == size(Bfull, 1)                  # sector dims partition the N-sector
        @test sort(kvals) ≈ ref                          # union of k-spectra = full spectrum
    end

    # eltype: real-phase sectors (k=0 and, for even L, k=L/2) are Float64;
    # generic momenta are ComplexF64.
    @test eltype(TranslationalFermionBasis(L=6, N=3, k=0)) <: Real
    @test eltype(TranslationalFermionBasis(L=6, N=3, k=3)) <: Real
    @test eltype(TranslationalFermionBasis(L=6, N=3, k=1)) <: Complex

    # JW operators are allowed on TranslationalFermionBasis but still rejected on
    # other permutation bases.
    Bk = TranslationalFermionBasis(L=6, N=3, k=0)
    @test trans_inv_fermion_operator("+-", [1, 2], Bk) isa EDKit.Operator
    @test_throws Exception fermion_operator("+", [1], TranslationalBasis(L=6, k=0, N=3))

    # base != 2 is rejected.
    @test_throws Exception TranslationalFermionBasis(L=4, N=2, base=3)

    # copy and order.
    @test copy(Bk).I == Bk.I
    @test EDKit.order(Bk) == 6
end
