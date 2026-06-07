# Regression tests for the 0.6.0 review-remediation roadmap.
# Self-runnable (`julia --project test/remediation_tests.jl`) and included by runtests.jl.
using EDKit, Test, LinearAlgebra

@testset "abelian-1: fixed N + inversion off half-filling errors" begin
    # Global spin-flip maps the N sector to the L-N sector. Off half-filling the
    # requested symmetry sector is empty/ill-defined and must error, not silently
    # return an incomplete basis.
    @test_throws ErrorException basis(L=6, N=2, z=1)
    @test_throws ErrorException basis(L=6, N=2, z=-1)
    # A custom inversion generator (flip all sites) off half-filling also errors.
    @test_throws ErrorException basis(L=6, N=2, symmetries=[(collect(1:6), 0, trues(6))])
    # Half-filling remains valid and builds a non-empty basis.
    B = basis(L=6, N=3, z=1)
    @test size(B, 1) > 0
end

@testset "productstate: correct on reduced bases, errors out of sector" begin
    L = 4
    v = [0, 1, 0, 1]

    # Reduced (Abelian) basis: must keep the orbit coefficient and place it at
    # the right component, with the basis element type.
    B = basis(; L, k=0)
    s = productstate(v, B)
    @test length(s) == size(B, 1)
    Texp = eltype(B) <: Integer ? Float64 : eltype(B)
    @test eltype(s) == Texp
    c, I = EDKit.index(B, collect(v))
    @test count(!iszero, s) == 1
    @test s[I] == c

    # Out-of-sector config must error, not silently write onto vector 1.
    Bp = basis(; L, N=2)                 # sum(dgt)==2 sector
    @test_throws ErrorException productstate([1, 1, 1, 1], Bp)  # sum==4, not in sector

    # Onsite basis still returns a Float64 unit vector (no behavior regression).
    Bt = TensorBasis(L=L, base=2)
    st = productstate(v, Bt)
    @test eltype(st) == Float64
    @test sum(st) == 1.0
    @test count(!iszero, st) == 1
end

# ---------------------------------------------------------------------------
# Phase 2 — integer-width sweep
# ---------------------------------------------------------------------------

@testset "bug-1/2/3: index capacity guard errors instead of silent overflow" begin
    # bug-1: Int32 + base>2 overflows the index type -> must error, not store wrapped reps.
    @test_throws ErrorException ProjectedBasis(Int32; L=20, N=0, base=3)
    @test_throws ErrorException basis(Int32; L=20, N=0, base=3)
    # bug-2: base=2 at L=63 overflowed the Gosper enumerator to an empty basis -> must error.
    @test_throws ErrorException ProjectedBasis(L=63, N=1, base=2)
    # bug-3: TensorBasis size overflowed to a negative value at L>=63 -> must error.
    @test_throws ErrorException TensorBasis(L=63, base=2)
    # Boundary positive: the largest valid base-2 Int64 size (L=62) still constructs.
    Bok = ProjectedBasis(L=62, N=1, base=2)
    @test size(Bok, 1) == 62
    # Ordinary small cases are unaffected.
    @test size(TensorBasis(L=10, base=2), 1) == 1024
    @test size(ProjectedBasis(L=8, N=4, base=2), 1) == binomial(8, 4)
end

@testset "abelian-2: AbelianBasis works with a non-default index dtype" begin
    # Int64 reference (half-filling, with translation symmetry).
    B64 = basis(L=6, N=3, k=0)
    # Same basis with a narrow index type must construct, not MethodError.
    B32 = basis(Int32; L=6, N=3, k=0)
    @test size(B32, 1) == size(B64, 1)
    @test eltype(B32.I) == Int32
    # Spectra agree: build the same Hamiltonian on both and compare sorted eigenvalues.
    mat = [1.0 0 0 0; 0 -1 2 0; 0 2 -1 0; 0 0 0 1]   # Heisenberg bond (XXZ-like) on 2 sites
    E64 = trans_inv_operator(mat, 2, B64) |> Array |> Hermitian |> eigvals
    E32 = trans_inv_operator(mat, 2, B32) |> Array |> Hermitian |> eigvals
    @test E32 ≈ E64
    # A non-default dtype also works on a base>2, no-symmetry-cap small case.
    B32b = basis(Int32; L=4, base=3, k=0)
    @test size(B32b, 1) > 0
end

@testset "parityflip-3 / sweep: symmetry bases guard index overflow" begin
    @test_throws ErrorException FlipBasis(L=63, p=1)
    @test_throws ErrorException ParityBasis(L=63, p=1)
    @test_throws ErrorException ParityFlipBasis(L=63, p=1, z=1)
    @test_throws ErrorException TranslationalBasis(L=63, k=0)
    # Small valid cases still construct.
    @test size(FlipBasis(L=6, p=1), 1) > 0
    @test size(TranslationalBasis(L=6, k=0), 1) > 0
end

# ---------------------------------------------------------------------------
# Phase 3 — eltype-correctness sweep
# ---------------------------------------------------------------------------

@testset "parityflip-1: Parity/Flip/ParityFlip have real eltype and type-stable index" begin
    L = 6
    @test eltype(ParityBasis(L=L, p=1)) == Float64
    @test eltype(FlipBasis(L=L, p=1)) == Float64
    @test eltype(ParityFlipBasis(L=L, p=1, z=1)) == Float64
    # index must be type-stable: returns a concrete Float64 tuple, not a Union with ComplexF64.
    for B in (ParityBasis(L=L, p=1), FlipBasis(L=L, p=1), ParityFlipBasis(L=L, p=1, z=1))
        r = @inferred EDKit.index(B, collect(zeros(Int, L)))
        @test r isa Tuple{Float64, <:Integer}
    end
end

@testset "abelian-5: real AbelianBasis sectors are Float64 and preserve spectrum" begin
    L = 6
    # Real-character sectors -> Float64; genuine momentum stays ComplexF64.
    @test eltype(basis(L=L, k=0)) == Float64
    @test eltype(basis(L=L, p=1)) == Float64
    @test eltype(basis(L=L, p=-1)) == Float64
    @test eltype(basis(L=L, z=-1)) == Float64
    @test eltype(basis(L=L, k=0, p=-1)) == Float64
    @test eltype(basis(L=L, k=1)) == ComplexF64

    # A real, reflection-symmetric Hamiltonian assembles real on a real sector,
    # and the parity sectors still partition the full spectrum (no value corruption).
    hmat = [1.0 0 0 0; 0 -1 2 0; 0 2 -1 0; 0 0 0 1]   # swap-symmetric real 2-site term
    Bp = basis(L=L, p=1)
    Hp = trans_inv_operator(hmat, 2, Bp) |> Array
    @test eltype(Hp) == Float64

    Efull = trans_inv_operator(hmat, 2, TensorBasis(L=L)) |> Array |> Hermitian |> eigvals
    Epar = Float64[]
    for s in (1, -1)
        append!(Epar, trans_inv_operator(hmat, 2, basis(L=L, p=s)) |> Array |> Hermitian |> eigvals)
    end
    @test sort(Epar) ≈ sort(Efull)
end

# ---------------------------------------------------------------------------
# Phase 4 — cached-sparse wiring
# ---------------------------------------------------------------------------

@testset "operator-2: matrix mul! consults the sparse cache" begin
    L = 6
    B = basis(L=L, N=3)
    H = trans_inv_operator([1.0 0 0 0; 0 -1 2 0; 0 2 -1 0; 0 0 0 1], 2, B)
    d = size(H, 1)
    T = promote_type(eltype(H), Float64)
    X = randn(T, d, 4)
    Yref = Array(H) * X

    clear_sparse_cache!()
    # Uncached path stays correct (regression guard): 3-arg accumulates, 5-arg scales.
    Y3 = zeros(T, d, 4); mul!(Y3, H, X);          @test Y3 ≈ Yref
    Y5 = randn(T, d, 4); mul!(Y5, H, X, 2.0, 0.0); @test Y5 ≈ 2 .* Yref

    # Cached path gives identical results.
    sparse!(H)
    Yc3 = zeros(T, d, 4); mul!(Yc3, H, X);          @test Yc3 ≈ Yref
    Yc5 = randn(T, d, 4); mul!(Yc5, H, X, 2.0, 0.0); @test Yc5 ≈ 2 .* Yref

    # Prove mul! actually READS the cache: poison the cached matrix and observe
    # the result change (matrix-free would ignore the poisoned entry).
    Sorig = EDKit._cached_sparse(H)
    EDKit._SPARSE_CACHE[objectid(H)] = 3 .* Sorig
    Yp = zeros(T, d, 4); mul!(Yp, H, X);  @test Yp ≈ 3 .* Yref
    clear_sparse_cache!()
end

# ---------------------------------------------------------------------------
# Phase 5 — long-tail correctness bugs
# ---------------------------------------------------------------------------

@testset "trans-1: order() returns the orbit size L/a" begin
    @test EDKit.order(TranslationalBasis(L=6, k=0, a=1)) == 6           # unchanged at a=1
    @test EDKit.order(TranslationalBasis(L=6, k=0, a=2)) == 3           # was 6
    @test EDKit.order(TranslationalBasis(L=6, k=0, a=3)) == 2           # was 6
    @test EDKit.order(TranslationParityBasis(L=6, k=0, p=1, a=2)) == 6  # 2*ncycle, was 12
    @test EDKit.order(TranslationFlipBasis(L=6, k=0, p=1, a=2)) == 6    # 2*ncycle, was 12
end

@testset "abelian-3: 2-arg index does not mutate the shared odometer" begin
    B = basis(L=8, k=1)                       # momentum sector with a non-trivial odometer
    s_before = copy(B.G.s)
    for _ in 1:50
        EDKit.index(B, rand(0:1, 8))
    end
    @test B.G.s == s_before                   # shared B.G must be left untouched (thread-safe)
    # And it still returns the same coefficient/index as the workspace-threaded internal call.
    dgt = rand(0:1, 8)
    @test EDKit.index(B, dgt) == EDKit.index(B, dgt, EDKit._shallow_workspace(B.G))
end

@testset "schmidt-1/2: Renyi entropy respects cutoff and normalizes" begin
    # schmidt-1: a noise eigenvalue below the cutoff is dropped (was ignored entirely).
    s = [0.25, 0.25, 0.25, 0.25, 1e-8]
    @test EDKit.entropy(s, α=0.5, cutoff=1e-6) ≈ log(4)
    # schmidt-2: unnormalized input is normalized before the Renyi formula.
    @test EDKit.entropy([0.5, 0.5, 0.5, 0.5], α=2) ≈ log(4)
end

@testset "expm-1: scaling-and-squaring is accurate at large norm" begin
    A = [0.0 5.0; -5.0 0.0]                    # ‖A‖ ≈ 5, where degree-10 Taylor diverges
    @test EDKit.expm(A) ≈ exp(A)
    M = randn(4, 4); B = 3 .* (M - M')         # antisymmetric, larger norm
    @test EDKit.expm(B) ≈ exp(B)
    v = randn(4)
    @test EDKit.expv(B, v) ≈ exp(B) * v
end

@testset "gapratio-1: sorting + degeneracy handling" begin
    @test gapratio([3.0, 1.0, 2.0]; sorted=false) == gapratio([1.0, 2.0, 3.0])
    r = gapratio([1.0, 1.0, 1.0, 2.0])               # 0/0 at the triple-degenerate gap
    @test isnan(r[1]) && r[2] == 0.0                 # was 1.0 / 0.0 (wrong)
    @test !isnan(meangapratio([1.0, 1.0, 1.0, 2.0])) # NaNs filtered out
end

@testset "operator-6: support dedup sums duplicate terms" begin
    X = [0.0 1.0; 1.0 0.0]
    op = operator([X, 2 .* X], [[1], [1]], TensorBasis(L=4))
    ref = operator([3 .* X], [[1]], TensorBasis(L=4))
    @test Array(op) ≈ Array(ref)
end

@testset "bug-4: change! accepts mismatched integer types for ind/base" begin
    dgt = zeros(Int32, 4)
    EDKit.change!(dgt, 6; base=2)     # ind::Int64, dgt::Int32 (was a MethodError)
    @test dgt == [0, 1, 0, 1]
    dgt3 = zeros(Int32, 3)
    EDKit.change!(dgt3, 9; base=3)    # base::Int64 (decimal 8 -> base-3 digits 0,2,2)
    @test dgt3 == [0, 2, 2]
end

@testset "parityflip-4: copy methods for parity/flip bases" begin
    for b in (ParityBasis(L=4, p=1), FlipBasis(L=4, p=1), ParityFlipBasis(L=4, p=1, z=1))
        c = copy(b)
        @test c.I === b.I            # representative list shared
        @test c.dgt !== b.dgt        # buffer independent
        @test c.dgt == b.dgt
    end
end
