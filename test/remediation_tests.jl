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
