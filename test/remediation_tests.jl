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
