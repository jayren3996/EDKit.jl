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
