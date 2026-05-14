using EDKit
using LinearAlgebra
using Test

const RUN_SLOW_MULTITHREAD_TESTS = get(ENV, "EDKIT_SLOW_TESTS", "0") == "1"

pxp_constraint(v::Vector{<:Integer}) =
    all(v[i] == 0 || v[mod1(i + 1, length(v))] == 0 for i in eachindex(v))

function assert_threaded_basis_matches_serial(serial, threaded)
    @test threaded.I == serial.I
    @test threaded.R ≈ serial.R
end

@testset "Multi-threaded PXP basis construction" begin
    L = RUN_SLOW_MULTITHREAD_TESTS ? 28 : 10

    @testset "TranslationalBasis" begin
        serial = TranslationalBasis(f = pxp_constraint, k = 0, L = L, threaded = false)
        threaded = TranslationalBasis(f = pxp_constraint, k = 0, L = L, threaded = true)
        assert_threaded_basis_matches_serial(serial, threaded)
    end

    @testset "TranslationParityBasis" begin
        serial = TranslationParityBasis(f = pxp_constraint, k = 0, p = 1, L = L, threaded = false)
        threaded = TranslationParityBasis(f = pxp_constraint, k = 0, p = 1, L = L, threaded = true)
        assert_threaded_basis_matches_serial(serial, threaded)
    end
end
