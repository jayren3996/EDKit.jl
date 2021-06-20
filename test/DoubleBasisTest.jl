using LinearAlgebra
include("../src/EDKit.jl")
using .EDKit
using Test

@testset "Basic" begin
    L = 10
    b1 = TranslationalBasis(x -> sum(x) == 2, 0, L)
    b2 = TranslationalBasis(x -> sum(x) == 0, 1, L)
    b = DoubleBasis(b1, b2)
    mat = rand(2,2) |> Hermitian
    H = trans_inv_operator(mat, 1, b) |> Array
    @test size(H) == (5, 0)
end

@testset "Spin-1 Random" begin
    L = 6
    Sp = begin
        p3 = [1]
        p2 = [2, 4, 10]
        p1 = [3, 5, 7, 11, 13, 19]
        z0 = [6, 8, 12, 14, 16, 20, 22]
        m1 = [9, 15, 17, 21, 23, 25]
        m2 = [18, 24, 26]
        m3 = [27]
        mat = zeros(ComplexF64, 27, 27)
        mat[m2, m3] .= rand(ComplexF64, 3, 1)
        mat[m1, m2] .= rand(ComplexF64, 6, 3)
        mat[z0, m1] .= rand(ComplexF64, 7, 6)
        mat[p1, z0] .= rand(ComplexF64, 6, 7)
        mat[p2, p1] .= rand(ComplexF64, 3, 6)
        mat[p3, p2] .= rand(ComplexF64, 1, 3)
        mat
    end
    XY = spin((1, "+-"), (1, "-+"), (1, "z1"), (1, "1z"), D=3)
    A = 0.0
    for n = 1:2L
        b1 = ProjectedBasis(x -> sum(x)==(n-1), L, base=3)
        b2 = ProjectedBasis(x -> sum(x)==n, L, base=3)
        basis = DoubleBasis(b1, b2)
        e1, v1 = trans_inv_operator(XY, 2, b1) |> Array |> Hermitian |> eigen
        e2, v2 = trans_inv_operator(XY, 2, b2) |> Array |> Hermitian |> eigen
        op = trans_inv_operator(Sp, 3, basis) |> Array
        A += abs2.(v1' * op * v2) |> sum
    end
    e, v = trans_inv_operator(XY, 2, L) |> Array |> Hermitian |> eigen
    op = trans_inv_operator(Sp, 3, L) |> Array
    B = abs2.(v' * op * v) |> sum
    @test A ≈ B atol = 1e-7
end
