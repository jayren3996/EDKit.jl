using LinearAlgebra
include("../src/EDKit.jl")
using .EDKit
using Test

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
    A = zeros(3^L, 3^L)
    p = 3^L
    for n = 1-L:L
        b1 = projectedbasis(x -> sum(x)==n, L, base=3)
        b2 = projectedbasis(x -> sum(x)==(n-1), L, base=3)
        l1, l2 = length(b1), length(b2)
        basis = doublebasis(b1, b2)
        op = trans_inv_operator(mat, 3, basis) |> Array
        A[p]
    end
    @test P == 3^L
    vals = trans_inv_operator(mat, 3, L) |> Array |> Hermitian |> eigvals
    @test vals ≈ sort!(E)
end