include("../src/EDKit.jl")
using LinearAlgebra
using .EDKit
using Test
#-------------------------------------------------------------------------------------------------------------------------
# Base-2
#-------------------------------------------------------------------------------------------------------------------------
@testset "Spin-1/2 XY" begin
    L = 10
    mat = spin((1, "xx"), (1, "yy"))
    E = zeros(2^L)
    P = 0
    for n = 0:L
        basis = projectedbasis(x->sum(x)==n, L)
        if (l = size(basis, 1)) > 0
            vals = trans_inv_operator(mat, 2, basis) |> Array |> Hermitian |> eigvals
            E[P+1:P+l] = vals
            P += l
        end
    end
    @test P == 2^L
    vals = trans_inv_operator(mat, 2, L) |> Array |> Hermitian |> eigvals
    @test vals ≈ sort!(E)
end
#-------------------------------------------------------------------------------------------------------------------------
@testset "Spin-1/2 Random" begin
    L = 10
    mat = begin
        M = zeros(ComplexF64, 8, 8)
        M[1,1] = rand()
        M[8,8] = rand()
        M[[2,3,5], [2,3,5]] = rand(ComplexF64, 3, 3) |> Hermitian
        M[[4,6,7], [4,6,7]] = rand(ComplexF64, 3, 3) |> Hermitian
        M
    end
    E = zeros(2^L)
    P = 0
    for n = 0:L
        basis = projectedbasis(x->sum(x)==n, L)
        if (l = size(basis, 1)) > 0
            vals = trans_inv_operator(mat, 3, basis) |> Array |> Hermitian |> eigvals
            E[P+1:P+l] = vals
            P += l
        end
    end
    @test P == 2^L
    vals = trans_inv_operator(mat, 3, L) |> Array |> Hermitian |> eigvals
    @test vals ≈ sort!(E)
end
#-------------------------------------------------------------------------------------------------------------------------
@testset "PXP_OBC" begin
    function fib(n::Integer)
        if n == 0
            return 1
        elseif n == 1
            return 2
        end
        a, b = 1, 2
        for m = 2:n
            a, b = b, a+b
        end
        b
    end
    mat = begin
        P = Diagonal([1, 1, 1, 0, 1, 1, 0, 0])
        X = [0 1; 1 0]
        P * kron(I(2), X, I(2)) * P
    end
    pxpf(v::Vector{Int}) = all(v[i]==0 || v[i+1]==0 for i=1:length(v)-1)
    for L = 3:20
        basis = projectedbasis(pxpf, L)
        @test size(basis, 1) == fib(L)
    end
    for L = 3:10
        basis = projectedbasis(pxpf, L)
        mats = fill(mat, L-2)
        inds = [[i,i+1, i+2] for i=1:L-2]
        E = operator(mats, inds, basis) |> Array |> Hermitian |> eigvals
        H = operator(mats, inds, L) |> Array
        vals = H[basis.I, basis.I] |> Hermitian |> eigvals
        @test norm(vals - E) ≈ 0.0
    end
end
#-------------------------------------------------------------------------------------------------------------------------
@testset "PXP_PBC" begin
    mat = begin
        P = Diagonal([1, 1, 1, 0, 1, 1, 0, 0])
        X = [0 1; 1 0]
        P * kron(I(2), X, I(2)) * P
    end
    pxpf(v::Vector{Int}) = all(v[i]==0 || v[mod(i, length(v))+1]==0 for i=1:length(v))
    for L = 3:14
        basis = projectedbasis(pxpf, L)
        E = trans_inv_operator(mat, 3, basis) |> Array |> Hermitian |> eigvals
        H = trans_inv_operator(mat, 3, L) |> Array
        vals = H[basis.I, basis.I] |> Hermitian |> eigvals
        @test norm(vals - E) ≈ 0.0
    end
end

#-------------------------------------------------------------------------------------------------------------------------
# Base-3
#-------------------------------------------------------------------------------------------------------------------------
@testset "Spin-1 Random" begin
    L = 6
    R3 = begin    
        m2 = [2, 4, 10]
        m1 = [3, 5, 7, 11, 13, 19]
        m0 = [6, 8, 12, 14, 16, 20, 22]
        p1 = [9, 15, 17, 21, 23, 25]
        p2 = [18, 24, 26]
        mat = zeros(ComplexF64, 27, 27)
        mat[1, 1]    = rand()
        mat[27, 27]  = rand()
        mat[m2, m2] .= rand(ComplexF64, 3, 3) |> Hermitian
        mat[m1, m1] .= rand(ComplexF64, 6, 6) |> Hermitian
        mat[m0, m0] .= rand(ComplexF64, 7, 7) |> Hermitian
        mat[p1, p1] .= rand(ComplexF64, 6, 6) |> Hermitian
        mat[p2, p2] .= rand(ComplexF64, 3, 3) |> Hermitian
        mat
    end
    E = zeros(3^L)
    P = 0
    for n = 0:2L
        basis = projectedbasis(x->sum(x)==n, L, base=3)
        if (l = size(basis, 1)) > 0
            vals = trans_inv_operator(mat, 3, basis) |> Array |> Hermitian |> eigvals
            E[P+1:P+l] = vals
            P += l
        end
    end
    @test P == 3^L
    vals = trans_inv_operator(mat, 3, L) |> Array |> Hermitian |> eigvals
    @test vals ≈ sort!(E)
end
