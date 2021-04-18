include("../src/EDKit.jl")
include("deprecate/ExactDiagonalization.jl")
using .EDKit
using LinearAlgebra
using Test
#-------------------------------------------------------------------------------------------------------------------------
directkron(mat, i, j, L, base=2) = kron(I(base^(i-1)), mat, I(base^(L-j)))
boundarykron(m1, m2, i, j, base=2) = kron(m1, I(base^(j-i+1)), m2)
#-------------------------------------------------------------------------------------------------------------------------
@testset "Base-2 Random Kron" begin
    L = 10
    for i=1:10
        mats = [rand(2^i, 2^i) for j=1:L-i+1]
        inds = [j:j+i-1 for j=1:L-i+1]
        M_test = operator(mats, inds, L) |> Array 
        M_targ = sum(directkron(mats[j], j, j+i-1, L) for j = 1:L-i+1)
        @test M_test ≈ M_targ
    end
end
#-------------------------------------------------------------------------------------------------------------------------
@testset "Base-3 Random Kron" begin
    L = 6
    for i=1:6
        mats = [rand(3^i, 3^i) for j=1:L-i+1]
        inds = [j:j+i-1 for j=1:L-i+1]
        M_test = operator(mats, inds, L) |> Array 
        M_targ = sum(directkron(mats[j], j, j+i-1, L, 3) for j = 1:L-i+1)
        @test M_test ≈ M_targ
    end
end
#-------------------------------------------------------------------------------------------------------------------------
@testset "Base-3 Boundary Kron" begin
    L = 6
    for i = 1:3, j = 1:3
        m1, m2 = rand(3^i, 3^i), rand(3^j, 3^j)
        mats = kron(m1, m2)
        inds = vcat(Vector(1:i), Vector(L-j+1:L))
        M_test = operator(mats, inds, L) |> Array 
        M_targ = boundarykron(m1, m2, i+1, L-j, 3)
        @test M_test ≈ M_targ
    end
end

#-------------------------------------------------------------------------------------------------------------------------
const X = [0 1; 1 0]
const Y = [0 -1; 1 0] * 1im
const Z = [1 0; 0 -1]
@testset "Operator Multiplication" begin
    L = 10
    mat = kron(I(2), X, Y, Z)
    v = rand(2^L)
    opt = trans_inv_operator(mat, 4, L)
    v2 = opt * v
    optm = Array(opt)
    v3 = optm * v
    @test v2 ≈ v3
end
