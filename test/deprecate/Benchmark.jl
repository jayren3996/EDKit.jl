include("../../src/EDKit.jl")
include("ExactDiagonalization.jl")
using .EDKit
import .OldED as OED
using LinearAlgebra
using Test
#-------------------------------------------------------------------------------------------------------------------------
const X = [0 1; 1 0]
const Y = [0 -1; 1 0] * 1im
const Z = [1 0; 0 -1]
#-------------------------------------------------------------------------------------------------------------------------
@testset "Random Real 4×4 Matrices" begin
    L = 14
    println("Random Real 4×4 Matrices, L = $L")
    mats = [rand(4, 4) |> Hermitian |> Array for i=1:L]
    inds = [mod.([i-1, i], L) .+ 1 for i=1:L]
    print("New :")
    @time op1 = operator(mats, inds, L) |> Array
    print("Old :")
    @time op2 = OED.operation(mats, inds, L) |> Array
end
#-------------------------------------------------------------------------------------------------------------------------
@testset "Random Complex 4×4 Matrices" begin
    L = 14
    println("Random Complex 4×4 Matrices, L = $L")
    mats = [rand(ComplexF64, 4, 4) |> Hermitian |> Array for i=1:L]
    inds = [mod.([i-1, i], L) .+ 1 for i=1:L]
    print("New :")
    @time op1 = operator(mats, inds, L) |> Array
    print("Old :")
    @time op2 = OED.operation(mats, inds, L) |> Array
end
#-------------------------------------------------------------------------------------------------------------------------
@testset "Sparse Matrix IXYZ" begin
    L = 14
    println("Sparse Matrix IXYZ, L = $L")
    mat = kron(I(2), X, Y, Z)
    print("New :")
    @time op1 = trans_inv_operator(mat, 4, L) |> Array
    print("Old :")
    @time op2 = OED.trans_inv_operation(mat, 1:4, L) |> Array
end
#-------------------------------------------------------------------------------------------------------------------------
@testset "Operator Multiplication" begin
    L = 14
    println("Operator Multiplication, L = $L")
    mat = kron(I(2), X, Y, Z)
    v = rand(2^L)
    opt = trans_inv_operator(mat, 4, L)
    print("Operator : ")
    @time v2 = opt * v
    print("Dense    : ")
    optm = Array(opt)
    @time v3 = optm * v
    @test v2 ≈ v3
end
