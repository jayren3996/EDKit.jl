using LinearAlgebra
using SparseArrays
include("../src/EDKit.jl")
include("deprecate/ExactDiagonalization.jl")
using .EDKit
import .OldED as OED
using Test
#-------------------------------------------------------------------------------------------------------------------------
# Pauli matrices
#-------------------------------------------------------------------------------------------------------------------------
const X = [0 1; 1 0]
const Y = [0 -1; 1 0] * 1im
const Z = [1 0; 0 -1]
@testset "PauliKron" begin
    XII = operator([X], [[1]], 3)
    @test Array(XII) ≈ kron(X, I(4))
    IYI = operator([Y], [[2]], 3)
    @test Array(IYI) ≈ kron(I(2), Y, I(2))
    IIZ = operator([Z], [[3]], 3)
    @test Array(IIZ) ≈ kron(I(4), Z)
    # Multiplication
    V = rand(8)
    M = rand(ComplexF64, 8, 2)
    @test (XII * V) ≈ (kron(X, I(4)) * V)
    @test (IYI * M) ≈ (kron(I(2), Y, I(2)) * M)
end

#-------------------------------------------------------------------------------------------------------------------------
# Randon Matrices
#-------------------------------------------------------------------------------------------------------------------------
@testset "Real Random Matrices" begin
    L = 12
    println("\nL = $L System:")
    mats = [rand(4, 4) |> Hermitian |> Array for i=1:L]
    inds = [mod.([i-1, i], L) .+ 1 for i=1:L]
    print("New :")
    @time op1 = operator(mats, inds, L) |> Array
    print("Old :")
    @time op2 = OED.operation(mats, inds, L) |> Array
    @test op1 ≈ op2
end

@testset "ComplexF64 Random Matrices" begin
    L = 13
    println("\nL = $L System:")
    mats = [rand(ComplexF64, 4, 4) |> Hermitian |> Array for i=1:L]
    inds = [mod.([i-1, i], L) .+ 1 for i=1:L]
    print("New :")
    @time op1 = operator(mats, inds, L) |> Array
    print("Old :")
    @time op2 = OED.operation(mats, inds, L) |> Array
    @test op1 ≈ op2
end

@testset "Benchmark Sparse Matrix" begin
    L = 14
    println("\nBenchmark Sparse Matrix for L = $L System:")
    mat = kron(I(2), X, Y, Z)
    print("New :")
    @time op1 = trans_inv_operator(mat, 4, L) |> Array
    print("Old :")
    @time op2 = OED.trans_inv_operation(mat, 1:4, L) |> Array
end

@testset "Benchmark Matrix Multiplication" begin
    L = 14
    println("\nBenchmark Matrix Multiplication for L = $L System:")
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