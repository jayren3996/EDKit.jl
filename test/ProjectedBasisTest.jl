using LinearAlgebra
include("../src/EDKit.jl")
using .EDKit
using Test
#-------------------------------------------------------------------------------------------------------------------------
# XY Model
#-------------------------------------------------------------------------------------------------------------------------
@testset "XY" begin
    L = 4
    println("------------------------------")
    println("XY Model with L = $L")
    println("------------------------------")
    basis = projectedbasis(x->sum(x)==2, L)
    println("Basis Vector:")
    for i=1:length(basis.I)
        change!(basis, i) 
        println("$i : |$(basis.dgt...)⟩") 
    end
    X = [0 1; 1 0]
    Y = [0 -1; 1 0] * 1im
    mat = div.(real(kron(X, X) + kron(Y, Y)), 2)
    H = trans_inv_operator(mat, 2, basis)
    Hmat = Array(H)
    println("\nMatrix:")
    display(Hmat)
    print("\n")
    HXY = [0  1  0  0  1  0; 
           1  0  1  1  0  1; 
           0  1  0  0  1  0; 
           0  1  0  0  1  0; 
           1  0  1  1  0  1; 
           0  1  0  0  1  0]
    @test Hmat == HXY
end

#-------------------------------------------------------------------------------------------------------------------------
# PXP Model
#-------------------------------------------------------------------------------------------------------------------------
@testset "PXP" begin
    L = 4
    println("\n------------------------------")
    println("PXP Model with L = $L")
    println("------------------------------")
    pxpf(v::Vector{Int}) = all(v[i]==0 || v[mod(i, length(v))+1]==0 for i=1:length(v))
    basis = projectedbasis(pxpf, L)
    println("Basis Vector:")
    for i=1:length(basis.I)
        change!(basis, i) 
        println("$i : |$(basis.dgt...)⟩") 
    end
    P = Diagonal([1, 1, 1, 0, 1, 1, 0, 0])
    X = [0 1; 1 0]
    mat = P * kron(I(2), X, I(2)) * P
    H = trans_inv_operator(mat, 3, basis)
    Hmat = Array(H)
    println("\nMatrix:")
    display(Hmat)
    print("\n")
    HPXP = [0  1  1  1  0  1  0; 
            1  0  0  0  1  0  0; 
            1  0  0  0  0  0  1; 
            1  0  0  0  1  0  0; 
            0  1  0  1  0  0  0; 
            1  0  0  0  0  0  1; 
            0  0  1  0  0  1  0]
    @test Hmat == HPXP
end
