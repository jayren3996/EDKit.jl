using LinearAlgebra
include("../src/EDKit.jl")
using .EDKit
using Test
#-------------------------------------------------------------------------------------------------------------------------
# XY Model
#-------------------------------------------------------------------------------------------------------------------------
L = 4
println("------------------------------")
println("XY Model with L = $L")
println("------------------------------")
println("Basis Vector:")


basis = translationalbasis(x->sum(x)==3, 0, L)

for i=1:length(basis.I)
    change!(basis, i)
    println("$i : |$(basis.dgt...)⟩")
end

H = begin
    X = [0 1; 1 0]
    Y = [0 -1; 1 0] * 1im
    mat = div.(real(kron(X, X) + kron(Y, Y)), 2)
    trans_inv_operator(mat, 2, basis)
end

Hmat = Array(H)
println("\nMatrix:")
display(Hmat)
println("\n")

#-------------------------------------------------------------------------------------------------------------------------
# PXP Model
#-------------------------------------------------------------------------------------------------------------------------
L = 26
println("------------------------------")
println("PXP Model with L = $L")
println("------------------------------")
println("Basis Vector:")

pxpf(v::Vector{Int}) = all(v[i]==0 || v[mod(i, length(v))+1]==0 for i=1:length(v))
basis = translationalbasis(pxpf, 0, L)

#for i=1:length(basis.I)
#    change!(basis, i) 
#    println("$i : |$(basis.dgt...)⟩") 
#end

H = begin
    P = Diagonal([1, 1, 1, 0, 1, 1, 0, 0])
    X = [0 1; 1 0]
    mat = P * kron(I(2), X, I(2)) * P
    trans_inv_operator(mat, 3, basis)
end

Hmat = Array(H)
println("\nMatrix:")
println(size(Hmat))
print("\n")
