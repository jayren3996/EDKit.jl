using LinearAlgebra
using SparseArrays

const DEV = true

if DEV
    include("../src/EDKit.jl")
    using .EDKit
else
    using EDKit
end

L = 6
delta = 1.1

mats = [
    fill(spin("XX"), L - 1);
    fill(spin("YY"), L - 1);
    fill(delta * spin("ZZ"), L - 1);
]

inds = vcat(
    [[i, i + 1] for i in 1:L-1],
    [[i, i + 1] for i in 1:L-1],
    [[i, i + 1] for i in 1:L-1],
)

H = operator(mats, inds, L)

psi = randn(ComplexF64, 2^L) |> normalize
Hpsi_linear = H * psi
H_dense = Array(H)
H_sparse = sparse(H)

println("Basic many-body operator construction")
println("  system size                 = $L")
println("  Hilbert-space dimension     = $(size(H, 1))")
println("  number of local terms       = $(length(H))")
println("  linear-vs-dense error       = $(norm(Hpsi_linear - H_dense * psi))")
println("  dense-vs-sparse error       = $(norm(H_dense - Matrix(H_sparse)))")
println("  lowest three eigenvalues    = $(eigvals(Hermitian(H_dense))[1:3])")
