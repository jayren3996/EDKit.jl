using LinearAlgebra

const DEV = true

if DEV
    include("../src/EDKit.jl")
    using .EDKit
else
    using EDKit
end

L = 8
bond = spin((1.0, "xx"), (1.0, "yy"))

full_basis = TensorBasis(L = L, base = 2)
nk_basis = basis(L = L, N = L ÷ 2, k = 0)
even_basis = basis(L = L, N = L ÷ 2, k = 0, p = 1)
odd_basis = basis(L = L, N = L ÷ 2, k = 0, p = -1)

full_H = trans_inv_operator(bond, 2, full_basis)
nk_H = trans_inv_operator(bond, 2, nk_basis)
even_H = trans_inv_operator(bond, 2, even_basis)
odd_H = trans_inv_operator(bond, 2, odd_basis)

nk_vals = eigvals(Hermitian(Array(nk_H))) |> sort
split_vals = vcat(
    eigvals(Hermitian(Array(even_H))),
    eigvals(Hermitian(Array(odd_H))),
) |> sort

println("Symmetry-sector reduction")
println("  full dimension              = $(size(full_basis, 1))")
println("  dim(N = L/2, k = 0)         = $(size(nk_basis, 1))")
println("  dim(N = L/2, k = 0, p = +1) = $(size(even_basis, 1))")
println("  dim(N = L/2, k = 0, p = -1) = $(size(odd_basis, 1))")
println("  parity recombination error  = $(norm(nk_vals - split_vals))")
println("  lowest five sector energies = $(nk_vals[1:5])")
