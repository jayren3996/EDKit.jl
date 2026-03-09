using LinearAlgebra
using Random
using SparseArrays
using Test
using ITensors
using ITensorMPS
using EDKit

Random.seed!(7)

orbit_order(B) = B isa EDKit.AbstractOnsiteBasis ? 1 : EDKit.order(B)

function sector_embedding(B::EDKit.AbstractBasis)
    full = TensorBasis(L = length(B), base = B.B)
    P = zeros(ComplexF64, size(full, 1), size(B, 1))
    for j in 1:size(full, 1)
        change!(full, j)
        B.dgt .= full.dgt
        C, pos = index(B)
        iszero(C) || (P[j, pos] = C / orbit_order(B))
    end
    P
end

function kron_embed(mat::AbstractMatrix, inds::AbstractVector{<:Integer}, L::Integer, base::Integer = 2)
    dims = fill(base, L)
    full = zeros(ComplexF64, prod(dims), prod(dims))
    all_sites = collect(1:L)
    rest = [i for i in all_sites if i ∉ inds]
    B = TensorBasis(L = L, base = base)
    d1 = zeros(Int, L)
    d2 = zeros(Int, L)
    for col in 1:size(full, 2)
        change!(d2, col, base = base)
        j = index(d2, inds, base = base)
        for row_local in 1:size(mat, 1)
            d1 .= d2
            change!(d1, inds, row_local, base = base)
            if all(d1[r] == d2[r] for r in rest)
                row = index(d1, base = base)
                full[row, col] += mat[row_local, j]
            end
        end
    end
    full
end
