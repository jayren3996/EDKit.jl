include("../../src/EDKit.jl")
using .EDKit
using LinearAlgebra
using DelimitedFiles
using Test
using PyCall
np = pyimport("numpy")

@testset "PXP-L-26" begin
    L, k, p = 26, 0, 1
    mat = begin
        P = Diagonal([1, 1, 1, 0, 1, 1, 0, 0])
        X = [0 1; 1 0]
        P * kron(I(2), X, I(2)) * P
    end
    pxpf(v::Vector{<:Integer}) = all(v[i]==0 || v[mod(i, length(v))+1]==0 for i=1:length(v))
    print("Searching time :")
    @time basis = TranslationParityBasis(pxpf, k, p, L, threaded=true)
    H = trans_inv_operator(mat, 2, basis) |> Array
    E = np.linalg.eigvalsh(H)
    PXPE = readdlm("evals_periodic_N26_k0_p0.npy.txt")
    @test norm(PXPE-E) ≈ 0.0 atol=1e-10
end

@testset "PXP-L-28" begin
    L, k, p = 28, 0, 1
    mat = begin
        P = Diagonal([1, 1, 1, 0, 1, 1, 0, 0])
        X = [0 1; 1 0]
        P * kron(I(2), X, I(2)) * P
    end
    pxpf(v::Vector{<:Integer}) = all(v[i]==0 || v[mod(i, length(v))+1]==0 for i=1:length(v))
    print("Searching time :")
    @time basis = TranslationParityBasis(pxpf, k, p, L, threaded=true)
    H = trans_inv_operator(mat, 2, basis) |> Array
    E = np.linalg.eigvalsh(H)
    PXPE = readdlm("evals_periodic_N28_k0_p0.npy.txt")
    @test norm(PXPE-E) ≈ 0.0 atol=1e-10
end