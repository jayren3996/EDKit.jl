using LinearAlgebra
include("../src/EDKit.jl")
using .EDKit
using Test
#-------------------------------------------------------------------------------------------------------------------------
# Spin-1/2 XY Model
#-------------------------------------------------------------------------------------------------------------------------
@testset "Spin-1/2 XY" begin
    mat = spin((1, "+-"), (1, "-+"), (1, "z1"), (1, "1z"))
    for L = 2:2:10
        for n = 0:L, k in [0, L ÷ 2]
            be = translationparitybasis(x->sum(x)==n, k, +1, L)
            bo = translationparitybasis(x->sum(x)==n, k, -1, L)
            ba = translationalbasis(x->sum(x)==n, k, L)
            ve = trans_inv_operator(mat, 2, be) |> Array |> Hermitian |> eigvals
            vo = trans_inv_operator(mat, 2, bo) |> Array |> Hermitian |> eigvals
            va = trans_inv_operator(mat, 2, ba) |> Array |> Hermitian |> eigvals
            E = sort(vcat(ve, vo))
            @test norm(E-va) ≈ 0.0 atol = 1e-12
        end
    end
end

#-------------------------------------------------------------------------------------------------------------------------
# Spin-1 XY Model
#-------------------------------------------------------------------------------------------------------------------------
@testset "Spin-1 XY" begin
    mat = spin((1, "+-"), (1, "-+"), (1, "z1"), (1, "1z"), D=3)
    for L = 2:2:6
        for n = 0:2L, k in [0, L ÷ 2]
            be = translationparitybasis(x->sum(x)==n, k, +1, L, base=3)
            bo = translationparitybasis(x->sum(x)==n, k, -1, L, base=3)
            ba = translationalbasis(x->sum(x)==n, k, L, base=3)
            ve = trans_inv_operator(mat, 2, be) |> Array |> Hermitian |> eigvals
            vo = trans_inv_operator(mat, 2, bo) |> Array |> Hermitian |> eigvals
            va = trans_inv_operator(mat, 2, ba) |> Array |> Hermitian |> eigvals
            E = sort(vcat(ve, vo))
            @test norm(E-va) ≈ 0.0 atol = 1e-12
        end
    end
end

#-------------------------------------------------------------------------------------------------------------------------
# PXP Model
#-------------------------------------------------------------------------------------------------------------------------
@testset "PXP" begin
    mat = begin
        P = Diagonal([1, 1, 1, 0, 1, 1, 0, 0])
        X = [0 1; 1 0]
        P * kron(I(2), X, I(2)) * P
    end
    pxpf(v::Vector{Int}) = all(v[i]==0 || v[mod(i, length(v))+1]==0 for i=1:length(v))
    for L = 4:2:20
        for k in [0, L ÷ 2]
            be = translationparitybasis(pxpf, k, +1, L)
            bo = translationparitybasis(pxpf, k, -1, L)
            ba = translationalbasis(pxpf, k, L)
            ve = trans_inv_operator(mat, 2, be) |> Array |> Hermitian |> eigvals
            vo = trans_inv_operator(mat, 2, bo) |> Array |> Hermitian |> eigvals
            va = trans_inv_operator(mat, 2, ba) |> Array |> Hermitian |> eigvals
            E = sort(vcat(ve, vo))
            @test norm(E-va) ≈ 0.0 atol = 1e-12
        end
    end
end
