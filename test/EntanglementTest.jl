include("../src/EDKit.jl")
using .EDKit
using LinearAlgebra, Test

#-------------------------------------------------------------------------------------------------------------------------
@testset "Random Hamiltonian" begin
    L = 6
    rm = randn(ComplexF64, 9, 9) |> Hermitian |> Array

    # Tensor Basis:
    B0 = TensorBasis(L=L, base=3)
    H0 = trans_inv_operator(rm, 2, B0)
    E, V = Array(H0) |> Hermitian |> eigen
    S = [ent_S(V[:, i], 1:L÷2, B0) for i in axes(V, 2)]
    
    # Momentum Resolved Bases:
    E2 = Float64[]
    S2 = Float64[]
    for i = 0:L-1
        B = TranslationalBasis(k=i, L=L, base=3)
        H = trans_inv_operator(rm, 2, B)
        e, v = Array(H) |> Hermitian |> eigen
        append!(E2, e)
        for j in axes(v, 2)
            push!(S2, ent_S(v[:, j], 1:L÷2, B))
        end
    end
    perm = sortperm(E2)
    @test E ≈ E2[perm]
    @test S ≈ S2[perm]
end
#-------------------------------------------------------------------------------------------------------------------------
@testset "Random: 2-site-cell" begin
    L = 6
    rm = randn(ComplexF64, 9, 9) |> Hermitian |> Array

    # Tensor Basis:
    B0 = TensorBasis(L=L, base=3)
    H0 = trans_inv_operator(rm, 2, B0)
    E, V = Array(H0) |> Hermitian |> eigen
    S = [ent_S(V[:, i], 1:L÷2, B0) for i in axes(V, 2)]
    
    # Momentum Resolved Bases:
    E2 = Float64[]
    S2 = Float64[]
    for i = 0:L÷2-1
        B = TranslationalBasis(k=i, L=L, base=3, a=2)
        H = trans_inv_operator(rm, 2, B)
        e, v = Array(H) |> Hermitian |> eigen
        append!(E2, e)
        for j in axes(v, 2)
            push!(S2, ent_S(v[:, j], 1:L÷2, B))
        end
    end
    perm = sortperm(E2)
    @test E ≈ E2[perm]
    @test S ≈ S2[perm]
end

#-------------------------------------------------------------------------------------------------------------------------
@testset "Qsymm: T+P" begin
    L = 6
    h = rand()
    perm = [1, 4, 7, 2, 5, 8, 3, 6, 9]
    rh = begin
        vecs = zeros(ComplexF64, 9, 4)
        vecs[8,1] = 1
        vecs[6,1] = -1
        vecs[7,2] = 1
        vecs[3,2] = -1
        vecs[4,3] = 1
        vecs[2,3] = -1
        vecs[3,4] = 1
        vecs[5,4] = -1
        vecs[7,4] = 1
        rm = rand(ComplexF64, 4,4) .- 0.5 |> Hermitian |> Array
        mat = vecs * rm * vecs'
        mat + mat[perm, perm]
    end

    function solve(basis)
        H  = trans_inv_operator(rh, 2, basis)
        H += trans_inv_operator(spin((h, "z"), D=3), 1, basis)
        e, v = Array(H) |> Hermitian |> eigen
        v
    end
    EE(v, inds, basis) = [ent_S(v[:, i], inds, basis) for i in axes(v, 2)]

    # K = 0
    ba = TranslationalBasis(k=0, L=L, base=3)
    be = TranslationParityBasis(k=0, p=+1, L=L, base=3)
    bo = TranslationParityBasis(k=0, p=-1, L=L, base=3)
    va, ve, vo = solve(ba), solve(be), solve(bo)
    for inds in Any[1:L÷2, [1,2], [3,5], [2,4,6]]
        eea = EE(va, inds, ba)
        eee = EE(ve, inds, be)
        eeo = EE(vo, inds, bo)
        @test norm(sort(eea) - sort(vcat(eee, eeo))) ≈ 0.0 atol = 1e-7
    end

    # K = π
    ba = TranslationalBasis(k=L÷2, L=L, base=3)
    be = TranslationParityBasis(k=L÷2, p=+1, L=L, base=3)
    bo = TranslationParityBasis(k=L÷2, p=-1, L=L, base=3)
    va, ve, vo = solve(ba), solve(be), solve(bo)
    for inds in Any[1:L÷2, [1,2], [3,5], [2,4,6]]
        eea = EE(va, inds, ba)
        eee = EE(ve, inds, be)
        eeo = EE(vo, inds, bo)
        @test norm(sort(eea) - sort(vcat(eee, eeo))) ≈ 0.0 atol = 1e-7
    end
end
#-------------------------------------------------------------------------------------------------------------------------
@testset "Qsymm: T+P 2-cell" begin
    L = 4
    h = rand()
    perm = [1, 4, 7, 2, 5, 8, 3, 6, 9]
    rh = begin
        vecs = zeros(ComplexF64, 9, 4)
        vecs[8,1] = 1
        vecs[6,1] = -1
        vecs[7,2] = 1
        vecs[3,2] = -1
        vecs[4,3] = 1
        vecs[2,3] = -1
        vecs[3,4] = 1
        vecs[5,4] = -1
        vecs[7,4] = 1
        rm = rand(ComplexF64, 4,4) .- 0.5 |> Hermitian |> Array
        mat = vecs * rm * vecs'
        mat + mat[perm, perm]
    end

    function solve(basis)
        H  = trans_inv_operator(rh, 2, basis)
        H += trans_inv_operator(spin((h, "z"), D=3), 1, basis)
        e, v = Array(H) |> Hermitian |> eigen
        v
    end
    EE(v, inds, basis) = [ent_S(v[:, i], inds, basis) for i in axes(v, 2)]

    # K = 0
    ba = TranslationalBasis(k=0, L=L, base=3, a=2)
    be = TranslationParityBasis(k=0, p=+1, L=L, base=3, a=2)
    bo = TranslationParityBasis(k=0, p=-1, L=L, base=3, a=2)
    va, ve, vo = solve(ba), solve(be), solve(bo)
    for inds in Any[[1,2], [1,3], [1,4,2]]
        eea = EE(va, inds, ba)
        eee = EE(ve, inds, be)
        eeo = EE(vo, inds, bo)
        @test norm(sort(eea) - sort(vcat(eee, eeo))) ≈ 0.0 atol = 1e-7
    end

    # K = π
    ba = TranslationalBasis(k=L÷2, L=L, base=3, a=2)
    be = TranslationParityBasis(k=L÷2, p=+1, L=L, base=3, a=2)
    bo = TranslationParityBasis(k=L÷2, p=-1, L=L, base=3, a=2)
    va, ve, vo = solve(ba), solve(be), solve(bo)
    for inds in Any[[1,4], [3,2], [2,4,3]]
        eea = EE(va, inds, ba)
        eee = EE(ve, inds, be)
        eeo = EE(vo, inds, bo)
        @test norm(sort(eea) - sort(vcat(eee, eeo))) ≈ 0.0 atol = 1e-7
    end
end

#-------------------------------------------------------------------------------------------------------------------------
# Random model with spin-flip symmetry
#-------------------------------------------------------------------------------------------------------------------------
@testset "Random spin-flip" begin
    L = 6
    perm = [9,8,7,6,5,4,3,2,1]
    rh = begin
        mat = rand(ComplexF64, 9, 9) |> Hermitian |> Array
        mat + mat[perm, perm]
    end

    function solve(basis)
        H  = trans_inv_operator(rh, 2, basis)
        e, v = Array(H) |> Hermitian |> eigen
        v
    end
    EE(v, inds, basis) = [ent_S(v[:, j], inds, basis) for j in axes(v, 2)]

    for i = 0:L-1
        ba = TranslationalBasis(k=i, L=L, base=3)
        be = TranslationFlipBasis(k=i, p=+1, L=L, base=3)
        bo = TranslationFlipBasis(k=i, p=-1, L=L, base=3)
        va, ve, vo = solve(ba), solve(be), solve(bo)
        for inds in Any[1:L÷2, [1,2], [3,5], [2,4,6]]
            eea = EE(va, inds, ba)
            eee = EE(ve, inds, be)
            eeo = EE(vo, inds, bo)
            @test norm(sort(eea) - sort(vcat(eee, eeo))) ≈ 0.0 atol = 1e-7
        end
    end
end

@testset "Random spin-flip 2-cell" begin
    L = 4
    perm = [9,8,7,6,5,4,3,2,1]
    rh = begin
        mat = rand(ComplexF64, 9, 9) |> Hermitian |> Array
        mat + mat[perm, perm]
    end

    function solve(basis)
        H  = trans_inv_operator(rh, 2, basis)
        e, v = Array(H) |> Hermitian |> eigen
        v
    end
    EE(v, inds, basis) = [ent_S(v[:, j], inds, basis) for j in axes(v, 2)]

    for i = 0:L÷2-1
        ba = TranslationalBasis(k=i, L=L, base=3, a=2)
        be = TranslationFlipBasis(k=i, p=+1, L=L, base=3, a=2)
        bo = TranslationFlipBasis(k=i, p=-1, L=L, base=3, a=2)
        va, ve, vo = solve(ba), solve(be), solve(bo)
        for inds in Any[[1,2], [3,1], [2,4,1]]
            eea = EE(va, inds, ba)
            eee = EE(ve, inds, be)
            eeo = EE(vo, inds, bo)
            @test norm(sort(eea) - sort(vcat(eee, eeo))) ≈ 0.0 atol = 1e-7
        end
    end
end

