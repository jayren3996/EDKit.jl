include("../src/EDKit.jl")
using .EDKit
using LinearAlgebra, Test

#-------------------------------------------------------------------------------------------------------------------------
# SO(3) Random Quasi-symmetric Model
#-------------------------------------------------------------------------------------------------------------------------
@testset "SO(3) qsymm translation" begin
    L = 6
    h = rand()
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
        vecs * rm * vecs'
    end
    H  = trans_inv_operator(rh, 2, L)
    H += trans_inv_operator(spin((h, "z"), D=3), 1, L)
    E, V = Array(H) |> Hermitian |> eigen

    for inds in Any[1:L÷2, [1,2], [3,5], [2,4,6]]
        ts = Float64[]
        for i = 0:L-1
            tb = TranslationalBasis(i, L, base = 3)
            H  = trans_inv_operator(rh, 2, tb)
            H += trans_inv_operator(spin((h, "z"), D=3), 1, tb)
            e, v = Array(H) |> Hermitian |> eigen
            ee = [ent_S(v[:, i], inds, tb) for i=1:size(v, 2)]
            append!(ts, ee)
        end
        EE = [ent_S(V[:, i], inds, TensorBasis(L, base=3)) for i=1:size(V, 2)]
        @test norm(sort(EE) - sort(ts)) ≈ 0.0 atol = 1e-7
    end
end

@testset "SO(3) qsymm translation + inversion" begin
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
    EE(v, inds, basis) = [ent_S(v[:, i], inds, basis) for i=1:size(v, 2)]

    # K = 0
    ba = TranslationalBasis(0, L, base=3)
    be = TranslationParityBasis(0, +1, L, base=3)
    bo = TranslationParityBasis(0, -1, L, base=3)
    va, ve, vo = solve(ba), solve(be), solve(bo)
    for inds in Any[1:L÷2, [1,2], [3,5], [2,4,6]]
        eea = EE(va, inds, ba)
        eee = EE(ve, inds, be)
        eeo = EE(vo, inds, bo)
        @test norm(sort(eea) - sort(vcat(eee, eeo))) ≈ 0.0 atol = 1e-7
    end

    # K = π
    ba = TranslationalBasis(L÷2, L, base=3)
    be = TranslationParityBasis(L÷2, +1, L, base=3)
    bo = TranslationParityBasis(L÷2, -1, L, base=3)
    va, ve, vo = solve(ba), solve(be), solve(bo)
    for inds in Any[1:L÷2, [1,2], [3,5], [2,4,6]]
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
    EE(v, inds, basis) = [ent_S(v[:, i], inds, basis) for i=1:size(v, 2)]

    for i = 0:L-1
        ba = TranslationalBasis(i, L, base=3)
        be = TranslationFlipBasis(i, +1, L, base=3)
        bo = TranslationFlipBasis(i, -1, L, base=3)
        va, ve, vo = solve(ba), solve(be), solve(bo)
        for inds in Any[1:L÷2, [1,2], [3,5], [2,4,6]]
            eea = EE(va, inds, ba)
            eee = EE(ve, inds, be)
            eeo = EE(vo, inds, bo)
            @test norm(sort(eea) - sort(vcat(eee, eeo))) ≈ 0.0 atol = 1e-7
        end
    end

end

