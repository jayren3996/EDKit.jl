include("../src/EDKit.jl")
using .EDKit
using LinearAlgebra, Test

#-------------------------------------------------------------------------------------------------------------------------
# SO(3) Random Quasi-symmetric Model
#-------------------------------------------------------------------------------------------------------------------------
@testset "SO(3) Quasisymmetry" begin
    L = 6
    h = 0.37
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
    ts = Float64[]
    for i = 0:L-1
        tb = translationalbasis(x->true, i, L, base = 3)
        H  = trans_inv_operator(rh, 2, tb)
        H += trans_inv_operator(spin((h, "z"), D=3), 1, tb)
        e, v = Array(H) |> Hermitian |> eigen
        ee = [ent_S(v[:, i], 1:L÷2, tb) for i=1:size(v, 2)]
        append!(ts, ee)
    end
    H  = trans_inv_operator(rh, 2, L)
    H += trans_inv_operator(spin((h, "z"), D=3), 1, L)
    e, v = Array(H) |> Hermitian |> eigen
    ee = [ent_S(v[:, i], 1:L÷2, tensorbasis(L, base=3)) for i=1:size(v, 2)]
    @test norm(sort(ee) - sort(ts)) ≈ 0.0 atol = 1e-7
end