using LinearAlgebra, Random, ITensors, ITensorMPS

DEV = true
if DEV
    include("../../src/EDKit.jl")
    using .EDKit
else
    using EDKit
end

Random.seed!(7)

orbit_order(B) = B isa EDKit.AbstractOnsiteBasis ? 1 : EDKit.order(B)

function sector_embedding(B::AbstractBasis)
    Bc = B
    full = TensorBasis(L=length(Bc), base=Bc.B)
    P = zeros(ComplexF64, size(full, 1), size(Bc, 1))
    for j in 1:size(full, 1)
        change!(full, j)
        Bc.dgt .= full.dgt
        C, pos = index(Bc)
        iszero(C) || (P[j, pos] = C / orbit_order(Bc))
    end
    P
end

function sector_to_full(v::AbstractVector, B::AbstractBasis)
    sector_embedding(B) * v
end

function roundtrip_error(B::AbstractBasis)
    P = sector_embedding(B)
    v_sector = randn(ComplexF64, size(B, 1)) |> normalize
    ψ_full = P * v_sector
    s = siteinds(B.B, length(B))
    ψ_mps = vec2mps(ψ_full, s)
    v_back = mps2vec(ψ_mps, B)
    norm(v_back - v_sector), norm(P' * ψ_full - v_sector)
end

function projection_mismatch(B::AbstractBasis)
    P = sector_embedding(B)
    ψ_full = randn(ComplexF64, size(P, 1)) |> normalize
    s = siteinds(B.B, length(B))
    ψ_mps = vec2mps(ψ_full, s)
    v_mps = mps2vec(ψ_mps, B)
    v_exact = P' * ψ_full
    norm(v_mps - v_exact)
end

L = 6
Nhalf = L ÷ 2

sectors = [
    ("N", basis(L=L, N=Nhalf)),
    ("k", basis(L=L, k=1)),
    ("p", basis(L=L, p=1)),
    ("z", basis(L=L, z=1)),
    ("N+k", basis(L=L, N=Nhalf, k=1)),
    ("N+p", basis(L=L, N=Nhalf, p=1)),
    ("N+z", basis(L=L, N=Nhalf, z=1)),
    ("k+p", basis(L=L, k=0, p=1)),
    ("k+z", basis(L=L, k=1, z=1)),
    ("p+z", basis(L=L, p=1, z=1)),
    ("N+k+p", basis(L=L, N=Nhalf, k=0, p=1)),
    ("N+k+z", basis(L=L, N=Nhalf, k=1, z=1)),
    ("N+p+z", basis(L=L, N=Nhalf, p=1, z=1)),
    ("k+p+z", basis(L=L, k=0, p=1, z=1)),
    ("N+k+p+z", basis(L=L, N=Nhalf, k=0, p=1, z=1)),
]

results = map(sectors) do (name, B)
    err_mps, err_exact = roundtrip_error(B)
    (name=name, dim=size(B, 1), mps2vec_roundtrip=err_mps, exact_projector_roundtrip=err_exact)
end

for row in results
    println(rpad(row.name, 10), " dim=", lpad(row.dim, 4),
            "   |mps2vec - v| = ", row.mps2vec_roundtrip,
            "   |P'Pv - v| = ", row.exact_projector_roundtrip)
end

Btest = basis(L=L, N=Nhalf, k=0, p=1)
projection_mismatch(Btest)

B = basis(L=L, N=Nhalf, k=0, p=1)
H = trans_inv_operator(spin((1.0, "xx"), (1.0, "yy"), (0.5, "zz")), 2, B)
vals, vecs = eigen(Hermitian(H))
v_sector = vecs[:, 1]

ψ_full = sector_to_full(v_sector, B)
ψ_mps = vec2mps(ψ_full, siteinds(B.B, length(B)))
v_back = mps2vec(ψ_mps, B)

println("ground-state recovery error = ", norm(v_back - v_sector))

