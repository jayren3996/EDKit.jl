using LinearAlgebra
using Random
using SparseArrays
using Test
using ITensors
using ITensorMPS
using EDKit

stable_rng(seed) = MersenneTwister(seed)

function random_state(B::EDKit.AbstractBasis; rng=stable_rng(7), T=ComplexF64)
    v = randn(rng, T, size(B, 1))
    isempty(v) ? v : normalize(v)
end

function random_vector(n::Integer; rng=stable_rng(7), T=ComplexF64)
    randn(rng, T, n)
end

function random_matrix(n::Integer, columns::Integer=3; rng=stable_rng(7), T=ComplexF64)
    randn(rng, T, n, columns)
end

orbit_order(B) = B isa EDKit.AbstractOnsiteBasis ? 1 : EDKit.order(B)

function sector_embedding(B::EDKit.AbstractBasis)
    full = TensorBasis(L = length(B), base = B.B)
    P = zeros(ComplexF64, size(full, 1), size(B, 1))
    full_dgt = similar(full.dgt)
    sector_dgt = similar(B.dgt)
    for j in 1:size(full, 1)
        change!(full, j, full_dgt)
        sector_dgt .= full_dgt
        C, pos = index(B, sector_dgt)
        iszero(C) || (P[j, pos] = C / orbit_order(B))
    end
    P
end

function assert_sector_embedding_isometry(B::EDKit.AbstractBasis; atol=1e-10, rtol=1e-10)
    P = sector_embedding(B)
    expected = Matrix{ComplexF64}(I, size(B, 1), size(B, 1))
    @test P' * P ≈ expected atol = atol rtol = rtol
    P
end

function dense_local_embedding(
    mat::AbstractMatrix,
    inds::AbstractVector{<:Integer},
    L::Integer;
    base::Integer=2,
)
    full_dim = base^L
    full = zeros(ComplexF64, full_dim, full_dim)
    row_dgt = zeros(Int, L)
    col_dgt = zeros(Int, L)

    for col in 1:full_dim
        change!(col_dgt, col, base = base)
        local_col = index(col_dgt, inds, base = base)
        for row_local in 1:size(mat, 1)
            val = mat[row_local, local_col]
            iszero(val) && continue
            row_dgt .= col_dgt
            change!(row_dgt, inds, row_local, base = base)
            row = index(row_dgt, base = base)
            full[row, col] += val
        end
    end

    full
end

kron_embed(mat::AbstractMatrix, inds::AbstractVector{<:Integer}, L::Integer, base::Integer = 2) =
    dense_local_embedding(mat, inds, L; base)

function assert_mul_variants_match_dense(H, x::AbstractVector; atol=1e-10, rtol=1e-10)
    assert_mul_variants_match_dense(H, x, Array(H); atol, rtol)
end

function assert_mul_variants_match_dense(H, x::AbstractVector, dense::AbstractMatrix; atol=1e-10, rtol=1e-10)
    expected = dense * x
    @test H * x ≈ expected atol = atol rtol = rtol
    @test EDKit.mul(H, x) ≈ expected atol = atol rtol = rtol

    accumulated = ones(ComplexF64, size(H, 1))
    mul!(accumulated, H, x)
    @test accumulated ≈ ones(ComplexF64, size(H, 1)) + expected atol = atol rtol = rtol

    scaled = ones(ComplexF64, size(H, 1))
    mul!(scaled, H, x, 2, 3)
    @test scaled ≈ 2 * expected .+ 3 atol = atol rtol = rtol
    expected
end

function assert_mul_variants_match_dense(H, x::AbstractMatrix; atol=1e-10, rtol=1e-10)
    assert_mul_variants_match_dense(H, x, Array(H); atol, rtol)
end

function assert_mul_variants_match_dense(H, x::AbstractMatrix, dense::AbstractMatrix; atol=1e-10, rtol=1e-10)
    expected = dense * x
    @test H * x ≈ expected atol = atol rtol = rtol
    @test EDKit.mul(H, x) ≈ expected atol = atol rtol = rtol

    accumulated = ones(ComplexF64, size(H, 1), size(x, 2))
    mul!(accumulated, H, x)
    @test accumulated ≈ ones(ComplexF64, size(H, 1), size(x, 2)) + expected atol = atol rtol = rtol

    scaled = ones(ComplexF64, size(H, 1), size(x, 2))
    mul!(scaled, H, x, 2, 3)
    @test scaled ≈ 2 * expected .+ 3 atol = atol rtol = rtol
    expected
end

function assert_operator_matches_dense(
    H;
    test_sparse::Bool=true,
    test_cached::Bool=true,
    rng=stable_rng(7),
    atol=1e-10,
    rtol=1e-10,
)
    dense = Array(H)
    @test size(dense) == size(H)

    if test_sparse
        @test sparse(H) ≈ sparse(dense) atol = atol rtol = rtol
    end

    v = random_vector(size(H, 2); rng, T=ComplexF64)
    assert_mul_variants_match_dense(H, v, dense; atol, rtol)

    M = random_matrix(size(H, 2), 3; rng, T=ComplexF64)
    assert_mul_variants_match_dense(H, M, dense; atol, rtol)

    if test_cached
        clear_sparse_cache!()
        @test sparse!(H) ≈ sparse(dense) atol = atol rtol = rtol
        @test H * M ≈ dense * M atol = atol rtol = rtol
        @test EDKit.mul(H, M) ≈ dense * M atol = atol rtol = rtol
        clear_sparse_cache!()
    end

    dense
end

function assert_threaded_basis_matches_serial(make_basis)
    serial = make_basis(false)
    threaded = make_basis(true)

    @test size(threaded) == size(serial)
    @test threaded.I == serial.I
    if hasproperty(threaded, :R)
        @test threaded.R ≈ serial.R
    end

    threaded
end
