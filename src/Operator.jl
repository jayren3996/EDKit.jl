"""
Operator

Construction and application of EDKit many-body operators.

This file defines the central `Operator` abstraction together with the helpers
that:
- build it from local terms,
- apply it matrix-free to vectors and matrices,
- convert it into dense or sparse explicit matrices,
- and create common local operators such as Pauli/spin products.
"""
#---------------------------------------------------------------------------------------------------
export operator, trans_inv_operator, sparse!, clear_sparse_cache!
"""
    Operator{Tv}

Many-body operator stored as a sum of local terms on a chosen basis.

Each term is represented by a sparse local matrix `M[n]` acting on lattice sites
`I[n]`, while `B` controls how those local actions are embedded into the full
Hilbert space.
"""
struct Operator{Tv<:Number, Tb<:AbstractBasis}
    M::Vector{SparseMatrixCSC{Tv, Int}}
    I::Vector{Vector{Int}}
    B::Tb
end
#---------------------------------------------------------------------------------------------------
eltype(opt::Operator{Tv, Tb}) where Tv where Tb = promote_type(Tv, eltype(opt.B))
length(opt::Operator) = length(opt.M)
size(opt::Operator) = size(opt.B)
size(opt::Operator, i::Integer) = size(opt.B, i)
function Base.display(opt::Operator)
    println("Operator of size $(size(opt)) with $(length(opt)) terms.")
end
#---------------------------------------------------------------------------------------------------
"""
    operator(mats::AbstractVector{<:AbstractMatrix}, inds::AbstractVector{<:AbstractVector}, B::AbstractBasis)

Construct an [`Operator`](@ref) from local matrices and the sites they act on.

Arguments:
- `mats`: local operator matrices.
- `inds`: site lists matching `mats`, using 1-based site labels.
- `B`: basis in which the operator should act.

Repeated site patterns are merged automatically by summing their local matrices.

Examples:
```julia
mat = spin((1.0, "xx"), (1.0, "yy"))
H = operator(fill(mat, L), [[i, mod(i, L) + 1] for i in 1:L], TensorBasis(L=L))
```
"""
function operator(mats::AbstractVector{<:AbstractMatrix}, inds::AbstractVector{<:AbstractVector}, B::AbstractBasis)
    num = length(mats)
    @assert num == length(inds) "Numbers mismatch: $num matrices and $(length(inds)) indices."
    dtype = promote_type(eltype.(mats)...)
    M = Vector{SparseMatrixCSC{dtype, Int64}}(undef, num)
    I = Vector{Vector{Int64}}(undef, num)
    N = 0
    for i = 1:num
        iszero(mats[i]) && continue
        ind = inds[i]
        pos = findfirst(x -> isequal(x, ind), view(I, 1:N))
        if isnothing(pos)
            N += 1
            I[N] = ind
            M[N] = sparse(mats[i])
        else
            M[pos] += mats[i]
        end
    end
    deleteat!(M, N+1:num)
    deleteat!(I, N+1:num)
    Operator(M, I, B)
end

function operator(mats::AbstractVector{<:AbstractMatrix}, inds::AbstractVector{<:AbstractVector}, L::Integer)
    base = find_base(size(mats[1], 1), length(inds[1]))
    basis = TensorBasis(L=L, base=base)
    operator(mats, inds, basis)
end
operator(mats::AbstractVector{<:AbstractMatrix}, inds::AbstractVector{<:Integer}, C) = operator(mats, [[i] for i in inds], C)
operator(mats::AbstractVector{<:AbstractMatrix}, L::Integer) = operator(mats, [[i] for i in 1:L], L)
operator(mats::AbstractVector{<:AbstractMatrix}, B::AbstractBasis) = operator(mats, [[i] for i in 1:length(B.dgt)], B)
operator(mat::AbstractMatrix, ind::AbstractVector{<:Integer}, C) = operator([mat], [ind], C)
#---------------------------------------------------------------------------------------------------
"""
    trans_inv_operator(mat::AbstractMatrix, ind, basis_or_length)

Construct a translationally invariant operator by translating one local term
around the lattice.

`ind` specifies the support of the seed term, for example `1:2` for a nearest-
neighbor bond term. The function then generates all translated copies with
periodic boundary conditions and sums them into a single [`Operator`](@ref).

This accepts either a concrete basis or just a system size `L`.
"""
function trans_inv_operator(mat::AbstractMatrix, ind::AbstractVector{<:Integer}, B::AbstractBasis)
    L = length(B.dgt)
    smat = sparse(mat)
    mats = fill(smat, L)
    inds = [mod.(ind .+ i, L) .+ 1 for i = -1:L-2]
    operator(mats, inds, B)
end

function trans_inv_operator(mat::AbstractMatrix, ind::AbstractVector{<:Integer}, L::Integer)
    base = find_base(size(mat, 1), length(ind))
    B = TensorBasis(L=L, base=base)
    trans_inv_operator(mat, ind, B)
end

trans_inv_operator(mat::AbstractMatrix, M::Integer, C) = trans_inv_operator(mat, 1:M, C)
#---------------------------------------------------------------------------------------------------
*(c::Number, opt::Operator) = operator(c .* opt.M, opt.I, opt.B)
*(opt::Operator, c::Number) = c * opt
/(opt::Operator, c::Number) = operator(opt.M ./ c, opt.I, opt.B)
function +(opt1::Operator, opt2::Operator)
    mats = vcat(opt1.M, opt2.M)
    inds = vcat(opt1.I, opt2.I)
    operator(mats, inds, opt1.B)
end
+(opt1::Operator, ::Nothing) = opt1
+(::Nothing, opt1::Operator) = opt1
-(opt::Operator) = Operator(-opt.M, opt.I, opt.B)
-(opt1::Operator, opt2::Operator) = opt1 + (-opt2)
#---------------------------------------------------------------------------------------------------
function LinearAlgebra.adjoint(opt::Operator)
    M = adjoint.(opt.M)
    operator(M, opt.I, opt.B)
end
#---------------------------------------------------------------------------------------------------
"""
    find_base(a, b)

Infer the local on-site dimension from a local operator acting on `b` sites.

Arguments:
- `a`: matrix dimension of the local operator.
- `b`: number of sites the local operator acts on.

Returns:
- The local base `D` satisfying `D^b == a`.

An error is thrown when no integer base is compatible with the given local
operator dimension.
"""
function find_base(a::Integer, b::Integer)
    isone(b) && return a
    for i in 2:a
        i^b == a && return i
    end
    error("Incompatible dimension: ($a, $b).")
end

#---------------------------------------------------------------------------------------------------
# Operator to matrices
#---------------------------------------------------------------------------------------------------
export addto!

mutable struct SparseTripletAccumulator{Tv, Ti<:Integer}
    rows::Vector{Ti}
    cols::Vector{Ti}
    vals::Vector{Tv}
    col::Ti
end

function SparseTripletAccumulator(::Type{Tv}, ::Type{Ti}=Int; sizehint::Integer=0) where {Tv, Ti<:Integer}
    rows = Vector{Ti}()
    cols = Vector{Ti}()
    vals = Vector{Tv}()
    if sizehint > 0
        sizehint!(rows, sizehint)
        sizehint!(cols, sizehint)
        sizehint!(vals, sizehint)
    end
    SparseTripletAccumulator{Tv, Ti}(rows, cols, vals, zero(Ti))
end

function _sparse_triplet_sizehint(opt::Operator)
    nterms = 0
    for M in opt.M
        nterms = min(nterms + nnz(M), 1_000_000)
    end
    min(size(opt, 2) * nterms, 1_000_000)
end

"""
    addto!(M::AbstractMatrix, opt::Operator)

Accumulate the matrix representation of `opt` into the preallocated array `M`.

This is useful when you want control over the output array type or want to reuse
existing storage instead of calling `Array(opt)` or `sparse(opt)`.
"""
function addto!(M::AbstractMatrix, opt::Operator)
    if _has_tensor_base2_kernel(opt)
        for j = 1:size(opt.B, 2)
            colmn!(view(M, :, j), opt, j, 1, 1)
        end
    else
        dgt = similar(opt.B.dgt)
        basis_workspace = _basis_index_workspace(opt.B)
        for j = 1:size(opt.B, 2)
            colmn!(view(M, :, j), opt, j, dgt, 1, 1, basis_workspace)
        end
    end
    M
end
#---------------------------------------------------------------------------------------------------
"""
    Array(opt::Operator)

Materialize `opt` as a dense matrix.

Returns:
- A newly allocated dense matrix of size `size(opt)`.

This is convenient for exact diagonalization, but may be prohibitively
expensive for large Hilbert spaces.
"""
function Array(opt::Operator)
    M = zeros(eltype(opt), size(opt))
    if size(M, 1) > 0 && size(M, 2) > 0
        addto!(M, opt)
    end
    M
end

Matrix(opt::Operator) = Array(opt)
#---------------------------------------------------------------------------------------------------
"""
    sparse(opt::Operator)

Materialize `opt` as an explicit sparse matrix.

Returns:
- A newly allocated sparse matrix of size `size(opt)`.

Compared with `Array(opt)`, this preserves sparsity but still constructs the
full operator explicitly.
"""
function SparseArrays.sparse(opt::Operator)
    Tv = eltype(opt)
    Ti = Int
    m, n = size(opt)
    (m == 0 || n == 0) && return SparseArrays.sparse(Ti[], Ti[], Tv[], m, n)

    # Per-column work is much cheaper on the TensorBasis base-2 fast kernel
    # than on reduced bases, so the parallelization break-even is much higher
    # there. The min-cols-per-thread thresholds below were measured on a
    # Julia 1.12 / 8-thread workstation.
    min_cols_per_thread = _has_tensor_base2_kernel(opt) ? 2048 : 256
    nt = min(Threads.nthreads(), max(1, n ÷ min_cols_per_thread))
    if nt == 1
        acc = SparseTripletAccumulator(Tv, Ti; sizehint=_sparse_triplet_sizehint(opt))
        _accumulate_columns!(acc, opt, 1:n)
        return SparseArrays.sparse(acc.rows, acc.cols, acc.vals, m, n)
    end

    ni = dividerange(n, nt)
    per_thread_hint = max(1, _sparse_triplet_sizehint(opt) ÷ nt)
    accs = [SparseTripletAccumulator(Tv, Ti; sizehint=per_thread_hint) for _ in 1:nt]
    Threads.@threads for i in 1:nt
        _accumulate_columns!(accs[i], opt, ni[i])
    end
    rows, cols, vals = _merge_accumulators(accs)
    SparseArrays.sparse(rows, cols, vals, m, n)
end

function _accumulate_columns!(acc::SparseTripletAccumulator, opt::Operator, range)
    if _has_tensor_base2_kernel(opt)
        for j in range
            colmn!(acc, opt, j, 1, 1)
        end
    else
        dgt = similar(opt.B.dgt)
        ws = _basis_index_workspace(opt.B)
        for j in range
            colmn!(acc, opt, j, dgt, 1, 1, ws)
        end
    end
end

function _merge_accumulators(accs::Vector{SparseTripletAccumulator{Tv, Ti}}) where {Tv, Ti}
    total = sum(length(a.rows) for a in accs)
    rows = Vector{Ti}(undef, total)
    cols = Vector{Ti}(undef, total)
    vals = Vector{Tv}(undef, total)
    offset = 0
    @inbounds for a in accs
        k = length(a.rows)
        copyto!(rows, offset + 1, a.rows, 1, k)
        copyto!(cols, offset + 1, a.cols, 1, k)
        copyto!(vals, offset + 1, a.vals, 1, k)
        offset += k
    end
    rows, cols, vals
end
#---------------------------------------------------------------------------------------------------
# Sparse-matrix cache for accelerated matrix multiplication
#---------------------------------------------------------------------------------------------------
const _SPARSE_CACHE = LRU{UInt, Any}(maxsize=8)

"""
    sparse!(opt::Operator) -> SparseMatrixCSC

Pre-compute and cache the sparse matrix representation of `opt`.

After calling `sparse!(opt)`, subsequent `opt * m` and `mul(opt, m)` calls
for **matrix** inputs will automatically use the cached sparse matrix via
SparseArrays SpMM, which is typically **10–1000× faster** than the default
matrix-free path.  Vector multiplication (`opt * v`) is unaffected and
always uses the matrix-free kernel.

The cache holds up to 8 operators (LRU eviction).  Call
[`clear_sparse_cache!`](@ref) when you no longer need the cached matrices
and want to reclaim memory.

!!! note "Snapshot semantics"
    `sparse!` snapshots the operator's contents at call time and keys the cache
    by object identity. Every public transformation (`c*opt`, `opt1+opt2`,
    `adjoint`, …) returns a fresh `Operator`, so this is transparent. The one
    exception is mutating an operator's stored matrices *in place* after
    `sparse!`: that is not tracked, and subsequent `*`/`mul`/`mul!` will use the
    stale cached matrix. Call `sparse!(opt)` again to refresh it.

# When to use

Call `sparse!` when you plan to multiply the same operator by a matrix more
than once (e.g. applying H to a block of eigenvectors, or inside an
iterative solver).  The one-time cost of building the sparse matrix is
amortized over all subsequent multiplications.

# Memory cost

The cached sparse matrix stores one `ComplexF64` and one `Int` per
structural nonzero, plus one `Int` per column.  For a typical Heisenberg
chain at half filling:

| System   | Basis dim | nnz     | Cache size |
|----------|-----------|---------|------------|
| L = 16   | 12 870    | 122 694 | ≈ 3 MiB    |
| L = 20   | 184 756   | 2.1 M   | ≈ 50 MiB   |

For very large systems where this overhead is prohibitive, skip `sparse!`
and rely on the matrix-free path instead.

# Example

```julia
H = trans_inv_operator(spin("xx","yy","zz"), 1:2, basis)

sparse!(H)             # one-time build + cache
result = H * states    # uses fast SpMM automatically
clear_sparse_cache!()  # free memory when done
```
"""
function sparse!(opt::Operator)
    S = sparse(opt)
    _SPARSE_CACHE[objectid(opt)] = S
    S
end

"""
    clear_sparse_cache!()

Release all cached sparse matrices created by [`sparse!`](@ref).

Call this after you are done with cached operator multiplication to free the
memory occupied by the sparse representations.
"""
clear_sparse_cache!() = empty!(_SPARSE_CACHE)

_cached_sparse(opt::Operator) = get(_SPARSE_CACHE, objectid(opt), nothing)
#---------------------------------------------------------------------------------------------------
LinearAlgebra.Hermitian(opt::Operator) = Array(opt) |> Hermitian
LinearAlgebra.Symmetric(opt::Operator) = Array(opt) |> Symmetric
LinearAlgebra.eigen(opt::Operator) = Array(opt) |> eigen
LinearAlgebra.eigvals(opt::Operator) = Array(opt) |>eigvals
LinearAlgebra.svd(opt::Operator) = Array(opt) |> svd
LinearAlgebra.svdvals(opt::Operator) = Array(opt) |> svdvals
#---------------------------------------------------------------------------------------------------
"""
    mul!(target, opt::Operator, v::AbstractVector)
    mul!(target, opt::Operator, v::AbstractVector, α, β)

Accumulate `opt * v` into the preallocated vector `target`.

The three-argument EDKit method intentionally accumulates into `target` because
the matrix-free column kernels are additive. The five-argument method follows
the standard `LinearAlgebra.mul!` convention `target = α * opt * v + β * target`.

Arguments:
- `target`: output buffer of length `size(opt, 1)`.
- `opt`: operator to apply.
- `v`: input coordinate vector.

Returns:
- The mutated `target`.

This is the single-threaded in-place application path underlying `opt * v`.
"""
@inline _apply_column_base2!(target, opt::Operator, j, rhs::AbstractVector, α) =
    colmn!(target, opt, j, α * rhs[j], 1)

@inline _apply_column_generic!(target, opt::Operator, j, rhs::AbstractVector, α, dgt, ws) =
    colmn!(target, opt, j, dgt, α * rhs[j], 1, ws)

@inline _apply_column_base2!(target, opt::Operator, j, rhs::AbstractMatrix, α) =
    _colmn_row!(target, opt, j, rhs, j, α)

@inline _apply_column_generic!(target, opt::Operator, j, rhs::AbstractMatrix, α, dgt, ws) =
    _colmn_row!(target, opt, j, dgt, rhs, j, α, ws)

function _apply_columns!(target, opt::Operator, range, rhs, α::Number=1)
    if _has_tensor_base2_kernel(opt)
        for j in range
            _apply_column_base2!(target, opt, j, rhs, α)
        end
    else
        dgt = similar(opt.B.dgt)
        ws = _basis_index_workspace(opt.B)
        for j in range
            _apply_column_generic!(target, opt, j, rhs, α, dgt, ws)
        end
    end
    target
end

mul!(target::AbstractVector, opt::Operator, v::AbstractVector) =
    _apply_columns!(target, opt, eachindex(v), v)

function mul!(target::AbstractVector, opt::Operator, v::AbstractVector, α::Number, β::Number)
    iszero(β) ? fill!(target, zero(eltype(target))) : (target .*= β)
    iszero(α) && return target
    _apply_columns!(target, opt, eachindex(v), v, α)
end

"""
    mul!(target, opt::Operator, m::AbstractMatrix)
    mul!(target, opt::Operator, m::AbstractMatrix, α, β)

Accumulate `opt * m` into the preallocated matrix `target`.

As for vectors, the three-argument method accumulates while the five-argument
method follows `target = α * opt * m + β * target`.
"""
function mul!(target::AbstractMatrix, opt::Operator, m::AbstractMatrix)
    S = _cached_sparse(opt)
    S !== nothing && return mul!(target, S, m, true, true)   # cached SpMM, accumulating (operator-2)
    _apply_columns!(target, opt, axes(m, 1), m)
end

function mul!(target::AbstractMatrix, opt::Operator, m::AbstractMatrix, α::Number, β::Number)
    S = _cached_sparse(opt)
    S !== nothing && return mul!(target, S, m, α, β)         # cached SpMM (operator-2)
    iszero(β) ? fill!(target, zero(eltype(target))) : (target .*= β)
    iszero(α) && return target
    _apply_columns!(target, opt, axes(m, 1), m, α)
end

export mul
"""
    mul(opt::Operator, v)
    mul(opt::Operator, m)

Multi-threaded multiplication of an [`Operator`](@ref) with a vector or matrix.

Unlike `opt * v`, this version parallelizes over input columns / basis states
using `Threads.@threads`.
"""
function mul(opt::Operator, v::AbstractVector)
    ctype = promote_type(eltype(opt), eltype(v))
    nt = min(Threads.nthreads(), max(1, length(v)))
    ni = dividerange(length(v), nt)
    Ms = [zeros(ctype, size(opt, 1)) for i in 1:nt]
    Threads.@threads for i in 1:nt
        _apply_columns!(Ms[i], opt, ni[i], v)
    end
    target = Ms[1]
    @inbounds for i in 2:nt
        axpy!(one(ctype), Ms[i], target)
    end
    target
end

function mul(opt::Operator, m::AbstractMatrix)
    ctype = promote_type(eltype(opt), eltype(m))
    S = _cached_sparse(opt)
    if S !== nothing
        return convert(Matrix{ctype}, S * m)
    end
    nt = min(Threads.nthreads(), max(1, size(m, 1)))
    ni = dividerange(size(m,1), nt)
    Ms = [zeros(ctype, size(opt, 1), size(m, 2)) for i in 1:nt]
    Threads.@threads for i in 1:nt
        _apply_columns!(Ms[i], opt, ni[i], m)
    end
    target = Ms[1]
    @inbounds for i in 2:nt
        axpy!(one(ctype), Ms[i], target)
    end
    target
end

function *(opt::Operator, v::AbstractVector)
    ctype = promote_type(eltype(opt), eltype(v))
    target = zeros(ctype, size(opt, 1))
    _apply_columns!(target, opt, eachindex(v), v)
end

function *(opt::Operator, m::AbstractMatrix)
    ctype = promote_type(eltype(opt), eltype(m))
    S = _cached_sparse(opt)
    if S !== nothing
        return convert(Matrix{ctype}, S * m)
    end
    target = zeros(ctype, size(opt, 1), size(m, 2))
    _apply_columns!(target, opt, axes(m, 1), m)
end

#---------------------------------------------------------------------------------------------------
# Helper functions
#---------------------------------------------------------------------------------------------------
@inline _accumulate!(target::AbstractVector, pos, C, val, coeff) =
    _accumulate!(target, pos, C, val, coeff, 1)

@inline _accumulate!(target::AbstractVector, pos, C, val, coeff, scale) =
    target[pos] += scale * coeff * C * val

@inline function _accumulate!(target::AbstractMatrix, pos, C, val, coeff)
    _accumulate!(target, pos, C, val, coeff, 1)
end

@inline function _accumulate!(target::AbstractMatrix, pos, C, val, coeff, scale)
    cv = scale * C * val
    @inbounds for k in axes(target, 2)
        target[pos, k] += cv * coeff[k]
    end
end

@inline function _accumulate!(target::SparseTripletAccumulator, pos, C, val, coeff, scale)
    iszero(C) && return nothing
    final_val = scale * coeff * C * val
    iszero(final_val) && return nothing
    push!(target.rows, pos)
    push!(target.cols, target.col)
    push!(target.vals, final_val)
    nothing
end

@inline function _accumulate_row!(target::AbstractMatrix, pos, C, val, m::AbstractMatrix, row, scale)
    cv = scale * C * val
    @inbounds for k in axes(target, 2)
        target[pos, k] += cv * m[row, k]
    end
end

@inline _basis_index_workspace(::AbstractBasis) = nothing
@inline _basis_index_workspace(b::AbelianBasis) = _shallow_workspace(b.G)

@inline _index_in_basis(b::AbstractBasis, dgt::AbstractVector, workspace) = index(b, dgt)
@inline _index_in_basis(b::AbelianBasis, dgt::AbstractVector, workspace::AbelianOperator) =
    index(b, dgt, workspace)

@inline _has_tensor_base2_kernel(opt::Operator) = _has_tensor_base2_kernel(opt.B)
@inline _has_tensor_base2_kernel(::AbstractBasis) = false
@inline _has_tensor_base2_kernel(b::TensorBasis) = b.B == 2

@inline function _tensor_base2_local_index(state::Integer, I::Vector{Int}, L::Integer)
    local_index = 0
    @inbounds for site in I
        local_index = (local_index << 1) | ((state >> (L - site)) & 1)
    end
    local_index + 1
end

@inline function _tensor_base2_replace(state::Integer, I::Vector{Int}, row::Integer, L::Integer)
    row_state = state
    local_row = row - 1
    @inbounds for p = length(I):-1:1
        site = I[p]
        mask = 1 << (L - site)
        if iszero(local_row & 1)
            row_state &= ~mask
        else
            row_state |= mask
        end
        local_row >>= 1
    end
    row_state + 1
end
#---------------------------------------------------------------------------------------------------
"""
    colmn!(target::AbstractVecOrMat, M::SparseMatrixCSC, I::Vector{Int}, b::AbstractBasis, coeff=1)

Internal helper used to apply one local term to the current basis state in `b`.

`b.dgt` is interpreted as the current many-body basis configuration. The result
of acting with the local matrix `M` on sites `I` is accumulated in `target`.
"""
function colmn!(target::AbstractVecOrMat, M::SparseMatrixCSC, I::Vector{Int}, b::AbstractBasis, coeff=1)
    colmn!(target, M, I, b, b.dgt, coeff)
end
"""
    colmn!(target, M, I, b, dgt, coeff)

Thread-safe variant of the local-term application that uses the supplied digit
buffer `dgt` instead of `b.dgt`.
"""
function colmn!(
    target::AbstractVecOrMat,
    M::SparseMatrixCSC,
    I::Vector{Int},
    b::AbstractBasis,
    dgt::AbstractVector,
    coeff=1,
    scale=1,
    basis_workspace=_basis_index_workspace(b),
)
    rows, vals = rowvals(M), nonzeros(M)
    j = index(dgt, I, base=b.B)
    change = false
    @inbounds for i in nzrange(M, j)
        row, val = rows[i], vals[i]
        change!(dgt, I, row, base=b.B)
        C, pos = _index_in_basis(b, dgt, basis_workspace)
        _accumulate!(target, pos, C, val, coeff, scale)
        change = true
    end
    change && change!(dgt, I, j, base=b.B)
    nothing
end

function colmn!(
    target::SparseTripletAccumulator,
    M::SparseMatrixCSC,
    I::Vector{Int},
    b::AbstractBasis,
    dgt::AbstractVector,
    coeff=1,
    scale=1,
    basis_workspace=_basis_index_workspace(b),
)
    rows, vals = rowvals(M), nonzeros(M)
    j = index(dgt, I, base=b.B)
    change = false
    @inbounds for i in nzrange(M, j)
        row, val = rows[i], vals[i]
        change!(dgt, I, row, base=b.B)
        C, pos = _index_in_basis(b, dgt, basis_workspace)
        _accumulate!(target, pos, C, val, coeff, scale)
        change = true
    end
    change && change!(dgt, I, j, base=b.B)
    nothing
end

function _colmn_row!(
    target::AbstractMatrix,
    M::SparseMatrixCSC,
    I::Vector{Int},
    b::AbstractBasis,
    dgt::AbstractVector,
    m::AbstractMatrix,
    row,
    scale=1,
    basis_workspace=_basis_index_workspace(b),
)
    rows, vals = rowvals(M), nonzeros(M)
    j = index(dgt, I, base=b.B)
    change = false
    @inbounds for i in nzrange(M, j)
        local_row, val = rows[i], vals[i]
        change!(dgt, I, local_row, base=b.B)
        C, pos = _index_in_basis(b, dgt, basis_workspace)
        _accumulate_row!(target, pos, C, val, m, row, scale)
        change = true
    end
    change && change!(dgt, I, j, base=b.B)
    nothing
end

function _colmn_row!(
    target::AbstractMatrix,
    M::SparseMatrixCSC,
    I::Vector{Int},
    b::TensorBasis,
    state::Integer,
    m::AbstractMatrix,
    row,
    scale=1,
)
    rows, vals = rowvals(M), nonzeros(M)
    L = length(b.dgt)
    j = _tensor_base2_local_index(state, I, L)
    @inbounds for i in nzrange(M, j)
        pos = _tensor_base2_replace(state, I, rows[i], L)
        _accumulate_row!(target, pos, 1, vals[i], m, row, scale)
    end
    nothing
end

function colmn!(
    target::Union{AbstractVector,SparseTripletAccumulator},
    M::SparseMatrixCSC,
    I::Vector{Int},
    b::TensorBasis,
    state::Integer,
    coeff=1,
    scale=1,
)
    rows, vals = rowvals(M), nonzeros(M)
    L = length(b.dgt)
    j = _tensor_base2_local_index(state, I, L)
    @inbounds for i in nzrange(M, j)
        pos = _tensor_base2_replace(state, I, rows[i], L)
        _accumulate!(target, pos, 1, vals[i], coeff, scale)
    end
    nothing
end
#---------------------------------------------------------------------------------------------------
"""
    colmn!(target, opt::Operator, j, coeff=1)

Apply column `j` of `opt` to `target`.

This is an internal helper used to build explicit matrix columns and to apply an
operator matrix-free. It assumes `target` has already been allocated and will
mutate the basis working buffer stored in `opt.B`.
"""
function colmn!(target::AbstractVecOrMat, opt::Operator, j::Integer, coeff=1)
    colmn!(target, opt, j, opt.B.dgt, coeff)
end
"""
    colmn!(target, opt, j, dgt, coeff)

Thread-safe variant that uses the supplied digit buffer `dgt`.
"""
function colmn!(
    target::AbstractVecOrMat,
    opt::Operator,
    j::Integer,
    dgt::AbstractVector,
    coeff=1,
    scale=1,
    basis_workspace=_basis_index_workspace(opt.B),
)
    b, M, I = opt.B, opt.M, opt.I
    r = change!(b, j, dgt)
    scaled = isone(r) ? scale : scale / r
    for i = 1:length(M)
        colmn!(target, M[i], I[i], b, dgt, coeff, scaled, basis_workspace)
    end
end

function colmn!(
    target::SparseTripletAccumulator,
    opt::Operator,
    j::Integer,
    dgt::AbstractVector,
    coeff=1,
    scale=1,
    basis_workspace=_basis_index_workspace(opt.B),
)
    b, M, I = opt.B, opt.M, opt.I
    target.col = j
    r = change!(b, j, dgt)
    scaled = isone(r) ? scale : scale / r
    for i = 1:length(M)
        colmn!(target, M[i], I[i], b, dgt, coeff, scaled, basis_workspace)
    end
end

function colmn!(
    target::Union{AbstractVector,SparseTripletAccumulator},
    opt::Operator{Tv,<:TensorBasis},
    j::Integer,
    coeff=1,
    scale=1,
) where {Tv}
    b, M, I = opt.B, opt.M, opt.I
    if b.B != 2
        dgt = similar(b.dgt)
        return colmn!(target, opt, j, dgt, coeff, scale, nothing)
    end
    target isa SparseTripletAccumulator && (target.col = j)
    state = j - 1
    for i = 1:length(M)
        colmn!(target, M[i], I[i], b, state, coeff, scale)
    end
    nothing
end

function _colmn_row!(
    target::AbstractMatrix,
    opt::Operator,
    j::Integer,
    dgt::AbstractVector,
    m::AbstractMatrix,
    row,
    scale=1,
    basis_workspace=_basis_index_workspace(opt.B),
)
    b, M, I = opt.B, opt.M, opt.I
    r = change!(b, j, dgt)
    scaled = isone(r) ? scale : scale / r
    for i = 1:length(M)
        _colmn_row!(target, M[i], I[i], b, dgt, m, row, scaled, basis_workspace)
    end
end

function _colmn_row!(
    target::AbstractMatrix,
    opt::Operator{Tv,<:TensorBasis},
    j::Integer,
    m::AbstractMatrix,
    row,
    scale=1,
) where {Tv}
    b, M, I = opt.B, opt.M, opt.I
    if b.B != 2
        dgt = similar(b.dgt)
        return _colmn_row!(target, opt, j, dgt, m, row, scale, nothing)
    end
    state = j - 1
    for i = 1:length(M)
        _colmn_row!(target, M[i], I[i], b, state, m, row, scale)
    end
    nothing
end

#---------------------------------------------------------------------------------------------------
# Spin Matrices
#---------------------------------------------------------------------------------------------------
const PAULI = sparse.(
    [ I(2), [0 1; 1 0], [0 -1im; 1im 0], [1 0; 0 -1] ]
)
const GELLMANN = sparse.([
    I(3),
    [0 1 0; 1 0 0; 0 0 0], 
    [0 -1im 0; 1im 0 0; 0 0 0], 
    [1 0 0; 0 -1 0; 0 0 0], 
    [0 0 1; 0 0 0; 1 0 0], 
    [0 0 -1im; 0 0 0; 1im 0 0], 
    [0 0 0; 0 0 1; 0 1 0],
    [0 0 0; 0 0 -1im; 0 1im 0], 
    [1 0 0; 0 1 0; 0 0 -2] / sqrt(3),
    [1 0 0; 0 1 0; 0 0 1]
])

"""
    σ(i::Integer)
    σ(Is::Integer...)
    σ(Is::Union{<:Tuple, <:AbstractVector})

Return Pauli-basis matrices indexed by integers.

Conventions:
- `σ(0)` is the `2 × 2` identity,
- `σ(1)`, `σ(2)`, `σ(3)` are the usual Pauli `X`, `Y`, `Z` matrices.

When multiple indices are supplied, the result is the Kronecker product in the
same order. This helper is mainly useful for compact low-level construction and
testing; `spin("X")`, `spin("XX")`, and related string interfaces are usually
more convenient for user-facing code.
"""
σ(i::Integer) = PAULI[i+1]
σ(Is::Integer...) = kron((σ(i) for i in Is)...)
σ(Is::Union{<:Tuple, <:AbstractVector}) = σ(Is...)

"""
    λ(i::Integer)
    λ(Is::Integer...)
    λ(Is::Union{<:Tuple, <:AbstractVector})

Return generalized Gell-Mann matrices indexed by integers.

Conventions:
- `λ(0)` is the `3 × 3` identity,
- `λ(1)` through `λ(8)` are the standard traceless Gell-Mann generators.

When several indices are passed, EDKit returns their Kronecker product in order.
This is the qutrit analogue of [`σ`](@ref).
"""
λ(i::Integer) = GELLMANN[i+1]
λ(Is::Integer...) = kron((λ(i) for i in Is)...)
λ(Is::Union{<:Tuple, <:AbstractVector}) = λ(Is...)
#---------------------------------------------------------------------------------------------------
export spin
"""
    spin_coeff(D::Integer)

Return the ladder-operator coefficients for a spin representation of local
dimension `D`.

The returned vector has length `D - 1` and stores the coefficients
`sqrt(i * (D - i))` used to assemble `S+`, `S-`, `Sx`, and `iSy` in EDKit's
local-spin convention. This helper is internal to the `spin_*` constructors but
is documented because it captures the normalization shared by the whole family.
"""
spin_coeff(D::Integer) = [sqrt(i*(D-i)) for i = 1:D-1]

"""
    spin_Sp(D::Integer)

Return the local raising operator `S+` for a `D`-dimensional spin space.

The matrix is sparse and uses the same coefficient convention as
[`spin_coeff`](@ref). For spin-`1/2`, this is the usual two-level raising
operator.
"""
spin_Sp(D::Integer) = sparse(1:D-1, 2:D, spin_coeff(D), D, D)

"""
    spin_Sm(D::Integer)

Return the local lowering operator `S-` for a `D`-dimensional spin space.

This is the Hermitian-adjoint partner of [`spin_Sp`](@ref) in EDKit's local
spin convention.
"""
spin_Sm(D::Integer) = sparse(2:D, 1:D-1, spin_coeff(D), D, D)

"""
    spin_Sx(D::Integer)

Return the local `Sx` operator for a `D`-dimensional spin space.

EDKit assembles this as `(S+ + S-) / 2`, using the coefficient convention from
[`spin_coeff`](@ref).
"""
function spin_Sx(D::Integer)
    coeff = spin_coeff(D) / 2
    sp = sparse(1:D-1, 2:D, coeff, D, D)
    sp + sp'
end

"""
    spin_iSy(D::Integer)

Return the antisymmetric matrix corresponding to `i Sy` for a local spin space
of dimension `D`.

This is intentionally `i Sy` rather than `Sy` itself, so that lowercase `"y"`
spin strings can remain real-valued in EDKit's many-body construction
convention. The final complex phase, when needed, is restored by [`spin`](@ref)
through `spin_ny`.
"""
function spin_iSy(D::Integer)
    coeff = spin_coeff(D) / 2
    sp = sparse(1:D-1, 2:D, coeff, D, D)
    sp - sp'
end

"""
    spin_Sz(D::Integer)

Return the diagonal `Sz` operator for a `D`-dimensional spin space.

The diagonal entries run from `J` down to `-J` with `J = (D - 1) / 2`, matching
the standard spin-`J` representation.
"""
function spin_Sz(D::Integer)
    J = (D-1) / 2
    sparse(1:D, 1:D, J:-1:-J)
end
#---------------------------------------------------------------------------------------------------
"""
    spin_dict(c::Char, D::Integer)

Translate one character from EDKit's spin-string syntax into a local sparse
operator matrix.

Supported symbols include:
- lowercase `"x"`, `"y"`, `"z"`, `"+"`, `"-"` for generic spin representations,
- `"1"` or `"I"` for the local identity,
- uppercase `"X"`, `"Y"`, `"Z"` for Pauli-style spin-`1/2` operators only.

This is an internal parser used by [`spin_product`](@ref). It throws when a
symbol is invalid for the requested local dimension.
"""
function spin_dict(c::Char, D::Integer)
    isequal(c, '+') && return spin_Sp(D)
    isequal(c, '-') && return spin_Sm(D)
    isequal(c, 'x') && return spin_Sx(D)
    isequal(c, 'y') && return spin_iSy(D)
    isequal(c, 'z') && return spin_Sz(D)
    isequal(c, '1') && return spdiagm(ones(D))
    isequal(c, 'I') && return spdiagm(ones(D))

    @assert isequal(D, 2) "Invalid spin symbol: $c."
    isequal(c, 'X') && return sparse([0 1; 1  0])
    isequal(c, 'Y') && return sparse([0 1;-1  0])
    isequal(c, 'Z') && return sparse([1 0; 0 -1])
    error("Invalid spin symbol: $c.")
end
#---------------------------------------------------------------------------------------------------
"""
    spin_ny(s::AbstractString, D::Integer)

Count how many `y`-type factors appear in the symbolic spin string `s`.

For `D = 2`, both lowercase `"y"` and uppercase `"Y"` count. For higher local
dimension only lowercase `"y"` is meaningful. This count is used by
[`spin`](@ref) to restore the correct phase relating the real-valued internal
`i Sy` convention to the requested physical operator.
"""
@inline function spin_ny(s::AbstractString, D::Integer)
    n = 0
    if isequal(D, 2)
        for c in s
            n += isequal(c, 'y') || isequal(c, 'Y')
        end
    else
        for c in s
            n += isequal(c, 'y')
        end
    end
    n
end

"""
    spin_product(s::AbstractString, D::Integer)

Form the raw Kronecker product associated with the symbolic spin string `s`.

Each character is converted by [`spin_dict`](@ref), then combined left-to-right
with `kron`. The result does not yet include the phase correction associated
with repeated `y` factors; [`spin`](@ref) applies that correction afterwards and
caches the final matrix.
"""
function spin_product(s::AbstractString, D::Integer)
    isempty(s) && error("Empty spin string.")
    mat = spin_dict(first(s), D)
    for c in Iterators.drop(s, 1)
        mat = kron(mat, spin_dict(c, D))
    end
    mat
end

const SPIN_CACHE = LRU{Tuple{Int, String}, SparseMatrixCSC{ComplexF64, Int}}(maxsize=256)
#---------------------------------------------------------------------------------------------------
"""
    spin(s::String; D::Integer=2)

Return the local operator encoded by the string `s`.

For `D = 2`, uppercase Pauli-style labels such as `"X"`, `"Y"`, `"Z"` are
supported, as well as lowercase spin-operator labels such as `"x"`, `"y"`,
`"z"`, `"+"`, `"-"`, and `"1"`. Multi-site strings like `"xx"` or `"x1z"`
build Kronecker products in the given order.

**Case convention** (D = 2):
- Uppercase `"X"`, `"Y"`, `"Z"` are the Pauli matrices, e.g. `spin("X") = [0 1; 1 0]`.
- Lowercase `"x"`, `"y"`, `"z"` are the spin-1/2 operators, e.g. `spin("x") = [0 1/2; 1/2 0]`.
  Equivalently, `spin("X") == 2 * spin("x")`.

For higher local dimension `D > 2` only the lowercase spin-operator labels are
meaningful (uppercase Pauli-style labels are rejected with an error).

Returns:
- A sparse `SparseMatrixCSC{ComplexF64, Int}` acting on `length(s)` sites.
"""
function spin(s::String; D::Integer=2)
    key = (Int(D), s)
    mat = get!(SPIN_CACHE, key) do
        ny = spin_ny(s, D)
        raw = spin_product(s, D)
        sign = iszero(mod(ny, 2)) ? (-1)^(ny÷2) : (-1im)^ny
        convert(SparseMatrixCSC{ComplexF64, Int}, sign * raw)
    end
    copy(mat)
end
spin(c::Number, s::String; D::Integer=2) = c * spin(s, D=D)
#---------------------------------------------------------------------------------------------------
"""
    spin(spins::Tuple{<:Number, String}...; D::Integer=2)

Build a linear combination of local operators.

Each argument should be a pair `(coefficient, "operator_string")`, for example
`spin((1.0, "xx"), (1.0, "yy"), (0.5, "zz"))`.

Returns:
- The sparse linear combination of all requested local operator strings.
"""
function spin(spins::Tuple{<:Number, String}...; D::Integer=2)
    sum(ci * spin(si, D=D) for (ci, si) in spins)
end

spin(spins::AbstractVector{<:Tuple{<:Number, String}}; D::Integer=2) = spin(spins..., D=D)
#---------------------------------------------------------------------------------------------------
"""
    operator(s::String, inds, basis)

Convenience wrapper around [`spin`](@ref) and [`operator`](@ref) that builds a
local operator matrix from a symbolic string before embedding it in the many-body
basis.
"""
function operator(s::String, inds::AbstractVector{<:Integer}, basis::AbstractBasis)
    mat = spin(s, D=basis.B)
    operator(mat, inds, basis)
end
#---------------------------------------------------------------------------------------------------
"""
    operator(s::String, inds, L; base=2)

Construct an [`Operator`](@ref) from a symbolic local operator string acting on
the full tensor-product basis of a length-`L` system.
"""
function operator(s::String, inds::AbstractVector{<:Integer}, L::Integer; base::Integer=2)
    basis = TensorBasis(L=L, base=base)
    operator(s, inds, basis)
end
#---------------------------------------------------------------------------------------------------
"""
    trans_inv_operator(s::String, inds, basis)

Convenience wrapper around [`spin`](@ref) and [`trans_inv_operator`](@ref) for
translation-invariant symbolic local terms.
"""
function trans_inv_operator(s::String, inds::AbstractVector{<:Integer}, basis::AbstractBasis)
    mat = spin(s, D=basis.B)
    trans_inv_operator(mat, inds, basis)
end
#---------------------------------------------------------------------------------------------------
"""
    trans_inv_operator(s::String, basis)

Construct a translation-invariant many-body operator from a symbolic local
operator string whose support size is inferred from `length(s)`.
"""
function trans_inv_operator(s::String, basis::AbstractBasis)
    mat = spin(s, D=basis.B)
    trans_inv_operator(mat, length(s), basis)
end
#---------------------------------------------------------------------------------------------------
"""
    trans_inv_operator(s::String, L; base=2)

Construct a translation-invariant symbolic operator on the full tensor-product
basis of a length-`L` system.
"""
function trans_inv_operator(s::String, L::Integer; base::Integer=2)
    basis = TensorBasis(L=L, base=base)
    trans_inv_operator(s, basis)
end
