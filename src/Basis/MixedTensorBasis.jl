export MixedTensorBasis
"""
    MixedTensorBasis <: AbstractOnsiteBasis

Full tensor-product basis with **per-site local dimensions**. Site `i` carries a
local Hilbert space of dimension `dims[i]`, so the global index is the
mixed-radix (big-endian) encoding

    I(i₁, …, i_L) = i₁·∏_{m>1} dims[m] + ⋯ + i_L + 1.

This generalizes [`TensorBasis`](@ref) (which fixes a single `base` for every
site) to e.g. alternating spin-½/spin-1 chains, spin-boson models, or bosons
with different per-site cutoffs. The scalar `TensorBasis` and its base-2 fast
kernel are untouched; mixed-radix `index`/`change!` are selected by passing the
per-site dimension vector as `base`.

Operators are built exactly as for any onsite basis: `operator(mat, sites, B)`
where `mat` acts on the local space of dimension `∏ dims[sites]`.
"""
struct MixedTensorBasis <: AbstractOnsiteBasis
    dgt::Vector{Int64}
    B::Vector{Int64}
end

"""
    MixedTensorBasis(; dims)

Construct a mixed-dimension tensor-product basis from the per-site local
dimensions `dims` (e.g. `dims=[2, 3, 2, 3]`). Each dimension must be `≥ 1`, and
the total Hilbert-space dimension `∏ dims` must fit a 64-bit index.
"""
function MixedTensorBasis(; dims::AbstractVector{<:Integer})
    isempty(dims) && error("`dims` must be non-empty.")
    all(≥(1), dims) || error("Every local dimension must be ≥ 1 (got $dims).")
    total = prod(big.(dims))
    total < typemax(Int64) || error("Total Hilbert dimension $total exceeds the Int64 index capacity.")
    B = Vector{Int64}(dims)
    MixedTensorBasis(zeros(Int64, length(dims)), B)
end

content(::MixedTensorBasis, i::Integer) = i
norm(::MixedTensorBasis, ::Integer) = 1
eltype(::MixedTensorBasis) = Int64
int_type(::MixedTensorBasis) = Int64
size(b::MixedTensorBasis, i::Integer) = isone(i) || isequal(i, 2) ? prod(b.B) : 1
size(b::MixedTensorBasis) = (l = prod(b.B); (l, l))
index(b::MixedTensorBasis) = 1, index(b.dgt, base=b.B)
index(b::MixedTensorBasis, dgt::AbstractVector) = 1, index(dgt, base=b.B)
copy(b::MixedTensorBasis) = MixedTensorBasis(deepcopy(b.dgt), copy(b.B))
