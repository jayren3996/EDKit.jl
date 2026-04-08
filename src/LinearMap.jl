#-----------------------------------------------------------------------------------------------------
# Helpers
#-----------------------------------------------------------------------------------------------------
"""
    index_nocheck(b)

Variant of [`index`](@ref) that never throws on absent states.

For `ProjectedBasis` this explicitly routes through `index(...; check=false)`,
which is important in map-building code where off-sector states should
contribute zero rather than abort the calculation.
"""
index_nocheck(b::AbstractBasis) = index(b)
index_nocheck(b::ProjectedBasis) = index(b, check=false)
index_nocheck(b::AbstractBasis, dgt::AbstractVector) = index(b, dgt)
index_nocheck(b::ProjectedBasis, dgt::AbstractVector) = index(b, dgt; check=false)
"""
    orbit_order(B::AbstractBasis)

Return the orbit size associated with basis `B` for embedding purposes.

Onsite bases have orbit size `1`; reduced bases use their `order(B)` value.
This helper lets embedding code treat both cases uniformly.
"""
orbit_order(B::AbstractBasis) = B isa AbstractOnsiteBasis ? 1 : order(B)

"""
    check_doublebasis_compatible(B1, B2)

Validate that two bases can be combined into a [`DoubleBasis`](@ref).

Compatibility means:
- same physical system size,
- same local on-site base.

An `ArgumentError` is thrown if either condition fails.
"""
function check_doublebasis_compatible(B1::AbstractBasis, B2::AbstractBasis)
    length(B1) == length(B2) || throw(ArgumentError("DoubleBasis requires bases with the same length."))
    B1.B == B2.B || throw(ArgumentError("DoubleBasis requires bases with the same local base."))
    nothing
end

"""
    basis_embedding(B::AbstractBasis)

Construct the explicit embedding matrix from coordinates in basis `B` to the
full tensor-product basis.

Returns:
- A dense matrix of size `dim(full) × dim(B)`.

This helper is the bridge between abstract symmetry-reduced coordinates and
ordinary full-space amplitudes. It is mainly used to build [`symmetrizer`](@ref)
matrices.
"""
function basis_embedding(B::AbstractBasis)
    full = TensorBasis(L = length(B), base = B.B)
    embed = zeros(ComplexF64, size(full, 1), size(B, 1))
    ord = orbit_order(B)
    dgt = similar(B.dgt)
    full_dgt = similar(full.dgt)
    for j in 1:size(full, 1)
        change!(full, j, full_dgt)
        dgt .= full_dgt
        coeff, pos = index_nocheck(B, dgt)
        iszero(coeff) || (embed[j, pos] = coeff / ord)
    end
    embed
end

#-----------------------------------------------------------------------------------------------------
# DoubleBasis
#-----------------------------------------------------------------------------------------------------
export DoubleBasis
"""
    DoubleBasis{Tb1<:AbstractBasis, Tb2<:AbstractBasis}

Basis for constructing transition matrix from one symmetry sector to another.
Note that DoubleBasis can be used as the projector, meaning that it will ignore the symmetry violation.
Therefore, extra care is needed when working with this basis.

Interpretation:
- `B1` is the target basis,
- `B2` is the source basis,
- `dgt` follows the current source-side representative as iteration proceeds.
"""
struct DoubleBasis{Tb1<:AbstractBasis, Tb2<:AbstractBasis} <: AbstractBasis
    dgt::Vector{Int64}
    B1::Tb1
    B2::Tb2
    B::Int64
end

"""
    DoubleBasis(B1, B2)

Construct a basis-to-basis map descriptor from source basis `B2` to target
basis `B1`.

Arguments:
- `B1`: target basis.
- `B2`: source basis.

Returns:
- A [`DoubleBasis`](@ref) object describing maps of size `size(B1, 1) × size(B2, 2)`.
"""
function DoubleBasis(B1::AbstractBasis, B2::AbstractBasis)
    check_doublebasis_compatible(B1, B2)
    DoubleBasis(B1.dgt, B1, B2, B2.B)
end
eltype(b::DoubleBasis) = promote_type(eltype(b.B1), eltype(b.B2))
function change!(b::DoubleBasis, i::Integer)
    N = change!(b.B2, i)
    b.dgt .= b.B2.dgt
    N
end
function change!(b::DoubleBasis, i::Integer, dgt::AbstractVector)
    change!(dgt, content(b.B2, i), base=b.B2.B)
    norm(b.B2, i)
end
index(b::DoubleBasis; check::Bool=false) = check ? index(b.B1) : index_nocheck(b.B1)
index(b::DoubleBasis, dgt::AbstractVector) = index(b.B1, dgt)
size(b::DoubleBasis) = size(b.B1, 1), size(b.B2, 2)
size(b::DoubleBasis, i::Integer) = isone(i) ? size(b.B1, 1) : isequal(i, 2) ? size(b.B2, 2) : 1
copy(b::DoubleBasis) = DoubleBasis(deepcopy(b.B1), deepcopy(b.B2))

"""
    (B::DoubleBasis)(v)

Apply the basis-overlap map encoded by `B` directly to a coordinate vector `v`
expressed in `B.B2`.

Returns:
- A vector of coordinates in `B.B1`.

The action is computed by summing overlaps through the full tensor-product
embedding, not by assuming the two bases share the same representative set.
"""
function (B::DoubleBasis)(v::AbstractVector)
    length(v) == size(B, 2) || throw(DimensionMismatch("Expected a vector of length $(size(B, 2)); got $(length(v))."))

    T = promote_type(ComplexF64, eltype(v), eltype(B.B1), eltype(B.B2))
    out = zeros(T, size(B, 1))
    full = TensorBasis(L = length(B), base = B.B)
    ord1 = orbit_order(B.B1)
    ord2 = orbit_order(B.B2)
    full_dgt = similar(full.dgt)
    dgt1 = similar(B.B1.dgt)
    dgt2 = similar(B.B2.dgt)

    for j in 1:size(full, 1)
        change!(full, j, full_dgt)
        dgt2 .= full_dgt
        coeff2, pos2 = index_nocheck(B.B2, dgt2)
        iszero(coeff2) && continue

        dgt1 .= full_dgt
        coeff1, pos1 = index_nocheck(B.B1, dgt1)
        iszero(coeff1) && continue

        out[pos1] += conj(coeff1 / ord1) * (coeff2 / ord2) * v[pos2]
    end

    out
end

export symmetrizer
"""
    symmetrizer(B::DoubleBasis)

Return the linear map from coordinates in `B.B2` to coordinates in `B.B1`.

This is the special basis-only map associated with a [`DoubleBasis`](@ref):

- if `B.B2` is a less symmetric basis, it acts as a symmetrization / projection
  into the target sector `B.B1`;
- in the general case it returns the overlap map between the two basis
  embeddings in the full tensor-product basis.

The returned matrix has size `size(B, 1) × size(B, 2)`.
"""
function symmetrizer(B::DoubleBasis)
    adjoint(basis_embedding(B.B1)) * basis_embedding(B.B2)
end
