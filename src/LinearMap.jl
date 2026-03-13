#-----------------------------------------------------------------------------------------------------
# Helpers
#-----------------------------------------------------------------------------------------------------
index_nocheck(b::AbstractBasis) = index(b)
index_nocheck(b::ProjectedBasis) = index(b, check=false)
orbit_order(B::AbstractBasis) = B isa AbstractOnsiteBasis ? 1 : order(B)

function check_doublebasis_compatible(B1::AbstractBasis, B2::AbstractBasis)
    length(B1) == length(B2) || throw(ArgumentError("DoubleBasis requires bases with the same length."))
    B1.B == B2.B || throw(ArgumentError("DoubleBasis requires bases with the same local base."))
    nothing
end

function basis_embedding(B::AbstractBasis)
    full = TensorBasis(L = length(B), base = B.B)
    embed = zeros(ComplexF64, size(full, 1), size(B, 1))
    ord = orbit_order(B)
    Bc = deepcopy(B)
    for j in 1:size(full, 1)
        change!(full, j)
        Bc.dgt .= full.dgt
        coeff, pos = index_nocheck(Bc)
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

Properties:
-----------
- `dgt`: Digits.
- `B1` : Basis of the target symmetry sector.
- `B2` : Basis of the starting symmetry sector.
- `B`  : Base.
"""
struct DoubleBasis{Tb1<:AbstractBasis, Tb2<:AbstractBasis} <: AbstractBasis
    dgt::Vector{Int64}
    B1::Tb1
    B2::Tb2
    B::Int64
end

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
index(b::DoubleBasis; check::Bool=false) = check ? index(b.B1) : index_nocheck(b.B1)
size(b::DoubleBasis) = size(b.B1, 1), size(b.B2, 2)
size(b::DoubleBasis, i::Integer) = isone(i) ? size(b.B1, 1) : isequal(i, 2) ? size(b.B2, 2) : 1
copy(b::DoubleBasis) = DoubleBasis(deepcopy(b.B1), deepcopy(b.B2))

function (B::DoubleBasis)(v::AbstractVector)
    length(v) == size(B, 2) || throw(DimensionMismatch("Expected a vector of length $(size(B, 2)); got $(length(v))."))

    T = promote_type(ComplexF64, eltype(v), eltype(B.B1), eltype(B.B2))
    out = zeros(T, size(B, 1))
    full = TensorBasis(L = length(B), base = B.B)
    B1c = deepcopy(B.B1)
    B2c = deepcopy(B.B2)
    ord1 = orbit_order(B.B1)
    ord2 = orbit_order(B.B2)

    for j in 1:size(full, 1)
        change!(full, j)
        B2c.dgt .= full.dgt
        coeff2, pos2 = index_nocheck(B2c)
        iszero(coeff2) && continue

        B1c.dgt .= full.dgt
        coeff1, pos1 = index_nocheck(B1c)
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
