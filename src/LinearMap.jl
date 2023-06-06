#-----------------------------------------------------------------------------------------------------
# Helpers
#-----------------------------------------------------------------------------------------------------
index_nocheck(b::AbstractBasis) = index(b)
index_nocheck(b::ProjectedBasis) = index(b, check=false)
cyclelength(b::TranslationalBasis) = length(b)
cyclelength(b::AbstractTranslationalParityBasis) = 2 * length(b)

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

DoubleBasis(B1::AbstractBasis, B2::AbstractBasis) = DoubleBasis(B1.dgt, B1, B2, B2.B)
eltype(b::DoubleBasis) = promote_type(eltype(b.B1), eltype(b.B2))
function change!(b::DoubleBasis, i::Integer)
    N = change!(b.B2, i)
    b.dgt .= b.B2.dgt
    N
end
index(b::DoubleBasis; check::Bool=false) = check ? index(b.B1) : index_nocheck(b.B1)
size(b::DoubleBasis) = size(b.B1, 1), size(b.B2, 2)
size(b::DoubleBasis, i::Integer) = isone(i) ? size(b.B1, 1) : isequal(i, 2) ? size(b.B2, 2) : 1
copy(b::DoubleBasis) = DoubleBasis(copy(b.B1), copy(b.B2))


#---------------------------------------------------------------------------------------------------
# Transform a vector from basis `B1` to basis `B2`.
#---------------------------------------------------------------------------------------------------
function (B::DoubleBasis{<:AbstractOnsiteBasis, <:AbstractOnsiteBasis})(v::AbstractVecOrMat)
    dtype = promote_type(eltype(B), eltype(v))
    out = isa(v, AbstractVector) ? zeros(dtype, size(B, 1)) : zeros(dtype, size(B, 1), size(v, 2))
    for i in axes(v, 1)
        change!(B, i)
        N, j = index(B)
        isa(v, AbstractVector) ? out[j] = N * v[i] : out[j, :] = N * v[i, :]
    end
    out
end

function (B::DoubleBasis{<:AbstractOnsiteBasis, <:AbstractPermuteBasis})(v::AbstractVecOrMat)
    dtype = promote_type(eltype(B), eltype(v))
    out = isa(v, AbstractVector) ? zeros(dtype, size(B, 1)) : zeros(dtype, size(B, 1), size(v, 2))
    B1, B2 = B.B1, B.B2
    L = cyclelength(B2)
    for i in 1:size(B,1)
        change!(B1, i)
        B2.dgt .= B1.dgt 
        N, j = index(B2)
        iszero(N) && continue
        isa(v, AbstractVector) ? out[i] = conj(N) * v[j] / L : out[i, :] = conj(N) * v[j, :] / L
    end
    out
end

function (B::DoubleBasis{<:AbstractPermuteBasis, <:AbstractOnsiteBasis})(v::AbstractVecOrMat)
    dtype = promote_type(eltype(B), eltype(v))
    out = isa(v, AbstractVector) ? zeros(dtype, size(B, 1)) : zeros(dtype, size(B, 1), size(v, 2))
    L = cyclelength(B.B1)
    for i in axes(v, 1)
        change!(B, i)
        N, j = index(B)
        iszero(N) && continue
        isa(v, AbstractVector) ? out[j] = L * v[i] / N : out[j, :] = L * v[i, :] / N
    end
    out
end

function (B::DoubleBasis{<:TranslationalBasis, <:AbstractTranslationalParityBasis})(v::AbstractVecOrMat)
    dtype = promote_type(eltype(B), eltype(v))
    out = isa(v, AbstractVector) ? zeros(dtype, size(B, 1)) : zeros(dtype, size(B, 1), size(v, 2))
    B1, B2 = B.B1, B.B2
    L = 2
    for i in 1:size(B,1)
        R1 = change!(B1, i)
        B2.dgt .= B1.dgt 
        N, j = index(B2)
        iszero(N) && continue
        isa(v, AbstractVector) ? out[i] = N * v[j] / R1 / L : out[i, :] = N * v[j, :] / R1 / L
    end
    out
end