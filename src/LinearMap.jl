#-----------------------------------------------------------------------------------------------------
# Double Basis
#-----------------------------------------------------------------------------------------------------
export DoubleBasis
"""
    DoubleBasis{Tb1<:AbstractBasis, Tb2<:AbstractBasis}

Basis for constructing transition matrix from one symmetry sector to another.

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
    DoubleBasis(B1::AbstractBasis, B2::AbstractBasis) = new{typeof(B1), typeof(B2)}(B1.dgt, B1, B2, B2.B)
end

@inline eltype(b::DoubleBasis) = promote_type(eltype(b.B1), eltype(b.B2))
@inline function change!(b::DoubleBasis, i::Integer)
    N = change!(b.B2, i)
    b.dgt .= b.B2.dgt
    N
end
@inline index(b::DoubleBasis) = index(b.B1)
@inline size(b::DoubleBasis) = size(b.B1, 1), size(b.B2, 2)
@inline size(b::DoubleBasis, i::Integer) = isone(i) ? size(b.B1, 1) : isequal(i, 2) ? size(b.B2, 2) : 1



