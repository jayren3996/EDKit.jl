#-----------------------------------------------------------------------------------------------------
# Double Basis
#-----------------------------------------------------------------------------------------------------
"""
    DoubleBasis{Tb1<:AbstractBasis, Tb2<:AbstractBasis}

Basis for constructing transition matrix from one symmetry sector to another.

Properties:
-----------
B1 : Basis of the target symmetry sector.
B2 : Basis of the starting symmetry sector.
"""
struct DoubleBasis{Tb1<:AbstractBasis, Tb2<:AbstractBasis} <: AbstractBasis
    B1::Tb1
    B2::Tb2
end
#-----------------------------------------------------------------------------------------------------
eltype(b::DoubleBasis) = promote_type(eltype(b.B1), eltype(b.B2))
change!(b::DoubleBasis, i::Integer) = change!(b.B2, i)
index(b::DoubleBasis) = index(b.B1)
size(b::DoubleBasis) = size(b.B1, 2), size(b.B2, 2)
size(b::DoubleBasis, i::Integer) = (i == 1) ? size(b.B1, 2) : (i == 2) : size(b.B2, 2) : 1
