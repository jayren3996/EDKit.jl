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
struct DoubleBasis{Tb<:AbstractBasis} <: AbstractBasis
    dgt::Vector{Int}
    B1::Tb
    B2::Tb
    B::Int
end
#-----------------------------------------------------------------------------------------------------
eltype(b::DoubleBasis) = promote_type(eltype(b.B1), eltype(b.B2))
function change!(b::DoubleBasis, i::Integer)
    N = change!(b.B2, i)
    b.dgt .= b.B2.dgt
    N
end
index(b::DoubleBasis) = index(b.B1)
size(b::DoubleBasis) = size(b.B1, 2), size(b.B2, 2)
size(b::DoubleBasis, i::Integer) = isone(i) ? size(b.B1, 2) : isequal(i, 2) ? size(b.B2, 2) : 1
#-----------------------------------------------------------------------------------------------------
export doublebasis
function doublebasis(B1::Tb, B2::Tb) where Tb <: AbstractBasis
    DoubleBasis(B1.dgt, B1, B2, B2.B)
end
