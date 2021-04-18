
struct DoubleBasis{Tb1<:AbstractBasis, Tb2<:AbstractBasis} <: AbstractBasis
    B1::Tb1
    B2::Tb2
end
#-----------------------------------------------------------------------------------------------------
eltype(b::DoubleBasis) = promote_type(eltype(b.B1), eltype(b.B2))
change!(b::DoubleBasis, i::Integer) = change!(b.B2, i)
index(b::DoubleBasis) = index(b.B1)
size(b::DoubleBasis) = size(b.B1, 2), size(b.B2, 2)
function size(b::DoubleBasis, i::Integer)
    if i == 1
        size(b.B1, 2)
    elseif i == 2
        size(b.B2, 2)
    else
        1
    end
end
#-----------------------------------------------------------------------------------------------------
