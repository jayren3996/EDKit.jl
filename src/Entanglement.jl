#-------------------------------------------------------------------------------------------------------------------------
# Shape
#-------------------------------------------------------------------------------------------------------------------------
function Base.Vector(v::AbstractVector, b::AbstractBasis)
    B, L = b.B, length(b.dgt)
    ctype = promote_type(eltype(v), eltype(b))
    target = zeros(ctype, B^L, 1)
    schmidt!(target, v, 1:L, b)
    reshape(target, :)
end
#-------------------------------------------------------------------------------------------------------------------------
export schmidt
function schmidt(v::AbstractVector, Aind::AbstractVector{<:Integer}, b::AbstractBasis)
    B, L, LA = b.B, length(b.dgt), length(Aind)
    ctype = promote_type(eltype(v), eltype(b))
    m = zeros(ctype, B^LA, B^(L-LA))
    schmidt!(m, v, Aind, b)
    m
end

#-------------------------------------------------------------------------------------------------------------------------
# Entaglement
#-------------------------------------------------------------------------------------------------------------------------
export entropy
function entropy(s::AbstractVector{<:Real}; cutoff::Real=1e-20)
    ent::Float64 = 0.0
    for si in s
        if si > cutoff
            ent -= si * log(si)
        else
            break
        end
    end
    ent
end
ent_spec(v::AbstractVector, Aind::AbstractVector{<:Integer}, b::AbstractBasis) = schmidt(v, Aind, b) |> svdvals
#-------------------------------------------------------------------------------------------------------------------------
export ent_S
function ent_S(v::AbstractVector, Aind::AbstractVector{<:Integer}, b::AbstractBasis; cutoff::Real=1e-20)
    s = ent_spec(v, Aind, b) .^ 2
    entropy(s, cutoff=cutoff)
end
