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
export ent_spec
function ent_spec(v::AbstractVector, Aind::AbstractVector{<:Integer}, b::AbstractBasis)
    m = schmidt(v, Aind, b)
    svdvals(m)
end
#-------------------------------------------------------------------------------------------------------------------------
export entropy
function entropy(s::AbstractVector{<:Real}; α::Real=1, cutoff::Real=1e-20)
    if isone(α)
        shannon_entropy(s, cutoff=cutoff)
    else
        renyi_entropy(s, α)
    end
end

function shannon_entropy(s::AbstractVector{<:Real}; cutoff::Real=1e-20)
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

renyi_entropy(s::AbstractVector{<:Real}, α::Real) = log(sum(s.^α)) / (1-α)
#-------------------------------------------------------------------------------------------------------------------------
export ent_S
function ent_S(v::AbstractVector, Aind::AbstractVector{<:Integer}, b::AbstractBasis; α::Real=1, cutoff::Real=1e-20)
    s = ent_spec(v, Aind, b) .^ 2
    entropy(s, α=α, cutoff=cutoff)
end
