export ent_spec, ent_S, entropy
"""
SchmidtMatrix{Tm <: Number, Ta <: SubArray, Tb <: SubArray, Ti <: Integer}

Struct that is used to construct the Schmidt matrix:

`|Ψ⟩ = ∑ᵢⱼ Mᵢᵢ|ψᵢ⟩ ⊗ |ψⱼ⟩`

Properties:
-----------
- `M`   : Schidt matrix.
- `A`   : View of region A.
- `B`   : View of region B.
- `base`: Base.
"""
struct SchmidtMatrix{Tm <: Number, Ta <: SubArray, Tb <: SubArray, Ti <: Integer}
    M::Matrix{Tm}
    A::Ta
    B::Tb
    base::Ti
end

Matrix(S::SchmidtMatrix) = S.M

"""
    SchmidtMatrix(T::DataType, b::AbstractBasis, Ainds::AbstractVector{Ta}) where Ta <: Integer

Construction of `SchmidtMatrix`.
"""
function SchmidtMatrix(T::DataType, b::AbstractBasis, Ainds::AbstractVector{Ta}) where Ta <: Integer
    L = length(b)
    B = base(b)
    Binds = Vector{Ta}(undef, L-length(Ainds))
    P = 1
    for i in range(one(Ta), stop=convert(Ta, L))
        if !in(i, Ainds)
            Binds[P] = i
            P += 1
        end
    end
    M = zeros(T, B^length(Ainds), B^length(Binds))
    SchmidtMatrix(M, view(digits(b), Ainds), view(digits(b), Binds), B)
end

function addto!(S::SchmidtMatrix, val::Number)
    ia = index(S.A, base=S.base)
    ib = index(S.B, base=S.base)
    S.M[ia, ib] += val
end

"""
    Vector(v::AbstractVector, b::AbstractBasis)

Convert a state in the symmetry sector to the original state in tensor-product Hilbert space.
"""
Vector(v::AbstractVector, b::AbstractBasis) = reshape(schmidt(v, 1:L, b), :)

"""
    ent_spec(v::AbstractVector, Aind::AbstractVector{<:Integer}, b::AbstractBasis)

Conpute the entanglement spectrum.
"""
ent_spec(v::AbstractVector, Aind::AbstractVector{<:Integer}, b::AbstractBasis) = svdvals(schmidt(v, Aind, b))

"""
    entropy(s::AbstractVector{<:Real}; α::Real=1, cutoff::Real=1e-20)

Compute the entropy of Schmidt values.  

Inputs:
-------
- `s`     : Schmidt values.
- `α`     : Renyi index.
- `cutoff`: Cutoff of the Schmidt values. 

Outputs:
--------
- `S`: Entanglement entropy.
"""
function entropy(s::AbstractVector{<:Real}; α::Real=1, cutoff::Real=1e-20)
    if isone(α)
        shannon_entropy(s, cutoff=cutoff)
    elseif iszero(α)
        renyi_zero_entropy(s, cutoff=cutoff)
    else
        renyi_entropy(s, α)
    end
end

"""
Compute Shannon entropy
"""
function shannon_entropy(s::AbstractVector{<:Real}; cutoff::Real=1e-20)
    ent = 0.0
    for si in s
        if si > cutoff
            ent -= si * log(si)
        else
            break
        end
    end
    ent
end

"""
Compute Renyi-0 entropy, which is thenumber of non-zero Schmidt values.
"""
function renyi_zero_entropy(s::AbstractVector{<:Real}; cutoff::Real=1e-20)
    N = 0
    for si in s
        if si > cutoff
            N += 1
        else
            break
        end
    end
    N
end

"""
Compute general Renyi entropy
"""
renyi_entropy(s::AbstractVector{<:Real}, α::Real) = log(sum(s.^α)) / (1-α)

"""
    ent_S(v::AbstractVector, Aind::AbstractVector{<:Integer}, b::AbstractBasis; α::Real=1, cutoff::Real=1e-20)

Conpute the entanglement entropy of a state.
"""
function ent_S(v::AbstractVector, Aind::AbstractVector{<:Integer}, b::AbstractBasis; α::Real=1, cutoff::Real=1e-20)
    s = ent_spec(v, Aind, b) .^ 2
    entropy(s, α=α, cutoff=cutoff)
end

function ent_S(v::AbstractVector, Aind::AbstractVector{<:Integer}, L::Integer; α::Real=1, cutoff::Real=1e-20)
    b = TensorBasis(L, base=round(Int, length(v)^(1/L) ) )
    s = ent_spec(v, Aind, b) .^ 2
    entropy(s, α=α, cutoff=cutoff)
end
