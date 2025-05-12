#-------------------------------------------------------------------------------------------------------------------------
"""
    SchmidtMatrix{Tm, Ta<:SubArray, Tb<:SubArray, Ti<:Integer}

Struct that is used to construct the Schmidt matrix:

`|Ψ⟩ = ∑ᵢⱼ Mᵢᵢ|ψᵢ⟩ ⊗ |ψⱼ⟩`

Properties:
-----------
- `M`   : Schidt matrix.
- `A`   : View of region A.
- `B`   : View of region B.
- `base`: Base.
"""
struct SchmidtMatrix{Tm <: Number, Ta <: SubArray, Tb <: SubArray, TB1 <: AbstractBasis, TB2 <: AbstractBasis}
    M::Matrix{Tm}
    A::Ta
    B::Tb
    B1::TB1
    B2::TB2
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    schmidtmatrix(T, b::AbstractBasis, Ainds::AbstractVector)

Construction of `SchmidtMatrix`.
"""
function schmidtmatrix(
    T::DataType, b::AbstractBasis, Ainds::AbstractVector{Ta},
    B1=nothing, B2=nothing
) where Ta <: Integer
    L = length(b)
    Binds = Vector{Ta}(undef, L-length(Ainds))
    P = 1
    for i in range(one(Ta), stop=convert(Ta, L))
        if !in(i, Ainds)
            Binds[P] = i
            P += 1
        end
    end
    B1 = isnothing(B1) ? TensorBasis(L=length(Ainds), base=b.B) : B1 
    B2 = isnothing(B2) ? TensorBasis(L=length(Binds), base=b.B) : B2 
    M = zeros(T, size(B1, 1), size(B2, 1))
    SchmidtMatrix(M, view(b.dgt, Ainds), view(b.dgt, Binds), B1, B2)
end
#-------------------------------------------------------------------------------------------------------------------------
function addto!(S::SchmidtMatrix, val::Number)
    S.B1.dgt .= S.A
    S.B2.dgt .= S.B
    _, ia = index(S.B1)
    _, ib = index(S.B2)
    S.M[ia, ib] += val
end

#-------------------------------------------------------------------------------------------------------------------------
# Entanglement Entropy
#-------------------------------------------------------------------------------------------------------------------------
"""
ent_spec(v::AbstractVector, Aind::AbstractVector{<:Integer}, b::AbstractBasis)

Conpute the entanglement spectrum.
"""
ent_spec(v::AbstractVector, Aind::AbstractVector{<:Integer}, b::AbstractBasis) = svdvals(schmidt(v, Aind, b))
#-------------------------------------------------------------------------------------------------------------------------
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
#-------------------------------------------------------------------------------------------------------------------------
"""
Compute Shannon entropy
"""
function shannon_entropy(s::AbstractVector{<:Real}; cutoff::Real=1e-20)
    ent = 0.0
    for si in s
        si > cutoff || break
        ent -= si * log(si)
    end
    ent
end
#-------------------------------------------------------------------------------------------------------------------------
"""
Compute Renyi-0 entropy, which is thenumber of non-zero Schmidt values.
"""
function renyi_zero_entropy(s::AbstractVector{<:Real}; cutoff::Real=1e-20)
    N = 0
    for si in s
        si > cutoff || break
        N += 1
    end
    N
end
#-------------------------------------------------------------------------------------------------------------------------
"""
Compute general Renyi entropy
"""
renyi_entropy(s::AbstractVector{<:Real}, α::Real) = log(sum(s.^α)) / (1-α)
#-------------------------------------------------------------------------------------------------------------------------
export ent_S
"""
ent_S(v::AbstractVector, Aind::AbstractVector{<:Integer}, b::AbstractBasis; α::Real=1, cutoff::Real=1e-20)

Conpute the entanglement entropy of a state.
"""
function ent_S(v::AbstractVector, Aind::AbstractVector{<:Integer}, b::AbstractBasis; α::Real=1, cutoff::Real=1e-20)
    s = ent_spec(v, Aind, b) .^ 2
    entropy(s, α=α, cutoff=cutoff)
end
#-------------------------------------------------------------------------------------------------------------------------
function ent_S(v::AbstractVector, Aind::AbstractVector{<:Integer}, L::Integer; α::Real=1, cutoff::Real=1e-20)
    b = TensorBasis(L, base=round(Int, length(v)^(1/L) ) )
    s = ent_spec(v, Aind, b) .^ 2
    entropy(s, α=α, cutoff=cutoff)
end


#-------------------------------------------------------------------------------------------------------------------------
# Specific Bases
#-------------------------------------------------------------------------------------------------------------------------
"""
schmidt(v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::AbstractOnsiteBasis)

Schmidt decomposition of state `v`, with respect to given lattice bipartition.

Inputs:
-------
- `v`    : State represented by a (abstract) vector. 
- `Ainds`: List of indices in subsystem `A`, the remaining indices are regarded as subsystem `B`.
- `b`    : Basis.

Outputs:
--------
- `S`: Matrix S in the decomposition: |v⟩ = Sᵢⱼ |Aᵢ⟩|Bⱼ⟩.
"""
function schmidt(v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::AbstractOnsiteBasis; B1=nothing, B2=nothing)
    S = schmidtmatrix(eltype(v), b, Ainds, B1, B2)
    for i = 1:length(v)
        change!(b, i)
        addto!(S, v[i])
    end
    S.M
end
#-------------------------------------------------------------------------------------------------------------------------
"""
Schmidt decomposition for `TranslationalBasis`.
"""
function schmidt(v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::TranslationalBasis; B1=nothing, B2=nothing)
    dgt, R, phase = b.dgt, b.R, b.C[2]
    S = schmidtmatrix(promote_type(eltype(v), eltype(b)), b, Ainds, B1, B2)
    for i = 1:length(v)
        change!(b, i)
        val = v[i] / R[i]
        for j in 1:length(dgt)÷b.A
            addto!(S, val)
            circshift!(dgt, b.A)
            val *= phase
        end
    end
    S.M
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    spinflip(v::AbstractVector{<:Integer}, base::Integer)

Flip spins Sz on each site.
"""
function spinflip(v::AbstractVector{<:Integer}, base::Integer)
    vf = Vector{eltype(v)}(undef, length(v))
    base -= 1
    for i = 1:length(vf)
        vf[i] = base - v[i]
    end
    vf
end

function spinflip!(v::AbstractVector{<:Integer}, base::Integer)
    base -= 1
    for i in eachindex(v)
        v[i] = base - v[i]
    end
    v
end
#-------------------------------------------------------------------------------------------------------------------------
"""
Helper function for `schmidt` on parity basis.
"""
function parity_schmidt(parity, v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::AbstractTranslationalParityBasis;B1=nothing, B2=nothing)
    dgt, R, phase = b.dgt, b.R, b.C[2]
    S = schmidtmatrix(promote_type(eltype(v), eltype(b)), b, Ainds, B1, B2)
    for i = 1:length(v)
        change!(b, i)
        val = v[i] / R[i]
        for j in 1:length(dgt)÷b.A
            addto!(S, val)
            circshift!(dgt, b.A)
            val *= phase
        end
        dgt .= parity(dgt)
        val *= b.P
        for j in 1:length(dgt)÷b.A
            addto!(S, val)
            circshift!(dgt, b.A)
            val *= phase
        end
    end
    S.M
end
#-------------------------------------------------------------------------------------------------------------------------
schmidt(v, Ainds, b::TranslationParityBasis;B1=nothing, B2=nothing) = parity_schmidt(reverse, v, Ainds, b; B1, B2)
schmidt(v, Ainds, b::TranslationFlipBasis;B1=nothing, B2=nothing) = parity_schmidt(x -> spinflip(x, b.B), v, Ainds, b; B1, B2)

#-------------------------------------------------------------------------------------------------------------------------
"""
Schmidt decomposition for `FlipBasis`.
"""
function schmidt(v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::FlipBasis; B1=nothing, B2=nothing)
    dgt, R, phase = b.dgt, b.R, b.P
    S = schmidtmatrix(promote_type(eltype(v), eltype(b)), b, Ainds, B1, B2)
    for i = 1:length(v)
        change!(b, i)
        val = v[i] / R[i]
        addto!(S, val)
        dgt .= spinflip(dgt, b.B)
        addto!(S, phase * val)
    end
    S.M
end
#-------------------------------------------------------------------------------------------------------------------------
"""
Schmidt decomposition for `ParityBasis`.
"""
function schmidt(v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::ParityBasis; B1=nothing, B2=nothing)
    dgt, R, phase = b.dgt, b.R, b.P
    S = schmidtmatrix(promote_type(eltype(v), eltype(b)), b, Ainds, B1, B2)
    for i = 1:length(v)
        change!(b, i)
        val = v[i] / R[i]
        addto!(S, val)
        reverse!(dgt)
        addto!(S, phase * val)
    end
    S.M
end
#-------------------------------------------------------------------------------------------------------------------------
"""
Schmidt decomposition for `ParityBasis`.
"""
function schmidt(v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::ParityFlipBasis; B1=nothing, B2=nothing)
    dgt, R, p1, p2 = b.dgt, b.R, b.P, b.Z
    S = schmidtmatrix(promote_type(eltype(v), eltype(b)), b, Ainds, B1, B2)
    for i = 1:length(v)
        # (P,Z) = (0,0)
        change!(b, i)
        val = v[i] / R[i]
        addto!(S, val)
        # (P,Z) = (1,0)
        reverse!(dgt)
        val *= p1
        addto!(S, val) 
        # (P,Z) = (1,1)
        dgt .= spinflip(dgt, b.B)
        val *= p2
        addto!(S, val) 
        # (P,Z) = (0,1)
        reverse!(dgt)
        val *= p1
        addto!(S, val) 
    end
    S.M
end