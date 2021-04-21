#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Translationnal + Parity Bases
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
    TranslationParityBasis

Basis with translational and reflection symmetries.

Properties:
-----------
- dgt::Vector{Int}    : Digits
- I::Vector{Int}      : Representing states
- R::Vector{Float64}  : Normalization
- K::Int              : {±1}, momentum 0/π
- P::Int              : {±1}, parity
- B::Int              : Base
"""
struct TranslationParityBasis <: AbstractBasis
    dgt::Vector{Int}        # Digits
    I::Vector{Int}          # Representing states
    R::Vector{Float64}      # Normalization
    K::Int                  # {±1}, momentum 0/π
    P::Int                  # {±1}, parity
    B::Int                  # Base
end
eltype(::TranslationParityBasis) = Float64
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function index(b::TranslationParityBasis)::Tuple{Float64, Int}
    reflect, i, t = parity_index(b.dgt, b.B)
    ind = binary_search(b.I, i)
    if ind == 0
        0.0, 1
    else
        N = reflect ? (b.P * b.K ^ t) : (b.K ^ t)
        N * b.R[ind], ind
    end
end

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
struct TranslationParityInfo
    F
    K::Int
    P::Int
    B::Int
    S::Vector{Float64}
end
translationparityinfo(F, K::Int, P::Int, B::Integer, S::Vector{Float64}) = TranslationParityInfo(F, K, P, Int(B), S)

function (info::TranslationParityInfo)(dgt::Vector{Int}, i::Integer)::Tuple{Bool, Float64}
    if info.F(dgt)
        isrep, reflect, r, m = parity_check(dgt, i, info.B)
        if isrep && (info.K ^ r == 1)
            Q, N = parity_norm(reflect, info.K, info.P, r, m, info.S)
            Q ? (Q, N) : (Q, 0.0)
        else
            (false, 0.0)
        end
    else
        (false, 0.0)
    end
end

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Construction
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
export translationparitybasis
function translationparitybasis(
    f, k::Integer, p::Integer, L::Integer; 
    base::Integer=2, alloc::Integer=1000, threaded=false
)
    dgt = zeros(Int, L)
    K = (k == 0) ? 1 : (2k == L) ? -1 : error("Momentum incompatible with parity.")
    P = (p == 1) ? 1 : (p == -1) ? -1 : error("Invalid parity $p.")
    S = append!([L/sqrt(i/2) for i = 1:L], sqrt(2))
    info = translationparityinfo(f, K, P, base, S)
    I, R = if threaded
        selectindexnorm_threaded(info, L, base=base, alloc=alloc)
    else
        selectindexnorm(info, L, 1:base^L, base=base, alloc=alloc)
    end
    TranslationParityBasis(dgt, I, R, K, P, base)
end

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Helper Functions
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
    parity_check(reverse!, dgt::AbstractVector{<:Integer}, I0::Integer, base::Integer)

Check whether a state is a valid representing state, and whether it is reflect-translation-invariant.

Inputs:
-------
- dgt     : Digits vector.
- I0      : Index of the index.
- base    : Base.

Outputs:
--------
Tuple{Bool, Bool, Int, Int} : 
    1. Whether `dgt` is a representing vector.
    2. Whether `dgt` is reflect-translation-invariant.
    3. Minimum R that Tᴿ⋅|dgt⟩ = |dgt⟩.
    4. Minimum M that Tᴹ⋅P⋅|dgt⟩ = |dgt⟩, 0 if it is not reflect-translation-invariant.
"""
function parity_check(dgt::AbstractVector{<:Integer}, I0::Integer, base::Integer)::Tuple{Bool, Bool, Int, Int}
    isrep, r = translation_check(dgt, I0, base)
    if isrep
        reverse!(dgt)
        for i = 0:r-1
            In = index(dgt, base=base)
            if In < I0
                change!(dgt, I0, base=base)
                return (false, false, r, 0)
            elseif In == I0
                return (true, true, r, i)
            end
            cyclebits!(dgt)
        end
        reverse!(dgt)
        (true, false, r, 0)
    else
        (false, false, r, 0) 
    end
end
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function parity_index(dgt::AbstractVector{<:Integer}, base::Integer)
    Im1, T1 = translation_index(dgt, base)
    reverse!(dgt)
    Im2, T2 = translation_index(dgt, base)
    reverse!(dgt)
    Im2 < Im1 ? (true, Im2, T2) : (false, Im1, T1)
end
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function parity_norm(reflect::Bool, k::Int, p::Int, r::Int, m::Int, S::Vector{Float64})
    if reflect
        p * k ^ m == 1 ? (true, S[end]*S[r]) : (false, 0.0)
    else
        (true, S[r])
    end
end
