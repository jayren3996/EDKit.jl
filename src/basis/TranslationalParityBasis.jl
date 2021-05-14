#-------------------------------------------------------------------------------------------------------------------------
# Abstract Translational + Parity Basis
#-------------------------------------------------------------------------------------------------------------------------
abstract type AbstractTranslationalParityBasis <: AbstractBasis end

#-------------------------------------------------------------------------------------------------------------------------
# Translational + Parity Basis
#-------------------------------------------------------------------------------------------------------------------------
"""
    TranslationParityBasis

Basis with translational and reflection symmetries.

Properties:
-----------
- dgt::Vector{BITTYPE}    : Digits
- I::Vector{Int}      : Representing states
- R::Vector{Float64}  : Normalization
- K::Int              : {±1}, momentum phase
- P::Int              : {±1}, parity
- B::Int              : Base
"""
struct TranslationParityBasis <: AbstractTranslationalParityBasis
    dgt::Vector{BITTYPE}    # Digits
    I::Vector{Int}          # Representing states
    R::Vector{Float64}      # Normalization
    C::Vector{Int}          # Momentum phase exp(1im * 0/π) = {±1}
    P::Int                  # {±1}, parity
    B::Int                  # Base
end

"""
    TranslationFlipBasis

Basis with translational and spin-flip symmetries.

Properties:
-----------
- dgt::Vector{BITTYPE}    : Digits
- I::Vector{Int}      : Representing states
- R::Vector{Float64}  : Normalization
- C::Vector{T}        : Momentum phase
- P::Int              : {±1}, parity
- B::Int              : Base
"""
struct TranslationFlipBasis{T} <: AbstractTranslationalParityBasis
    dgt::Vector{BITTYPE}    # Digits
    I::Vector{Int}          # Representing states
    R::Vector{Float64}      # Normalization
    C::Vector{T}            # Momentum phase exp(-1im * 2π/L * K)
    P::Int                  # {±1}, parity
    B::Int                  # Base
end

eltype(::TranslationParityBasis) = Float64
eltype(::TranslationFlipBasis{T}) where T = T

#-------------------------------------------------------------------------------------------------------------------------
# Indexing
#-------------------------------------------------------------------------------------------------------------------------
function parity_index(parity!, b::AbstractTranslationalParityBasis)
    reflect, i, t = parity_index(parity!, b.dgt, b.B)
    ind = binary_search(b.I, i)
    N, I = if iszero(ind)
        0.0, 1
    else
        n = reflect ? b.P * b.C[t+1] : b.C[t+1]
        n * b.R[ind], ind
    end
    N, I
end

index(b::TranslationParityBasis) = parity_index(reverse!, b)
index(b::TranslationFlipBasis) = parity_index(x -> spinflip!(x, b.B), b)

#-------------------------------------------------------------------------------------------------------------------------
# Judge
#-------------------------------------------------------------------------------------------------------------------------
struct TranslationParityJudge
    F                   # Projective selection
    parity!             # Parity operation acting on digits
    K::Int              # Momentum
    P::Int              # Parity eigenvalue {±1}
    L::Int              # Length
    B::Int              # Base
    C::Vector{Float64}  # Normalization coefficients
end

function (judge::TranslationParityJudge)(dgt::AbstractVector{<:Integer}, i::Integer)
    Q, N = if judge.F(dgt)
        isrep, reflect, r, m = parity_check(judge.parity!, dgt, i, judge.B)
        if isrep && iszero(mod(judge.K * r, judge.L))
            if reflect
                trans = mod(judge.K * m, judge.L)
                if ( isequal(2 * trans, judge.L) && isequal(judge.P, -1) ) || ( iszero(trans) && isone(judge.P) )
                    (true, judge.C[judge.L+r])
                else
                    (false, 0.0)
                end
            else
                (true, judge.C[r])
            end
        else
            (false, 0.0)
        end
    else
        (false, 0.0)
    end
    Q, N
end

#-------------------------------------------------------------------------------------------------------------------------
# Construction
#-------------------------------------------------------------------------------------------------------------------------
export translationparitybasis, translationflipbasis
function translationparitybasis(
    f, K::Integer, P::Integer, L::Integer; 
    base::Integer=2, alloc::Integer=1000, threaded=false
)
    @assert iszero(K) || isequal(2K, L) "Momentum incompatible with parity."
    @assert isone(P) || isequal(P, -1) "Invalid parity"
    judge = begin
        N2 = [2L / sqrt(i) for i = 1:L]
        N1 = N2 ./ sqrt(2)
        TranslationParityJudge(f, reverse!, K, P, L, base, vcat(N1, N2))
    end

    I, R = if threaded
        selectindexnorm_threaded(judge, L, base=base, alloc=alloc)
    else
        selectindexnorm(judge, L, 1:base^L, base=base, alloc=alloc)
    end

    C = if iszero(K)
        fill(1.0, L)
    else
        [iseven(i) ? 1.0 : -1.0 for i=0:L-1]
    end

    TranslationParityBasis(zeros(Int, L), I, R, C, P, base)
end

translationparitybasis(
    K::Integer, P::Integer, L::Integer; 
    base::Integer=2, alloc::Integer=1000, threaded=false
) = translationparitybasis(x -> true, K, P, L, base=base, alloc=alloc, threaded=threaded)

function translationflipbasis(
    f, K::Integer, P::Integer, L::Integer; 
    base::Integer=2, alloc::Integer=1000, threaded=false
)
    @assert isone(P) || isequal(P, -1) "Invalid parity"
    judge = begin
        N2 = [2L / sqrt(i) for i = 1:L]
        N1 = N2 ./ sqrt(2)
        TranslationParityJudge(f, x -> spinflip!(x, base), K, P, L, base, vcat(N1, N2))
    end

    I, R = if threaded
        selectindexnorm_threaded(judge, L, base=base, alloc=alloc)
    else
        selectindexnorm(judge, L, 1:base^L, base=base, alloc=alloc)
    end

    C = if iszero(K)
        fill(1.0, L)
    elseif isequal(2K, L)
        [iseven(i) ? 1.0 : -1.0 for i=0:L-1]
    else
        [exp(-1im * 2π/L * K * i) for i=0:L-1]
    end

    TranslationFlipBasis(zeros(BITTYPE, L), I, R, C, P, base)
end

translationflipbasis(
    K::Integer, P::Integer, L::Integer; 
    base::Integer=2, alloc::Integer=1000, threaded=false
) = translationflipbasis(x -> true, K, P, L, base=base, alloc=alloc, threaded=threaded)

#-------------------------------------------------------------------------------------------------------------------------
# Schmidt Form
#-------------------------------------------------------------------------------------------------------------------------
function parity_schmidt!(
    parity!, target::AbstractMatrix, v::AbstractVector, 
    Ainds::AbstractVector{<:Integer}, b::AbstractTranslationalParityBasis
)
    dgt, R, phase = b.dgt, b.R, b.C[2]
    Binds = complement(Ainds, length(dgt))
    for i = 1:length(v)
        val = v[i] / R[i]
        change!(dgt, b.I[i], base=b.B)
        cycle_write!(target, dgt, val, phase, Ainds, Binds, b.B)
        cyclebits!(dgt)
        parity!(dgt)
        cycle_write!(target, dgt, b.P * val, phase, Ainds, Binds, b.B)
    end
    target
end

function schmidt!(target::AbstractMatrix, v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::TranslationParityBasis)
    parity_schmidt!(reverse!, target, v, Ainds, b)
end

function schmidt!(target::AbstractMatrix, v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::TranslationFlipBasis)
    parity_schmidt!(x -> spinflip!(x, b.B), target, v, Ainds, b)
end
