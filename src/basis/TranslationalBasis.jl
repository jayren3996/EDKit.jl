#-------------------------------------------------------------------------------------------------------------------------
# TranslationalBasis
#-------------------------------------------------------------------------------------------------------------------------
"""
    TranslationalBasis

Basis for subspace that is spanned by momentum states, which can also incorporate projected restriction.

Properties:
-----------
- dgt : Vector{BITTYPE}
- I   : Vector{Int}, list of indicies.
- R   : Vector{Float64}, list of normalization.
- C   : Vector{Float64/ComplexF64}, unit phase factor.
- B   : Int, levels on each site.
"""
struct TranslationalBasis{T <: Number} <: AbstractBasis
    dgt::Vector{BITTYPE}
    I::Vector{Int}
    R::Vector{Float64}
    C::Vector{T}
    B::UInt8
    TranslationalBasis(dgt::Vector{BITTYPE}, I, R, C::Vector{T}, B::Integer) where T = new{T}(dgt, I, R, C, UInt8(B))
end
eltype(::TranslationalBasis{T}) where T = T
#-------------------------------------------------------------------------------------------------------------------------
function index(b::TranslationalBasis{Tc})::Tuple{Tc, Int} where Tc
    Im, T = translation_index(b.dgt, b.B)
    ind = binary_search(b.I, Im)
    if iszero(ind)
        zero(Tc), 1
    else
        N::Tc = b.C[T+1] * b.R[ind]
        N, ind
    end
end

#-------------------------------------------------------------------------------------------------------------------------
# Judge
#-------------------------------------------------------------------------------------------------------------------------
struct TranslationJudge
    F                   # Projective selection
    K::Int              # Momentum
    B::Int              # Base
    C::Vector{Float64}  # Normalization coefficient
end

function (judge::TranslationJudge)(dgt::AbstractVector{<:Integer}, i::Integer)
    Q, N = if judge.F(dgt)
        c, r = translation_check(dgt, i, judge.B)
        if c && iszero( mod(r * judge.K, length(dgt) ) )
            (true, judge.C[r])
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
export translationalbasis
function translationalbasis(f, k::Integer, L::Integer; base::Integer=2, alloc::Integer=1000, threaded::Bool=false)
    dgt = zeros(BITTYPE, L)
    judge = TranslationJudge(f, k, base, [L/sqrt(i) for i = 1:L])
    I, R = if threaded
        selectindexnorm_threaded(judge, L, base=base, alloc=alloc)
    else
        selectindexnorm(judge, L, 1:Int(base)^L, base=base, alloc=alloc)
    end
    C = if iszero(k)
        fill(1.0, L)
    elseif isequal(2k, L)
        [iseven(i) ? 1.0 : -1.0 for i=0:L-1] 
    else
        [exp(-1im*2π/L*k*i) for i=0:L-1]
    end
    TranslationalBasis(dgt, I, R, C, base)
end
#-------------------------------------------------------------------------------------------------------------------------
function translationalbasis(k::Integer, L::Integer; base::Integer=2, alloc::Integer=1000, threaded::Bool=false)
    translationalbasis(x -> true, k, L, base=base, alloc=alloc, threaded=threaded)
end

#-------------------------------------------------------------------------------------------------------------------------
# Reshape
#-------------------------------------------------------------------------------------------------------------------------
function schmidt!(target::AbstractMatrix, v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::TranslationalBasis)
    dgt, R, phase = b.dgt, b.R, b.C[2]
    Binds = complement(Ainds, length(dgt))
    for i = 1:length(v)
        change!(dgt, b.I[i], base=b.B)
        cycle_write!(target, dgt, v[i] / R[i], phase, Ainds, Binds, b.B)
    end
    target
end

function cycle_write!(
    target::AbstractMatrix, dgt::AbstractVector{<:Integer}, val::Number, phase::Number,
    Ainds::AbstractVector{<:Integer}, Binds::AbstractVector{<:Integer}, base::Integer
)
    perm_element!(target, dgt, val, Ainds, Binds, base)
    for _ = 1:length(dgt)-1
        val *= phase
        cyclebits!(dgt)
        perm_element!(target, dgt, val, Ainds, Binds, base)
    end
end
