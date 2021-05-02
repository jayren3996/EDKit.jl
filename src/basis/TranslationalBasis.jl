#-------------------------------------------------------------------------------------------------------------------------
# TranslationalBasis
#-------------------------------------------------------------------------------------------------------------------------
"""
    TranslationalBasis

Basis for subspace that is spanned by momentum states, which can also incorporate projected restriction.

Properties:
-----------
- dgt : Vector{Int}
- I   : Vector{Int}, list of indicies.
- R   : Vector{Float64}, list of normalization.
- C   : Vector{Float64/ComplexF64}, unit phase factor.
- B   : Int, levels on each site.
"""
struct TranslationalBasis{T <: Number} <: AbstractBasis
    dgt::Vector{Int}
    I::Vector{Int}
    R::Vector{Float64}
    C::Vector{T}
    B::Int
end
eltype(::TranslationalBasis{T}) where T = T
#-------------------------------------------------------------------------------------------------------------------------
function index(b::TranslationalBasis{Tc})::Tuple{Tc, Int} where Tc
    Im, T = translation_index(b.dgt, b.B)
    ind = binary_search(b.I, Im)
    if ind == 0
        zero(Tc), 1
    else
        N::Tc = (T == 0) ? b.R[ind] : b.C[T] * b.R[ind]
        N, ind
    end
end

#-------------------------------------------------------------------------------------------------------------------------
# Judge
#-------------------------------------------------------------------------------------------------------------------------
struct TranslationInfo
    F
    K::Int
    B::Int
    C::Vector{Float64}
    TranslationInfo(F, K::Integer, B::Integer, C::AbstractVector{<:Real}) = new(F, Int(K), Int(B), Vector{Float64}(C))
end
function (info::TranslationInfo)(dgt::Vector{Int}, i::Integer)::Tuple{Bool, Float64}
    if info.F(dgt)
        c, r = translation_check(dgt, i, info.B)
        (Q = c && (info.K * r % length(dgt) == 0)) ? (Q, info.C[r]) : (Q, 0.0)
    else
        false, 0.0
    end
end

#-------------------------------------------------------------------------------------------------------------------------
# Construction
#-------------------------------------------------------------------------------------------------------------------------
export translationalbasis
function translationalbasis(f, k::Integer, L::Integer; base::Integer=2, alloc::Integer=1000, threaded::Bool=false)
    dgt = zeros(Int, L)
    info = TranslationInfo(f, k, base, [L/sqrt(i) for i = 1:L])
    I, R = threaded ? selectindexnorm_threaded(info, L, base=base, alloc=alloc) : selectindexnorm(info, L, 1:base^L, base=base, alloc=alloc)
    C = (k == 0) ? fill(1.0, L-1) : (2k == L) ? [iseven(i) ? 1.0 : -1.0 for i=1:L-1] : [exp(-1im*2π/L*k*i) for i=1:L-1]
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
    dgt, R, C = b.dgt, b.R, b.C[1]
    Binds = complement(Ainds, length(dgt))
    for i = 1:length(v)
        change!(dgt, b.I[i], base=b.B)
        cycle_write!(target, dgt, v[i] / R[i], C, Ainds, Binds, b.B)
    end
    target
end

function cycle_write!(
    target::AbstractMatrix, dgt::AbstractVector{<:Integer}, val::Number, phase::Number,
    Ainds::AbstractVector{<:Integer}, Binds::AbstractVector{<:Integer}, base::Integer
)
    perm_element!(target, dgt, val, Ainds, Binds, base)
    for j = 1:length(dgt)-1
        val *= phase
        cyclebits!(dgt)
        perm_element!(target, dgt, val, Ainds, Binds, base)
    end
end