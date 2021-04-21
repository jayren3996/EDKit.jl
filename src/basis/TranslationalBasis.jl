#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Translational Basis
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
struct TranslationalBasis{T <: Number} <: AbstractBasis
    dgt::Vector{Int}
    I::Vector{Int}
    R::Vector{Float64}
    C::T
    B::Int
end
eltype(::TranslationalBasis{T}) where T = T
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function index(b::TranslationalBasis)
    Im, T = translation_index(b.dgt, b.B)
    ind = binary_search(b.I, Im)
    if ind == 0
        return parse(eltype(b), "0"), 1
    else
        N = b.C ^ T * b.R[ind]
        return N, ind
    end
end

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Construct
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
export translationalbasis
function translationalbasis(
    f, k::Integer, L::Integer; 
    base::Integer=2, alloc::Integer=1000, threaded::Bool=false
)
    dgt = zeros(Int, L)
    C = [L/sqrt(i) for i = 1:L]
    info = TranslationInfo(f, k, base, C)
    I, R = if threaded
        selectindexnorm_threaded(info, L, base=base, alloc=alloc)
    else
        selectindexnorm(info, L, 1:base^L, base=base, alloc=alloc)
    end
    expk = (k == 0) ? 1.0 : (2k == L) ? -1.0 : exp(-1im*2π/L*k)
    TranslationalBasis(dgt, I, R, expk, base)
end
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function translationalbasis(
    k::Integer, L::Integer; 
    base::Integer=2, alloc::Integer=1000, threaded::Bool=false
)
    translationalbasis(x->true, k, L, base=base, alloc=alloc, threaded=threaded)
end

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Helper Functions
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function cyclebits!(dgt::AbstractVector{<:Integer})
    p = dgt[end]
    for i = 1:length(dgt)
        p, dgt[i] = dgt[i], p
    end
    dgt
end
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function translation_check(dgt::AbstractVector{<:Integer}, I0::Integer, base::Integer)
    cyclebits!(dgt)
    for i=1:length(dgt)-1
        if (In = index(dgt, base=base)) < I0
            change!(dgt, I0, base=base)
            return false, 0
        elseif In == I0
            return true, i
        end
        cyclebits!(dgt)
    end
    true, length(dgt)
end
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function translation_index(dgt::AbstractVector{<:Integer}, base::Integer)
    I0 = index(dgt, base=base)
    Im, T = I0, 0
    cyclebits!(dgt)
    for i=1:length(dgt)-1
        if (In = index(dgt, base=base)) == I0
            break
        elseif In < Im
            Im, T = In, i
        end
        cyclebits!(dgt)
    end
    Im, T
end

