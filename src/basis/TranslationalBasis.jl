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
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
eltype(::TranslationalBasis{T}) where T = T
change!(b::TranslationalBasis, i::Integer) = (change!(b.dgt, b.I[i], base=b.B); b.R[i])
function index(b::TranslationalBasis)
    Im, T = indmin(b.dgt, b.B)
    ind = binary_search(b.I, Im)
    if ind == 0
        return parse(eltype(b), "0"), 1
    else
        N = b.C ^ T * b.R[ind]
        return N, ind
    end
end
size(b::TranslationalBasis, i::Integer) = (i == 2 || i == 1) ? length(b.I) : 1
size(b::TranslationalBasis) = (l=length(b.I); (l,l))

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Construct
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
export translationalbasis
function translationalbasis(f, k::Integer, L::Integer; base::Integer=2, alloc::Integer=1000)
    dgt = zeros(Int, L)
    cl = [L/sqrt(i) for i = 1:L]
    I, R = allocvector(Int, alloc), allocvector(Float64, alloc)
    for i = 1:base^L
        change!(dgt, i, base=base)
        if f(dgt)
            c, r = checkstate(dgt, i, base)
            if c && (k * r % L == 0)
                append!(I, i)
                append!(R, cl[r])
            end
        end
    end
    expk = (k == 0) ? 1.0 : (2k == L) ? -1.0 : exp(-1im*2π/L*k)
    TranslationalBasis(dgt, Vector(I), Vector(R), expk, base)
end

translationalbasis(k::Integer, L::Integer; base::Integer=2, alloc::Integer=1000) = translationalbasis(x->true, k, L, base=base, alloc=alloc)
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
function checkstate(dgt::AbstractVector{<:Integer}, I0::Integer, base::Integer)
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
function indmin(dgt::AbstractVector{<:Integer}, base::Integer)
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
