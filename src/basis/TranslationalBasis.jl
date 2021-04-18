#-----------------------------------------------------------------------------------------------------
# Translational Basis
#-----------------------------------------------------------------------------------------------------
struct TranslationalBasis{T <: Number} <: AbstractBasis
    dgt::Vector{Int}
    I::Vector{Int}
    R::Vector{Float64}
    C::T
    B::Int
end
eltype(::TranslationalBasis{T}) where T = T
#-----------------------------------------------------------------------------------------------------
function change!(b::TranslationalBasis, i::Integer)::Float64
    change!(b.dgt, b.I[i], base=b.B)
    1 / b.R[i]
end
#-----------------------------------------------------------------------------------------------------
function index(b::TranslationalBasis)::Tuple{ComplexF64, Int}
    Im, T = indmin(b.dgt, b.B)
    ind = binary_search(b.I, Im)
    N = b.C ^ T * b.R[ind]
    N, ind
end
#-----------------------------------------------------------------------------------------------------
size(b::TranslationalBasis, i::Integer) = (i == 2 || i == 1) ? length(b.I) : 1
size(b::TranslationalBasis) = (l=length(b.I); (l,l))
#-----------------------------------------------------------------------------------------------------
# Construct
#-----------------------------------------------------------------------------------------------------
export translationalbasis
function translationalbasis(f, k::Integer, L::Integer; base::Integer=2)
    dgt = zeros(Int, L)
    I, R = Int[], Float64[]
    for i = 1:base^L
        change!(dgt, i)
        if f(dgt)
            c, r = checkstate(dgt, i, base)
            if c && (k * r % L == 0)
                append!(I, i)
                append!(R, L/sqrt(r))
            end
        end
    end
    expk = exp(-1im*2π/L*k)
    C = imag(expk) == 0 ? real(expk) : expk
    TranslationalBasis(dgt, I, R, C, base)
end

#-----------------------------------------------------------------------------------------------------
# Helper Functions
#-----------------------------------------------------------------------------------------------------
function movebits!(dgt::AbstractVector{<:Integer})
    p = dgt[end]
    for i = 1:length(dgt)
        p, dgt[i] = dgt[i], p
    end
    dgt
end
#-----------------------------------------------------------------------------------------------------
function checkstate(dgt::AbstractVector{<:Integer}, I0::Integer, base::Integer)
    for i=1:length(dgt)-1
        movebits!(dgt)
        In = index(dgt, base=base)
        if In < I0
            change!(dgt, I0, base=base)
            return false, 0
        elseif In == I0
            return true, i
        end
    end
    movebits!(dgt)
    true, length(dgt)
end
#-----------------------------------------------------------------------------------------------------
function indmin(dgt::AbstractVector{<:Integer}, base::Integer)
    I0 = index(dgt, base=base)
    Im, T = I0, 0
    for i=1:length(dgt)-1
        movebits!(dgt)
        In = index(dgt, base=base)
        if In == Im
            return Im, T
        elseif In < Im
            Im, T = In, i
        end
    end
    movebits!(dgt)
    Im, T
end
