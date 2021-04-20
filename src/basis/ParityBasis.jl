#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Parity Bases
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
struct TranslationParityBasis <: AbstractBasis
    dgt::Vector{Int}
    I::Vector{Int}
    R::Vector{Float64}
    K::Bool
    P::Bool
    B::Int
end
eltype(::TranslationParityBasis) = Float64
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function index(b::TranslationParityBasis)
    p, i, t = indparity(reverse!, b.dgt, b.B)
    ind = binary_search(b.I, i)
    if ind == 0
        0.0, 1
    else
        N = if p
            K = 2 * b.K - 1
            K ^ t
        else
            K, P = 2 * b.K - 1, 2 * b.P - 1
            P * K ^ t
        end
        N * b.R[ind], ind
    end
end
change!(b::TranslationParityBasis, i::Integer) = (change!(b.dgt, b.I[i], base=b.B); b.R[i])
size(b::TranslationParityBasis, i::Integer) = (i == 2 || i == 1) ? length(b.I) : 1
size(b::TranslationParityBasis) = (l=length(b.I); (l,l))

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Construction
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
export translationparitybasis
function translationparitybasis(f, k::Int, p::Int, L::Integer; base::Integer=2, alloc::Integer=1000, threaded=false)
    K = (k == 0) ? true : (2k == L) ? false : error("Momentum incompatible with parity.")
    P = (p == 1) ? true : (p == -1) ? false : error("Invalid parity $p.")
    dgt = zeros(Int, L)
    sqrtl = [L/sqrt(i/2) for i = 1:L]
    cospl = [sqrt(1+cos(2π*i/L)) for i = 0:L-1]
    cosml = [sqrt(1-cos(2π*i/L)) for i = 0:L-1]
    I, R = if threaded
        paritysearch_threaded(reverse!, f, k, P, L, sqrtl, cospl, cosml, base=base, alloc=alloc)
    else
        paritysearch(reverse!, f, k, P, L, sqrtl, cospl, cosml, 1:base^L, base=base, alloc=alloc)
    end
    TranslationParityBasis(dgt, I, R, K, P, base)
end

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Helper Functions
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function checkparity(parity!, dgt::AbstractVector{<:Integer}, I0::Integer, base::Integer)::Tuple{Bool, Bool, Int, Int}
    b1, r = checkstate(dgt, I0, base)
    if !b1 
        return (false, false, r, 0) 
    end
    parity!(dgt)
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
    parity!(dgt)
    true, false, r, 0
end
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function indparity(parity!, dgt::AbstractVector{<:Integer}, base::Integer)
    Im1, T1 = indmin(dgt, base)
    parity!(dgt)
    Im2, T2 = indmin(dgt, base)
    parity!(dgt)
    if Im2 < Im1
        false, Im2, T2
    else
        true, Im1, T1
    end
end
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function paritynorm(b::Bool, k::Int, p::Bool, r::Int, m::Int, L::Integer, sqrtl::Vector{Float64}, cospl::Vector{Float64}, cosml::Vector{Float64})
    Q, N = if b
        km = k * m
        if p
            (2 * rem(km, L) == L) ? (false, 0.0) : (true, cospl[km % L + 1] * sqrtl[r])
        else
            (km % L == 0) ? (false, 0.0) : (true, cosml[km % L + 1] * sqrtl[r])
        end
    else
        true, sqrtl[r]
    end
end

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Search Basis
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function paritysearch(
    parity!, f, k::Int, p::Bool, L::Integer, sqrtl::Vector, cospl::Vector, cosml::Vector, rg::UnitRange; 
    base::Integer=2, alloc::Integer=1000
)
    dgt = zeros(Int, L)
    I, R = allocvector(Int, alloc), allocvector(Float64, alloc)
    for i in rg
        change!(dgt, i, base=base)
        if f(dgt)
            b1, b2, r, m = checkparity(parity!, dgt, i, base)
            if b1 && (k * r % L == 0)
                Q, N = paritynorm(b2, k, p, r, m, L, sqrtl, cospl, cosml)
                if Q
                    append!(I, i)
                    append!(R, N)
                end
            end
        end
    end
    Vector(I), Vector(R)
end
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Multi-threads
function paritysearch_threaded(
    parity!, f, k::Int, p::Bool, L::Integer, sqrtl::Vector, cospl::Vector, cosml::Vector;
    base::Integer=2, alloc::Integer=1000
)
    nt = Threads.nthreads()
    ni = begin
        list = Vector{UnitRange{Int64}}(undef, nt)
        N = base^L
        eacht = N ÷ nt
        for i = 1:nt-1
            list[i] = (i-1)*eacht+1:i*eacht
        end
        list[nt] = eacht*(nt-1)+1:base^L
        list
    end
    nI, nR = Vector{Vector{Int}}(undef, nt), Vector{Vector{Float64}}(undef, nt)
    Threads.@threads for ti in 1:nt
        nI[ti], nR[ti] = paritysearch(parity!, f, k, p, L, sqrtl, cospl, cosml, ni[ti], base=base, alloc=alloc)
    end
    vcat(nI...), vcat(nR...)
end
