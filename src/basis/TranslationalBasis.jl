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
    Im, T = indmin(b.dgt, b.B)
    ind = binary_search(b.I, Im)
    if ind == 0
        return parse(eltype(b), "0"), 1
    else
        N = b.C ^ T * b.R[ind]
        return N, ind
    end
end
change!(b::TranslationalBasis, i::Integer) = (change!(b.dgt, b.I[i], base=b.B); b.R[i])
size(b::TranslationalBasis, i::Integer) = (i == 2 || i == 1) ? length(b.I) : 1
size(b::TranslationalBasis) = (l=length(b.I); (l,l))

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Construct
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
export translationalbasis
function translationalbasis(
    f, k::Integer, L::Integer; 
    base::Integer=2, alloc::Integer=1000, threaded::Bool=false
)
    dgt = zeros(Int, L)
    cl = [L/sqrt(i) for i = 1:L]
    I, R = if threaded
        translationsearch_threaded(f, k, L, cl, base=base, alloc=alloc)
    else
        translationsearch(f, k, L, cl, 1:base^L, base=base, alloc=alloc)
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

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Search Basis
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function translationsearch(
    f, k::Integer, L::Integer, cl::Vector, rg::UnitRange; 
    base::Integer=2, alloc::Integer=1000
)
    dgt = zeros(Int, L)
    I, R = allocvector(Int, alloc), allocvector(Float64, alloc)
    for i in rg
        change!(dgt, i, base=base)
        if f(dgt)
            c, r = checkstate(dgt, i, base)
            if c && (k * r % L == 0)
                append!(I, i)
                append!(R, cl[r])
            end
        end
    end
    Vector(I), Vector(R)
end
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Multi-threads
function translationsearch_threaded(
    f, k::Integer, L::Integer, cl::Vector; 
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
        nI[ti], nR[ti] = translationsearch(f, k, L, cl, ni[ti], base=base, alloc=alloc)
    end
    vcat(nI...), vcat(nR...)
end
