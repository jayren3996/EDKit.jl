#-------------------------------------------------------------------------------------------------------------------------------------------------------
# Abstract Basis Type
#-------------------------------------------------------------------------------------------------------------------------------------------------------
abstract type AbstractBasis end
eltype(::AbstractBasis) = ComplexF64

#-------------------------------------------------------------------------------------------------------------------------------------------------------
change!(b::AbstractBasis, i::Integer) = (change!(b.dgt, b.I[i], base=b.B); b.R[i])
size(b::AbstractBasis, i::Integer) = (i == 1 || i == 2) ? length(b.I) : 1
size(b::AbstractBasis) = (l=length(b.I); (l,l))
length(b::AbstractBasis) = length(b.dgt)

#-------------------------------------------------------------------------------------------------------------------------------------------------------
# Digits ⟷ Index
#-------------------------------------------------------------------------------------------------------------------------------------------------------
"""
    index(dgt::AbstractVector{<:Integer})

Digit ⟹ index. 

The method compute the polynomial

number = bits[i] * base^(L-i) + 1

in the most efficient way.
"""
function index(dgt::AbstractVector{<:Integer}; base::Integer=2)::Int
    N = 0
    for i = 1:length(dgt)
        N = muladd(base, N, dgt[i])
    end
    N + 1
end
#-------------------------------------------------------------------------------------------------------------------------------------------------------
function index(dgt::AbstractVector{<:Integer}, sites::AbstractVector{<:Integer}; base::Integer=2)::Int
    N = 0
    for i in sites
        N = muladd(base, N, dgt[i])
    end
    N + 1
end
#-------------------------------------------------------------------------------------------------------------------------------------------------------
"""
    change!(dgt::Digit, index::Integer)

Index ⟹ Digit. 

The method compute the bits vector and write to bits.
"""
function change!(dgt::AbstractVector{<:Integer}, ind::Integer; base::Integer=2)
    N = ind - 1
    for i = length(dgt):-1:1
        N, dgt[i] = divrem(N, base)
    end
end
#-------------------------------------------------------------------------------------------------------------------------------------------------------
function change!(dgt::AbstractVector{<:Integer}, sites::AbstractVector{<:Integer}, ind::Integer; base::Integer=2)
    N = ind - 1
    for i = length(sites):-1:1
        N, dgt[sites[i]] = divrem(N, base)
    end
end

#-------------------------------------------------------------------------------------------------------------------------------------------------------
# Search index
#-------------------------------------------------------------------------------------------------------------------------------------------------------
"""
    binary_search(list::AbstractVector{<:Integer}, i<:Integer)

Return the position of i in a sorted list using binary search algorithm.
"""
function binary_search(list::AbstractVector{<:Integer}, i::Integer)
    l, r = 1, length(list)
    c = (l+r) ÷ 2
    while true
        t = list[c]
        (i < t) ? (r = c - 1) : (i > t) ? (l = c + 1) : break
        (l > r) ? (c = 0; break) : (c = (l + r) ÷ 2)
    end
    c
end

#-------------------------------------------------------------------------------------------------------------------------------------------------------
# Select basis & norm
#-------------------------------------------------------------------------------------------------------------------------------------------------------
function selectindex(f, L::Integer, rg::UnitRange; base::Integer=2, alloc::Integer=1000)
    dgt = zeros(Int, L)
    I = Int[]
    sizehint!(I, alloc)
    for i in rg
        change!(dgt, i, base=base)
        if f(dgt)
            append!(I, i)
        end
    end
    I
end
#-------------------------------------------------------------------------------------------------------------------------------------------------------
function selectindexnorm(f, L::Integer, rg::UnitRange; base::Integer=2, alloc::Integer=1000)
    dgt = zeros(Int, L)
    I, R = Int[], Float64[]
    sizehint!(I, alloc)
    sizehint!(R, alloc)
    for i in rg
        change!(dgt, i, base=base)
        Q, N = f(dgt, i)
        if Q
            append!(I, i)
            append!(R, N)
        end
    end
    I, R
end
#-------------------------------------------------------------------------------------------------------------------------------------------------------
# Multi-threads
function dividerange(maxnum::Integer, nthreads::Integer)
    list = Vector{UnitRange{Int64}}(undef, nthreads)
    eachthreads, left = divrem(maxnum, nthreads)
    start = 1
    for i = 1:left
        stop  = start + eachthreads
        list[i] = start:stop
        start = stop+1
        
    end
    for i = left+1:nthreads-1
        stop  = start + eachthreads - 1
        list[i] = start:stop
        start = stop+1
    end
    list[nthreads] = start:maxnum
    list
end
#-------------------------------------------------------------------------------------------------------------------------------------------------------
function selectindex_threaded(f, L::Integer; base::Integer=2, alloc::Integer=1000)
    nt = Threads.nthreads()
    ni = dividerange(base^L, nt)
    nI = Vector{Vector{Int}}(undef, nt)
    Threads.@threads for ti in 1:nt
        nI[ti] = selectindex(f, L, ni[ti], base=base, alloc=alloc)
    end
    vcat(nI...)
end
#-------------------------------------------------------------------------------------------------------------------------------------------------------
function selectindexnorm_threaded(f, L::Integer; base::Integer=2, alloc::Integer=1000)
    nt = Threads.nthreads()
    ni = dividerange(base^L, nt)
    nI = Vector{Vector{Int}}(undef, nt)
    nR = Vector{Vector{Float64}}(undef, nt)
    Threads.@threads for ti in 1:nt
        nI[ti], nR[ti] = selectindexnorm(f, L, ni[ti], base=base, alloc=alloc)
    end
    vcat(nI...), vcat(nR...)
end

