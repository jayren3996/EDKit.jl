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

