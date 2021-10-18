"""
    index(dgt::AbstractVector{T}; base::Integer=2) where T <: Integer

Convert a digits to the integer index using the relation: 

`number = ∑ᵢ bits[i] * base^(L-i) + 1`

The summaton is evaluate using the efficieint polynomial evaluation method.
"""
@inline function index(dgt::AbstractVector{T}; base::Integer=2) where T <: Integer
    N = zero(T)
    for i = 1:length(dgt)
        N *= base
        N += dgt[i]
    end
    N + one(T)
end

"""
    index(dgt::AbstractVector{T}, sites::AbstractVector{<:Integer}; base::Integer=2) where T <: Integer

Convert a sub-digits (subarray of `dgt`) to the integer index.
"""
@inline function index(dgt::AbstractVector{T}, sites::AbstractVector{<:Integer}; base::Integer=2) where T <: Integer
    N = zero(T)
    for i in sites
        N *= base
        N += dgt[i]
    end
    N + one(T)
end

"""
    change!(dgt::AbstractVector{T}, ind::Integer; base::Integer=2) where T <: Integer

Change the digits to that with the target index.  
This method is the inverse of `index`.
"""
@inline function change!(dgt::AbstractVector{T}, ind::Integer; base::Integer=2) where T <: Integer
    N = convert(T, ind) - one(T)
    for i = length(dgt):-1:1
        N, dgt[i] = divrem(N, base)
    end
end

"""
    change!(dgt::AbstractVector{T}, sites::AbstractVector{<:Integer}, ind::Integer; base::Integer=2) where T <: Integer

Change the sub-digits to that with the target index.  
This method is the inverse of `index`.
"""
@inline function change!(dgt::AbstractVector{T}, sites::AbstractVector{<:Integer}, ind::Integer; base::Integer=2) where T <: Integer
    N = convert(T, ind) - one(T)
    for i = length(sites):-1:1
        N, dgt[sites[i]] = divrem(N, base)
    end
end

"""
    binary_search(list::AbstractVector{<:Integer}, i<:Integer)

Return the position of i in a sorted list using binary search algorithm.
"""
function binary_search(list::AbstractVector{<:Integer}, i::Integer)
    l::Int = 1
    r::Int = length(list)
    c::Int = (l + r) ÷ 2
    while true
        t = list[c]
        (i < t) ? (r = c - 1) : (i > t) ? (l = c + 1) : break
        (l > r) ? (c = 0; break) : (c = (l + r) ÷ 2)
    end
    c
end

"""
    selectindex(f, L, rg; base=2, alloc=1000)

Select the legal indices for a basis.

Inputs:
-------
- `f`    : Function for digits that tells whether a digits is valid.
- `L`    : Length of the system.
- `rg`   : Range of interation.
- `base` : (Optional) Base, `base=2` by default.
- `alloc`: (Optional) Pre-allocation memory for the list of indices, `alloc=1000` by default.

Outputs:
--------
- `I`: List of indices in a basis.
"""
function selectindex(
    f, 
    L::Integer, 
    rg::UnitRange; 
    base::Integer=2, 
    alloc::Integer=1000
)
    dgt = zeros(Int64, L)
    I = Int[]
    sizehint!(I, alloc)
    for i in rg
        change!(dgt, i, base=base)
        f(dgt) && append!(I, i)
    end
    I
end

"""
    selectindexnorm(f, L, rg; base=2, alloc=1000)

Select the legal indices and corresponding norm for a basis.

Inputs:
-------
- `f`    : Function for digits that tells whether a digits is valid as well as its normalization.
- `L`    : Length of the system.
- `rg`   : Range of interation.
- `base` : (Optional) Base, `base=2` by default.
- `alloc`: (Optional) Pre-allocation memory for the list of indices, `alloc=1000` by default.

Outputs:
--------
- `I`: List of indices in a basis.
- `R`: List of normalization for each states.
"""
function selectindexnorm(
    f, 
    L::Integer, 
    rg::UnitRange; 
    base::Integer=2, 
    alloc::Integer=1000
)
    dgt = zeros(Int64, L)
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

"""
    dividerange(maxnum::Integer, nthreads::Integer)

Divide the interation range equally according to the number of threads.

Inputs:
-------
- `maxnum`  : Maximum number of range.
- `nthreads`: Number of threads.

Outputs:
--------
- `list`: List of ranges.
"""
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

"""
    selectindex_threaded(f, L, rg; base=2, alloc=1000)

Select the legal indices with multi-threads.

Inputs:
-------
- `f`    : Function for digits that tells whether a digits is valid as well as its normalization.
- `L`    : Length of the system.
- `base` : (Optional) Base, `base=2` by default.
- `alloc`: (Optional) Pre-allocation memory for the list of indices, `alloc=1000` by default.

Outputs:
--------
- `I`: List of indices in a basis.
"""
function selectindex_threaded(
    f, 
    L::Integer; 
    base::Integer=2, 
    alloc::Integer=1000
)
    nt = Threads.nthreads()
    ni = dividerange(base^L, nt)
    nI = Vector{Vector{Int}}(undef, nt)
    Threads.@threads for ti in 1:nt
        nI[ti] = selectindex(f, L, ni[ti], base=base, alloc=alloc)
    end
    I = vcat(nI...)
    I
end

"""
    selectindexnorm_threaded(f, L, rg; base=2, alloc=1000)

Select the legal indices and corresponding norm for a basis with multi-threads.

Inputs:
-------
- `f`    : Function for digits that tells whether a digits is valid as well as its normalization.
- `L`    : Length of the system.
- `base` : (Optional) Base, `base=2` by default.
- `alloc`: (Optional) Pre-allocation memory for the list of indices, `alloc=1000` by default.

Outputs:
--------
- `I`: List of indices in a basis.
- `R`: List of normalization for each states.
"""
function selectindexnorm_threaded(
    f, 
    L::Integer; 
    base::Integer=2, 
    alloc::Integer=1000
)
    nt = Threads.nthreads()
    ni = dividerange(base^L, nt)
    nI = Vector{Vector{Int}}(undef, nt)
    nR = Vector{Vector{Float64}}(undef, nt)
    Threads.@threads for ti in 1:nt
        nI[ti], nR[ti] = selectindexnorm(f, L, ni[ti], base=base, alloc=alloc)
    end
    I, R = vcat(nI...), vcat(nR...)
    I, R
end


# Level statistics
export gapratio, meangapratio
function gapratio(E::AbstractVector{<:Real})
    dE = diff(E)
    r = zeros(length(dE)-1)
    for i = 1:length(r)
        if dE[i] < dE[i+1]
            r[i] = dE[i]/dE[i+1]
        elseif dE[i] > dE[i+1]
            r[i] = dE[i+1]/dE[i]
        else
            r[i] = 1.0
        end
    end
    r
end

meangapratio(E::AbstractVector{<:Real}) = sum(gapratio(E)) / (length(E) - 2)