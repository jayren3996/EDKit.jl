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
