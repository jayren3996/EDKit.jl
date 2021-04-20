#-----------------------------------------------------------------------------------------------------
# Projected Basis
#-----------------------------------------------------------------------------------------------------
"""
    ProjectedBasis

Basis for subspace that is spanned only by product states.

Properties:
-----------
- dgt : Vector{Int}
- I   : Vector{Int}, list of indicies.
- B   : Int, levels on each site.
"""
struct ProjectedBasis <: AbstractBasis
    dgt::Vector{Int}
    I::Vector{Int}
    B::Int
end
#-----------------------------------------------------------------------------------------------------
eltype(::ProjectedBasis) = Int
change!(b::ProjectedBasis, i::Integer) = (change!(b.dgt, b.I[i], base=b.B); 1)
index(b::ProjectedBasis) = 1, binary_search(b.I, index(b.dgt, base=b.B))
size(b::ProjectedBasis, i::Integer) = (i == 2 || i == 1) ? length(b.I) : 1
size(b::ProjectedBasis) = (l=length(b.I); (l,l))

#-----------------------------------------------------------------------------------------------------
# Pre-allocated vector
#-----------------------------------------------------------------------------------------------------
mutable struct AllocVector{Tv}
    V::Vector{Tv}
    P::Int
    A::Int
end
allocvector(T::DataType, alloc::Int) = AllocVector(Vector{T}(undef, alloc), 0, alloc)
function append!(av::AllocVector{Tv}, e::Tv) where Tv
    if av.P == length(av.V)
        av.V = vcat(av.V, Vector{Tv}(undef, av.A))
    end
    av.P += 1
    av.V[av.P] = e 
end
Vector(av::AllocVector) = deleteat!(av.V, av.P+1 : length(av.V))

#-----------------------------------------------------------------------------------------------------
# Construction
#-----------------------------------------------------------------------------------------------------
export projectedbasis
function projectedbasis(f, L::Integer; base::Integer=2, alloc::Integer=1000)
    dgt = zeros(Int, L)
    I = allocvector(Int, alloc)
    for i = 1:base^L
        change!(dgt, i, base=base)
        if f(dgt)
            append!(I, i)
        end
    end
    ProjectedBasis(dgt, Vector(I), base)
end

#-----------------------------------------------------------------------------------------------------
# Helper functions
#-----------------------------------------------------------------------------------------------------
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
