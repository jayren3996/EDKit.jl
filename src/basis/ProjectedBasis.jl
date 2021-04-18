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
eltype(::ProjectedBasis) = Int
#-----------------------------------------------------------------------------------------------------
function change!(b::ProjectedBasis, i::Integer)::Int
    change!(b.dgt, b.I[i], base=b.B)
    1
end
#-----------------------------------------------------------------------------------------------------
function index(b::ProjectedBasis)::Tuple{Int, Int}
    ind = binary_search(b.I, index(b.dgt, base=b.B))
    1, ind
end
#-----------------------------------------------------------------------------------------------------
size(b::ProjectedBasis, i::Integer) = (i == 2 || i == 1) ? length(b.I) : 1
size(b::ProjectedBasis) = (l=length(b.I); (l,l))

#-----------------------------------------------------------------------------------------------------
# Construct basis
#-----------------------------------------------------------------------------------------------------
export projectedbasis
function projectedbasis(f, L::Integer; base::Integer=2)
    dgt = zeros(Int, L)
    I = []
    for i = 1:base^L
        change!(dgt, i, base=base)
        if f(dgt)
            append!(I, i)
        end
    end
    ProjectedBasis(dgt, I, base)
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
        if i < t
            r = c - 1
        elseif i > t
            l = c + 1
        else
            break
        end
        if l > r
            return 0
        else
            c = (l+r) ÷ 2
        end
    end
    c
end
