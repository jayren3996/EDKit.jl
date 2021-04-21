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
function index(b::ProjectedBasis)::Tuple{Int, Int}
    i = binary_search(b.I, b.dgt)
    if i > 0
        1, i
    else
        error("No such symmetry.")
    end
end
#-----------------------------------------------------------------------------------------------------
# Construction
#-----------------------------------------------------------------------------------------------------
export projectedbasis
function projectedbasis(f, L::Integer; base::Integer=2, alloc::Integer=1000, threaded::Bool=false)
    dgt = zeros(Int, L)
    I = if threaded
        selectindex_threaded(f, L, base=base, alloc=alloc)
    else
        selectindex(f, L, 1:base^L, base=base, alloc=alloc)
    end
    ProjectedBasis(dgt, I, base)
end
