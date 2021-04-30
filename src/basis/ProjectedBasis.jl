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
change!(b::ProjectedBasis, i::Integer) = (change!(b.dgt, b.I[i], base=b.B); 1)
#-----------------------------------------------------------------------------------------------------
function index(b::ProjectedBasis)::Tuple{Int, Int}
    i = index(b.dgt, base=b.B)
    ind = binary_search(b.I, i)
    ind > 0 ? (1, ind) : error("No such symmetry.")
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

#-----------------------------------------------------------------------------------------------------
# To Vector
#-----------------------------------------------------------------------------------------------------
function schmidt!(target::AbstractMatrix, v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::ProjectedBasis)
    Binds = Int[i for i = 1:length(b.dgt) if !in(i, Ainds)]
    dgt = b.dgt
    for i = 1:length(v)
        change!(dgt, b.I[i], base=b.B)
        ia = index(dgt, Ainds, base=b.B)
        ib = index(dgt, Binds, base=b.B)
        target[ia, ib] = v[i]
    end
    target
end
