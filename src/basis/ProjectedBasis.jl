#-------------------------------------------------------------------------------------------------------------------------
# Projected Basis
#-------------------------------------------------------------------------------------------------------------------------
"""
    ProjectedBasis

Basis for subspace that is spanned only by product states.

Properties:
-----------
- dgt : Vector{BITTYPE}
- I   : Vector{Int}, list of indicies.
- B   : Int, levels on each site.
"""
struct ProjectedBasis <: AbstractBasis
    dgt::Vector{BITTYPE}
    I::Vector{Int}
    B::UInt64
    ProjectedBasis(dgt::Vector{BITTYPE}, I::Vector{Int}, B::Integer) = new(dgt, I, UInt64(B))
end
eltype(::ProjectedBasis) = Int
change!(b::ProjectedBasis, i::Integer) = (change!(b.dgt, b.I[i], base=b.B); 1)
#-------------------------------------------------------------------------------------------------------------------------
function index(b::ProjectedBasis)::Tuple{Int, Int}
    i = index(b.dgt, base=b.B)
    ind = binary_search(b.I, i)
    ind > 0 ? (1, ind) : error("No such symmetry.")
end

#-------------------------------------------------------------------------------------------------------------------------
# Construction
#-------------------------------------------------------------------------------------------------------------------------
export projectedbasis
function projectedbasis(f, L::Integer; base::Integer=2, alloc::Integer=1000, threaded::Bool=false)
    dgt = zeros(BITTYPE, L)
    I = threaded ? selectindex_threaded(f, L, base=base, alloc=alloc) : selectindex(f, L, 1:base^L, base=base, alloc=alloc)
    ProjectedBasis(dgt, I, base)
end

#-------------------------------------------------------------------------------------------------------------------------
# Schmidt form
#-------------------------------------------------------------------------------------------------------------------------
function schmidt!(target::AbstractMatrix, v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::ProjectedBasis)
    dgt = b.dgt
    Binds = complement(Ainds, length(dgt))
    for i = 1:length(v)
        change!(dgt, b.I[i], base=b.B)
        perm_element!(target, dgt, v[i], Ainds, Binds, b.B)
    end
    target
end
