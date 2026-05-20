export SpinlessFermionBasis
"""
    SpinlessFermionBasis{T<:Integer}

Occupation-number basis for `L` spinless fermions, optionally restricted to a
fixed particle-number sector.

Internally identical to a `base=2` [`ProjectedBasis`](@ref): each digit `dgt[i]`
is `0` for an empty site and `1` for an occupied site. The distinct type
exists so that fermion-aware operator constructors (see
[`fermion`](@ref) and [`fermion_operator`](@ref)) can dispatch on it.
"""
struct SpinlessFermionBasis{T <: Integer} <: AbstractOnsiteBasis
    dgt::Vector{T}
    I::Vector{T}
    B::T
end

"""
    SpinlessFermionBasis(dtype::DataType=Int64; L, N=nothing,
                        f=nothing, alloc=1000, threaded=true, small_N=true)

Construct a spinless-fermion occupation basis on `L` sites.

Arguments:
- `dtype` : Integer type for stored representatives.
- `L`     : Number of sites.
- `N`     : (Optional) fixed particle-number sector.
- `f`     : (Optional) predicate `f(dgt) -> Bool` further restricting which
            occupation strings are kept.
- `alloc`, `threaded`, `small_N` : delegated to the same options as
  [`ProjectedBasis`](@ref).
"""
function SpinlessFermionBasis(dtype::DataType=Int64;
    L::Integer, N::Union{Nothing, Integer}=nothing, f=nothing,
    alloc::Integer=1000, threaded::Bool=true, small_N::Bool=true,
)
    base = convert(dtype, 2)
    # `selectindex`/`selectindex_threaded` call `f(dgt)` unconditionally, so when
    # no extra predicate is supplied we substitute a trivial true filter.
    f_full = isnothing(f) ? (_ -> true) : f
    I = if isnothing(N)
        threaded ? selectindex_threaded(f_full, L, base=base, alloc=alloc) :
                   selectindex(f_full, L, 1:base^L, base=base, alloc=alloc)
    elseif small_N
        # selectindex_N(_, L, M) returns indices with Hamming weight L - M
        # (that matches ProjectedBasis's spin convention where N counts
        # dgt=0 entries). For fermions we want N occupations, i.e. Hamming
        # weight N, so call it with L - N.
        selectindex_N(f, L, L - N, base=base)
    else
        g = isnothing(f) ? x -> sum(x) == N : x -> (sum(x) == N && f(x))
        threaded ? selectindex_threaded(g, L, base=base, alloc=alloc) :
                   selectindex(g, L, 1:base^L, base=base, alloc=alloc)
    end
    SpinlessFermionBasis(zeros(dtype, L), I, base)
end

copy(b::SpinlessFermionBasis) = SpinlessFermionBasis(deepcopy(b.dgt), b.I, b.B)

function index(b::SpinlessFermionBasis; check::Bool=false)
    index(b, b.dgt; check)
end

function index(b::SpinlessFermionBasis, dgt::AbstractVector; check::Bool=false)
    i = index(dgt, base=b.B)
    ind = binary_search(b.I, i)
    ind > 0 && return 1, ind
    check ? error("No such basis state.") : return zero(eltype(b)), one(ind)
end
