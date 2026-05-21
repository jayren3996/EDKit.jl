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
    SpinlessFermionBasis(dtype::DataType=Int64; L, N=nothing, nf=nothing,
                        f=nothing, alloc=1000, threaded=true, small_N=true)

Construct a spinless-fermion occupation basis on `L` sites.

Arguments:
- `dtype` : Integer type for stored representatives.
- `L`     : Number of sites.
- `N`     : (Optional) particle-number sector. Either a single `Integer`
            (one sector) or an `AbstractVector{<:Integer}` (the union of
            several fixed-N sectors, with sorted representatives).
- `nf`    : (Optional) filling fraction. Equivalent to `N = nf * L`; errors
            if `nf * L` is not an integer or if combined with `N`.
- `f`     : (Optional) predicate `f(dgt) -> Bool` further restricting which
            occupation strings are kept.
- `alloc`, `threaded`, `small_N` : delegated to the same options as
  [`ProjectedBasis`](@ref).
"""
function SpinlessFermionBasis(dtype::DataType=Int64;
    L::Integer, N=nothing, nf::Union{Nothing, Real}=nothing, f=nothing,
    alloc::Integer=1000, threaded::Bool=true, small_N::Bool=true,
)
    base = convert(dtype, 2)
    if !isnothing(nf)
        isnothing(N) || error("Specify either N or nf, not both.")
        scaled = nf * L
        Nint = round(Int, scaled)
        # Tight tolerance: accepts float drift from rationals like 1/3 (≈2e-16
        # of slop on L*1/3) but rejects clearly non-integer fillings like
        # 0.5 + 1e-9. Use `nf = a // b` for exact rational fillings.
        isapprox(scaled, Nint; atol=1e-10) ||
            error("nf=$nf does not yield an integer particle number for L=$L (nf*L=$scaled).")
        N = Nint
    end

    if N isa Integer
        0 <= N <= L ||
            error("N=$N is out of range 0:$L for a $L-site spinless fermion basis.")
    elseif N isa AbstractVector
        isempty(N) && error("Multi-sector N must contain at least one particle number (got empty vector).")
        all(0 .<= N .<= L) ||
            error("Every entry of N must lie in 0:$L (got $N).")
    end

    I = if isnothing(N)
        _fermion_select(dtype, L, nothing, f, base, alloc, threaded, small_N)
    elseif N isa AbstractVector
        sectors = sort!(unique(N))
        parts = [_fermion_select(dtype, L, n, f, base, alloc, threaded, small_N) for n in sectors]
        sort!(vcat(parts...))
    else
        _fermion_select(dtype, L, N, f, base, alloc, threaded, small_N)
    end
    SpinlessFermionBasis(zeros(dtype, L), I, base)
end

function _fermion_select(dtype, L, N, f, base, alloc, threaded, small_N)
    # `selectindex`/`selectindex_threaded` call `f(dgt)` unconditionally, so when
    # no extra predicate is supplied we substitute a trivial true filter.
    f_full = isnothing(f) ? (_ -> true) : f
    if isnothing(N)
        threaded ? selectindex_threaded(f_full, L, base=base, alloc=alloc) :
                   selectindex(f_full, L, 1:base^L, base=base, alloc=alloc)
    elseif small_N
        # selectindex_N(_, L, M) returns indices with Hamming weight L - M
        # (matches ProjectedBasis's spin convention where N counts dgt=0
        # entries). For fermions N counts dgt=1 entries, so pass L - N.
        selectindex_N(f, L, L - N, base=base)
    else
        g = isnothing(f) ? x -> sum(x) == N : x -> (sum(x) == N && f(x))
        threaded ? selectindex_threaded(g, L, base=base, alloc=alloc) :
                   selectindex(g, L, 1:base^L, base=base, alloc=alloc)
    end
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
