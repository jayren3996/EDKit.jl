export ProjectedBasis
"""
    ProjectedBasis{Ti}

Basis obtained by selecting a subset of product states from the full tensor
product Hilbert space.

This behaves like [`TensorBasis`](@ref), but only retains basis vectors whose
digit representation satisfies a user-defined predicate or a fixed charge `N`.

This is the standard choice for constrained Hilbert spaces where states remain
product-state-labelled, but not every product state is allowed.
"""
struct ProjectedBasis{T <: Integer} <: AbstractOnsiteBasis
    dgt::Vector{T}
    I::Vector{T}
    B::T
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    copy(b::ProjectedBasis)

Return a copy of `b` with an independent mutable digit buffer.

The stored representative list `I` is shared, while `dgt` is duplicated so that
iterative algorithms can mutate the copy safely.
"""
copy(b::ProjectedBasis) = ProjectedBasis(deepcopy(b.dgt), b.I, b.B)
#-------------------------------------------------------------------------------------------------------------------------
"""
    index(b::ProjectedBasis; check=false)

Interpret the current digit buffer `b.dgt` as coordinates in the projected
basis.

Arguments:
- `b`: projected basis whose working digit buffer has already been set.
- `check`: if `true`, throw an error when the current digits are not contained
  in the basis; otherwise return a zero coefficient.

Returns:
- `(1, i)` when the state is present as the `i`th projected basis vector.
- `(0, 1)` when the state is not present and `check == false`.

The zero-coefficient branch is important during operator assembly, where an
off-sector state should silently contribute nothing instead of aborting the
matrix construction.
"""
function index(b::ProjectedBasis; check::Bool=false)
    index(b, b.dgt; check)
end
function index(b::ProjectedBasis, dgt::AbstractVector; check::Bool=false)
    i = index(dgt, base=b.B)
    ind = binary_search(b.I, i)
    ind > 0 && return 1, ind
    check ? error("No such symmetry.") : return zero(eltype(b)), one(ind)
end


#-------------------------------------------------------------------------------------------------------------------------
# Select basis vectors
#-------------------------------------------------------------------------------------------------------------------------
"""
    binary_search(list::AbstractVector{<:Integer}, i::Integer)

Return the position of `i` inside a sorted integer list, or `0` if it is absent.

This helper is used pervasively by reduced bases to locate canonical
representatives inside their stored index arrays without a linear scan.
"""
function binary_search(list::AbstractVector{<:Integer}, i::Integer)
    isempty(list) && return 0

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

function _binary_fixed_weight_indices(::Type{T}, L::Integer, n::Integer) where T <: Integer
    (0 <= n <= L) || return T[]
    if iszero(n)
        return T[one(T)]
    elseif n == L
        return T[(one(T) << L)]
    end
    state = (one(T) << n) - one(T)
    limit = one(T) << L
    out = T[]
    sizehint!(out, binomial(L, n))
    while state < limit
        push!(out, state + one(T))
        c = state & -state
        r = state + c
        state = (((r ⊻ state) >> 2) ÷ c) | r
    end
    out
end

function _complement_digits!(dgt::AbstractVector, fdgt::AbstractVector, base)
    bm1 = base - 1
    @inbounds for i in eachindex(dgt, fdgt)
        dgt[i] = bm1 - fdgt[i]
    end
    dgt
end

#-------------------------------------------------------------------------------------------------------------------------
# Construction
#-------------------------------------------------------------------------------------------------------------------------
"""
    dividerange(maxnum::Integer, nthreads::Integer)

Split the integer range `1:maxnum` into approximately equal chunks for
multi-threaded iteration.

Inputs:
-------
- `maxnum`  : Maximum number of range.
- `nthreads`: Number of threads.

Returns:
- `list`: List of ranges.

This helper centralizes EDKit's convention for static work partitioning across
threads.
"""
function dividerange(maxnum::T, nthreads::Integer) where T <: Integer
    list = Vector{UnitRange{T}}(undef, nthreads)
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
#-------------------------------------------------------------------------------------------------------------------------
"""
    selectindex(f, L, rg; base=2, alloc=1000)

Enumerate all product-state indices in `rg` whose digit strings satisfy `f`.

Inputs:
-------
- `f`    : Function for digits that tells whether a digits is valid.
- `L`    : Length of the system.
- `rg`   : Range of interation.
- `base` : (Optional) Base, `base=2` by default.
- `alloc`: (Optional) Pre-allocation memory for the list of indices, `alloc=1000` by default.

Returns:
- A vector of 1-based product-state indices whose decoded digits satisfy `f`.

This helper performs the actual scan used by [`ProjectedBasis`](@ref) and other
reduced-basis constructors.
"""
function selectindex(f, L::Integer, rg::UnitRange{T}; base::Integer=2, alloc::Integer=1000) where T <: Integer
    dgt = zeros(T, L)
    I = T[]
    sizehint!(I, alloc)
    for i in rg
        change!(dgt, i, base=base)
        f(dgt) && push!(I, i)
    end
    I
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    selectindex_threaded(f, L, rg; base=2, alloc=1000)

Threaded version of [`selectindex`](@ref) over the full Hilbert space.

Inputs:
-------
- `f`    : Function for digits that tells whether a digits is valid as well as its normalization.
- `L`    : Length of the system.
- `base` : (Optional) Base, `base=2` by default.
- `alloc`: (Optional) Pre-allocation memory for the list of indices, `alloc=1000` by default.

Returns:
- A concatenated vector of accepted product-state indices.

The result ordering matches the natural increasing index order because each
thread scans a disjoint monotone block and the blocks are concatenated in order.
"""
function selectindex_threaded(f, L::Integer; base::T=2, alloc::Integer=1000) where T <: Integer
    maxnum = base^L
    nt = maxnum < Threads.nthreads() ? Int(maxnum) : Threads.nthreads()
    ni = dividerange(maxnum, nt)
    nI = Vector{Vector{T}}(undef, nt)
    Threads.@threads for ti in 1:nt
        nI[ti] = selectindex(f, L, ni[ti], base=base, alloc=alloc)
    end
    I = vcat(nI...)
    I
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    selectindex_N(f, L, N; base=2, alloc=1000, sorted=true)

Enumerate a fixed-charge subset of product states without scanning the entire
Hilbert space.

Arguments:
- `f`: optional predicate applied after the fixed-`N` filter.
- `L`: system size.
- `N`: target charge in EDKit's convention.
- `base`: local base.
- `alloc`: initial capacity hint.
- `sorted`: whether to sort the resulting representative indices.

Returns:
- A vector of basis indices compatible with the requested fixed-charge sector.

This helper is used when direct combinatorial enumeration is cheaper than a full
Hilbert-space scan.
"""
function selectindex_N(f, L::Integer, N::Integer; base::T=2, alloc::Integer=1000, sorted::Bool=true) where T <: Integer
    if base == 2 && sorted
        candidates = _binary_fixed_weight_indices(T, L, L - N)
        isnothing(f) && return candidates
        I = T[]
        sizehint!(I, length(candidates))
        dgt = zeros(T, L)
        for ind in candidates
            change!(dgt, ind; base=base)
            f(dgt) && push!(I, ind)
        end
        return I
    end
    I = T[]
    sizehint!(I, alloc)
    dgt = Vector{T}(undef, L)
    for fdgt in multiexponents(L, N)
        all(b < base for b in fdgt) || continue
        _complement_digits!(dgt, fdgt, base)
        isnothing(f) || f(dgt) || continue
        ind = index(dgt, base=base)
        push!(I, ind)
    end
    sorted ? sort!(I) : I
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    ProjectedBasis(dtype::DataType=Int64; L, f=nothing, N=nothing, base=2, alloc=1000, threaded=true, small_N=true)

Construct a projected basis on `L` sites.

Inputs:
-------
- `dtype`   : Data type for indices.
- `L`       : Length of the system.
- `f`       : Predicate on digit vectors. Only states with `f(dgt) == true` are kept.
- `N`       : Fixed charge sector. In EDKit's digit convention this corresponds to
              `sum(dgt) == L*(base-1) - N`.
- `base`    : Base, default = 2.
- `alloc`   : Size of the prealloc memory for the basis content, used only in multithreading, default = 1000.
- `threaded`: Whether use the multithreading, default = true.
- `small_N` : Whether to enumerate fixed-`N` states directly instead of scanning the full space.

Returns:
- `b` : `ProjectedBasis`.

This constructor is useful for constrained Hilbert spaces such as PXP/Rydberg
constraints or fixed-magnetization sectors.

Notes:
- If both `f` and `N` are provided, the state must satisfy both filters.
- The charge convention follows EDKit's digit encoding, not necessarily the most
  obvious particle-count convention, so docstrings elsewhere often state it
  explicitly as `sum(dgt) == L*(base-1) - N`.
"""
function ProjectedBasis(dtype::DataType=Int64;
    L::Integer, f=nothing, N::Union{Nothing, Integer}=nothing,
    base::Integer=2, alloc::Integer=1000, 
    threaded::Bool=true, small_N::Bool=true
)
    base = convert(dtype, base)
    I = if isnothing(N)
        threaded ? selectindex_threaded(f, L, base=base, alloc=alloc) : selectindex(f, L, 1:base^L, base=base, alloc=alloc)
    elseif small_N
        selectindex_N(f, L, N, base=base)
    else
        num = L * (base-1) - N
        g = isnothing(f) ? x -> sum(x) == num : x -> (sum(x) == num && f(x))
        threaded ? selectindex_threaded(g, L, base=base, alloc=alloc) : selectindex(g, L, 1:base^L, base=base, alloc=alloc)
    end
    ProjectedBasis(zeros(dtype, L), I, base)
end
