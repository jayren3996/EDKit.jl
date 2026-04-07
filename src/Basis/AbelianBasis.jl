#-------------------------------------------------------------------------------------------------------------------------
# Benes Network for O(log L) bit permutation (base=2 only)
#-------------------------------------------------------------------------------------------------------------------------
"""
    BenesNetwork

Compiled bit-permutation network for O(log L) application of a site permutation
to an integer state.  Only used when `base=2`.

Each layer is a pair `(mask, shift)` encoding a set of independent swaps: every
bit pair whose positions differ by `shift` and whose lower bit is marked in
`mask` is swapped simultaneously.
"""
struct BenesNetwork
    layers::Vector{Tuple{UInt64, Int}}   # (mask, shift) pairs
end

"""
    compile_benes(perm::Vector{Int}, L::Int)

Compile a site permutation into a [`BenesNetwork`](@ref).

`perm[i] = j` means site `i` maps to site `j` (1-based).  In integer
representation site `i` occupies bit `L-i` (big-endian, matching EDKit's
`index()`), so the bit at position `L-i` moves to position `L-j`.

The implementation decomposes the permutation into transpositions and groups
compatible swaps (same bit-distance) into shared butterfly layers.
"""
function compile_benes(perm::Vector{Int}, L::Int)
    # Build bit-level permutation: bit_perm[bit_src+1] = bit_dst (0-indexed bits)
    bit_perm = Vector{Int}(undef, L)
    for i in 1:L
        bit_perm[L - i + 1] = L - perm[i]  # 0-indexed: bit (L-i) -> bit (L-perm[i])
    end

    # Decompose into transpositions using selection-sort on the bit permutation
    swaps = Tuple{Int,Int}[]  # ordered list of (lo, hi) transpositions
    current = collect(0:L-1)  # current[pos+1] = which original bit is at position pos
    where_is = collect(0:L-1)  # where_is[bit+1] = current position of original bit

    # bit_perm[src+1] = dst: the bit originally at position src should end at position dst
    for dst in 0:L-1
        # Which original bit should end up at position dst?
        # We need inverse: inv_perm[dst+1] = src such that bit_perm[src+1] = dst
        src_bit = -1
        for s in 0:L-1
            if bit_perm[s+1] == dst
                src_bit = s
                break
            end
        end
        cur_pos = where_is[src_bit+1]
        if cur_pos != dst
            push!(swaps, (min(dst, cur_pos), max(dst, cur_pos)))
            other_bit = current[dst+1]
            current[dst+1], current[cur_pos+1] = current[cur_pos+1], current[dst+1]
            where_is[src_bit+1] = dst
            where_is[other_bit+1] = cur_pos
        end
    end

    # Each swap is a transposition that must be applied in order.
    # Convert each into its own delta-swap layer.
    layers = Tuple{UInt64, Int}[]
    for (lo, hi) in swaps
        shift = hi - lo
        mask = UInt64(1) << lo
        push!(layers, (mask, shift))
    end

    BenesNetwork(layers)
end

"""
    apply_benes(bn::BenesNetwork, state::UInt64, L::Int)

Apply a compiled Benes network to an integer state, returning the permuted state.

Uses the delta-swap technique: for each layer `(mask, shift)`, all bit pairs at
distance `shift` whose lower position is marked in `mask` are swapped.
"""
@inline function apply_benes(bn::BenesNetwork, state::UInt64, L::Int)
    for (mask, shift) in bn.layers
        # Delta swap: swap bits at positions marked in mask with those shifted up by shift
        delta = ((state >> shift) ⊻ state) & mask
        state = state ⊻ delta ⊻ (delta << shift)
    end
    state
end

"""
    apply_perm_int(state::UInt64, perm::Vector{Int}, L::Int)

Apply a site permutation to an integer state using direct bit extraction.
Fallback when Benes network is not available or for verification.
"""
@inline function apply_perm_int(state::UInt64, perm::Vector{Int}, L::Int)
    result = UInt64(0)
    @inbounds for i in 1:L
        bit = (state >> (L - i)) & UInt64(1)
        result |= bit << (L - perm[i])
    end
    result
end


#-------------------------------------------------------------------------------------------------------------------------
# Abelian Operator
#-------------------------------------------------------------------------------------------------------------------------
"""
    AbelianOperator{Tp <: Number}

Internal representation of a product of commuting cyclic symmetry generators.

An `AbelianOperator` stores enough state to iterate through all elements of a
finite Abelian group action on a digit string while tracking the accumulated
phase associated with the current group element.

This is the symmetry backend used by [`AbelianBasis`](@ref).
"""
struct AbelianOperator{Tp <: Number}
    s::Vector{Int}              # state counters (mutable odometer)
    g::Vector{Int}              # orders of each cyclic factor
    c::Vector{Vector{Tp}}       # phase tables
    perm::Vector{Vector{Int}}   # permutation arrays (site i -> site perm[i])
    inv::Vector{BitVector}      # spin-inversion flags per generator
    benes::Vector{Union{Nothing, BenesNetwork}}  # compiled Benes networks
    inv_masks::Vector{UInt64}   # precomputed XOR masks for spin inversion (base=2)
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    AbelianOperator(order::Int, k::Integer, perm::Vector{Int}; inv=falses(length(perm)))

Construct one elementary cyclic generator from a permutation array.

Arguments:
- `order`: order of the cyclic group.
- `k`: momentum or character label selecting the phase representation.
- `perm`: permutation array where `perm[i] = j` means site `i` maps to site `j`.
- `inv`: BitVector indicating which sites get spin-flipped after permutation.

Returns:
- An [`AbelianOperator`](@ref) initialized at the identity group element.
"""
function AbelianOperator(order::Int, k::Integer, perm::Vector{Int}; inv=falses(length(perm)))
    L = length(perm)
    c = if iszero(k)
        ones(order)
    elseif 2k == order
        [iseven(j) ? 1 : -1 for j in 0:order-1]
    else
        phase = 1im * 2π * k / order
        [exp(phase * j) for j in 0:order-1]
    end
    inv_bv = BitVector(inv)
    # Try to compile a Benes network if no inversion is needed
    bn = try
        compile_benes(perm, L)
    catch
        nothing
    end
    # Precompute inversion mask for base=2: XOR mask with 1s at inverted sites
    inv_mask = UInt64(0)
    for i in 1:L
        if inv_bv[i]
            inv_mask |= UInt64(1) << (L - i)
        end
    end
    benes_vec = Union{Nothing, BenesNetwork}[bn]
    AbelianOperator([1], [order], [c], [perm], [inv_bv], benes_vec, [inv_mask])
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    +(g1::AbelianOperator, g2::AbelianOperator)

Form the direct-product combination of two commuting Abelian generators.
"""
function +(g1::AbelianOperator, g2::AbelianOperator)
    AbelianOperator(
        [g1.s; g2.s], [g1.g; g2.g], [g1.c; g2.c],
        [g1.perm; g2.perm], [g1.inv; g2.inv],
        [g1.benes; g2.benes], [g1.inv_masks; g2.inv_masks]
    )
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    order(g::AbelianOperator)

Return the total number of group elements represented by `g`.
"""
function order(g::AbelianOperator)
    prod(g.g)
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    phase(g::AbelianOperator)

Return the current character phase associated with the internal group state of
`g`.

theta(s) = prod_j exp(i*k*s_j)
"""
function phase(g::AbelianOperator)
    prod(g.c[i][g.s[i]] for i in eachindex(g.c))
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    init!(g::AbelianOperator)

Reset the internal group-element counters of `g` to the identity element.
"""
function init!(g::AbelianOperator)
    for i in eachindex(g.s)
        g.s[i] = 1
    end
    g
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    apply_perm!(dgt::Vector, perm::Vector{Int})

Apply a site permutation to a digit buffer in-place.
`perm[i] = j` means the value at site `i` moves to site `j`.
"""
function apply_perm!(dgt::Vector, perm::Vector{Int})
    tmp = copy(dgt)
    @inbounds for i in eachindex(perm)
        dgt[perm[i]] = tmp[i]
    end
    dgt
end

"""
    apply_inv!(dgt::Vector, inv::BitVector, base::Integer)

Apply spin inversion (complement) to sites marked in `inv`.
"""
function apply_inv!(dgt::Vector, inv::BitVector, base::Integer)
    bm1 = base - 1
    @inbounds for i in eachindex(inv)
        if inv[i]
            dgt[i] = bm1 - dgt[i]
        end
    end
    dgt
end

#-------------------------------------------------------------------------------------------------------------------------
"""
    (ag::AbelianOperator)(dgt, base)

Apply the next group action to `dgt` and advance the internal group state.

This mutates both `dgt` and the internal counters stored in `ag`.
"""
function (ag::AbelianOperator)(dgt::Vector, base::Integer)
    for i in eachindex(ag.s)
        apply_perm!(dgt, ag.perm[i])
        any(ag.inv[i]) && apply_inv!(dgt, ag.inv[i], base)
        (ag.s[i] = ag.s[i] + 1) > ag.g[i] ? (ag.s[i] = 1) : break
    end
    dgt
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    check_min(dgt, g::AbelianOperator; base=2)

Check whether `dgt` is the canonical representative of its full Abelian orbit.

Returns:
- `(true, n)` when `dgt` is canonical and has stabilizer size `n`,
- `(false, 0)` otherwise.

This helper is used during basis construction to decide whether a product state
should be stored as a representative.
"""
function check_min(dgt, g::AbelianOperator; base=2)
    init!(g)
    I0 = index(dgt; base)
    N = 1
    for _ in 2:order(g)
        g(dgt, base)
        In = index(dgt; base)
        if In < I0
            return false, 0
        elseif In == I0
            abs(phase(g)-1) < 1e-14 || return false, 0
            N += 1
        end
    end
    true, N
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    shift_canonical!(dgt, g::AbelianOperator; base=2)

Find the canonical representative of the Abelian orbit of `dgt` and record the
group element that maps the current state to that representative.

|d_tilde> = prod_j g_j^s_j |d>

Returns:
- `Im`: the representative index,
- `g`: the same operator with internal counters set to the canonicalizing group
  element.
"""
function shift_canonical!(dgt, g::AbelianOperator; base=2)
    init!(g)
    Im = index(dgt; base)
    ms = g.s[:]
    for _ in 1:order(g)
        g(dgt, base)
        In = index(dgt; base)
        if In < Im
            Im = In
            ms .= g.s
        end
    end
    g.s .= ms
    Im, g
end

#-------------------------------------------------------------------------------------------------------------------------
# Integer-state orbit search for base=2 with Benes networks
#-------------------------------------------------------------------------------------------------------------------------
"""
    _apply_generator_int(state::UInt64, g::AbelianOperator, gen_idx::Int, L::Int)

Apply generator `gen_idx` of an AbelianOperator to an integer state (0-based).
Uses Benes network if available, otherwise falls back to direct bit permutation.
Also applies spin inversion mask if present.
"""
@inline function _apply_generator_int(state::UInt64, g::AbelianOperator, gen_idx::Int, L::Int)
    bn = g.benes[gen_idx]
    if bn !== nothing
        state = apply_benes(bn, state, L)
    else
        state = apply_perm_int(state, g.perm[gen_idx], L)
    end
    state = state ⊻ g.inv_masks[gen_idx]
    state
end

"""
    check_min_int(state::UInt64, g::AbelianOperator, L::Int)

Integer-state version of `check_min` for base=2.
Operates on 0-indexed UInt64 states, avoiding digit buffer allocation.

Returns `(is_canonical::Bool, stabilizer_size::Int)`.
"""
function check_min_int(state::UInt64, g::AbelianOperator, L::Int)
    init!(g)
    s0 = state
    cur = s0
    N = 1
    for _ in 2:order(g)
        # Apply next group element (odometer increment)
        for i in eachindex(g.s)
            cur = _apply_generator_int(cur, g, i, L)
            (g.s[i] = g.s[i] + 1) > g.g[i] ? (g.s[i] = 1) : break
        end
        if cur < s0
            return false, 0
        elseif cur == s0
            abs(phase(g) - 1) < 1e-14 || return false, 0
            N += 1
        end
    end
    true, N
end

"""
    shift_canonical_int(state::UInt64, g::AbelianOperator, L::Int)

Integer-state version of `shift_canonical!` for base=2.
Operates on 0-indexed UInt64 states.

Returns `(min_state_1based::Int, g)` with `g.s` set to the canonicalizing element.
"""
function shift_canonical_int(state::UInt64, g::AbelianOperator, L::Int)
    init!(g)
    best = state
    ms = g.s[:]
    cur = state
    for _ in 1:order(g)
        for i in eachindex(g.s)
            cur = _apply_generator_int(cur, g, i, L)
            (g.s[i] = g.s[i] + 1) > g.g[i] ? (g.s[i] = 1) : break
        end
        if cur < best
            best = cur
            ms .= g.s
        end
    end
    g.s .= ms
    Int(best) + 1, g  # Return 1-based index
end

#-------------------------------------------------------------------------------------------------------------------------
# Gosper's Hack: enumerate integers with exactly N set bits
#-------------------------------------------------------------------------------------------------------------------------
"""
    _gosper_next(state::UInt64)

Return the next integer with the same number of set bits (Gosper's hack).
Returns `UInt64(0)` with overflow when no more exist (caller must detect).
"""
@inline function _gosper_next(state::UInt64)
    c = state & (-state)       # lowest set bit
    r = state + c              # carry into next group
    # (((r ^ state) >> 2) / c) | r  but using bit ops
    diff = r ⊻ state
    diff = (diff >> 2) ÷ c
    r | diff
end

"""
    _gosper_enumerate(L::Int, N::Int)

Return a sorted vector of all 1-based indices whose 0-based representation
has exactly `N` set bits among `L` bits.
"""
function _gosper_enumerate(L::Int, N::Int)
    N == 0 && return [1]  # only state 0 -> 1-based index 1
    N == L && return [Int(UInt64(1) << L)]  # all bits set -> 2^L (1-based: 2^L - 1 + 1)

    result = Int[]
    sizehint!(result, binomial(L, N))

    state = (UInt64(1) << N) - UInt64(1)  # smallest with N bits set
    max_state = UInt64(1) << L
    while state < max_state
        push!(result, Int(state) + 1)  # 1-based
        state = _gosper_next(state)
    end
    result
end

#-------------------------------------------------------------------------------------------------------------------------
# Abelian Basis
#-------------------------------------------------------------------------------------------------------------------------
"""
    AbelianBasis

Basis built from a collection of commuting discrete symmetries represented as
Abelian actions on digit strings.

This is the general symmetry-reduction backend used by the high-level
[`basis`](@ref) constructor when one or more symmetry quantum numbers such as
translation momentum, reflection parity, or spin-flip parity are requested.
"""
struct AbelianBasis{Ti <: Integer, Tg <: Number} <: AbstractPermuteBasis
    dgt::Vector{Ti}         # Digits
    I::Vector{Ti}           # Representing states
    R::Vector{Float64}      # Normalization
    G::AbelianOperator{Tg}  # Generator
    B::Ti                   # Base
end
order(b::AbelianBasis) = order(b.G)
#-------------------------------------------------------------------------------------------------------------------------
"""
    _has_all_benes(G::AbelianOperator)

Check if all generators have compiled Benes networks (required for integer path).
"""
_has_all_benes(G::AbelianOperator) = all(bn !== nothing for bn in G.benes)

function AbelianBasis(
    dtype::DataType=Int64;
    L::Integer, G::AbelianOperator, base::Integer=2, f=x->true,
    N::Union{Nothing, Integer}=nothing,
    alloc=1000, threaded::Bool=true
)
    Ng = order(G)

    C = zeros(Ng)
    for i in eachindex(C)
        iszero(mod(Ng, i)) && (C[i] = sqrt(Ng * i))
    end

    # Dispatch between construction paths
    use_gosper = (base == 2 && N !== nothing && 0 <= N <= L)
    use_int_path = (base == 2 && _has_all_benes(G))

    I, R = if use_gosper && use_int_path
        # Path 1: Gosper + integer orbit search
        candidates = _gosper_enumerate(L, N)
        _abelian_select_int(G, L, candidates, C; threaded)
    elseif use_gosper
        # Path 2: Gosper + digit path
        candidates = _gosper_enumerate(L, N)
        _abelian_select_gosper(G, L, candidates, C; base, threaded)
    elseif threaded
        # Path 3: Full scan + digit path (threaded)
        nt = Threads.nthreads()
        ni = dividerange(base^L, nt)
        nI = Vector{Vector{dtype}}(undef, nt)
        nR = Vector{Vector{Float64}}(undef, nt)
        Threads.@threads for ti in 1:nt
            nI[ti], nR[ti] = _abelian_select(f, deepcopy(G), L, ni[ti], C; base, alloc)
        end
        vcat(nI...), vcat(nR...)
    else
        _abelian_select(f, G, L, 1:base^L, C; base, alloc)
    end
    AbelianBasis(zeros(dtype, L), I, R, G, base)
end

"""
    _abelian_select_int(G, L, candidates, C; threaded=false)

Integer-path basis selection using Benes networks.
`candidates` is a sorted vector of 1-based indices.
"""
function _abelian_select_int(G::AbelianOperator, L::Int, candidates::Vector{Int}, C; threaded::Bool=false)
    if threaded && length(candidates) > 1000
        nt = Threads.nthreads()
        chunk_size = cld(length(candidates), nt)
        nI = Vector{Vector{Int}}(undef, nt)
        nR = Vector{Vector{Float64}}(undef, nt)
        Threads.@threads for ti in 1:nt
            lo = (ti - 1) * chunk_size + 1
            hi = min(ti * chunk_size, length(candidates))
            lo > hi && (nI[ti] = Int[]; nR[ti] = Float64[]; continue)
            g_local = deepcopy(G)
            Is = Int[]
            Rs = Float64[]
            for idx in lo:hi
                state = UInt64(candidates[idx] - 1)
                Q, n = check_min_int(state, g_local, L)
                Q || continue
                push!(Is, candidates[idx])
                push!(Rs, C[n])
            end
            nI[ti] = Is
            nR[ti] = Rs
        end
        vcat(nI...), vcat(nR...)
    else
        Is = Int[]
        Rs = Float64[]
        g_local = deepcopy(G)
        for idx in candidates
            state = UInt64(idx - 1)
            Q, n = check_min_int(state, g_local, L)
            Q || continue
            push!(Is, idx)
            push!(Rs, C[n])
        end
        Is, Rs
    end
end

"""
    _abelian_select_gosper(G, L, candidates, C; base=2, threaded=false)

Gosper-based candidate list with digit-buffer orbit search.
"""
function _abelian_select_gosper(G::AbelianOperator, L::Int, candidates::Vector{Int}, C; base::Integer=2, threaded::Bool=false)
    Is = Int[]
    Rs = Float64[]
    dgt = zeros(Int, L)
    g_local = deepcopy(G)
    for idx in candidates
        change!(dgt, idx; base=Int(base))
        Q, n = check_min(dgt, g_local; base=Int(base))
        Q || continue
        push!(Is, idx)
        push!(Rs, C[n])
    end
    Is, Rs
end

"""
    _abelian_select(f, G, L, rg, C; base=2, alloc=1000)

Low-level worker that scans a range of product-state indices and selects
canonical Abelian-orbit representatives.

This helper is used by both the threaded and non-threaded [`AbelianBasis`](@ref)
constructor paths.
"""
function _abelian_select(
    f,
    G,
    L::Integer,
    rg::UnitRange{T},
    C;
    base::Integer=2,
    alloc::Integer=1000
) where T <: Integer
    dgt = zeros(T, L)
    Is = T[]
    Rs = Float64[]
    sizehint!(Is, alloc)
    sizehint!(Rs, alloc)
    for i in rg
        change!(dgt, i; base)
        f(dgt) || continue
        Q, n = check_min(dgt, G; base)
        Q || continue
        push!(Is, i)
        push!(Rs, C[n])
    end
    Is, Rs
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    index(B::AbelianBasis)

Interpret the current digit buffer as coordinates in the Abelian-reduced basis.

The returned coefficient includes both the stored normalization and the phase
associated with the group element that maps the current digits to the canonical
representative.
"""
function index(B::AbelianBasis)
    if B.B == 2 && _has_all_benes(B.G)
        # Integer path for base=2
        state = UInt64(0)
        L = length(B.dgt)
        @inbounds for i in 1:L
            state = (state << 1) | UInt64(B.dgt[i])
        end
        Im, g = shift_canonical_int(state, B.G, L)
        ind = binary_search(B.I, Im)
        if iszero(ind)
            return zero(eltype(B)), one(eltype(B.I))
        else
            return phase(g) * B.R[ind], ind
        end
    else
        Im, g = shift_canonical!(B.dgt, B.G; base=B.B)
        ind = binary_search(B.I, Im)
        if iszero(ind)
            return zero(eltype(B)), one(eltype(B.I))
        else
            return phase(g) * B.R[ind], ind
        end
    end
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    schmidt(v, Ainds, b::AbelianBasis; B1=nothing, B2=nothing)

Construct the Schmidt matrix of a state represented in an [`AbelianBasis`](@ref).

Unlike onsite bases, each basis coefficient must be expanded over the entire
Abelian orbit with the correct symmetry phase before it contributes to the
bipartite decomposition.
"""
function schmidt(v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::AbelianBasis; B1=nothing, B2=nothing)
    dgt, R, g = b.dgt, b.R, b.G
    S = schmidtmatrix(promote_type(eltype(v), eltype(b)), b, Ainds, B1, B2)
    for i in eachindex(v)
        init!(g)
        change!(b, i)
        val = v[i] / R[i]
        addto!(S, val)
        for _ in 2:order(g)
            g(dgt, b.B)
            addto!(S, phase(g) * val)
        end
    end
    S.M
end
#-------------------------------------------------------------------------------------------------------------------------
export basis
"""
    basis(dtype::DataType=Int64; L, f=nothing, base=2, N=nothing, k=nothing, a=1, p=nothing, z=nothing, symmetries=nothing, threaded=base^L>3000)

High-level basis constructor for the common symmetry combinations supported by
EDKit.

Keywords:
- `L`: system size.
- `f`: optional predicate used to impose additional local constraints.
- `base`: local Hilbert-space dimension.
- `N`: fixed U(1)-like charge sector.
- `k`: translation quantum number.
- `a`: unit-cell length used together with `k`.
- `p`: reflection-parity eigenvalue +/-1.
- `z`: spin-flip eigenvalue +/-1.
- `symmetries`: vector of `(order, k, perm)` or `(order, k, perm, inv)` tuples
  for arbitrary generators (e.g., 2D/3D lattice symmetries).

Return value:
- `TensorBasis` if no symmetry or constraint is requested.
- `ProjectedBasis` if only `f` and/or `N` is requested.
- `AbelianBasis` when one or more discrete symmetries (`k`, `p`, `z`) are used.

Notes:
- `k` and `p` are only simultaneously valid in the symmetry-compatible momentum
  sectors handled by the underlying basis implementation.
- `N` and `z` are only compatible in the half-filling sector for spin-1/2
  systems, mirroring the restrictions of the dedicated basis types.
- This is the most convenient user-facing entry point when you want to combine
  several commuting symmetries without manually choosing a concrete basis type.
"""
function basis(
    dtype::DataType=Int64;
    L::Integer,
    f=nothing,
    base::Integer=2,
    N::Union{Nothing, Integer}=nothing,
    k::Union{Nothing, Integer}=nothing, a::Integer=1,
    p::Union{Nothing, Integer}=nothing,
    z::Union{Nothing, Integer}=nothing,
    symmetries=nothing,
    threaded::Bool=base^L>3000
)
    gs = AbelianOperator[]

    # Translation symmetry: cyclic shift by `a` sites
    if !isnothing(k)
        perm_t = Vector{Int}(undef, L)
        for i in 1:L
            # circshift!(dgt, a) moves dgt[i] to dgt[mod1(i-a, L)]
            # So site i's content goes to site mod1(i-a, L)
            perm_t[i] = mod1(i - a, L)
        end
        push!(gs, AbelianOperator(L÷a, k, perm_t))
    end

    # Parity (reflection) symmetry: reverse!
    if !isnothing(p)
        perm_p = Vector{Int}(undef, L)
        for i in 1:L
            perm_p[i] = L - i + 1
        end
        push!(gs, AbelianOperator(2, isone(-p) ? 1 : 0, perm_p))
    end

    # Spin-flip symmetry
    if !isnothing(z)
        perm_z = collect(1:L)  # identity permutation
        inv_z = trues(L)       # invert all sites
        push!(gs, AbelianOperator(2, isone(-z) ? 1 : 0, perm_z; inv=inv_z))
    end

    # Custom symmetries
    if !isnothing(symmetries)
        for sym in symmetries
            if length(sym) == 3
                ord, kk, prm = sym
                push!(gs, AbelianOperator(ord, kk, prm))
            elseif length(sym) == 4
                ord, kk, prm, inv_flag = sym
                push!(gs, AbelianOperator(ord, kk, prm; inv=inv_flag))
            end
        end
    end

    if isempty(gs)
        isnothing(f) && isnothing(N) && return TensorBasis(;L, base)
        return ProjectedBasis(dtype; L, f, N, base, threaded)
    else
        g = if isnothing(N)
            isnothing(f) ? x -> true : f
        else
            num = L*(base-1)-N
            isnothing(f) ? x -> (sum(x) == num) : x -> (sum(x) == num && f(x))
        end
        return AbelianBasis(dtype; L, G=sum(gs), base, f=g, N, threaded)
    end
end
