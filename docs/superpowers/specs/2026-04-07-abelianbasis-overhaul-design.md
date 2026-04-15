# AbelianBasis Overhaul: Array-Based Permutations, Benes Networks, and 2D/3D Support

**Date:** 2026-04-07
**Status:** Design
**Motivation:** QuSpin's `spin_basis_general` demonstrates that storing symmetries as integer permutation arrays (instead of opaque closures) enables (a) Benes-network O(log L) bit permutation for base=2, (b) trivial user-defined 2D/3D lattice symmetries via array comprehensions, and (c) Gosper's hack for particle-conservation enumeration. EDKit's `AbelianBasis` already has the correct group-theoretic structure (direct product of cyclic groups) but uses closures, missing these optimizations.

## Scope

Modify `AbelianBasis` and `AbelianOperator` in-place. Break the existing closure-based API. No new types. The high-level `basis(...)` constructor gains a new `symmetries` keyword for arbitrary user-defined symmetry generators.

## 1. AbelianOperator: Array-Based Permutations

### Current State

```julia
struct AbelianOperator{Tp <: Number}
    s::Vector{Int}              # state counters (mutable odometer)
    g::Vector{Int}              # orders of each cyclic factor
    c::Vector{Vector{Tp}}       # phase tables
    f::Vector                   # permutation closures (opaque)
end
```

The `f` field stores arbitrary functions. This prevents any algebraic analysis of the permutations.

### New Design

```julia
struct AbelianOperator{Tp <: Number}
    s::Vector{Int}              # state counters (mutable odometer) [unchanged]
    g::Vector{Int}              # orders of each cyclic factor [unchanged]
    c::Vector{Vector{Tp}}       # phase tables [unchanged]
    perm::Vector{Vector{Int}}   # permutation arrays (1-indexed site maps)
    inv::Vector{BitVector}      # spin-inversion masks per generator
    benes::Vector{Any}          # compiled Benes networks (base=2 only, nothing otherwise)
end
```

**Permutation arrays:** `perm[k]` is a length-L vector where `perm[k][i] = j` means site `i` maps to site `j` under generator `k`. 1-indexed (Julia convention).

**Inversion masks:** `inv[k]` is a length-L BitVector. If `inv[k][i]` is true, the digit at site `i` is complemented (`base-1-d`) after the permutation is applied.

**Benes networks:** `benes[k]` stores the pre-compiled Benes network for generator `k` when `base=2`, or `nothing` when `base>2` or when the Benes compilation is not applicable. See Section 3.

### Constructor

```julia
# New API:
AbelianOperator(order::Int, k::Integer, perm::Vector{Int}; inv=falses(length(perm)))
```

- `order`: cyclic group order. Validated against the actual period of the permutation (must divide; error if not).
- `k`: character label (momentum quantum number).
- `perm`: 1-indexed permutation array of length L.
- `inv`: optional spin-inversion mask.

**Auto-period computation:** The period of a permutation+inversion is computed by repeatedly applying it until identity is recovered. The constructor validates that `order` equals or is a multiple of this period.

**Phase table construction:** Unchanged from current code (special-cased for k=0, 2k=g, general complex).

### Group Action Application

**For base > 2 (digit-buffer path):**
```julia
function apply_generator!(dgt::Vector, perm::Vector{Int}, inv::BitVector, base::Integer)
    tmp = copy(dgt)  # or use a pre-allocated buffer
    @inbounds for i in eachindex(perm)
        d = tmp[i]
        dgt[perm[i]] = inv[i] ? (base - 1 - d) : d
    end
end
```

**For base=2 (integer path):**
```julia
function apply_generator(state::Integer, benes::BenesNetwork, inv_mask::Integer)
    apply_benes(benes, state ^ inv_mask)  # XOR for spin inversion, then Benes permute
end
```

### Odometer Iteration

The odometer logic in `(ag::AbelianOperator)(dgt)` changes from calling `ag.f[i](dgt)` to calling `apply_generator!(dgt, ag.perm[i], ag.inv[i], base)`. For the integer path, a separate method operates on integer states directly.

## 2. Integer-State Orbit Search

### Current Bottleneck

`shift_canonical!` and `check_min` currently work on digit buffers:
1. Apply group action to digit buffer: O(L)
2. Convert digit buffer to integer index: O(L)
3. Compare integers: O(1)

Each step is O(L), and there are O(|G|) iterations, giving O(|G| * L) per lookup.

### New: Direct Integer-State Operations

For base=2, operate directly on integer states without digit buffers:

```julia
function shift_canonical_int(state::Integer, ag::AbelianOperator)
    init!(ag)
    best = state
    best_s = copy(ag.s)
    current = state
    for _ in 1:order(ag)
        current = apply_next_int(current, ag)  # advance odometer + apply on integer
        if current < best
            best = current
            best_s .= ag.s
        end
    end
    ag.s .= best_s
    best, ag
end
```

With Benes networks, each `apply_next_int` call costs O(log L) instead of O(L), giving total cost O(|G| * log L) vs O(|G| * L).

For base > 2, a similar integer path using `_int_translate`/`_int_reverse`/`_int_spinflip` (already in AbstractBasis.jl) avoids the digit-buffer overhead. However, for general permutations (not just translate/reverse/flip), the digit-buffer path is needed. We can still avoid repeated `index(dgt)` calls by computing the integer directly from the permuted digits in a single pass.

### Dispatch Strategy

- `base == 2` and all generators have Benes networks: integer path with Benes
- `base > 2` or Benes not available: digit-buffer path (improved with pre-allocated temp buffer)

## 3. Benes Network for Base=2

### What is a Benes Network?

A Benes network is a circuit of `2*ceil(log2(L))-1` layers that can implement any permutation of L elements using only local swap operations. Each layer is defined by:
- A **mask**: which positions participate in swaps
- A **shift**: how far apart swap partners are

Applying one layer costs O(1) bitwise operations (mask, shift, XOR). The full network costs O(log L) layers.

### Data Structure

```julia
struct BenesNetwork
    masks::Vector{UInt64}    # bit masks for each layer
    shifts::Vector{Int}      # shift amounts for each layer
    nlayers::Int
end
```

### Compilation

At `AbelianOperator` construction time (when base=2):

1. Pad L to the next power of 2: `n = nextpow(2, L)`
2. Build the Benes network topology for n elements
3. Route the specific permutation through the network using a recursive algorithm:
   - Split the permutation into "upper" and "lower" sub-problems
   - Assign each element to pass through the upper or lower sub-network
   - Recurse on sub-networks
4. Extract the masks and shifts for each layer

This is a well-known algorithm (Knuth TAOCP Vol 4, Waksman 1968).

### Application

```julia
function apply_benes(bn::BenesNetwork, state::UInt64)
    s = state
    @inbounds for i in 1:bn.nlayers
        # Butterfly swap: exchange bits at positions differing by shifts[i]
        # where masks[i] indicates which positions to swap
        t = ((s >> bn.shifts[i]) ^ s) & bn.masks[i]
        s = s ^ t ^ (t << bn.shifts[i])
    end
    s
end
```

Each layer is 4 bitwise operations. For L=16: ~7 layers = ~28 operations, vs L=16 digit manipulations in the current code.

### Spin Inversion Combination

For base=2, spin inversion at site `i` means flipping bit `i`. This is a single XOR with a precomputed bitmask:

```julia
inv_mask = sum(1 << (L-i) for i where inv[i] is true)
state_after = apply_benes(bn, state ^ inv_mask)
```

The XOR is folded into the Benes application: `apply_benes(bn, state ^ inv_mask)`.

## 4. Gosper's Hack for Particle Conservation

### When It Applies

When `base == 2` and a particle number constraint `N` is active (i.e., only states with exactly `N` set bits are valid).

### Algorithm

```julia
function gosper_next(state::T) where T <: Integer
    c = state & (-state)       # lowest set bit
    r = state + c              # carry propagation
    (((r ^ state) >> 2) div c) | r  # fill from right
end
```

Starting from `(1 << N) - 1` (smallest integer with N set bits), repeatedly calling `gosper_next` generates all integers with exactly N set bits in ascending order.

### Integration

In `_abelian_select`, when `base == 2` and the predicate `f` is a pure particle-number filter:

```julia
function _abelian_select_pcon(G, L, N, C; alloc=1000)
    state = (one(Int64) << N) - one(Int64)  # first state with N bits
    max_state = state << (L - N)             # last state with N bits
    Is = Int64[]
    Rs = Float64[]
    sizehint!(Is, alloc)
    sizehint!(Rs, alloc)
    while state <= max_state
        # state is 0-indexed; EDKit uses 1-indexed
        Q, n = check_min_int(state, G)
        if Q
            push!(Is, state + 1)
            push!(Rs, C[n])
        end
        state = gosper_next(state)
    end
    Is, Rs
end
```

**Speedup:** For L=20, N=10: scans 184,756 states instead of 1,048,576 (5.7x). For L=24, N=12: scans 2,704,156 instead of 16,777,216 (6.2x).

### Detection

The `basis(...)` constructor already wraps particle-number constraints as `x -> sum(x) == num`. We detect this case at construction time:
- If `base == 2` and `N` is specified and no other predicate `f` is given: use Gosper path
- Otherwise: fall back to full enumeration with predicate filtering

## 5. Updated `basis(...)` Constructor

### New `symmetries` Keyword

```julia
function basis(
    dtype::DataType=Int64;
    L::Integer,
    f=nothing,
    base::Integer=2,
    N::Union{Nothing, Integer}=nothing,
    k::Union{Nothing, Integer}=nothing, a::Integer=1,
    p::Union{Nothing, Integer}=nothing,
    z::Union{Nothing, Integer}=nothing,
    symmetries=nothing,  # NEW: Vector of (perm, quantum_number) tuples
    threaded::Bool=base^L>3000
)
```

**Behavior:**
- If `symmetries` is provided, it takes precedence over `k`, `p`, `z` (which must be `nothing`)
- Each element of `symmetries` is a tuple `(perm::Vector{Int}, q::Integer)` or `(perm::Vector{Int}, q::Integer, inv::BitVector)`
- The `AbelianOperator` is constructed from these generators
- `N` and `f` still work as before (particle conservation + custom predicate)

**Existing 1D API remains unchanged:**
```julia
basis(; L=12, N=6, k=0, p=1, z=1)  # still works, generates permutation arrays internally
```

**New 2D example:**
```julia
Lx, Ly = 4, 3; L = Lx * Ly
sites = [(x, y) for y in 0:Ly-1 for x in 0:Lx-1]

T_x = [mod(x+1, Lx) + Lx*y + 1 for (x,y) in sites]
T_y = [x + Lx*mod(y+1, Ly) + 1 for (x,y) in sites]
P_x = [Lx-1-x + Lx*y + 1 for (x,y) in sites]

B = basis(; L, N=L/2, base=2,
    symmetries=[(T_x, 0), (T_y, 0), (P_x, 0)])
```

### Internal Changes to 1D Path

The existing `k`, `p`, `z` keywords now generate permutation arrays instead of closures:

```julia
# Translation: circshift by a sites
push!(gs, AbelianOperator(L/a, k, [mod1(i+a, L) for i in 1:L]))

# Parity (reflection): reverse site order
push!(gs, AbelianOperator(2, isone(-p) ? 1 : 0, [L+1-i for i in 1:L]))

# Spin-flip: identity permutation with full inversion mask
push!(gs, AbelianOperator(2, isone(-z) ? 1 : 0, collect(1:L); inv=trues(L)))
```

## 6. Thread Safety

The current `AbelianOperator` stores mutable state (`s` counters) that makes it thread-unsafe for concurrent `index()` calls. This is already handled by `deepcopy(G)` in the threaded basis constructor (line 202).

**No change needed for the overhaul.** The digit-buffer path requires per-thread copies of both `dgt` and `G`, which is already the pattern. The integer path avoids digit buffers entirely but still needs per-thread `G` copies for the odometer state.

For operator application (`colmn!` in Operator.jl), each basis object has its own `dgt` and `G`, so single-threaded operator application is safe. The multi-threaded `mul()` already creates per-thread basis copies.

## 7. Schmidt Decomposition Compatibility

The `schmidt` method for `AbelianBasis` (line 273-287 in AbelianBasis.jl) iterates over orbit elements by calling `g(dgt)`. This must be updated to use the new permutation-array group action. The logic is the same; only the mechanism of applying the group element changes.

## 8. Files Modified

| File | Changes |
|------|---------|
| `src/Basis/AbelianBasis.jl` | Complete rewrite of `AbelianOperator` struct, constructors, `check_min`, `shift_canonical!`, `index`, `schmidt`. Add Benes network, integer orbit search, Gosper's hack. |
| `src/Basis/AbstractBasis.jl` | No changes needed (integer helpers already exist). |
| `test/basis_tests.jl` | Update tests for new API. Add 2D lattice tests. |
| `test/core_tests.jl` | Update any `AbelianOperator` construction calls. |

## 9. Testing Strategy

1. **Correctness regression:** For all existing 1D test cases (translation, parity, flip, combinations), verify that the new array-based `AbelianBasis` produces identical basis states and identical eigenvalues as the old closure-based version.

2. **2D lattice test:** Construct a 4x3 Heisenberg model on a square lattice with (T_x, T_y, P_x) symmetries. Verify eigenvalues match the full-space (`TensorBasis`) diagonalization projected onto the corresponding symmetry sector.

3. **Benes network test:** For several random permutations of lengths 4, 8, 16, 32: verify that the Benes network produces the same result as direct digit-buffer permutation on all 2^L states.

4. **Gosper's hack test:** Verify that the particle-conserving enumeration produces the same basis as the full-scan path for L=12, N=6.

5. **Performance benchmark:** Compare basis construction and operator application times for L=16, base=2, translation+parity+flip, before and after the overhaul.

## 10. Non-Goals

- **Non-Abelian symmetries:** Out of scope. The direct-product-of-cyclic-groups structure is preserved.
- **Fermion/boson bases:** Out of scope. This overhaul is spin-only.
- **Lattice helper objects:** Users provide raw permutation arrays (QuSpin style). No `SquareLattice` abstraction.
- **Time evolution / dynamics:** Separate feature, not part of this overhaul.
