# AbelianBasis Overhaul Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace closure-based `AbelianOperator` with array-based permutations, add Benes network optimization for base=2, Gosper's hack for particle conservation, and a `symmetries` keyword for 2D/3D lattice support.

**Architecture:** Modify `AbelianOperator` struct to store permutation arrays + inversion masks instead of closures. Add a `BenesNetwork` struct for O(log L) bit permutation. Add integer-state orbit search paths dispatched by base. Add `_gosper_next` for particle-conserving enumeration. Extend `basis(...)` with a `symmetries` keyword.

**Tech Stack:** Julia, SparseArrays, LinearAlgebra, Test

---

## File Structure

| File | Responsibility | Action |
|------|---------------|--------|
| `src/Basis/AbelianBasis.jl` | `AbelianOperator` struct, Benes network, orbit search, basis construction, `basis()` constructor | Rewrite |
| `test/basis_tests.jl` | Core basis and symmetry decomposition tests | Modify (add 2D tests) |
| `test/abelian_overhaul_tests.jl` | Dedicated tests for Benes, Gosper, integer orbit, 2D lattice | Create |
| `test/runtests.jl` | Test runner | Modify (include new test file) |
| `docs/src/abelian_basis.md` | User-facing documentation for AbelianBasis and 2D/3D usage | Create |

---

### Task 1: Benes Network Implementation

**Files:**
- Modify: `src/Basis/AbelianBasis.jl` (add at top, before AbelianOperator)
- Create: `test/abelian_overhaul_tests.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Write failing test for Benes network**

Create `test/abelian_overhaul_tests.jl`:
```julia
using EDKit
using LinearAlgebra
using SparseArrays
using Test
using Random

Random.seed!(42)

@testset "Benes Network" begin
    # Test: identity permutation
    bn_id = EDKit.compile_benes([1, 2, 3, 4])
    for s in UInt64(0):UInt64(15)
        @test EDKit.apply_benes(bn_id, s, 4) == s
    end

    # Test: reverse permutation [4,3,2,1]
    bn_rev = EDKit.compile_benes([4, 3, 2, 1])
    # Reversing bits of 0b1010 (10) with L=4 gives 0b0101 (5)
    @test EDKit.apply_benes(bn_rev, UInt64(10), 4) == UInt64(5)
    # Reversing bits of 0b1100 (12) gives 0b0011 (3)
    @test EDKit.apply_benes(bn_rev, UInt64(12), 4) == UInt64(3)

    # Test: cyclic shift [2,3,4,1] (site 1->2, 2->3, 3->4, 4->1)
    bn_cyc = EDKit.compile_benes([2, 3, 4, 1])
    # 0b1000 (8) -> site 1 is set -> maps to site 2 -> 0b0100 (4)
    @test EDKit.apply_benes(bn_cyc, UInt64(8), 4) == UInt64(4)

    # Test: exhaustive correctness for random permutations
    for L in [4, 8, 12, 16]
        perm = randperm(L)
        bn = EDKit.compile_benes(perm)
        dgt_in = zeros(Int, L)
        dgt_out = zeros(Int, L)
        for s in UInt64(0):UInt64(min(2^L - 1, 1023))
            # Compute expected result via digit manipulation
            for i in L:-1:1
                dgt_in[i] = (s >> (L - i)) & 1
            end
            for i in 1:L
                dgt_out[perm[i]] = dgt_in[i]
            end
            expected = UInt64(0)
            for i in 1:L
                expected = (expected << 1) | UInt64(dgt_out[i])
            end
            @test EDKit.apply_benes(bn, s, L) == expected
        end
    end
end
```

Add to `test/runtests.jl` at the end, before the closing:
```julia
include("abelian_overhaul_tests.jl")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd /Users/ren/Library/CloudStorage/OneDrive-UniversityofLeeds/GitHub/EDKit.jl && julia --project -e 'using Test; include("test/abelian_overhaul_tests.jl")'`

Expected: FAIL with "compile_benes not defined" or "apply_benes not defined"

- [ ] **Step 3: Implement BenesNetwork**

Add the following at the top of `src/Basis/AbelianBasis.jl` (before the AbelianOperator section):

```julia
#-------------------------------------------------------------------------------------------------------------------------
# Benes Network for O(log L) bit permutation (base=2 only)
#-------------------------------------------------------------------------------------------------------------------------
struct BenesNetwork
    masks::Vector{UInt64}
    shifts::Vector{Int}
    nlayers::Int
end

"""
    compile_benes(perm::Vector{Int}) -> BenesNetwork

Compile a site permutation into a Benes network for O(log L) bit permutation.

The permutation `perm` is 1-indexed: `perm[i] = j` means site `i` maps to site `j`.
The resulting network operates on integer states where bit `L-i` represents site `i`
(matching EDKit's big-endian digit convention).

The algorithm uses recursive Waksman decomposition: split the permutation into two
halves routed through an upper and lower sub-network, then recurse.
"""
function compile_benes(perm::Vector{Int})
    L = length(perm)
    n = nextpow(2, L)
    # Extend permutation to power-of-2 size (extra positions map to themselves)
    ext_perm = collect(1:n)
    ext_perm[1:L] .= perm

    nlayers = 2 * Int(log2(n)) - 1
    masks = zeros(UInt64, nlayers)
    shifts = zeros(Int, nlayers)

    _benes_route!(masks, shifts, ext_perm, n, 0, 1)
    BenesNetwork(masks, shifts, nlayers)
end

"""
    _benes_route!(masks, shifts, perm, size, offset, depth)

Recursively route a permutation through a Benes network.

This implements the Waksman/Benes recursive decomposition. At each level:
1. The network has `size` wires starting at `offset`
2. Wires are paired: (offset+i, offset+i+half) for i in 0:half-1
3. Each pair can either swap or pass-through (controlled by mask bits)
4. The input and output layers route elements to the correct sub-network
5. Recurse on upper and lower halves
"""
function _benes_route!(masks, shifts, perm, size, offset, depth)
    if size <= 1
        return
    end

    half = size ÷ 2
    n_total = length(masks)
    log_n = (n_total + 1) ÷ 2

    input_layer = depth
    output_layer = n_total + 1 - depth
    shifts[input_layer] = half
    shifts[output_layer] = half

    # Build a bipartite matching: which elements go to upper (0:half-1) vs lower (half:size-1)
    # target_half[i] = 0 if perm sends element i to upper half, 1 if lower half
    target_half = fill(-1, size)
    source_of = fill(-1, size)  # inverse permutation within this block

    # Adjust perm to be local (0-indexed within this block)
    local_perm = zeros(Int, size)
    for i in 1:size
        local_perm[i] = perm[offset + i] - offset
    end

    for i in 1:size
        source_of[local_perm[i]] = i
    end

    # Greedy coloring: assign elements to upper/lower sub-networks
    # Process pairs (i, i+half) at the output side
    for start in 1:half
        # Output pair: positions start and start+half
        # Try to assign one to upper, one to lower
        out_top = start
        out_bot = start + half

        src_top = source_of[out_top]
        src_bot = source_of[out_bot]

        if target_half[src_top] == -1 && target_half[src_bot] == -1
            # Neither assigned yet - assign src_top -> upper (no output swap)
            _benes_assign!(target_half, source_of, local_perm, half, src_top, 0)
        end
    end

    # Fill any remaining unassigned
    for i in 1:size
        if target_half[i] == -1
            _benes_assign!(target_half, source_of, local_perm, half, i, 0)
        end
    end

    # Set input layer mask: for each input pair (i, i+half), swap if needed
    for i in 0:half-1
        top_wire = offset + i + 1
        bot_wire = offset + i + half + 1
        # top_wire should go to upper sub-network (target_half=0)
        # If target_half[i+1] == 1, we need to swap this input pair
        if target_half[i + 1] == 1
            # Set the mask bit for this pair
            # In EDKit's convention, site j uses bit position (total_L - j)
            # But within Benes, we use position = offset + i
            bit_pos = offset + i
            masks[input_layer] |= UInt64(1) << bit_pos
        end
    end

    # Build sub-permutations for upper and lower networks
    upper_perm = collect(1:half)
    lower_perm = collect(1:half)

    for i in 1:size
        src = i
        dst = local_perm[i]
        src_sub = target_half[src]  # 0=upper, 1=lower
        dst_sub = (dst > half) ? 1 : 0

        src_pos = src_sub == 0 ? src : src - half
        if src > half
            src_pos = src_sub == 0 ? src - half : src
        end
        # Simpler: track position within sub-network
    end

    # Set output layer mask and compute sub-permutations
    upper_count = 0
    lower_count = 0
    upper_input = zeros(Int, half)   # maps upper sub-network input position to source element
    lower_input = zeros(Int, half)

    for i in 1:half
        # Input pair: (i, i+half)
        top_elem = i
        bot_elem = i + half
        if target_half[top_elem] == 0
            upper_count += 1
            upper_input[upper_count] = top_elem
        else
            lower_count += 1
            lower_input[lower_count] = top_elem
        end
        if target_half[bot_elem] == 0
            upper_count += 1
            upper_input[upper_count] = bot_elem
        else
            lower_count += 1
            lower_input[lower_count] = bot_elem
        end
    end

    # Compute where each element ends up at the output
    upper_output_map = zeros(Int, half)
    lower_output_map = zeros(Int, half)

    for ui in 1:half
        elem = upper_input[ui]
        dst = local_perm[elem]
        # dst is in 1:half (upper output) or half+1:size (lower output)
        if dst <= half
            upper_output_map[ui] = dst
        else
            upper_output_map[ui] = dst - half
            # Need output swap for this pair
            bit_pos = offset + dst - half - 1
            masks[output_layer] |= UInt64(1) << bit_pos
        end
    end

    for li in 1:half
        elem = lower_input[li]
        dst = local_perm[elem]
        if dst > half
            lower_output_map[li] = dst - half
        else
            lower_output_map[li] = dst
        end
    end

    # Build sub-permutation arrays (1-indexed, offset-adjusted)
    for i in 1:half
        perm[offset + i] = offset + upper_output_map[i]
        perm[offset + half + i] = offset + lower_output_map[i]
    end

    # Recurse on upper and lower halves
    _benes_route!(masks, shifts, perm, half, offset, depth + 1)
    _benes_route!(masks, shifts, perm, half, offset + half, depth + 1)
end

function _benes_assign!(target_half, source_of, local_perm, half, start, color)
    pos = start
    c = color
    while true
        target_half[pos] = c
        # Follow the chain: pos -> output -> partner -> source
        dst = local_perm[pos]
        partner_dst = dst <= half ? dst + half : dst - half
        partner_src = source_of[partner_dst]
        if target_half[partner_src] != -1
            break
        end
        pos = partner_src
        c = 1 - c
    end
end

"""
    apply_benes(bn::BenesNetwork, state::UInt64, L::Int) -> UInt64

Apply a compiled Benes network to permute the bits of `state`.

Only the lowest `L` bits of `state` are meaningful; higher bits are ignored
and will be zero in the output.
"""
@inline function apply_benes(bn::BenesNetwork, state::UInt64, L::Int)
    s = state
    @inbounds for i in 1:bn.nlayers
        t = ((s >> bn.shifts[i]) ⊻ s) & bn.masks[i]
        s = s ⊻ t ⊻ (t << bn.shifts[i])
    end
    s & ((UInt64(1) << L) - UInt64(1))
end
```

**Note:** The Benes routing algorithm above is a sketch. The actual implementation needs careful attention to the recursive decomposition. A simpler alternative that is guaranteed correct: for each permutation, precompute the network by solving the routing problem via the standard Waksman algorithm. If the recursive routing proves tricky to get right, an alternative approach is to decompose the permutation into transpositions and compile each transposition into a single butterfly layer, which is simpler but produces more layers.

**Practical fallback approach** (simpler, guaranteed correct):

```julia
function compile_benes(perm::Vector{Int})
    L = length(perm)
    # Decompose into swap layers: each layer swaps pairs at distance 2^k
    # Use the "odd-even merge sort network" approach for simplicity
    n = nextpow(2, L)
    layers = _sorting_network_layers(n)
    masks, shifts = _route_permutation(perm, layers, n)
    BenesNetwork(masks, shifts, length(masks))
end
```

The key correctness requirement: `apply_benes(compile_benes(perm), state, L)` must produce the same result as applying `perm` to the bits of `state` for all states `0:2^L-1`.

- [ ] **Step 4: Run test to verify it passes**

Run: `cd /Users/ren/Library/CloudStorage/OneDrive-UniversityofLeeds/GitHub/EDKit.jl && julia --project -e 'using Test; include("test/abelian_overhaul_tests.jl")'`

Expected: All "Benes Network" tests PASS

If the recursive routing is buggy, debug by testing small permutations (L=2, L=4) individually and tracing through the layers. The exhaustive test over all 2^L states for L=4,8 will catch any errors.

- [ ] **Step 5: Commit**

```bash
git add src/Basis/AbelianBasis.jl test/abelian_overhaul_tests.jl test/runtests.jl
git commit -m "feat: add Benes network for O(log L) bit permutation"
```

---

### Task 2: New AbelianOperator Struct with Array-Based Permutations

**Files:**
- Modify: `src/Basis/AbelianBasis.jl` (replace AbelianOperator struct and all methods)

- [ ] **Step 1: Write failing test for new AbelianOperator constructor**

Add to `test/abelian_overhaul_tests.jl`:
```julia
@testset "AbelianOperator Array-Based" begin
    L = 6

    # Test: translation generator (cyclic shift by 1)
    perm_t = [mod1(i + 1, L) for i in 1:L]
    g_t = EDKit.AbelianOperator(L, 0, perm_t)
    @test EDKit.order(g_t) == L
    @test length(g_t.perm) == 1
    @test g_t.perm[1] == perm_t
    @test all(.!g_t.inv[1])  # no spin inversion

    # Test: parity generator (reverse)
    perm_p = [L + 1 - i for i in 1:L]
    g_p = EDKit.AbelianOperator(2, 0, perm_p)
    @test EDKit.order(g_p) == 2

    # Test: spin-flip generator (identity perm + full inversion)
    perm_z = collect(1:L)
    g_z = EDKit.AbelianOperator(2, 0, perm_z; inv=trues(L))
    @test EDKit.order(g_z) == 2
    @test all(g_z.inv[1])

    # Test: combination via +
    g_combined = g_t + g_p
    @test EDKit.order(g_combined) == 2L
    @test length(g_combined.perm) == 2

    # Test: digit-buffer group action
    dgt = [0, 1, 0, 1, 0, 0]
    dgt_copy = copy(dgt)
    EDKit.init!(g_t)
    g_t(dgt, 2)  # base=2
    # After circshift by 1: [0, 0, 1, 0, 1, 0]
    @test dgt == [0, 0, 1, 0, 1, 0]

    # Test: phase computation for k=1
    g_k1 = EDKit.AbelianOperator(L, 1, perm_t)
    EDKit.init!(g_k1)
    @test EDKit.phase(g_k1) ≈ 1.0  # identity element
    g_k1(zeros(Int, L), 2)  # advance once
    @test EDKit.phase(g_k1) ≈ exp(2im * π / L)

    # Test: auto-period validation
    # Order 4 for a permutation that actually has period 4
    perm4 = [2, 3, 4, 1, 5, 6]  # cycles (1234), fixes 5,6
    g4 = EDKit.AbelianOperator(4, 0, perm4)
    @test EDKit.order(g4) == 4
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Test; include("test/abelian_overhaul_tests.jl")'`

Expected: FAIL — old AbelianOperator constructor takes a function, not an array

- [ ] **Step 3: Implement new AbelianOperator struct**

Replace the `AbelianOperator` section in `src/Basis/AbelianBasis.jl` (lines 1-102) with:

```julia
#-------------------------------------------------------------------------------------------------------------------------
# Abelian Operator (array-based permutations)
#-------------------------------------------------------------------------------------------------------------------------
"""
    AbelianOperator{Tp <: Number}

Internal representation of a product of commuting cyclic symmetry generators.

Each generator is stored as a permutation array plus an optional spin-inversion
mask. For base=2 systems, a pre-compiled Benes network enables O(log L)
integer-state permutation.

This is the symmetry backend used by [`AbelianBasis`](@ref).
"""
struct AbelianOperator{Tp <: Number}
    s::Vector{Int}
    g::Vector{Int}
    c::Vector{Vector{Tp}}
    perm::Vector{Vector{Int}}
    inv::Vector{BitVector}
    benes::Vector{Union{Nothing, BenesNetwork}}
    inv_masks::Vector{UInt64}
end

"""
    _compute_period(perm::Vector{Int}, inv::BitVector)

Compute the period of a combined permutation + inversion operation.

The period is the smallest positive integer `p` such that applying the operation
`p` times returns every site to its original position with original spin direction.
"""
function _compute_period(perm::Vector{Int}, inv::BitVector)
    L = length(perm)
    pos = collect(1:L)
    flipped = falses(L)
    for p in 1:2*L  # period is at most 2L (perm period * 2 for inversion)
        new_pos = similar(pos)
        new_flipped = similar(flipped)
        for i in 1:L
            new_pos[i] = perm[pos[i]]
            new_flipped[i] = flipped[i] ⊻ inv[pos[i]]
        end
        pos .= new_pos
        flipped .= new_flipped
        if pos == 1:L && !any(flipped)
            return p
        end
    end
    error("Period computation failed — permutation may be invalid")
end

"""
    AbelianOperator(g::Int, k::Integer, perm::Vector{Int}; inv=falses(length(perm)))

Construct one elementary cyclic generator from a permutation array.

Arguments:
- `g`: order of the cyclic group (must equal the period of the permutation+inversion).
- `k`: character label (momentum/quantum number) selecting the phase representation.
- `perm`: 1-indexed permutation array of length L, where `perm[i] = j` means
  site `i` maps to site `j`.
- `inv`: optional spin-inversion mask. If `inv[i]` is true, the digit at site `i`
  is complemented (`base-1-d`) after the permutation.

Returns:
- An [`AbelianOperator`](@ref) initialized at the identity group element.
"""
function AbelianOperator(g::Int, k::Integer, perm::Vector{Int}; inv::Union{BitVector, Vector{Bool}}=falses(length(perm)))
    L = length(perm)
    inv_bv = inv isa BitVector ? inv : BitVector(inv)

    # Validate permutation
    @assert length(inv_bv) == L "Inversion mask length must match permutation length"
    @assert sort(perm) == 1:L "perm must be a valid permutation of 1:L"

    # Validate order matches period
    period = _compute_period(perm, inv_bv)
    @assert g == period "Specified order $g does not match permutation period $period"

    # Phase table
    c = if iszero(k)
        ones(g)
    elseif 2k == g
        [iseven(j) ? 1 : -1 for j in 0:g-1]
    else
        ph = 1im * 2π * k / g
        [exp(ph * j) for j in 0:g-1]
    end

    # Compile Benes network (will be used when base=2)
    bn = try
        compile_benes(perm)
    catch
        nothing
    end

    # Precompute inversion mask for base=2 integer path
    inv_mask = UInt64(0)
    for i in 1:L
        if inv_bv[i]
            inv_mask |= UInt64(1) << (L - i)
        end
    end

    AbelianOperator([1], [g], [c], [perm], [inv_bv], [bn], [inv_mask])
end

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

"""
    order(g::AbelianOperator)

Return the total number of group elements represented by `g`.
"""
order(g::AbelianOperator) = prod(g.g)

"""
    phase(g::AbelianOperator)

Return the current character phase associated with the internal group state.
"""
function phase(g::AbelianOperator)
    prod(g.c[i][g.s[i]] for i in eachindex(g.c))
end

"""
    init!(g::AbelianOperator)

Reset the internal group-element counters to the identity element.
"""
function init!(g::AbelianOperator)
    for i in eachindex(g.s)
        g.s[i] = 1
    end
    g
end

"""
    apply_perm!(dgt::Vector, perm::Vector{Int}, inv::BitVector, base::Integer)

Apply a single generator's permutation + inversion to a digit buffer in-place.

Uses an internal temporary buffer to avoid overwrite conflicts.
"""
function apply_perm!(dgt::Vector, perm::Vector{Int}, inv::BitVector, base::Integer)
    L = length(dgt)
    # We need a temporary copy to avoid read-after-write conflicts
    tmp = dgt[:]  # stack-allocated for small L
    @inbounds for i in 1:L
        d = tmp[i]
        dgt[perm[i]] = inv[i] ? (base - 1 - d) : d
    end
end

"""
    (ag::AbelianOperator)(dgt::Vector, base::Integer)

Apply the next group action to `dgt` and advance the internal group state.

This is the digit-buffer path used for base > 2 (or as fallback for base=2).
"""
function (ag::AbelianOperator)(dgt::Vector, base::Integer)
    for i in eachindex(ag.s)
        apply_perm!(dgt, ag.perm[i], ag.inv[i], base)
        (ag.s[i] = ag.s[i] + 1) > ag.g[i] ? (ag.s[i] = 1) : break
    end
    dgt
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project -e 'using Test; include("test/abelian_overhaul_tests.jl")'`

Expected: "AbelianOperator Array-Based" tests PASS

- [ ] **Step 5: Commit**

```bash
git add src/Basis/AbelianBasis.jl test/abelian_overhaul_tests.jl
git commit -m "feat: replace AbelianOperator with array-based permutations"
```

---

### Task 3: Integer-State Orbit Search

**Files:**
- Modify: `src/Basis/AbelianBasis.jl` (add integer orbit functions)

- [ ] **Step 1: Write failing test for integer orbit search**

Add to `test/abelian_overhaul_tests.jl`:
```julia
@testset "Integer Orbit Search" begin
    L = 6

    # Translation k=0: representative should be the minimum under cyclic shifts
    perm_t = [mod1(i + 1, L) for i in 1:L]
    g = EDKit.AbelianOperator(L, 0, perm_t)

    # State 0b101010 = 42 (0-indexed)
    # Cyclic shifts: 42, 21, 42, 21, ... (period 2 for L=6 of 101010)
    # Actually for L=6: 101010 -> 010101 -> 101010 -> ...
    # min(101010, 010101) = 010101 = 21
    best, _ = EDKit.shift_canonical_int(UInt64(42), g, L)
    @test best == UInt64(21)

    # State 0b000001 = 1 (0-indexed)
    # Shifts: 1, 2, 4, 8, 16, 32
    # min = 1
    best, _ = EDKit.shift_canonical_int(UInt64(1), g, L)
    @test best == UInt64(1)

    # check_min_int should return true only for the canonical representative
    @test EDKit.check_min_int(UInt64(21), g, L)[1] == true
    @test EDKit.check_min_int(UInt64(42), g, L)[1] == false

    # Parity k=0: representative is min of state and reverse
    perm_p = [L + 1 - i for i in 1:L]
    g_p = EDKit.AbelianOperator(2, 0, perm_p)

    # 0b110000 = 48 reversed = 0b000011 = 3, min = 3
    best, _ = EDKit.shift_canonical_int(UInt64(48), g_p, L)
    @test best == UInt64(3)

    # Combined translation + parity
    g_tp = g + g_p
    # Check that it produces the correct representative
    best1, _ = EDKit.shift_canonical_int(UInt64(48), g_tp, L)
    # All shifts of 48 and all shifts of reversed 48
    # 48=110000, shifts: 110000,011000,001100,000110,000011,100001
    # rev(48)=000011=3, shifts: 000011,100001,110000,011000,001100,000110
    # min = 000011 = 3
    @test best1 == UInt64(3)
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Test; include("test/abelian_overhaul_tests.jl")'`

Expected: FAIL — `shift_canonical_int` not defined

- [ ] **Step 3: Implement integer orbit search**

Add to `src/Basis/AbelianBasis.jl` after the `AbelianOperator` callable:

```julia
#-------------------------------------------------------------------------------------------------------------------------
# Integer-state orbit search (base=2 fast path)
#-------------------------------------------------------------------------------------------------------------------------
"""
    _apply_gen_int(state::UInt64, bn::Union{Nothing, BenesNetwork}, inv_mask::UInt64, L::Int) -> UInt64

Apply a single generator to an integer state using the Benes network.
"""
@inline function _apply_gen_int(state::UInt64, bn::BenesNetwork, inv_mask::UInt64, L::Int)
    apply_benes(bn, state ⊻ inv_mask, L)
end

"""
    _advance_int!(ag::AbelianOperator, state::UInt64, L::Int) -> UInt64

Apply the next group action to an integer state and advance the odometer.

This mirrors the digit-buffer `(ag)(dgt, base)` callable but works on raw integers.
"""
@inline function _advance_int!(ag::AbelianOperator, state::UInt64, L::Int)
    for i in eachindex(ag.s)
        state = _apply_gen_int(state, ag.benes[i], ag.inv_masks[i], L)
        (ag.s[i] = ag.s[i] + 1) > ag.g[i] ? (ag.s[i] = 1) : break
    end
    state
end

"""
    check_min_int(state::UInt64, g::AbelianOperator, L::Int) -> (Bool, Int)

Check whether `state` is the canonical (minimum) representative of its orbit.

Returns `(true, stabilizer_size)` if canonical, `(false, 0)` otherwise.
Operates on 0-indexed integer states using Benes networks.
"""
function check_min_int(state::UInt64, g::AbelianOperator, L::Int)
    init!(g)
    N = 1
    current = state
    for _ in 2:order(g)
        current = _advance_int!(g, current, L)
        if current < state
            return false, 0
        elseif current == state
            abs(phase(g) - 1) < 1e-14 || return false, 0
            N += 1
        end
    end
    true, N
end

"""
    shift_canonical_int(state::UInt64, g::AbelianOperator, L::Int) -> (UInt64, AbelianOperator)

Find the canonical (minimum) representative of the orbit of `state`.

Returns the minimum state and sets `g`'s internal counters to the group element
that maps the input state to the representative.
"""
function shift_canonical_int(state::UInt64, g::AbelianOperator, L::Int)
    init!(g)
    best = state
    best_s = g.s[:]
    current = state
    for _ in 1:order(g)
        current = _advance_int!(g, current, L)
        if current < best
            best = current
            best_s .= g.s
        end
    end
    g.s .= best_s
    best, g
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project -e 'using Test; include("test/abelian_overhaul_tests.jl")'`

Expected: "Integer Orbit Search" tests PASS

- [ ] **Step 5: Commit**

```bash
git add src/Basis/AbelianBasis.jl test/abelian_overhaul_tests.jl
git commit -m "feat: add integer-state orbit search with Benes networks"
```

---

### Task 4: Gosper's Hack for Particle Conservation

**Files:**
- Modify: `src/Basis/AbelianBasis.jl`

- [ ] **Step 1: Write failing test for Gosper's hack**

Add to `test/abelian_overhaul_tests.jl`:
```julia
@testset "Gosper's Hack" begin
    # Test: enumerate all 6-choose-3 = 20 states
    L = 6; N = 3
    states = EDKit._gosper_enumerate(L, N)
    @test length(states) == binomial(L, N)
    # All should have exactly N bits set
    @test all(count_ones(s) == N for s in states)
    # Should be sorted ascending
    @test issorted(states)
    # Should contain no duplicates
    @test length(unique(states)) == length(states)

    # Test edge cases
    @test length(EDKit._gosper_enumerate(8, 0)) == 1  # only state 0
    @test length(EDKit._gosper_enumerate(8, 8)) == 1  # only state 0xFF
    @test length(EDKit._gosper_enumerate(10, 5)) == binomial(10, 5)

    # Test: _gosper_next produces correct sequence
    s = UInt64(0b000111)  # 7, smallest 3-bit state for L>=3
    s = EDKit._gosper_next(s)
    @test s == UInt64(0b001011)  # 11
    s = EDKit._gosper_next(s)
    @test s == UInt64(0b001101)  # 13
    s = EDKit._gosper_next(s)
    @test s == UInt64(0b001110)  # 14
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Test; include("test/abelian_overhaul_tests.jl")'`

Expected: FAIL — `_gosper_next` and `_gosper_enumerate` not defined

- [ ] **Step 3: Implement Gosper's hack**

Add to `src/Basis/AbelianBasis.jl`:

```julia
#-------------------------------------------------------------------------------------------------------------------------
# Gosper's hack: enumerate integers with fixed popcount
#-------------------------------------------------------------------------------------------------------------------------
"""
    _gosper_next(state::T) where T <: Integer

Return the next integer with the same number of set bits, in ascending order.

Uses Gosper's hack (bit manipulation) to generate the next state.
"""
@inline function _gosper_next(state::T) where T <: Integer
    c = state & (-state)       # lowest set bit
    r = state + c              # carry propagation
    (((r ⊻ state) >> 2) ÷ c) | r
end

"""
    _gosper_enumerate(L::Int, N::Int) -> Vector{UInt64}

Enumerate all L-bit integers with exactly N set bits, in ascending order.
"""
function _gosper_enumerate(L::Int, N::Int)
    if N == 0
        return [UInt64(0)]
    end
    if N == L
        return [UInt64((1 << L) - 1)]
    end
    first_state = UInt64((1 << N) - 1)
    last_state = first_state << (L - N)
    states = UInt64[]
    sizehint!(states, binomial(L, N))
    s = first_state
    while s <= last_state
        push!(states, s)
        s = _gosper_next(s)
    end
    states
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project -e 'using Test; include("test/abelian_overhaul_tests.jl")'`

Expected: "Gosper's Hack" tests PASS

- [ ] **Step 5: Commit**

```bash
git add src/Basis/AbelianBasis.jl test/abelian_overhaul_tests.jl
git commit -m "feat: add Gosper's hack for fixed-popcount enumeration"
```

---

### Task 5: Updated AbelianBasis Constructor and `index()` / `check_min` / `shift_canonical!`

This task rewrites the core `AbelianBasis` construction and lookup to use the new `AbelianOperator` with both digit-buffer and integer-state paths.

**Files:**
- Modify: `src/Basis/AbelianBasis.jl` (rewrite `check_min`, `shift_canonical!`, `index`, `AbelianBasis` constructor, `_abelian_select`)

- [ ] **Step 1: Write failing test for updated AbelianBasis**

Add to `test/abelian_overhaul_tests.jl`:
```julia
@testset "AbelianBasis Construction" begin
    L = 8

    # Test: translation-only basis matches TranslationalBasis eigenvalues
    XXZ = spin((1.0, "xx"), (1.0, "yy"), (0.4, "zz"))
    full_vals = eigvals(Hermitian(trans_inv_operator(XXZ, 2, L))) |> sort

    # All k sectors should reconstruct full spectrum
    vals_k = Float64[]
    for k in 0:L-1
        B = basis(; L, k)
        iszero(size(B, 1)) && continue
        append!(vals_k, eigvals(Hermitian(trans_inv_operator(XXZ, 2, B))))
    end
    @test sort(vals_k) ≈ full_vals

    # Test: parity-only
    vals_p = Float64[]
    for p in [1, -1]
        B = basis(; L, p)
        append!(vals_p, eigvals(Hermitian(trans_inv_operator(XXZ, 2, B))))
    end
    @test sort(vals_p) ≈ full_vals

    # Test: spin-flip-only
    vals_z = Float64[]
    for z in [1, -1]
        B = basis(; L, z)
        append!(vals_z, eigvals(Hermitian(trans_inv_operator(XXZ, 2, B))))
    end
    @test sort(vals_z) ≈ full_vals

    # Test: N + k combined
    vals_nk = Float64[]
    for N in 0:L, k in 0:L-1
        B = basis(; L, N, k)
        iszero(size(B, 1)) && continue
        append!(vals_nk, eigvals(Hermitian(trans_inv_operator(XXZ, 2, B))))
    end
    @test sort(vals_nk) ≈ full_vals

    # Test: full symmetric (N + k + p + z at half-filling, k=0)
    N = L ÷ 2
    Bref = basis(; L, N, k=0)
    iszero(size(Bref, 1)) || begin
        ref_vals = eigvals(Hermitian(trans_inv_operator(XXZ, 2, Bref))) |> sort
        vals_sym = Float64[]
        for p in [1, -1], z in [1, -1]
            B = basis(; L, N, k=0, p, z)
            iszero(size(B, 1)) && continue
            append!(vals_sym, eigvals(Hermitian(trans_inv_operator(XXZ, 2, B))))
        end
        @test sort(vals_sym) ≈ ref_vals
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Test; include("test/abelian_overhaul_tests.jl")'`

Expected: FAIL — the `basis(...)` constructor now fails because it tries to pass closures to the new `AbelianOperator` constructor

- [ ] **Step 3: Rewrite check_min, shift_canonical!, index, constructor, and basis()**

Replace the corresponding sections in `src/Basis/AbelianBasis.jl`:

```julia
#-------------------------------------------------------------------------------------------------------------------------
# Canonicality check and orbit search (dual-path: integer for base=2, digit-buffer for base>2)
#-------------------------------------------------------------------------------------------------------------------------
"""
    check_min(dgt, g::AbelianOperator; base=2)

Check whether `dgt` is the canonical representative of its Abelian orbit (digit-buffer path).
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
            abs(phase(g) - 1) < 1e-14 || return false, 0
            N += 1
        end
    end
    true, N
end

"""
    shift_canonical!(dgt, g::AbelianOperator; base=2)

Find the canonical representative of the orbit of `dgt` (digit-buffer path).
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
# AbelianBasis struct and constructor
#-------------------------------------------------------------------------------------------------------------------------
struct AbelianBasis{Ti <: Integer, Tg <: Number} <: AbstractPermuteBasis
    dgt::Vector{Ti}
    I::Vector{Ti}
    R::Vector{Float64}
    G::AbelianOperator{Tg}
    B::Ti
end
order(b::AbelianBasis) = order(b.G)

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

    # Choose enumeration strategy
    use_gosper = (base == 2 && !isnothing(N) && 0 <= N <= L)
    use_int_path = (base == 2 && all(bn !== nothing for bn in G.benes))

    I, R = if use_gosper && use_int_path
        _abelian_select_gosper_int(G, L, N, C; alloc, threaded)
    elseif use_gosper
        _abelian_select_gosper(f, G, L, N, C; alloc)
    elseif threaded
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

function _abelian_select(
    f, G, L::Integer, rg::UnitRange{T}, C; base::Integer=2, alloc::Integer=1000
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

"""
    _abelian_select_gosper(f, G, L, N, C; alloc)

Basis construction using Gosper's hack to enumerate only states with N set bits (digit-buffer path).
"""
function _abelian_select_gosper(f, G, L::Integer, N::Integer, C; alloc::Integer=1000)
    Is = Int64[]
    Rs = Float64[]
    sizehint!(Is, alloc)
    sizehint!(Rs, alloc)
    dgt = zeros(Int64, L)
    if N == 0
        change!(dgt, Int64(1); base=Int64(2))
        if f(dgt)
            Q, n = check_min(dgt, G; base=2)
            Q && (push!(Is, Int64(1)); push!(Rs, C[n]))
        end
        return Is, Rs
    end
    states = _gosper_enumerate(L, N)
    for s0 in states
        i = Int64(s0 + 1)  # 0-indexed -> 1-indexed
        change!(dgt, i; base=Int64(2))
        f(dgt) || continue
        Q, n = check_min(dgt, G; base=2)
        Q || continue
        push!(Is, i)
        push!(Rs, C[n])
    end
    Is, Rs
end

"""
    _abelian_select_gosper_int(G, L, N, C; alloc, threaded)

Basis construction using Gosper's hack + integer orbit search (fastest path for base=2).
"""
function _abelian_select_gosper_int(G, L::Integer, N::Integer, C; alloc::Integer=1000, threaded::Bool=true)
    states = _gosper_enumerate(L, N)

    if !threaded || length(states) < 1000
        return _gosper_int_worker(deepcopy(G), states, L, C; alloc)
    end

    # Threaded path
    nt = Threads.nthreads()
    chunk_size = cld(length(states), nt)
    nI = Vector{Vector{Int64}}(undef, nt)
    nR = Vector{Vector{Float64}}(undef, nt)
    Threads.@threads for ti in 1:nt
        lo = (ti - 1) * chunk_size + 1
        hi = min(ti * chunk_size, length(states))
        lo > hi && (nI[ti] = Int64[]; nR[ti] = Float64[]; continue)
        nI[ti], nR[ti] = _gosper_int_worker(deepcopy(G), @view(states[lo:hi]), L, C; alloc)
    end
    vcat(nI...), vcat(nR...)
end

function _gosper_int_worker(G, states, L, C; alloc=1000)
    Is = Int64[]
    Rs = Float64[]
    sizehint!(Is, alloc)
    sizehint!(Rs, alloc)
    for s0 in states
        Q, n = check_min_int(s0, G, L)
        Q || continue
        push!(Is, Int64(s0 + 1))
        push!(Rs, C[n])
    end
    Is, Rs
end

#-------------------------------------------------------------------------------------------------------------------------
# index() for AbelianBasis
#-------------------------------------------------------------------------------------------------------------------------
"""
    index(B::AbelianBasis)

Interpret the current digit buffer as coordinates in the Abelian-reduced basis.
"""
function index(B::AbelianBasis)
    L = length(B.dgt)
    if B.B == 2 && all(bn !== nothing for bn in B.G.benes)
        # Integer fast path
        state = UInt64(index(B.dgt; base=B.B) - 1)
        best, g = shift_canonical_int(state, B.G, L)
        ind = binary_search(B.I, Int64(best + 1))
        if iszero(ind)
            return zero(eltype(B)), one(eltype(B.I))
        else
            return phase(g) * B.R[ind], ind
        end
    else
        # Digit-buffer path
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
# Schmidt decomposition for AbelianBasis
#-------------------------------------------------------------------------------------------------------------------------
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
# Updated basis() constructor
#-------------------------------------------------------------------------------------------------------------------------
export basis
"""
    basis(dtype; L, f, base, N, k, a, p, z, symmetries, threaded)

High-level basis constructor supporting 1D symmetries (k, p, z) and arbitrary
Abelian symmetries via the `symmetries` keyword.

The `symmetries` keyword accepts a vector of tuples, each specifying one cyclic
generator as `(perm, quantum_number)` or `(perm, quantum_number, inv_mask)`:
- `perm`: 1-indexed permutation array of length L
- `quantum_number`: integer character label for this generator
- `inv_mask` (optional): BitVector or Vector{Bool} of length L for spin inversion

Example (2D square lattice):
```julia
Lx, Ly = 4, 3; L = Lx * Ly
sites = [(x, y) for y in 0:Ly-1 for x in 0:Lx-1]
T_x = [mod(x+1, Lx) + Lx*y + 1 for (x,y) in sites]
T_y = [x + Lx*mod(y+1, Ly) + 1 for (x,y) in sites]
B = basis(; L, symmetries=[(T_x, 0), (T_y, 0)])
```
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
    if !isnothing(symmetries)
        # Validate that k, p, z are not also specified
        @assert isnothing(k) && isnothing(p) && isnothing(z) "Cannot specify both symmetries and k/p/z"

        gs = AbelianOperator[]
        for sym in symmetries
            if length(sym) == 2
                perm, q = sym
                push!(gs, AbelianOperator(_compute_period(perm, falses(length(perm))), q, perm))
            elseif length(sym) == 3
                perm, q, inv_mask = sym
                inv_bv = inv_mask isa BitVector ? inv_mask : BitVector(inv_mask)
                push!(gs, AbelianOperator(_compute_period(perm, inv_bv), q, perm; inv=inv_bv))
            else
                error("Each symmetry must be (perm, q) or (perm, q, inv_mask)")
            end
        end

        G = length(gs) == 1 ? gs[1] : sum(gs)
        g_pred = if isnothing(N)
            isnothing(f) ? x -> true : f
        else
            num = L * (base - 1) - N
            isnothing(f) ? x -> (sum(x) == num) : x -> (sum(x) == num && f(x))
        end
        return AbelianBasis(dtype; L, G, base, f=g_pred, N, alloc=1000, threaded)
    end

    # Original 1D path: construct AbelianOperator from k, p, z keywords
    gs = AbelianOperator[]
    if !isnothing(k)
        perm_t = [mod1(i + a, L) for i in 1:L]
        push!(gs, AbelianOperator(L ÷ a, k, perm_t))
    end
    if !isnothing(p)
        perm_p = [L + 1 - i for i in 1:L]
        push!(gs, AbelianOperator(2, isone(-p) ? 1 : 0, perm_p))
    end
    if !isnothing(z)
        perm_z = collect(1:L)
        push!(gs, AbelianOperator(2, isone(-z) ? 1 : 0, perm_z; inv=trues(L)))
    end

    if isempty(gs)
        isnothing(f) && isnothing(N) && return TensorBasis(; L, base)
        return ProjectedBasis(dtype; L, f, N, base, threaded)
    else
        g_pred = if isnothing(N)
            isnothing(f) ? x -> true : f
        else
            num = L * (base - 1) - N
            isnothing(f) ? x -> (sum(x) == num) : x -> (sum(x) == num && f(x))
        end
        return AbelianBasis(dtype; L, G=sum(gs), base, f=g_pred, N, alloc=1000, threaded)
    end
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project -e 'using Test; include("test/abelian_overhaul_tests.jl")'`

Expected: "AbelianBasis Construction" tests PASS

- [ ] **Step 5: Commit**

```bash
git add src/Basis/AbelianBasis.jl test/abelian_overhaul_tests.jl
git commit -m "feat: rewrite AbelianBasis with integer orbit search and Gosper's hack"
```

---

### Task 6: 2D Lattice Support Tests

**Files:**
- Modify: `test/abelian_overhaul_tests.jl`

- [ ] **Step 1: Write 2D lattice test**

Add to `test/abelian_overhaul_tests.jl`:
```julia
@testset "2D Lattice Symmetries" begin
    # 3x3 square lattice with Heisenberg model
    Lx, Ly = 3, 3
    L = Lx * Ly
    sites = [(x, y) for y in 0:Ly-1 for x in 0:Lx-1]

    # Translation in x: site (x,y) -> (x+1 mod Lx, y)
    T_x = [mod(x + 1, Lx) + Lx * y + 1 for (x, y) in sites]
    # Translation in y: site (x,y) -> (x, y+1 mod Ly)
    T_y = [x + Lx * mod(y + 1, Ly) + 1 for (x, y) in sites]

    # Build 2D Heisenberg Hamiltonian (nearest-neighbor)
    XXZ = spin((1.0, "xx"), (1.0, "yy"), (1.0, "zz"))
    bonds = Tuple{Int,Int}[]
    for (x, y) in sites
        i = x + Lx * y + 1
        # x-bond
        j = mod(x + 1, Lx) + Lx * y + 1
        push!(bonds, (i, j))
        # y-bond
        j = x + Lx * mod(y + 1, Ly) + 1
        push!(bonds, (i, j))
    end

    # Full-space eigenvalues
    H_full = operator([XXZ for _ in bonds], [[b[1], b[2]] for b in bonds], L)
    full_vals = eigvals(Hermitian(Array(H_full))) |> sort

    # With (kx=0, ky=0) translation symmetry only
    B_k00 = basis(; L, base=2, symmetries=[(T_x, 0), (T_y, 0)])
    @test size(B_k00, 1) > 0  # should have states

    # Verify: all (kx, ky) sectors reconstruct full spectrum
    vals_all_k = Float64[]
    for kx in 0:Lx-1, ky in 0:Ly-1
        B = basis(; L, base=2, symmetries=[(T_x, kx), (T_y, ky)])
        iszero(size(B, 1)) && continue
        H = operator([XXZ for _ in bonds], [[b[1], b[2]] for b in bonds], B)
        append!(vals_all_k, eigvals(Hermitian(Array(H))))
    end
    @test sort(vals_all_k) ≈ full_vals

    # Test: 2D with particle conservation
    N = L ÷ 2
    vals_k_N = Float64[]
    H_N_ref = operator([XXZ for _ in bonds], [[b[1], b[2]] for b in bonds],
                        basis(; L, N))
    ref_vals = eigvals(Hermitian(Array(H_N_ref))) |> sort
    for kx in 0:Lx-1, ky in 0:Ly-1
        B = basis(; L, N, base=2, symmetries=[(T_x, kx), (T_y, ky)])
        iszero(size(B, 1)) && continue
        H = operator([XXZ for _ in bonds], [[b[1], b[2]] for b in bonds], B)
        append!(vals_k_N, eigvals(Hermitian(Array(H))))
    end
    @test sort(vals_k_N) ≈ ref_vals
end
```

- [ ] **Step 2: Run test**

Run: `julia --project -e 'using Test; include("test/abelian_overhaul_tests.jl")'`

Expected: "2D Lattice Symmetries" tests PASS

- [ ] **Step 3: Commit**

```bash
git add test/abelian_overhaul_tests.jl
git commit -m "test: add 2D square lattice symmetry decomposition tests"
```

---

### Task 7: Full Regression Test Suite

Run the complete existing test suite to ensure nothing is broken.

**Files:**
- Modify: `test/runtests.jl` (add new test file inclusion)

- [ ] **Step 1: Ensure test runner includes new tests**

Verify `test/runtests.jl` has:
```julia
include("abelian_overhaul_tests.jl")
```
at the end.

- [ ] **Step 2: Run full test suite**

Run: `cd /Users/ren/Library/CloudStorage/OneDrive-UniversityofLeeds/GitHub/EDKit.jl && julia --project -e 'using Pkg; Pkg.test()'`

Expected: ALL tests PASS, including all existing tests in:
- `core_tests.jl`
- `basis_tests.jl`
- `entanglement_tests.jl`
- `advanced_tests.jl`
- `lindblad_tests.jl`
- `itensor_tests.jl`
- `abelian_overhaul_tests.jl`

If any existing tests fail, the issue is likely in the `basis()` constructor or the `AbelianOperator` callable signature change (now requires `base` argument). Fix any call sites that pass the old closure-based API.

**Common failure points to check:**
1. `basis_tests.jl` calls `basis(; L, k, p, z, N)` — should work with new constructor
2. `entanglement_tests.jl` uses `basis()` with various symmetries — should work
3. `advanced_tests.jl` uses `DoubleBasis` with `AbelianBasis` — `index(B)` interface unchanged
4. `AbelianBasisTest.jl` is NOT included in `runtests.jl` (it's a standalone file) — skip for now

- [ ] **Step 3: Fix any failures**

If tests fail, read the error messages carefully. The most likely issues are:
- The `(ag::AbelianOperator)(dgt)` callable now requires a `base` argument: `ag(dgt, base)`. The `schmidt` function in AbelianBasis.jl already passes `b.B`. Check if any other call sites need updating.
- Any direct construction of `AbelianOperator` with closures must be updated to use permutation arrays.

- [ ] **Step 4: Commit fixes**

```bash
git add -A
git commit -m "fix: resolve regression failures from AbelianBasis overhaul"
```

---

### Task 8: Entanglement Tests with 2D Basis

Verify entanglement calculations work correctly with the new 2D AbelianBasis.

**Files:**
- Modify: `test/abelian_overhaul_tests.jl`

- [ ] **Step 1: Write entanglement test for 2D basis**

Add to `test/abelian_overhaul_tests.jl`:
```julia
@testset "2D Entanglement" begin
    Lx, Ly = 3, 2
    L = Lx * Ly
    sites = [(x, y) for y in 0:Ly-1 for x in 0:Lx-1]

    T_x = [mod(x + 1, Lx) + Lx * y + 1 for (x, y) in sites]
    T_y = [x + Lx * mod(y + 1, Ly) + 1 for (x, y) in sites]

    XXZ = spin((1.0, "xx"), (1.0, "yy"), (0.5, "zz"))
    bonds = Tuple{Int,Int}[]
    for (x, y) in sites
        i = x + Lx * y + 1
        push!(bonds, (i, mod(x + 1, Lx) + Lx * y + 1))
        push!(bonds, (i, x + Lx * mod(y + 1, Ly) + 1))
    end

    # Full space calculation
    H_full = operator([XXZ for _ in bonds], [[b[1], b[2]] for b in bonds], L)
    E_full, V_full = eigen(Hermitian(Array(H_full)))
    gs_full = V_full[:, 1]  # ground state
    S_full = ent_S(gs_full, 1:L÷2, TensorBasis(; L, base=2))

    # Symmetry-reduced: reconstruct ground state entropy
    # At (kx=0, ky=0), the ground state should be in this sector
    B_sym = basis(; L, base=2, symmetries=[(T_x, 0), (T_y, 0)])
    if size(B_sym, 1) > 0
        H_sym = operator([XXZ for _ in bonds], [[b[1], b[2]] for b in bonds], B_sym)
        E_sym, V_sym = eigen(Hermitian(Array(H_sym)))
        gs_sym = V_sym[:, 1]
        S_sym = ent_S(gs_sym, 1:L÷2, B_sym)
        # Ground state energy should match
        @test E_sym[1] ≈ E_full[1] atol=1e-10
        # Entanglement entropy should match
        @test S_sym ≈ S_full atol=1e-10
    end
end
```

- [ ] **Step 2: Run test**

Run: `julia --project -e 'using Test; include("test/abelian_overhaul_tests.jl")'`

Expected: "2D Entanglement" tests PASS

- [ ] **Step 3: Commit**

```bash
git add test/abelian_overhaul_tests.jl
git commit -m "test: add 2D entanglement entropy verification"
```

---

### Task 9: Performance Benchmark

**Files:**
- Modify: `test/abelian_overhaul_tests.jl`

- [ ] **Step 1: Write performance comparison**

Add to `test/abelian_overhaul_tests.jl`:
```julia
@testset "Performance Sanity Check" begin
    # Verify basis construction completes in reasonable time for L=16
    L = 16
    t = @elapsed begin
        B = basis(; L, N=L÷2, k=0)
    end
    @test size(B, 1) > 0
    @info "L=$L, N=$(L÷2), k=0 basis construction: $(round(t, digits=3))s, dim=$(size(B, 1))"

    # Verify operator multiplication works
    XXZ = spin((1.0, "xx"), (1.0, "yy"), (0.5, "zz"))
    H = trans_inv_operator(XXZ, 2, B)
    v = randn(ComplexF64, size(B, 1))
    v ./= norm(v)
    t_mul = @elapsed (w = H * v)
    @test norm(w) > 0
    @info "Matrix-free multiply: $(round(t_mul, digits=4))s"

    # 2D performance check
    Lx, Ly = 4, 3
    L2 = Lx * Ly
    sites = [(x, y) for y in 0:Ly-1 for x in 0:Lx-1]
    T_x = [mod(x + 1, Lx) + Lx * y + 1 for (x, y) in sites]
    T_y = [x + Lx * mod(y + 1, Ly) + 1 for (x, y) in sites]
    t2 = @elapsed begin
        B2 = basis(; L=L2, N=L2÷2, base=2, symmetries=[(T_x, 0), (T_y, 0)])
    end
    @test size(B2, 1) > 0
    @info "2D 4x3 basis (N=$(L2÷2), kx=0, ky=0): $(round(t2, digits=3))s, dim=$(size(B2, 1))"
end
```

- [ ] **Step 2: Run test**

Run: `julia --project -e 'using Test; include("test/abelian_overhaul_tests.jl")'`

Expected: PASS with timing info printed. No specific timing thresholds — this is a sanity check.

- [ ] **Step 3: Commit**

```bash
git add test/abelian_overhaul_tests.jl
git commit -m "test: add performance sanity benchmarks for AbelianBasis"
```

---

### Task 10: Documentation

**Files:**
- Create: `docs/src/abelian_basis.md`

- [ ] **Step 1: Write documentation**

Create `docs/src/abelian_basis.md`:
```markdown
# AbelianBasis: General Symmetry-Reduced Bases

`AbelianBasis` is EDKit's most general symmetry-reduction mechanism. It constructs
a basis of states that are simultaneous eigenstates of an arbitrary collection of
commuting discrete symmetries, each specified as a permutation of lattice sites
(optionally with spin inversion).

## Quick Start

### 1D Chain (via convenience keywords)

For 1D periodic chains, the `basis()` constructor provides shorthand keywords:

```julia
# Heisenberg chain, L=12, half-filling, zero momentum, even parity
B = basis(; L=12, N=6, k=0, p=1)

# With spin-flip symmetry
B = basis(; L=12, N=6, k=0, p=1, z=1)
```

### 2D Lattice (via `symmetries` keyword)

For arbitrary lattice geometries, pass symmetry generators as permutation arrays:

```julia
Lx, Ly = 4, 3
L = Lx * Ly
sites = [(x, y) for y in 0:Ly-1 for x in 0:Lx-1]

# Define symmetry generators as permutation arrays
T_x = [mod(x+1, Lx) + Lx*y + 1 for (x,y) in sites]  # translation in x
T_y = [x + Lx*mod(y+1, Ly) + 1 for (x,y) in sites]   # translation in y
P_x = [Lx-1-x + Lx*y + 1 for (x,y) in sites]          # reflection in x

# Build basis with (kx=0, ky=0, px=+1) symmetry sector
B = basis(; L, N=L÷2, base=2,
    symmetries=[(T_x, 0), (T_y, 0), (P_x, 0)])
```

### 3D Lattice

The same approach works for any dimension:

```julia
Lx, Ly, Lz = 2, 2, 3
L = Lx * Ly * Lz
sites = [(x,y,z) for z in 0:Lz-1 for y in 0:Ly-1 for x in 0:Lx-1]

T_x = [mod(x+1,Lx) + Lx*y + Lx*Ly*z + 1 for (x,y,z) in sites]
T_y = [x + Lx*mod(y+1,Ly) + Lx*Ly*z + 1 for (x,y,z) in sites]
T_z = [x + Lx*y + Lx*Ly*mod(z+1,Lz) + 1 for (x,y,z) in sites]

B = basis(; L, symmetries=[(T_x, 0), (T_y, 0), (T_z, 0)])
```

## Permutation Array Format

A permutation array `perm` of length `L` defines a site mapping: `perm[i] = j`
means site `i` maps to site `j`. Arrays are 1-indexed (Julia convention).

Common patterns:
- **Translation**: `[mod1(i+a, L) for i in 1:L]` (shift by `a` sites)
- **Reflection**: `[L+1-i for i in 1:L]` (reverse site order)
- **Point-group**: any valid permutation of `1:L`

### Spin Inversion

Some symmetries also flip spins (e.g., spin-flip symmetry where `|up> <-> |down>`).
Specify this with a third element in the symmetry tuple:

```julia
# Spin-flip: identity permutation + full inversion
z_flip = (collect(1:L), 0, trues(L))

# Sublattice inversion: flip spins on even sites only
even_sites = BitVector([iseven(i) for i in 1:L])
z_sub = (collect(1:L), 0, even_sites)

B = basis(; L, symmetries=[z_flip])
```

## How It Works

### Group Structure

Each symmetry generator defines a cyclic group. The `AbelianOperator` combines
multiple generators into the direct product G = Z_{p1} x Z_{p2} x ... x Z_{pn},
where p_i is the period of generator i.

### Basis Construction

1. Enumerate candidate states (all `base^L` states, or only those with fixed
   particle number via Gosper's hack for base=2).
2. For each candidate, check if it is the canonical (minimum-index) representative
   of its symmetry orbit.
3. Store representatives with their normalization factors.

### Performance Optimizations

- **Benes networks (base=2):** Permutations are pre-compiled into O(log L) bit
  circuits for fast integer-state manipulation.
- **Integer orbit search (base=2):** Orbit representatives are found by operating
  directly on integer states, avoiding digit-buffer overhead.
- **Gosper's hack (base=2, fixed N):** Only states with exactly N particles are
  enumerated, reducing the search space by binomial(L,N)/2^L.
- **Multi-threading:** Basis construction is parallelized across threads.

## API Reference

### `basis(; L, symmetries, N, f, base, threaded)`

High-level constructor. Returns `AbelianBasis` when symmetries are specified.

**Arguments:**
- `L::Integer`: system size (number of sites)
- `symmetries`: vector of `(perm, q)` or `(perm, q, inv)` tuples
- `N::Integer` (optional): particle number constraint
- `f` (optional): custom predicate on digit vectors
- `base::Integer` (default 2): local Hilbert space dimension
- `threaded::Bool` (default `base^L > 3000`): enable multi-threading

### `AbelianOperator(order, k, perm; inv)`

Low-level generator constructor. Normally called internally by `basis()`.

**Arguments:**
- `order::Int`: period of the cyclic group
- `k::Integer`: character label (quantum number)
- `perm::Vector{Int}`: 1-indexed permutation array
- `inv::BitVector` (optional): spin-inversion mask

## Examples

### 2D Heisenberg Model

```julia
using EDKit, LinearAlgebra

Lx, Ly = 4, 3
L = Lx * Ly
sites = [(x, y) for y in 0:Ly-1 for x in 0:Lx-1]

# Symmetry generators
T_x = [mod(x+1, Lx) + Lx*y + 1 for (x,y) in sites]
T_y = [x + Lx*mod(y+1, Ly) + 1 for (x,y) in sites]

# Nearest-neighbor bonds
XXZ = spin((1.0, "xx"), (1.0, "yy"), (1.0, "zz"))
bonds = Tuple{Int,Int}[]
for (x, y) in sites
    i = x + Lx * y + 1
    push!(bonds, (i, mod(x+1, Lx) + Lx*y + 1))
    push!(bonds, (i, x + Lx*mod(y+1, Ly) + 1))
end

# Ground state in (kx=0, ky=0, N=L/2) sector
B = basis(; L, N=L÷2, base=2, symmetries=[(T_x, 0), (T_y, 0)])
H = operator([XXZ for _ in bonds], [[b[1], b[2]] for b in bonds], B)
E, V = eigen(Hermitian(Array(H)))
println("Ground state energy per site: ", E[1] / L)
```

### Triangular Lattice

```julia
Lx, Ly = 3, 3
L = Lx * Ly
sites = [(x, y) for y in 0:Ly-1 for x in 0:Lx-1]

T_x = [mod(x+1, Lx) + Lx*y + 1 for (x,y) in sites]
T_y = [x + Lx*mod(y+1, Ly) + 1 for (x,y) in sites]

# Triangular lattice bonds: horizontal, vertical, and diagonal
bonds = Tuple{Int,Int}[]
for (x, y) in sites
    i = x + Lx * y + 1
    push!(bonds, (i, mod(x+1, Lx) + Lx*y + 1))           # horizontal
    push!(bonds, (i, x + Lx*mod(y+1, Ly) + 1))             # vertical
    push!(bonds, (i, mod(x+1, Lx) + Lx*mod(y+1, Ly) + 1)) # diagonal
end

B = basis(; L, N=L÷2, base=2, symmetries=[(T_x, 0), (T_y, 0)])
H = operator([XXZ for _ in bonds], [[b[1], b[2]] for b in bonds], B)
E = eigvals(Hermitian(Array(H)))
println("Triangular Heisenberg ground state energy: ", E[1])
```
```

- [ ] **Step 2: Commit documentation**

```bash
git add docs/src/abelian_basis.md
git commit -m "docs: add AbelianBasis documentation with 2D/3D examples"
```

---

### Task 11: Final Full Test Run and Cleanup

- [ ] **Step 1: Run full test suite one final time**

Run: `cd /Users/ren/Library/CloudStorage/OneDrive-UniversityofLeeds/GitHub/EDKit.jl && julia --project -e 'using Pkg; Pkg.test()'`

Expected: ALL tests PASS

- [ ] **Step 2: Run the standalone AbelianBasisTest.jl to verify compatibility**

Run: `cd /Users/ren/Library/CloudStorage/OneDrive-UniversityofLeeds/GitHub/EDKit.jl && julia --project test/AbelianBasisTest.jl`

Expected: ALL tests PASS. If this file constructs `AbelianOperator` with closures (it uses `basis()` which does it internally), it should work with the new array-based internals since the `basis()` constructor now generates permutation arrays.

- [ ] **Step 3: Update CLAUDE.md with new API information**

Add to `CLAUDE.md` after the "Conventions" section:

```markdown
## 2D/3D Lattice Support

`AbelianBasis` supports arbitrary lattice geometries via the `symmetries` keyword in `basis()`. Each symmetry generator is a permutation array:

```julia
# 2D square lattice example
basis(; L=12, N=6, symmetries=[(T_x, kx), (T_y, ky)])
```

For base=2 systems, Benes networks provide O(log L) bit permutation and Gosper's hack accelerates particle-conserving enumeration.
```

- [ ] **Step 4: Final commit**

```bash
git add CLAUDE.md
git commit -m "docs: update CLAUDE.md with 2D/3D lattice support info"
```
