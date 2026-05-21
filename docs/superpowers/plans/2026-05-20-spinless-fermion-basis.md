# SpinlessFermionBasis Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:executing-plans`. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `SpinlessFermionBasis` type plus a small fermion-operator helper layer that bakes the left-string Jordan-Wigner transformation into local operator matrices, so users can write fermion Hamiltonians without hand-rolling σ_z chains.

**Architecture:** `SpinlessFermionBasis` is a new `AbstractOnsiteBasis` subtype whose Hilbert space is identical to `ProjectedBasis` at `base=2` (occupation labels per site, optional fixed-N filter), but tags the type for fermion-aware operator construction. A `fermion(op, span)` helper builds local matrices (with JW chain) and `fermion_operator(op, sites, B)` wraps `operator(...)`. The JW convention is the **left** string: `c_j = (∏_{k<j} σ_z^k) σ⁺_j` in EDKit's spin-coordinate convention (where dgt=0 ↔ empty, dgt=1 ↔ occupied, so the fermionic creation `c†` maps to spin lowering `σ⁻`).

**Tech Stack:** Julia ≥ 1.9, SparseArrays, EDKit internals (`AbstractOnsiteBasis`, `operator`, `spin_Sp`, `spin_Sm`, `spin_Sz`).

**Scope (MUST):**
- Basis type with optional `N` filter
- `fermion(:n)`, `fermion(:c†c, span)`, `fermion(:cc, span)`, `fermion(:c†c†, span)` local matrices
- `fermion_operator(op, sites, B)` high-level constructor
- Tests against hand-written JW spin operators on `TensorBasis`

**Scope (DEFERRED to a separate plan):**
- Spinful fermion basis (separate ↑/↓ tracking)
- Symmetry-reduced fermion bases (translation, parity)
- 4-fermion interaction operators in one call
- Operator-string parsing API (`fermion("c†c", [i,j])`)

**File structure:**
- Create: `src/Basis/SpinlessFermionBasis.jl` — the basis type
- Create: `src/FermionOperator.jl` — fermion operator helpers
- Modify: `src/EDKit.jl` — register the two new files (after `AbelianBasis.jl` and after `Operator.jl` respectively)
- Create: `test/fermion_tests.jl` — tests
- Modify: `test/runtests.jl` — include the new test file
- Create: `docs/src/manual/fermions.md` — usage docs (after the MVP is green)

---

## Task 1: SpinlessFermionBasis skeleton

**Files:**
- Create: `src/Basis/SpinlessFermionBasis.jl`
- Modify: `src/EDKit.jl`
- Test: `test/fermion_tests.jl` (create), `test/runtests.jl`

- [ ] **Step 1.1: Create the test file with a failing skeleton test**

Create `test/fermion_tests.jl`:

```julia
@testset "Spinless Fermion Basis" begin
    @testset "Basis construction" begin
        B = SpinlessFermionBasis(L = 4)
        @test size(B, 1) == 2^4
        @test size(B, 2) == 2^4
        @test B.B == 2
        @test length(B.I) == 2^4

        Bn = SpinlessFermionBasis(L = 6, N = 3)
        @test size(Bn, 1) == binomial(6, 3)

        # N=1 must mean exactly one occupied site (dgt=1), not one empty site.
        # binomial(6,1) is symmetric so size alone wouldn't catch a swapped
        # convention — check digit sums directly.
        Bone = SpinlessFermionBasis(L = 6, N = 1)
        @test size(Bone, 1) == 6
        for i in 1:size(Bone, 1)
            change!(Bone, i)
            @test sum(Bone.dgt) == 1
        end
    end
end
```

Register it in `test/runtests.jl` by appending:

```julia
include("fermion_tests.jl")
```

- [ ] **Step 1.2: Run test to verify it fails**

```bash
cd /Users/ren/Library/CloudStorage/OneDrive-UniversityofLeeds/GitHub/EDKit.jl
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -20
```

Expected: `UndefVarError: SpinlessFermionBasis not defined`.

- [ ] **Step 1.3: Create the basis file with struct + constructor**

Create `src/Basis/SpinlessFermionBasis.jl`:

```julia
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
    I = if isnothing(N)
        threaded ? selectindex_threaded(f, L, base=base, alloc=alloc) :
                   selectindex(f, L, 1:base^L, base=base, alloc=alloc)
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
```

Note: `dgt=0` means empty, `dgt=1` means occupied. The `N` keyword counts occupations directly (number of `dgt=1` entries). This differs from `ProjectedBasis`, where `N` follows EDKit's spin convention (number of `dgt=0` entries) — the constructor compensates by passing `L - N` to `selectindex_N`.

- [ ] **Step 1.4: Register the file in `src/EDKit.jl`**

Open `src/EDKit.jl`. Find the line that includes `Basis/AbelianBasis.jl` (around line 34) and append after it:

```julia
include("Basis/SpinlessFermionBasis.jl")
```

- [ ] **Step 1.5: Run test to verify pass**

```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -10
```

Expected: PASS.

- [ ] **Step 1.6: Commit**

```bash
git add src/Basis/SpinlessFermionBasis.jl src/EDKit.jl test/fermion_tests.jl test/runtests.jl
git commit -m "$(cat <<'EOF'
Add SpinlessFermionBasis skeleton

A new AbstractOnsiteBasis subtype that mirrors ProjectedBasis at
base=2 but tags the type for fermion-aware operator construction.
Constructor supports optional fixed-N particle sectors and reuses
selectindex/_threaded/_N from the projected-basis machinery.

Operator dispatch (fermion(), fermion_operator()) lands in a
follow-up commit.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Basis methods (copy, index)

**Files:**
- Modify: `src/Basis/SpinlessFermionBasis.jl`
- Test: `test/fermion_tests.jl`

- [ ] **Step 2.1: Add failing tests**

Append to the `@testset "Basis construction"` block in `test/fermion_tests.jl`:

```julia
        # copy gives an independent digit buffer
        B = SpinlessFermionBasis(L = 4)
        Bc = copy(B)
        @test Bc isa SpinlessFermionBasis
        @test Bc.I === B.I    # share the representative list
        @test Bc.dgt !== B.dgt # but own dgt

        # index round-trip
        change!(B, 5)
        coeff, pos = index(B)
        @test coeff == 1
        @test pos == 5
```

- [ ] **Step 2.2: Run to verify failure**

```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | grep -A3 "Test Failed\|got exception" | head -20
```

Expected: MethodError on `copy(::SpinlessFermionBasis)` and `index(::SpinlessFermionBasis)`.

- [ ] **Step 2.3: Add methods to `src/Basis/SpinlessFermionBasis.jl`**

Append to the file:

```julia
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
```

- [ ] **Step 2.4: Also extend `index_nocheck` in `src/LinearMap.jl`**

Open `src/LinearMap.jl`. Find lines 13-16 (`index_nocheck` for ProjectedBasis). Append below them:

```julia
index_nocheck(b::SpinlessFermionBasis) = index(b, check=false)
index_nocheck(b::SpinlessFermionBasis, dgt::AbstractVector) = index(b, dgt; check=false)
```

- [ ] **Step 2.5: Run to verify pass**

```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -10
```

Expected: PASS.

- [ ] **Step 2.6: Commit**

```bash
git add src/Basis/SpinlessFermionBasis.jl src/LinearMap.jl test/fermion_tests.jl
git commit -m "$(cat <<'EOF'
Implement SpinlessFermionBasis copy / index methods

Mirror the ProjectedBasis implementations of copy() and index() so
SpinlessFermionBasis works with all generic algorithms that expect
the (coefficient, position) protocol. Also extend index_nocheck in
LinearMap.jl so DoubleBasis/symmetrizer machinery accepts the new
basis type.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Number operator + nearest-neighbor hop

**Files:**
- Create: `src/FermionOperator.jl`
- Modify: `src/EDKit.jl`
- Test: `test/fermion_tests.jl`

**Background:** In EDKit's spin convention, `dgt=0` is "spin up" (Sz=+½) and `dgt=1` is "spin down" (Sz=-½). For spinless fermions we identify `dgt=0 ↔ empty`, `dgt=1 ↔ occupied`, so:
- `c_j` (annihilation, takes occupied→empty) ↔ `σ⁺_j` (spin raising in EDKit's convention)
- `c†_j` (creation, takes empty→occupied) ↔ `σ⁻_j`
- `n_j = c†_j c_j = σ⁻ σ⁺ = diag([0, 1])` (no JW string — number is local)

Left-string JW: `c_j = (σ_z^1 σ_z^2 ⋯ σ_z^{j-1}) σ⁺_j`, `c†_j = (σ_z^1 ⋯ σ_z^{j-1}) σ⁻_j`.

For two-fermion products, the inter-site σ_z chains partially cancel. For nearest neighbors (j = i+1) there is no chain in between, so `c†_i c_{i+1}` is just `σ⁻_i σ⁺_{i+1}` (a `kron(σ⁻, σ⁺)` 4×4 matrix).

- [ ] **Step 3.1: Add failing tests**

Append to `test/fermion_tests.jl`:

```julia
    @testset "Number and nearest-neighbor hop" begin
        # Number operator: 1-site, diag([0, 1])
        @test fermion(:n) ≈ [0.0 0.0; 0.0 1.0]

        # Nearest-neighbor c†_1 c_2 on 2 sites: σ⁻ ⊗ σ⁺ in EDKit convention
        # (c† = σ⁻ at site 1, c = σ⁺ at site 2, no JW chain between adjacent sites)
        expected = kron(Array(EDKit.spin_Sm(2)), Array(EDKit.spin_Sp(2)))
        @test fermion(:c†c, 2) ≈ expected

        # 1-site n operator embedded as an Operator on a TensorBasis matches σ⁻σ⁺
        B = SpinlessFermionBasis(L = 4)
        Bt = TensorBasis(L = 4, base = 2)
        n_op_fermion = Array(fermion_operator(:n, [2], B))
        n_op_spin = Array(operator(diagm([0.0, 1.0]), [2], Bt))
        @test n_op_fermion ≈ n_op_spin
    end
```

- [ ] **Step 3.2: Run to verify failure**

```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | grep -A3 "Test Failed\|got exception\|UndefVarError" | head -10
```

Expected: `UndefVarError: fermion`, `fermion_operator`.

- [ ] **Step 3.3: Create `src/FermionOperator.jl`**

```julia
export fermion, fermion_operator

"""
    fermion(op::Symbol[, span::Integer]; convention::Symbol=:left)

Return the local matrix representation of a spinless-fermion operator on a
contiguous block of `span` sites.

Supported `op` values:
- `:n` — the single-site number operator (`span` not used, default 1).
- `:c†c` — two-fermion product `c†_1 c_span`. The Jordan-Wigner string of
  `σ_z` is inserted on sites `2, 3, …, span-1` whenever `span > 2`.
- `:cc` — pair annihilation `c_1 c_span` (anticommuting; includes the
  appropriate sign).
- `:c†c†` — pair creation `c†_1 c†_span`.

`convention` selects the Jordan-Wigner direction. Only `:left` is implemented
in this MVP; the docstring keeps the keyword so a future `:right` variant can
be added without breaking callers.

The returned matrix lives on `span` adjacent sites; embed it into a many-body
basis with [`operator`](@ref) or [`fermion_operator`](@ref).
"""
function fermion(op::Symbol, span::Integer=1; convention::Symbol=:left)
    convention === :left || error("Only the :left Jordan-Wigner convention is implemented (got $convention).")
    if op === :n
        span == 1 || error(":n is a single-site operator (got span=$span).")
        return [0.0 0.0; 0.0 1.0]
    end
    span >= 2 || error("Two-fermion operators require span >= 2 (got $span).")
    σ⁺ = Array(spin_Sp(2))
    σ⁻ = Array(spin_Sm(2))
    σ_z = Array(spin_Sz(2)) * 2  # spin_Sz returns Sz = σ_z / 2; we want σ_z = diag(1, -1)
    Id = Matrix{Float64}(I, 2, 2)

    # Build the chain of one-site factors for sites 1..span
    factors = Vector{Matrix{Float64}}(undef, span)
    if op === :c†c
        # c†_1 c_span = σ⁻_1 σ_z^2 ⋯ σ_z^{span-1} σ⁺_span
        factors[1] = σ⁻
        for k in 2:span-1; factors[k] = σ_z; end
        factors[span] = σ⁺
        sign = 1.0
    elseif op === :cc
        # c_1 c_span = -σ⁺_1 σ_z^2 ⋯ σ_z^{span-1} σ⁺_span
        # The − sign comes from the on-site product σ⁺ σ_z = −σ⁺ at site 1
        # (σ⁺ anticommutes with σ_z, whereas σ⁻ commutes with σ_z in matrix
        # product form: σ⁻ σ_z = +σ⁻ — see :c†c† below).
        factors[1] = σ⁺
        for k in 2:span-1; factors[k] = σ_z; end
        factors[span] = σ⁺
        sign = -1.0
    elseif op === :c†c†
        # c†_1 c†_span = +σ⁻_1 σ_z^2 ⋯ σ_z^{span-1} σ⁻_span
        # On-site product σ⁻ σ_z = +σ⁻ contributes no overall sign — this is
        # asymmetric with :cc above and is required for {c†_i, c†_j}=0 to hold
        # (cross-check: c†_i c†_j = -c†_j c†_i).
        factors[1] = σ⁻
        for k in 2:span-1; factors[k] = σ_z; end
        factors[span] = σ⁻
        sign = 1.0
    else
        error("Unsupported fermion operator: $op. Supported: :n, :c†c, :cc, :c†c†.")
    end

    mat = factors[1]
    for k in 2:span
        mat = kron(mat, factors[k])
    end
    sign * mat
end

"""
    fermion_operator(op::Symbol, sites::AbstractVector{<:Integer},
                     B::AbstractBasis; convention::Symbol=:left)

Build the full-system [`Operator`](@ref) corresponding to a fermion operator
`op` placed on the supplied `sites`.

`sites` must be a one- or two-element vector. The function determines the
contiguous support `min(sites):max(sites)`, builds the local matrix via
[`fermion`](@ref) (including any Jordan-Wigner `σ_z` chain on intermediate
sites), then dispatches to [`operator`](@ref).

If `sites` is unsorted (e.g. `[3, 1]` for `c†_3 c_1`), the operator is
adjusted by the appropriate fermionic sign so that the returned matrix
represents the operator literally written by the caller.
"""
function fermion_operator(op::Symbol, sites::AbstractVector{<:Integer},
        B::AbstractBasis; convention::Symbol=:left)
    convention === :left || error("Only the :left convention is implemented.")
    if op === :n
        length(sites) == 1 || error(":n takes exactly one site.")
        return operator(fermion(:n), [sites[1]], B)
    end
    length(sites) == 2 || error("Two-fermion operators take exactly two sites.")
    i, j = sites[1], sites[2]
    i == j && error(":$op on the same site is identically zero (Pauli exclusion).")
    swapped = i > j
    lo, hi = minmax(i, j)
    span = hi - lo + 1
    local_op = if swapped
        # Swapped-site reductions (lo < hi, sites[1] > sites[2]):
        #  c†_hi c_lo  = (c†_lo c_hi)†         (Hermitian conjugate, no sign)
        #  c_hi c_lo   = -c_lo c_hi             (fermionic anticommutation)
        #  c†_hi c†_lo = -c†_lo c†_hi          (fermionic anticommutation)
        if op === :c†c
            fermion(:c†c, span)'
        elseif op === :cc
            -fermion(:cc, span)
        elseif op === :c†c†
            -fermion(:c†c†, span)
        else
            error("Unsupported swap-mode fermion operator: $op")
        end
    else
        fermion(op, span)
    end
    operator(local_op, collect(lo:hi), B)
end
```

- [ ] **Step 3.4: Register the new file in `src/EDKit.jl`**

Open `src/EDKit.jl`. Find the `include("Operator.jl")` line (around line 39) and append after it:

```julia
include("FermionOperator.jl")
```

- [ ] **Step 3.5: Run to verify pass**

```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -10
```

Expected: PASS.

- [ ] **Step 3.6: Commit**

```bash
git add src/FermionOperator.jl src/EDKit.jl test/fermion_tests.jl
git commit -m "$(cat <<'EOF'
Add fermion() and fermion_operator() helpers (n, c†c, cc, c†c†)

fermion(op, span) returns local matrices for spinless fermion
operators on a contiguous block of sites, with the left-string
Jordan-Wigner σ_z chain baked in for spans > 2. fermion_operator()
wraps operator() to build the full-system operator on a given basis.

This MVP supports the four operators that cover free-fermion hopping
(c†c), BCS-style pairing (cc, c†c†), and chemical potential terms
(n). Higher-order products (4-fermion interactions in one call) are
left for a follow-up.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Non-adjacent hop with JW string verification

**Files:**
- Test: `test/fermion_tests.jl`

The implementation already supports `span > 2`. This task adds the regression test that proves it.

- [ ] **Step 4.1: Add a failing test that compares against hand-written JW**

Append to `test/fermion_tests.jl`:

```julia
    @testset "Non-adjacent hop carries the JW string" begin
        # c†_1 c_3 on a 3-site block: σ⁻_1 σ_z^2 σ⁺_3
        σ⁺ = Array(EDKit.spin_Sp(2))
        σ⁻ = Array(EDKit.spin_Sm(2))
        σ_z = [1.0 0.0; 0.0 -1.0]
        expected = kron(kron(σ⁻, σ_z), σ⁺)
        @test fermion(:c†c, 3) ≈ expected

        # c†_1 c_4 on a 4-site block: σ⁻_1 σ_z^2 σ_z^3 σ⁺_4
        expected4 = kron(kron(kron(σ⁻, σ_z), σ_z), σ⁺)
        @test fermion(:c†c, 4) ≈ expected4

        # Compare a many-body matrix built via fermion_operator against the
        # explicit JW spin construction on a TensorBasis.
        L = 5
        Bf = SpinlessFermionBasis(L = L)
        Bt = TensorBasis(L = L, base = 2)
        # c†_1 c_4
        H_fermion = Array(fermion_operator(:c†c, [1, 4], Bf))
        H_spin = Array(operator(kron(σ⁻, σ_z, σ_z, σ⁺), [1, 2, 3, 4], Bt))
        @test H_fermion ≈ H_spin
        # c†_2 c_5
        H_fermion2 = Array(fermion_operator(:c†c, [2, 5], Bf))
        H_spin2 = Array(operator(kron(σ⁻, σ_z, σ_z, σ⁺), [2, 3, 4, 5], Bt))
        @test H_fermion2 ≈ H_spin2

        # Pair operators: :cc has a leading − from σ⁺σ_z = −σ⁺, while :c†c† has
        # no overall sign because σ⁻σ_z = +σ⁻.
        @test fermion(:cc, 2) ≈ -kron(σ⁺, σ⁺)
        @test fermion(:c†c†, 2) ≈ kron(σ⁻, σ⁻)
        @test fermion(:cc, 4) ≈ -kron(kron(kron(σ⁺, σ_z), σ_z), σ⁺)
        @test fermion(:c†c†, 4) ≈ kron(kron(kron(σ⁻, σ_z), σ_z), σ⁻)

        # Swapped sites: c†_3 c_1 = (c†_1 c_3)† = +σ⁺_1 σ_z^2 σ⁻_3.
        # Verify both against the hand-built JW spin operator.
        L = 5
        Bf3 = SpinlessFermionBasis(L = L)
        Bt3 = TensorBasis(L = L, base = 2)
        forward = Array(fermion_operator(:c†c, [1, 3], Bf3))
        forward_spin = Array(operator(kron(σ⁻, σ_z, σ⁺), [1, 2, 3], Bt3))
        @test forward ≈ forward_spin
        reverse = Array(fermion_operator(:c†c, [3, 1], Bf3))
        reverse_spin = Array(operator(kron(σ⁺, σ_z, σ⁻), [1, 2, 3], Bt3))
        @test reverse ≈ reverse_spin

        # Fermionic anticommutation: {c†_1, c†_3} = 0 (cross-check on the signs
        # produced by both the natural and the swapped branch of :c†c†).
        c13 = Array(fermion_operator(:c†c†, [1, 3], Bf3))
        c31 = Array(fermion_operator(:c†c†, [3, 1], Bf3))
        @test c13 + c31 ≈ zeros(2^L, 2^L) atol = 1e-12
    end
```

- [ ] **Step 4.2: Run to verify pass**

```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -10
```

Expected: PASS (the implementation from Task 3 already covers this).

- [ ] **Step 4.3: Commit**

```bash
git add test/fermion_tests.jl
git commit -m "$(cat <<'EOF'
Add JW-string regression test for non-adjacent fermion hops

Verify that fermion(:c†c, span) for span > 2 inserts the σ_z chain
on intermediate sites and that the resulting full-system Operator on
a SpinlessFermionBasis matches the hand-written spin construction on
a TensorBasis.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Hermiticity and free-fermion spectrum sanity

**Files:**
- Test: `test/fermion_tests.jl`

Add two more physics-grounded checks: the hopping operator is Hermitian after symmetrization, and the free-fermion tight-binding chain has the expected single-particle spectrum.

- [ ] **Step 5.1: Add failing tests**

Append to `test/fermion_tests.jl`:

```julia
    @testset "Hermiticity and free-fermion spectrum" begin
        L = 6
        Bf = SpinlessFermionBasis(L = L, N = 1)  # single-particle sector

        # H = -Σ_i (c†_i c_{i+1} + h.c.) with PBC
        t = 1.0
        H = sum(
            -t * (
                fermion_operator(:c†c, [i, mod1(i + 1, L)], Bf) +
                fermion_operator(:c†c, [mod1(i + 1, L), i], Bf)
            )
            for i in 1:L
        )
        Hm = Array(H)
        @test Hm ≈ Hm'  # Hermitian

        # Single-particle spectrum: -2t cos(2π k / L) for k = 0, 1, …, L-1
        eigvals_h = eigvals(Hermitian(Hm))
        analytic = sort([-2 * t * cos(2π * k / L) for k in 0:L-1])
        @test eigvals_h ≈ analytic atol = 1e-10
    end
```

- [ ] **Step 5.2: Run to verify pass**

```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -10
```

Expected: PASS.

- [ ] **Step 5.3: Commit**

```bash
git add test/fermion_tests.jl
git commit -m "$(cat <<'EOF'
Add Hermiticity and free-fermion spectrum tests for fermion_operator

Build H = -t Σ_i (c†_i c_{i+1} + h.c.) on the single-particle sector
of a SpinlessFermionBasis and verify (a) the matrix is Hermitian and
(b) its eigenvalues match the analytic -2t cos(2π k/L) tight-binding
dispersion at L=6.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Documentation page

**Files:**
- Create: `docs/src/manual/fermions.md`

- [ ] **Step 6.1: Create the docs page**

Create `docs/src/manual/fermions.md`:

```markdown
# Spinless fermions

EDKit provides a `SpinlessFermionBasis` type and a small helper layer that
inserts the Jordan-Wigner string automatically when you build local fermion
operators.

## Basis

```julia
B = SpinlessFermionBasis(L = 8)               # full Hilbert space
Bn = SpinlessFermionBasis(L = 8, N = 4)       # fixed-N sector
```

Each digit `dgt[i]` is `0` for an empty site and `1` for an occupied site.
The basis is internally a `base=2` projected occupation basis; the distinct
type lets [`fermion_operator`](@ref) dispatch on it.

## Local operators

The `fermion(op, span)` helper returns the local matrix for an operator on a
contiguous block of `span` sites, with the Jordan-Wigner `σ_z` chain inserted
on the intermediate sites:

| `op` | meaning | span |
|------|---------|------|
| `:n` | number operator `n_i = c†_i c_i` | 1 |
| `:c†c` | hop `c†_1 c_span` | ≥ 2 |
| `:cc` | pair annihilation `c_1 c_span` | ≥ 2 |
| `:c†c†` | pair creation `c†_1 c†_span` | ≥ 2 |

The default Jordan-Wigner convention is left-string:
`c_j = (σ_z^1 σ_z^2 ⋯ σ_z^{j-1}) σ⁺_j`. EDKit identifies `dgt=0` with empty
and `dgt=1` with occupied, so `c†` maps to `σ⁻` and `c` maps to `σ⁺`.

## Full operators

`fermion_operator(op, sites, B)` embeds the local matrix into the many-body
basis:

```julia
B = SpinlessFermionBasis(L = 6)

# Number on site 3
N3 = fermion_operator(:n, [3], B)

# Nearest-neighbor hop c†_1 c_2
hop = fermion_operator(:c†c, [1, 2], B)

# Long-range hop c†_1 c_5 with σ_z chain on sites 2, 3, 4
long_hop = fermion_operator(:c†c, [1, 5], B)
```

A simple tight-binding chain with PBC:

```julia
t = 1.0
H = sum(
    -t * (
        fermion_operator(:c†c, [i, mod1(i+1, L)], B) +
        fermion_operator(:c†c, [mod1(i+1, L), i], B)
    )
    for i in 1:L
)
```

## Limitations of the MVP

The current implementation supports:

- Single-site `n`
- Two-fermion products: `c†c`, `cc`, `c†c†` (with the JW string)

Not yet supported:

- Spinful fermion basis (separate ↑/↓ tracking)
- Four-fermion interaction terms in a single call (build them as products of
  `n` operators or as sums of `c†c` calls)
- Symmetry-reduced fermion bases (translation, parity)
- Operator-string parsing API
```

- [ ] **Step 6.2: Verify the docs build (if Documenter is set up)**

If `docs/Project.toml` exists and `docs/make.jl` defines a build, run it:

```bash
cd docs && julia --project=. make.jl 2>&1 | tail -10 && cd ..
```

If no docs build is wired up, skip this verification — the markdown file is still useful as a reference.

- [ ] **Step 6.3: Commit**

```bash
git add docs/src/manual/fermions.md
git commit -m "$(cat <<'EOF'
Document SpinlessFermionBasis and fermion operator helpers

A short manual page covering basis construction, the supported
fermion local operators, and a tight-binding example. Also enumerates
the explicit scope limits of the MVP (no spinful, no 4-fermion calls,
no symmetry reduction) so users know where to expect rough edges.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Self-review checklist

After all tasks complete:

- [ ] `julia --project -e 'using Pkg; Pkg.test()'` — full suite green
- [ ] `JULIA_NUM_THREADS=8 julia --project -t 8 -e 'using Pkg; Pkg.test()'` — green at 8 threads
- [ ] `git log --oneline -8` shows six new commits, one per task
- [ ] No `@assert` introduced that would silently disable under `--check-bounds=no`
- [ ] `fermion_operator` does not export internal helper functions that aren't part of the public API
- [ ] Docstring of `SpinlessFermionBasis` clarifies the `dgt=0 ↔ empty, dgt=1 ↔ occupied` convention so users don't get confused by the spin↔fermion sign mapping

## Out of scope (separate plans)

- Spinful fermion basis with `(N↑, N↓)` sectors
- 4-fermion interaction calls (`fermion_operator(:c†c†cc, [i, j, k, l], B)`)
- Operator-string parsing (`fermion_operator("+-", [i, j], B)`)
- Symmetry-resolved fermion bases (momentum, parity, particle-hole)
- Fermionic Schmidt decomposition / entanglement-spectrum conventions
