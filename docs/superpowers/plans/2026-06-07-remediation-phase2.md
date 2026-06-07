# EDKit.jl Remediation — Phase 2: Integer-width sweep

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Convert every silent integer-overflow path in basis construction into a loud, clear error, and make `AbelianBasis` work with a non-default index `dtype`.

**Architecture:** One shared constructor-time capacity guard, `_check_index_capacity(dtype, base, L)`, computed in `BigInt` so it never suffers the overflow it detects. Apply it to every basis constructor. Separately, fix `AbelianBasis`'s non-default-`dtype` crash by converting the representative vector and base to `dtype` at the struct-construction boundary.

**Tech Stack:** Julia ≥1.10, `Test`, `Pkg.test()`. Branch `remediation-0.6.0`.

**Findings covered:** `bug-1`, `bug-2`, `bug-3` (capacity guard), `abelian-2` (dtype threading), `abelian-7` (subsumed by guard + existing Gosper guard), `parityflip-3` (silent path closed by guard; cosmetic `M::Ti` parked).

---

## Key design facts (verified against the code)

- The largest 1-based index for a base-`b`, `L`-site system is `b^L` (index = number + 1, number ∈ [0, b^L−1]). So the index type must satisfy `b^L < typemax(dtype)`. The `<` (strict, via `>=` in the guard) leaves +1 headroom for `FlipBasis`/`ParityFlipBasis`'s `M = base^L + 1`.
- For the default `Int64`, the guard caps base-2 at `L ≤ 62` (`2^63 ≥ typemax(Int64)`), exactly the bound `_gosper_enumerate` already enforces — so the `abelian-7` inconsistency disappears for the default dtype.
- `AbelianBasis`'s internal selection helpers use `Int` representatives; with the capacity guard guaranteeing values fit `dtype`, converting `I`/`base` to `dtype` only at line ~639 is sufficient and avoids threading `dtype` through every helper.
- Existing tests use small `L` (≤ ~16), far below any guard boundary, so the guard cannot regress them.

---

## Task 2.1: Capacity guard for `bug-1`, `bug-2`, `bug-3`

**Files:**
- Modify: `src/Basis/AbstractBasis.jl` (add helper after `_charge_predicate`, ~line 499)
- Modify: `src/Basis/AbstractBasis.jl` (TensorBasis constructor, ~line 246)
- Modify: `src/Basis/ProjectedBasis.jl` (constructor, ~line 242)
- Modify: `src/Basis/AbelianBasis.jl` (constructor top, ~line 606)
- Modify: `test/remediation_tests.jl` (append Phase 2 section)

- [ ] **Step 1: Append the failing tests to `test/remediation_tests.jl`**

```julia
@testset "bug-1/2/3: index capacity guard errors instead of silent overflow" begin
    # bug-1: Int32 + base>2 overflows the index type -> must error, not store wrapped reps.
    @test_throws ErrorException ProjectedBasis(Int32; L=20, N=0, base=3)
    @test_throws ErrorException basis(Int32; L=20, N=0, base=3)
    # bug-2: base=2 at L=63 overflowed the Gosper enumerator to an empty basis -> must error.
    @test_throws ErrorException ProjectedBasis(L=63, N=1, base=2)
    # bug-3: TensorBasis size overflowed to a negative value at L>=63 -> must error.
    @test_throws ErrorException TensorBasis(L=63, base=2)
    # Boundary positive: the largest valid base-2 Int64 size (L=62) still constructs.
    Bok = ProjectedBasis(L=62, N=1, base=2)
    @test size(Bok, 1) == 62
    # Ordinary small cases are unaffected.
    @test size(TensorBasis(L=10, base=2), 1) == 1024
    @test size(ProjectedBasis(L=8, N=4, base=2), 1) == binomial(8, 4)
end
```

- [ ] **Step 2: Run the new suite to verify it FAILS**

Run: `julia --project test/remediation_tests.jl`
Expected: FAIL — the four `@test_throws` currently do not throw (silent overflow / empty basis / negative size).

- [ ] **Step 3: Add the guard helper**

In `src/Basis/AbstractBasis.jl`, find:

```julia
@inline _charge_predicate(::Nothing, num::Integer) = x -> sum(x) == num
@inline _charge_predicate(f, num::Integer) = x -> (sum(x) == num && f(x))
```

and insert immediately after it:

```julia
#-------------------------------------------------------------------------------------------------------------------------
"""
    _check_index_capacity(dtype, base, L)

Error if a base-`base` system on `L` sites cannot be indexed in `dtype` without
overflow. The largest 1-based index is `base^L`; some bases also store
`base^L + 1`, so we require `base^L < typemax(dtype)`. Computed in `BigInt` so the
check itself never overflows.
"""
function _check_index_capacity(dtype::DataType, base::Integer, L::Integer)
    if big(base)^L >= typemax(dtype)
        error(
            "A base-$base system on $L sites needs 1-based indices up to $(big(base)^L), " *
            "which does not fit the index type $dtype (typemax = $(typemax(dtype))). " *
            "Construct the basis with a wider index dtype, e.g. Int128."
        )
    end
    nothing
end
```

- [ ] **Step 4: Apply the guard in the TensorBasis constructor**

In `src/Basis/AbstractBasis.jl`, find:

```julia
function TensorBasis(;L::Integer, base::Integer=2)
    dgt = zeros(Int64, L)
    B = Int64(base)
    TensorBasis(dgt, B)
end
```

replace with:

```julia
function TensorBasis(;L::Integer, base::Integer=2)
    _check_index_capacity(Int64, base, L)
    dgt = zeros(Int64, L)
    B = Int64(base)
    TensorBasis(dgt, B)
end
```

- [ ] **Step 5: Apply the guard in the ProjectedBasis constructor**

In `src/Basis/ProjectedBasis.jl`, find:

```julia
    base = convert(dtype, base)
    I = if isnothing(N)
```

replace with:

```julia
    _check_index_capacity(dtype, base, L)
    base = convert(dtype, base)
    I = if isnothing(N)
```

- [ ] **Step 6: Apply the guard in the AbelianBasis constructor**

In `src/Basis/AbelianBasis.jl`, find:

```julia
    Ng = order(G)

    C = zeros(Ng)
```

replace with:

```julia
    _check_index_capacity(dtype, base, L)
    Ng = order(G)

    C = zeros(Ng)
```

- [ ] **Step 7: Run the new suite to verify it PASSES**

Run: `julia --project test/remediation_tests.jl`
Expected: PASS — all Phase 2 capacity-guard assertions green; Phase 0/1 testsets still green.

- [ ] **Step 8: Run the full suite (no regressions)**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS.

- [ ] **Step 9: Update CHANGELOG**

In `CHANGELOG.md`, under `### Breaking`, append:

```markdown
- Basis constructors (`TensorBasis`, `ProjectedBasis`, `AbelianBasis`, and the
  symmetry bases) now error when `base^L` does not fit the chosen index `dtype`,
  instead of silently overflowing to wrong/negative/empty representatives
  (`bug-1`, `bug-2`, `bug-3`).
```

- [ ] **Step 10: Commit**

```bash
git add src/Basis/AbstractBasis.jl src/Basis/ProjectedBasis.jl src/Basis/AbelianBasis.jl test/remediation_tests.jl CHANGELOG.md
git commit -m "fix(basis): error on index-type overflow at construction (bug-1/2/3)"
```

---

## Task 2.2: `abelian-2` — non-default dtype no longer crashes

**Bug:** `basis(Int32; L=6, N=3, k=0)` throws a `MethodError`: the selection helpers return `Vector{Int}` representatives and the raw `Int` base, which cannot unify with the `AbelianBasis{Int32}` struct fields. Fix: convert the representative vector and base to `dtype` at the struct-construction boundary (safe because Task 2.1's guard guarantees the values fit).

**Files:**
- Modify: `src/Basis/AbelianBasis.jl` (~line 639)
- Modify: `test/remediation_tests.jl` (append)

- [ ] **Step 1: Append the failing test**

```julia
@testset "abelian-2: AbelianBasis works with a non-default index dtype" begin
    # Int64 reference (half-filling, with translation symmetry).
    B64 = basis(L=6, N=3, k=0)
    # Same basis with a narrow index type must construct, not MethodError.
    B32 = basis(Int32; L=6, N=3, k=0)
    @test size(B32, 1) == size(B64, 1)
    @test eltype(B32.I) == Int32
    # Spectra agree: build the same Hamiltonian on both and compare sorted eigenvalues.
    mat = [1.0 0 0 0; 0 -1 2 0; 0 2 -1 0; 0 0 0 1]   # Heisenberg bond (XXZ-like) on 2 sites
    E64 = trans_inv_operator(mat, 2, B64) |> Array |> Hermitian |> eigvals
    E32 = trans_inv_operator(mat, 2, B32) |> Array |> Hermitian |> eigvals
    @test E32 ≈ E64
    # A non-default dtype also works on a base>2, no-symmetry-cap small case.
    B32b = basis(Int32; L=4, base=3, k=0)
    @test size(B32b, 1) > 0
end
```

- [ ] **Step 2: Run the suite to verify the new testset FAILS**

Run: `julia --project test/remediation_tests.jl`
Expected: FAIL — `basis(Int32; L=6, N=3, k=0)` currently throws `MethodError` while constructing `AbelianBasis`.

- [ ] **Step 3: Apply the source fix**

In `src/Basis/AbelianBasis.jl`, find:

```julia
    AbelianBasis(zeros(dtype, L), I, R, G, base)
end
```

replace with:

```julia
    AbelianBasis(zeros(dtype, L), convert(Vector{dtype}, I), R, G, convert(dtype, base))
end
```

- [ ] **Step 4: Run the suite to verify it PASSES**

Run: `julia --project test/remediation_tests.jl`
Expected: PASS.

- [ ] **Step 5: Run the full suite (no regressions)**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS.

- [ ] **Step 6: Update CHANGELOG**

In `CHANGELOG.md`, under `### Fixed`, append:

```markdown
- `AbelianBasis` / `basis(...)` with a non-default index `dtype` (e.g. `Int32`)
  no longer throws a `MethodError`; representatives and base are converted to the
  requested type at construction (`abelian-2`). Note: the base-2 integer fast
  path remains bounded by 64-bit orbit machinery, so `L ≤ 62`/`64` still applies
  regardless of `dtype` (`abelian-7`; `_gosper_enumerate` already errors above
  this, and the new capacity guard reports it with a clearer message).
```

- [ ] **Step 7: Commit**

```bash
git add src/Basis/AbelianBasis.jl test/remediation_tests.jl CHANGELOG.md
git commit -m "fix(basis): AbelianBasis supports non-default index dtype (abelian-2)"
```

---

## Task 2.3: Extend the capacity guard across the symmetry bases

Apply the same guard to the remaining `dtype`-taking constructors so the integer-width sweep is complete. Each insertion is one line before the existing `base = convert(dtype, base)`. This closes `parityflip-3`'s silent-corruption path (`FlipBasis(L=63)` now errors before the `Int` `M = base^L + 1` field can wrap).

**Files (each: insert one line before `base = convert(dtype, base)`):**
- `src/Basis/FlipBasis.jl:78`
- `src/Basis/ParityBasis.jl:81`
- `src/Basis/ParityFlipBasis.jl:93`
- `src/Basis/TranslationalBasis.jl:265`
- `src/Basis/TranslationalParityBasis.jl:156`
- `src/Basis/TranslationalFlipBasis.jl:140`
- Modify: `test/remediation_tests.jl` (append)

- [ ] **Step 1: Append the failing test**

```julia
@testset "parityflip-3 / sweep: symmetry bases guard index overflow" begin
    @test_throws ErrorException FlipBasis(L=63, p=1)
    @test_throws ErrorException ParityBasis(L=63, p=1)
    @test_throws ErrorException ParityFlipBasis(L=63, p=1, z=1)
    @test_throws ErrorException TranslationalBasis(L=63, k=0)
    # Small valid cases still construct.
    @test size(FlipBasis(L=6, p=1), 1) > 0
    @test size(TranslationalBasis(L=6, k=0), 1) > 0
end
```

(`ParityFlipBasis` takes `p` and `z`; `ParityBasis`/`FlipBasis` take `p`; `TranslationalBasis` takes `k` — confirmed against the constructor signatures.)

- [ ] **Step 2: Run the suite to verify it FAILS**

Run: `julia --project test/remediation_tests.jl`
Expected: FAIL — these constructors currently overflow silently rather than throwing.

- [ ] **Step 3: Insert the guard in each symmetry constructor**

In each of the six files, find the line:

```julia
    base = convert(dtype, base)
```

and replace it with:

```julia
    _check_index_capacity(dtype, base, L)
    base = convert(dtype, base)
```

(One edit per file: `FlipBasis.jl`, `ParityBasis.jl`, `ParityFlipBasis.jl`, `TranslationalBasis.jl`, `TranslationalParityBasis.jl`, `TranslationalFlipBasis.jl`.)

- [ ] **Step 4: Run the suite to verify it PASSES**

Run: `julia --project test/remediation_tests.jl`
Expected: PASS.

- [ ] **Step 5: Run the full suite (no regressions)**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/Basis/FlipBasis.jl src/Basis/ParityBasis.jl src/Basis/ParityFlipBasis.jl src/Basis/TranslationalBasis.jl src/Basis/TranslationalParityBasis.jl src/Basis/TranslationalFlipBasis.jl test/remediation_tests.jl
git commit -m "fix(basis): extend index-capacity guard to symmetry bases (parityflip-3 sweep)"
```

---

## Self-Review

**Spec coverage (Phase 2 of the design spec):**
- `bug-1` (Int32/base>2 accumulator overflow) → Task 2.1 ✓
- `bug-2` (`_binary_fixed_weight_indices` empty basis at L=63) → Task 2.1 ✓ (guard rejects L=63 base-2)
- `bug-3` (TensorBasis size overflow) → Task 2.1 ✓
- `abelian-2` (non-default dtype crash) → Task 2.2 ✓
- `abelian-7` (Gosper/Benes cap) → subsumed: existing `_gosper_enumerate` guard + Task 2.1 guard; documented in Task 2.2 CHANGELOG ✓
- `parityflip-3` (silent corruption at L≥63) → Task 2.3 (guard closes the path); cosmetic `M::Ti` parked with rationale ✓

**Placeholder scan:** None. The one explicit "confirm the keyword" note in Task 2.3 Step 1 is a real verification instruction, not a code placeholder.

**Type consistency:** `_check_index_capacity(dtype::DataType, base::Integer, L::Integer)` is called identically everywhere with the constructor's own `dtype`/`base`/`L` (TensorBasis passes `Int64`). The guard uses `>=` against `typemax(dtype)` to leave the `+1` headroom that `FlipBasis`/`ParityFlipBasis` need for `M = base^L + 1`.

---

## Parked (with rationale)

- `parityflip-3` cosmetic field type (`M::Int` → `M::Ti` in `FlipBasis`/`ParityFlipBasis`): the capacity guard makes the silent-corruption path unreachable, so this is pure type-consistency with no functional effect under the supported `Int64` path. Changing it touches the flip index arithmetic for zero benefit; deferred to a future Int128-large-L extension (which would also need the index machinery made width-generic).
