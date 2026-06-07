# EDKit.jl Remediation — Phase 4: cached-sparse wiring

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Make `mul!` honor the `sparse!` cache so block / iterative-solver code actually gets the documented SpMM acceleration, and document the cache's snapshot semantics.

**Architecture:** The matrix `*(opt,m)` and `mul(opt,m)` paths already consult `_cached_sparse`; the matrix `mul!` methods did not. Add the same cache check to both matrix `mul!` methods and delegate to SparseArrays' in-place SpMM. The cached matrix `S = sparse(opt)` carries `eltype(opt)`, so `mul!(target, S, m, α, β)` has the same target-type requirement the matrix-free path already imposes.

**Tech Stack:** Julia ≥1.10, `Test`, `Pkg.test()`. Branch `remediation-0.6.0`.

**Findings covered:** `operator-2`, `operator-5` (Task 4.1); `operator-1` (Task 4.2, documentation). `operator-4` parked (see end).

---

## Key design facts (verified against the code)

- `Operator.jl:444-451`: the matrix `mul!` (3-arg accumulating, 5-arg standard) call `_apply_columns!` directly — never `_cached_sparse`. The vector `mul!`/`*` correctly stay matrix-free by design.
- `Operator.jl:478-511`: `mul(opt,m)` and `*(opt,m)` already do `S = _cached_sparse(opt); S !== nothing && ...`. Task 4.1 makes `mul!` consistent with them.
- `Operator.jl:354-356`: `sparse!` stores `sparse(opt)` under `objectid(opt)`. `sparse(opt)` uses `Tv = eltype(opt)` (line 247), so `S` is `Float64` for real operators and `ComplexF64` for complex — matching what the matrix-free path produces.
- The EDKit 3-arg `mul!(target, opt, m)` **accumulates** (`target += opt*m`); the SpMM equivalent is `mul!(target, S, m, true, true)`. The 5-arg follows `target = α·opt·m + β·target`, i.e. `mul!(target, S, m, α, β)`.
- `S` is a `SparseMatrixCSC`, so `mul!(target, S, m, …)` dispatches to SparseArrays' in-place SpMM, not back into EDKit's `Operator` method — no recursion.

---

## Task 4.1: `operator-2` / `operator-5` — matrix `mul!` uses the cache

**Files:**
- Modify: `src/Operator.jl:444-451`
- Modify: `test/remediation_tests.jl` (append Phase 4 section)

- [ ] **Step 1: Append the failing test**

```julia
# ---------------------------------------------------------------------------
# Phase 4 — cached-sparse wiring
# ---------------------------------------------------------------------------

@testset "operator-2: matrix mul! consults the sparse cache" begin
    L = 6
    B = basis(L=L, N=3)
    H = trans_inv_operator([1.0 0 0 0; 0 -1 2 0; 0 2 -1 0; 0 0 0 1], 2, B)
    d = size(H, 1)
    T = promote_type(eltype(H), Float64)
    X = randn(T, d, 4)
    Yref = Array(H) * X

    clear_sparse_cache!()
    # Uncached path stays correct (regression guard): 3-arg accumulates, 5-arg scales.
    Y3 = zeros(T, d, 4); mul!(Y3, H, X);          @test Y3 ≈ Yref
    Y5 = randn(T, d, 4); mul!(Y5, H, X, 2.0, 0.0); @test Y5 ≈ 2 .* Yref

    # Cached path gives identical results.
    sparse!(H)
    Yc3 = zeros(T, d, 4); mul!(Yc3, H, X);          @test Yc3 ≈ Yref
    Yc5 = randn(T, d, 4); mul!(Yc5, H, X, 2.0, 0.0); @test Yc5 ≈ 2 .* Yref

    # Prove mul! actually READS the cache: poison the cached matrix and observe
    # the result change (matrix-free would ignore the poisoned entry).
    Sorig = EDKit._cached_sparse(H)
    EDKit._SPARSE_CACHE[objectid(H)] = 3 .* Sorig
    Yp = zeros(T, d, 4); mul!(Yp, H, X);  @test Yp ≈ 3 .* Yref
    clear_sparse_cache!()
end
```

- [ ] **Step 2: Run the suite to verify it FAILS**

Run: `julia --project test/remediation_tests.jl`
Expected: FAIL — the poison assertion `Yp ≈ 3 .* Yref` fails because the current matrix `mul!` ignores the cache (it returns `Yref`, not `3·Yref`). The other assertions pass.

- [ ] **Step 3: Apply the source fix**

In `src/Operator.jl`, find:

```julia
mul!(target::AbstractMatrix, opt::Operator, m::AbstractMatrix) =
    _apply_columns!(target, opt, axes(m, 1), m)

function mul!(target::AbstractMatrix, opt::Operator, m::AbstractMatrix, α::Number, β::Number)
    iszero(β) ? fill!(target, zero(eltype(target))) : (target .*= β)
    iszero(α) && return target
    _apply_columns!(target, opt, axes(m, 1), m, α)
end
```

replace with:

```julia
function mul!(target::AbstractMatrix, opt::Operator, m::AbstractMatrix)
    S = _cached_sparse(opt)
    S !== nothing && return mul!(target, S, m, true, true)   # cached SpMM, accumulating (operator-2)
    _apply_columns!(target, opt, axes(m, 1), m)
end

function mul!(target::AbstractMatrix, opt::Operator, m::AbstractMatrix, α::Number, β::Number)
    S = _cached_sparse(opt)
    S !== nothing && return mul!(target, S, m, α, β)         # cached SpMM (operator-2)
    iszero(β) ? fill!(target, zero(eltype(target))) : (target .*= β)
    iszero(α) && return target
    _apply_columns!(target, opt, axes(m, 1), m, α)
end
```

- [ ] **Step 4: Run the suite to verify it PASSES**

Run: `julia --project test/remediation_tests.jl`
Expected: PASS — including the poison assertion (mul! now reads the cache).

- [ ] **Step 5: Run the full suite (no regressions)**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS. The `DoubleBasis`/symmetrizer and any `mul!`-based tests exercise this path.

- [ ] **Step 6: Update CHANGELOG**

In `CHANGELOG.md`, under `### Fixed`, append:

```markdown
- `mul!(target, opt::Operator, m::AbstractMatrix, …)` now uses the `sparse!`
  cache (SpMM) like `*` and `mul` already did, so block and iterative-solver code
  that calls the standard `mul!` gets the documented acceleration instead of
  silently falling back to the matrix-free path (`operator-2`, `operator-5`).
```

- [ ] **Step 7: Commit**

```bash
git add src/Operator.jl test/remediation_tests.jl CHANGELOG.md
git commit -m "fix(operator): matrix mul! consults the sparse cache (operator-2/5)"
```

---

## Task 4.2: `operator-1` — document the cache's snapshot semantics

The cache is keyed on `objectid(opt)`. Because `Operator` is immutable and every public transformation builds a fresh `Operator`, the only way to desync is mutating an operator's stored matrices in place after `sparse!` — which leaves `objectid` unchanged and returns a stale matrix. Per-lookup content hashing would tax the hot path for this internal-only case, so we document the contract instead.

**Files:**
- Modify: `src/Operator.jl` (the `sparse!` docstring, ~line 308-353)

- [ ] **Step 1: Add a snapshot-semantics note to the `sparse!` docstring**

In `src/Operator.jl`, find the line in the `sparse!` docstring:

```julia
The cache holds up to 8 operators (LRU eviction).  Call
[`clear_sparse_cache!`](@ref) when you no longer need the cached matrices
and want to reclaim memory.
```

replace with:

```julia
The cache holds up to 8 operators (LRU eviction).  Call
[`clear_sparse_cache!`](@ref) when you no longer need the cached matrices
and want to reclaim memory.

!!! note "Snapshot semantics"
    `sparse!` snapshots the operator's contents at call time and keys the cache
    by object identity. Every public transformation (`c*opt`, `opt1+opt2`,
    `adjoint`, …) returns a fresh `Operator`, so this is transparent. The one
    exception is mutating an operator's stored matrices *in place* after
    `sparse!`: that is not tracked, and subsequent `*`/`mul`/`mul!` will use the
    stale cached matrix. Call `sparse!(opt)` again to refresh it.
```

- [ ] **Step 2: Sanity-check the package still loads (docstring change only)**

Run: `julia --project -e 'using EDKit; println("ok")'`
Expected: `ok`.

- [ ] **Step 3: Update CHANGELOG**

In `CHANGELOG.md`, under `### Fixed`, append:

```markdown
- Documented `sparse!`'s snapshot semantics: the cache is keyed by object
  identity, so mutating an operator's stored matrices in place after `sparse!`
  requires a fresh `sparse!` call to refresh the cache (`operator-1`).
```

- [ ] **Step 4: Commit**

```bash
git add src/Operator.jl CHANGELOG.md
git commit -m "docs(operator): document sparse! cache snapshot semantics (operator-1)"
```

---

## Self-Review

**Spec coverage (Phase 4):**
- `operator-2` (mul! ignores cache) → Task 4.1 ✓
- `operator-5` (in-place path absent from cache) → Task 4.1 ✓ (both 3-arg and 5-arg matrix mul! covered)
- `operator-1` (objectid keying stale on mutation) → Task 4.2 (documented) ✓
- `operator-4` (threading in-place mul!) → parked with rationale ✓

**Placeholder scan:** None.

**Type consistency:** The cached `S` is `SparseMatrixCSC{eltype(opt)}`. The 3-arg delegates to `mul!(target, S, m, true, true)` (accumulate, matching the matrix-free 3-arg semantics); the 5-arg delegates to `mul!(target, S, m, α, β)`. Both impose the same `target` element-type requirement as the pre-existing matrix-free path (`target` must hold `opt*m`); the test allocates `target` with `promote_type(eltype(H), eltype(X))`.

---

## Parked (with rationale)

- `operator-4` (thread the in-place matrix `mul!`): the matrix-free `mul!` is intentionally allocation-free for repeated iterative-solver calls; threading it requires per-thread buffers + a reduction, adding an allocation per call that hurts that exact use case. Parallel throughput is already available two ways — `sparse!`+`mul!` (cached SpMM, Task 4.1) and the threaded `mul(opt, m)`. Revisit only if profiling shows the uncached threaded `mul!` is a real bottleneck for a concrete workload.
- `operator-1` content-hash cache key: rejected in favor of documentation because hashing all nonzeros on every `_cached_sparse` lookup would tax the multiply hot path for an internal-only mutation scenario.
