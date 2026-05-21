# Phase 1 — Critical Bug Fixes Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:executing-plans` to execute this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the five critical/suspected bugs identified in the 2026-05-20 multi-agent code review, with TDD discipline and regression tests for each.

**Architecture:** Each task is a self-contained TDD cycle: write failing test → run → minimal fix → run → commit. All tests live in the existing `test/` flat layout; no new test files needed except for thread-safety regression (added to `entanglement_tests.jl`).

**Tech Stack:** Julia ≥ 1.9, SparseArrays, LRU caches, `Test` stdlib. Tests are run via `julia --project -e 'using Pkg; Pkg.test()'` from the repo root.

**Source review:** See [agent reports captured in conversation 2026-05-20]. Counterparts: [src/Schmidt.jl](../../../src/Schmidt.jl), [src/ToolKit.jl](../../../src/ToolKit.jl), [src/algorithms/QIM.jl](../../../src/algorithms/QIM.jl), [src/Operator.jl](../../../src/Operator.jl), [src/Basis/TranslationalFlipBasis.jl](../../../src/Basis/TranslationalFlipBasis.jl).

---

## Task 1: `meangapratio` divide-by-zero for small spectra

**Files:**
- Modify: `src/ToolKit.jl:37`
- Test: `test/core_tests.jl` (append to existing `gapratio`/`meangapratio` block at line 166-168)

**Background:** `meangapratio(E) = sum(gapratio(E)) / (length(E) - 2)`. For `length(E) ≤ 2`, denominator is `≤ 0` and `gapratio` is empty, giving `0/0 = NaN`. We make this explicit and well-defined.

- [ ] **Step 1.1: Write the failing test**

Append to `test/core_tests.jl` *inside* the existing `@testset` that contains `meangapratio(E)` (around line 166-168, after the `gapratio` assertions):

```julia
    @test isnan(meangapratio(Float64[]))
    @test isnan(meangapratio([0.0]))
    @test isnan(meangapratio([0.0, 1.0]))
    @test meangapratio([0.0, 1.0, 3.0]) ≈ 0.5
```

- [ ] **Step 1.2: Run the test to verify it fails**

```bash
cd /Users/ren/Library/CloudStorage/OneDrive-UniversityofLeeds/GitHub/EDKit.jl
julia --project -e 'using Pkg; Pkg.test(test_args=["--core"])' 2>&1 | tail -30
```

Or run the file directly:
```bash
julia --project test/core_tests.jl 2>&1 | tail -20
```

Expected: tests for `Float64[]`, `[0.0]`, `[0.0, 1.0]` FAIL (DomainError or `NaN` mismatch / `BoundsError`).

- [ ] **Step 1.3: Fix `meangapratio`**

Edit `src/ToolKit.jl:37`:

```julia
function meangapratio(E::AbstractVector{<:Real})
    r = gapratio(E)
    isempty(r) ? NaN : sum(r) / length(r)
end
```

Also update `gapratio` at `src/ToolKit.jl:15-28` to handle short inputs cleanly: change line 17 from `r = zeros(length(dE)-1)` to handle `length(dE) ≤ 1` by returning `Float64[]`. New `gapratio`:

```julia
function gapratio(E::AbstractVector{<:Real})
    dE = diff(E)
    length(dE) < 2 && return Float64[]
    r = zeros(length(dE) - 1)
    for i = 1:length(r)
        if dE[i] < dE[i+1]
            r[i] = dE[i]/dE[i+1]
        elseif dE[i] > dE[i+1]
            r[i] = dE[i+1]/dE[i]
        else
            r[i] = 1.0
        end
    end
    r
end
```

- [ ] **Step 1.4: Run the test to verify it passes**

```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -20
```

Expected: all tests PASS, including the new edge cases.

- [ ] **Step 1.5: Commit**

```bash
git add src/ToolKit.jl test/core_tests.jl
git commit -m "$(cat <<'EOF'
Fix meangapratio NaN for spectra with fewer than 3 levels

Guard gapratio against empty diff and make meangapratio return NaN
rather than 0/0 when there are fewer than 3 ratios to average.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: `qimsolve`/`simplify` null-pivot crash

**Files:**
- Modify: `src/algorithms/QIM.jl:65-82`
- Test: `test/advanced_tests.jl` (append to existing `@testset` at line 89)

**Background:** `simplify(h)` calls `findfirst(x -> abs(x) > 1e-7, h[:,i])`, which can return `nothing` when the elimination produces an all-zero column. The next line then indexes with `nothing` and throws. Fix: skip columns whose pivot is missing.

- [ ] **Step 2.1: Write the failing test**

Append to `test/advanced_tests.jl` *inside* the existing `@testset` for qimsolve at lines 89-93:

```julia
    @testset "simplify handles columns that eliminate to zero" begin
        # Two linearly dependent columns: elimination zeroes the second one.
        h = Float64[1.0 1.0;
                    0.0 0.0]
        result = EDKit.simplify(copy(h))
        @test result[:, 1] ≈ [1.0, 0.0]
        @test all(iszero, result[:, 2])

        # qimsolve with linearly dependent ops should not crash.
        dep_ops = [Float64[1 0; 0 0], Float64[1 0; 0 0]]
        dep_sol = qimsolve(dep_ops, [1.0, 0.0])
        @test size(dep_sol, 1) == 2
    end
```

- [ ] **Step 2.2: Run the test to verify it fails**

```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -40
```

Expected: ArgumentError or MethodError on indexing with `nothing` inside `simplify`.

- [ ] **Step 2.3: Fix `simplify`**

Edit `src/algorithms/QIM.jl:65-82`. Replace the function with:

```julia
function simplify(h)
    n = size(h, 2)
    for i in 1:n
        j = findfirst(x -> abs(x) > 1e-7, h[:,i])
        isnothing(j) && continue
        for k in 1:n
            isequal(k, i) && continue
            h[:,k] .-= h[j,k] / h[j,i] * h[:,i]
        end
    end
    for i in 1:n
        j = findfirst(x -> abs(x) > 1e-7, h[:,i])
        isnothing(j) && continue
        h[:,i] ./= h[j,i]
    end
    for i in eachindex(h)
        abs(h[i]) < 1e-10 && (h[i] = 0.0)
    end
    h
end
```

- [ ] **Step 2.4: Run the test to verify it passes**

```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -30
```

Expected: all tests PASS.

- [ ] **Step 2.5: Commit**

```bash
git add src/algorithms/QIM.jl test/advanced_tests.jl
git commit -m "$(cat <<'EOF'
Guard qimsolve/simplify against null pivots

simplify now skips columns where elimination has zeroed the pivot
position, preventing crashes on linearly dependent operator lists.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Strongly type `SPIN_CACHE`

**Files:**
- Modify: `src/Operator.jl:1055` (cache declaration) and `:1070-1079` (the `spin` function)
- Test: `test/core_tests.jl` (add type-stability check)

**Background:** `LRU{Tuple{Int,String}, Any}` forces dynamic dispatch on every cache hit, which slows down hot-path operator construction. Pin to `SparseMatrixCSC{ComplexF64, Int}` and convert at insertion time. This is safe because (a) spin matrices are always sparse, and (b) downstream code already promotes to complex when combining real and complex spin operators.

- [ ] **Step 3.1: Write the type-stability test**

Append to `test/core_tests.jl` *inside the same `@testset` as the existing gapratio tests* (after the new Task 1 assertions):

```julia
    # SPIN_CACHE is concretely typed so cache hits are type-stable.
    @test valtype(EDKit.SPIN_CACHE) === SparseMatrixCSC{ComplexF64, Int}
    @test spin("X") isa SparseMatrixCSC{ComplexF64, Int}
    @test spin("xx") isa SparseMatrixCSC{ComplexF64, Int}
    # Values are unchanged numerically.
    @test Matrix(spin("X")) ≈ ComplexF64[0 1; 1 0]
    @test Matrix(spin("Y")) ≈ ComplexF64[0 -im; im 0]
    @test Matrix(spin("Z")) ≈ ComplexF64[1 0; 0 -1]
```

Add `using SparseArrays` at the top of `test/core_tests.jl` if not already imported. Check first:
```bash
grep -n "using SparseArrays\|^using" /Users/ren/Library/CloudStorage/OneDrive-UniversityofLeeds/GitHub/EDKit.jl/test/core_tests.jl | head -5
```
If missing, add `using SparseArrays` after the existing `using EDKit`-style imports. (Note: `runtests.jl` already imports `SparseArrays`, so this may already be available transitively — verify by running.)

- [ ] **Step 3.2: Run the test to verify it fails**

```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -30
```

Expected: `valtype` test FAILS (returns `Any`), `isa SparseMatrixCSC{ComplexF64,Int}` tests FAIL for `spin("X")` (it's `SparseMatrixCSC{Int,Int}`).

- [ ] **Step 3.3: Fix `SPIN_CACHE` typing**

Edit `src/Operator.jl:1055`. Change:

```julia
const SPIN_CACHE = LRU{Tuple{Int, String}, Any}(maxsize=256)
```

to:

```julia
const SPIN_CACHE = LRU{Tuple{Int, String}, SparseMatrixCSC{ComplexF64, Int}}(maxsize=256)
```

Edit `src/Operator.jl:1070-1079`. Change the `spin(s::String; D::Integer=2)` body to convert at insertion:

```julia
function spin(s::String; D::Integer=2)
    key = (Int(D), s)
    mat = get!(SPIN_CACHE, key) do
        ny = spin_ny(s, D)
        raw = spin_product(s, D)
        sign = iszero(mod(ny, 2)) ? (-1)^(ny÷2) : (-1im)^ny
        convert(SparseMatrixCSC{ComplexF64, Int}, sign * raw)
    end
    copy(mat)
end
```

- [ ] **Step 3.4: Run the full test suite to verify nothing else breaks**

```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -40
```

Expected: all tests PASS. If any downstream test fails because it expected `spin("X")` to be `Int`-typed (e.g. `eltype(spin("X")) === Int`), document and update.

- [ ] **Step 3.5: Commit**

```bash
git add src/Operator.jl test/core_tests.jl
git commit -m "$(cat <<'EOF'
Strongly type SPIN_CACHE to SparseMatrixCSC{ComplexF64,Int}

Replace the Any value type with a concrete sparse-complex type so
cached spin operators have a stable inferred type, eliminating
dynamic dispatch in the hot operator-construction path.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: `Schmidt.addto!` thread-safety

**Files:**
- Modify: `src/Schmidt.jl:16-22` (struct), `:39-57` (constructor), `:65-71` (`addto!`)
- Test: `test/entanglement_tests.jl` (add threaded regression test)

**Background:** `addto!` writes into `S.B1.dgt` and `S.B2.dgt` — fields of the subsystem-basis objects, which are *shared* across concurrent `schmidt` calls. Two threads calling `schmidt` simultaneously corrupt each other's index lookups. The fix: give `SchmidtMatrix` its own per-instance digit buffers so each Schmidt assembly is self-contained.

- [ ] **Step 4.1: Write the failing threaded regression test**

Append to `test/entanglement_tests.jl` (end of file or a new `@testset`):

```julia
@testset "Schmidt thread safety" begin
    L = 8
    B = TensorBasis(L=L, base=2)
    Ainds = collect(1:L÷2)

    # Make a batch of distinct random states.
    nstates = 32
    rng = Random.MersenneTwister(42)
    states = [normalize!(randn(rng, ComplexF64, size(B, 1))) for _ in 1:nstates]

    # Reference: serial schmidt SVD spectra.
    serial_spectra = [svdvals(EDKit.schmidt(v, Ainds, B)) for v in states]

    # Parallel: each thread independently computes schmidt on the same B.
    threaded_spectra = Vector{Vector{Float64}}(undef, nstates)
    Threads.@threads for i in 1:nstates
        threaded_spectra[i] = svdvals(EDKit.schmidt(states[i], Ainds, B))
    end

    for i in 1:nstates
        @test threaded_spectra[i] ≈ serial_spectra[i]
    end
end
```

- [ ] **Step 4.2: Run the threaded test to verify it fails**

```bash
julia --project -t 8 -e 'using Pkg; Pkg.test()' 2>&1 | tail -30
```

Expected: failures (mismatched spectra) when run with multiple threads. May be non-deterministic; if it passes by luck on the first run, run several times:

```bash
for i in 1 2 3 4 5; do
  julia --project -t 8 -e 'using Pkg; Pkg.test()' 2>&1 | grep -E "FAIL|PASS|Test Failed" | head -5
done
```

- [ ] **Step 4.3: Refactor `SchmidtMatrix` to own its own buffers**

Edit `src/Schmidt.jl:16-22`. Replace the struct definition with:

```julia
struct SchmidtMatrix{Tm <: Number, Ta <: SubArray, Tb <: SubArray, TB1 <: AbstractBasis, TB2 <: AbstractBasis, Td1 <: AbstractVector, Td2 <: AbstractVector}
    M::Matrix{Tm}
    A::Ta
    B::Tb
    B1::TB1
    B2::TB2
    dgt1::Td1
    dgt2::Td2
end
```

Edit `src/Schmidt.jl:39-57`. Replace `schmidtmatrix` with:

```julia
function schmidtmatrix(
    T::DataType, b::AbstractBasis, Ainds::AbstractVector{Ta},
    B1=nothing, B2=nothing;
    dgt::AbstractVector=b.dgt
) where Ta <: Integer
    L = length(b)
    Binds = Vector{Ta}(undef, L-length(Ainds))
    P = 1
    for i in range(one(Ta), stop=convert(Ta, L))
        if !in(i, Ainds)
            Binds[P] = i
            P += 1
        end
    end
    B1 = isnothing(B1) ? TensorBasis(L=length(Ainds), base=b.B) : B1
    B2 = isnothing(B2) ? TensorBasis(L=length(Binds), base=b.B) : B2
    M = zeros(T, size(B1, 1), size(B2, 1))
    dgt1 = similar(B1.dgt)
    dgt2 = similar(B2.dgt)
    SchmidtMatrix(M, view(dgt, Ainds), view(dgt, Binds), B1, B2, dgt1, dgt2)
end
```

Edit `src/Schmidt.jl:65-71`. Replace `addto!`:

```julia
function addto!(S::SchmidtMatrix, val::Number)
    S.dgt1 .= S.A
    S.dgt2 .= S.B
    _, ia = index(S.B1, S.dgt1)
    _, ib = index(S.B2, S.dgt2)
    S.M[ia, ib] += val
end
```

- [ ] **Step 4.4: Run the full test suite**

```bash
julia --project -t 8 -e 'using Pkg; Pkg.test()' 2>&1 | tail -40
```

Expected: all tests PASS, including the new threaded test.

- [ ] **Step 4.5: Verify no other call sites need updating**

```bash
grep -rn "SchmidtMatrix(" /Users/ren/Library/CloudStorage/OneDrive-UniversityofLeeds/GitHub/EDKit.jl/src/ /Users/ren/Library/CloudStorage/OneDrive-UniversityofLeeds/GitHub/EDKit.jl/test/
```

If any external constructor calls exist, they must pass `dgt1` and `dgt2` too. Expected: only `schmidtmatrix(...)` constructs `SchmidtMatrix`, so this should be empty besides the definition.

- [ ] **Step 4.6: Commit**

```bash
git add src/Schmidt.jl test/entanglement_tests.jl
git commit -m "$(cat <<'EOF'
Make SchmidtMatrix thread-safe with per-instance buffers

addto! previously wrote into the shared subsystem-basis dgt fields,
which corrupts concurrent schmidt() calls. Each SchmidtMatrix now
owns its own dgt1/dgt2 buffers, matching the buffer-passing contract
documented in CLAUDE.md.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: `TranslationFlipBasis` — input validation for odd `L`

**Files:**
- Modify: `src/Basis/TranslationalFlipBasis.jl:122-158` (constructor — add assertion)
- Test: `test/basis_tests.jl` (add input-validation test)

**Background:** The parity-check at line 79 (`iseven(div(n * judge.K, judge.L÷2))`) silently truncates `L÷2` for odd `L`. Without a containing assertion, callers with odd `L` get incorrect basis filtering rather than a clear error. We add an explicit precondition and document the constraint. The deeper "is the convention correct?" question goes to Phase 2 — for now we just refuse inputs we don't trust.

- [ ] **Step 5.1: Write the failing test**

Append to `test/basis_tests.jl`. First check the file structure:

```bash
grep -n "@testset\|TranslationFlipBasis\|^using\|^end" /Users/ren/Library/CloudStorage/OneDrive-UniversityofLeeds/GitHub/EDKit.jl/test/basis_tests.jl | head -30
```

Add a new `@testset` at the end of the file:

```julia
@testset "TranslationFlipBasis input validation" begin
    # Even L is supported.
    @test_nowarn TranslationFlipBasis(L=4, k=0, p=1)
    @test_nowarn TranslationFlipBasis(L=6, k=0, p=1)

    # Odd L should error (the internal parity check uses L÷2 which truncates).
    @test_throws AssertionError TranslationFlipBasis(L=5, k=0, p=1)
    @test_throws AssertionError TranslationFlipBasis(L=7, k=0, p=1)
end
```

- [ ] **Step 5.2: Run the test to verify it fails**

```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -30
```

Expected: `@test_throws AssertionError` for `L=5` FAILS (the constructor returns a basis rather than throwing).

- [ ] **Step 5.3: Add the precondition assertion**

Edit `src/Basis/TranslationalFlipBasis.jl:122-129`. Insert a new `@assert` after the existing `@assert isnothing(N) || isequal(...)` line (currently line 129):

```julia
function TranslationFlipBasis(
    dtype::DataType=Int64; f=nothing, k::Integer=0, p::Integer=1, L::Integer, N::Union{Nothing, Integer}=nothing,
    a::Integer=1, base::Integer=2, alloc::Integer=1000, threaded::Bool=false, small_N::Bool=false
)
    len, check_a = divrem(L, a)
    @assert iszero(check_a) "Length of unit-cell $a incompatible with L=$L"
    @assert isone(p) || isone(-p) "Invalid parity"
    @assert isnothing(N) || isequal(2N, L*(base-1)) "N = $N not compatible."
    @assert iseven(L) "TranslationFlipBasis requires even L (got L=$L); the flip-parity check at TranslationalFlipBasis.jl:79 truncates for odd L."
    ...
```

- [ ] **Step 5.4: Run the test to verify it passes**

```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -20
```

Expected: PASS.

- [ ] **Step 5.5: Update the docstring**

Edit `src/Basis/TranslationalFlipBasis.jl:99-121`. Update the `Notes:` section in the docstring to mention the even-`L` requirement:

```julia
"""
TranslationFlipBasis(f, k, p, L; base=2, alloc=1000, threaded=true)

Construct a translation-plus-spin-flip basis.

Arguments:
- `f`       : Selection function for the basis contents.
- `k`       : Momentum number from 0 to L-1.
- `p`       : Eigenvalue under spin flip, `+1` or `-1`.
- `L`       : Length of the system. Must be even.
- `base`    : Base, default = 2.
- `alloc`   : Size of the prealloc memory for the basis content, used only in multithreading, default = 1000.
- `threaded`: Whether use the multithreading, default = true.

Outputs:
--------
- `b`: TranslationFlipBasis.

Notes:
- If `N` is provided, it must be compatible with spin flip.
- `L` must be even — the flip-parity selector uses `L÷2` integer
  division and is not correct for odd `L`.
- The internal momentum convention follows the same sign choice as
  [`TranslationalBasis`](@ref).
"""
```

- [ ] **Step 5.6: Commit**

```bash
git add src/Basis/TranslationalFlipBasis.jl test/basis_tests.jl
git commit -m "$(cat <<'EOF'
Require even L in TranslationFlipBasis

The flip-parity selector uses div(n*K, L÷2), which silently
truncates and gives incorrect sector filtering when L is odd.
Refuse odd L at construction with a clear assertion until the
deeper convention question is resolved in Phase 2.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Self-review checklist

After all tasks complete:

- [ ] Run full suite single-threaded: `julia --project -e 'using Pkg; Pkg.test()'`
- [ ] Run full suite multi-threaded (8 threads): `julia --project -t 8 -e 'using Pkg; Pkg.test()'`
- [ ] Confirm `git log --oneline -5` shows five new commits, one per task.
- [ ] No leftover `@assert` patterns that depend on `--check-bounds=yes` (out of scope for Phase 1 but flag if encountered).
- [ ] No introduced backwards-incompatibility beyond `TranslationFlipBasis(L=odd)` (which was always broken).

## Out of scope for Phase 1

- `_SPARSE_CACHE` typing (still `Any`) — needs design discussion since the eltype varies per operator. Defer to Phase 2 alongside the broader Operator.jl refactor.
- The deeper `TranslationFlipBasis` correctness question (is `div(n*K, L÷2)` the right expression even for even `L`?) — needs cross-check vs `TranslationalBasis` half-spectra; deferred to Phase 2.
- Refactor of the 8 near-identical `mul!`/`mul`/`*` paths in Operator.jl — Phase 2.
- `apply_perm!` Schmidt allocation removal — Phase 3 (performance).
