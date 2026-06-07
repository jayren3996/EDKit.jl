# EDKit.jl Remediation — Phase 0 + Phase 1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stand up a CI/test safety net, then fix the two critical silent-wrong-answer bugs (`abelian-1`, `productstate`) from the review, each with a regression test.

**Architecture:** Work on branch `remediation-0.6.0`. Add a CI workflow + `CHANGELOG.md`, then a dedicated regression suite `test/remediation_tests.jl` (self-runnable and wired into `runtests.jl`). Each bug fix is TDD: failing test → source fix → green test → full-suite check → commit.

**Tech Stack:** Julia ≥1.10, `Test` stdlib, `Pkg.test()`, GitHub Actions (`julia-actions/*`).

**Scope note:** This plan covers only Phase 0 + Phase 1 of the roadmap in `docs/superpowers/specs/2026-06-07-review-remediation-design.md`. Phases 2–8 (integer-width sweep, eltype sweep, cached-sparse unification, long-tail bugs, optimizations, extensions, uncertain triage) each get their own just-in-time plan, because several depend on discoveries made while implementing earlier phases.

---

## File Structure

| File | Responsibility | Action |
|---|---|---|
| `.github/workflows/CI.yml` | Run `Pkg.test()` on PRs across Julia/OS matrix | Create |
| `CHANGELOG.md` | Human-readable record of 0.6.0 changes | Create |
| `test/remediation_tests.jl` | Regression tests for every remediation finding, organized by finding id; self-runnable + included by `runtests.jl` | Create |
| `test/runtests.jl` | Master test includes | Modify (add one `include`) |
| `test/itensor_tests.jl` | ITensor/Pauli suite (already in `runtests.jl`) | Modify (adopt orphaned Pauli round-trips, `ext-1`) |
| `src/Basis/AbelianBasis.jl` | AbelianBasis constructor | Modify (`abelian-1` guard, ~line 613) |
| `src/ToolKit.jl` | `productstate` | Modify (`productstate` fix, lines 107-113) |

---

## Phase 0 — Test & CI foundation

### Task 0.1: Confirm green baseline

**Files:** none (verification only)

- [ ] **Step 1: Run the full suite on the untouched branch**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: completes with `Testing EDKit tests passed` and exit code 0. (Already verified once at plan-writing time — re-confirm before changing anything.)

If it does NOT pass, stop and triage before proceeding — the safety net must be green first.

### Task 0.2: Add CI workflow

**Files:**
- Create: `.github/workflows/CI.yml`

- [ ] **Step 1: Create the workflow**

```yaml
name: CI

on:
  push:
    branches: [main]
  pull_request:
  workflow_dispatch:

concurrency:
  group: ci-${{ github.ref }}
  cancel-in-progress: true

jobs:
  test:
    name: Julia ${{ matrix.version }} - ${{ matrix.os }}
    runs-on: ${{ matrix.os }}
    strategy:
      fail-fast: false
      matrix:
        version: ['1.10', '1']
        os: [ubuntu-latest, macOS-latest]
    steps:
      - uses: actions/checkout@v4
      - uses: julia-actions/setup-julia@v2
        with:
          version: ${{ matrix.version }}
      - uses: julia-actions/cache@v2
      - uses: julia-actions/julia-buildpkg@v1
      - uses: julia-actions/julia-runtest@v1
```

- [ ] **Step 2: Validate YAML locally**

Run: `julia --project -e 'isfile(".github/workflows/CI.yml") && println("ok")'`
Expected: `ok` (full CI runs on GitHub once pushed; local check just confirms the file exists).

- [ ] **Step 3: Commit**

```bash
git add .github/workflows/CI.yml
git commit -m "ci: run Pkg.test() on Julia 1.10+latest, Ubuntu+macOS"
```

### Task 0.3: Add CHANGELOG

**Files:**
- Create: `CHANGELOG.md`

- [ ] **Step 1: Create the changelog**

```markdown
# Changelog

All notable changes to EDKit.jl are recorded here, following
[Keep a Changelog](https://keepachangelog.com/) and (loosely) SemVer.

## [Unreleased] — targeting 0.6.0

This release is the outcome of a full section-by-section review. It includes
breaking changes; see each entry below.

### Breaking
- _none yet_

### Fixed
- _none yet_

### Performance
- _none yet_

### Added
- Continuous-integration workflow running the test suite on Julia 1.10 and
  latest across Ubuntu and macOS.
```

- [ ] **Step 2: Commit**

```bash
git add CHANGELOG.md
git commit -m "docs: add CHANGELOG for 0.6.0 remediation"
```

### Task 0.4: Adopt orphaned Pauli round-trip tests (`ext-1`)

The Pauli MPS round-trip tests in `test/TensorTest.jl:198-221` are not run by `Pkg.test()`. Move equivalent testsets into `test/itensor_tests.jl`, which is in the master suite.

**Files:**
- Modify: `test/itensor_tests.jl` (append at end)

- [ ] **Step 1: Append the two testsets to `test/itensor_tests.jl`**

```julia
#----------------------------------------------------------------------------------------------------
# Pauli MPS / MPO round-trips (adopted from TensorTest.jl, ext-1)
#----------------------------------------------------------------------------------------------------
@testset "Pauli PMPS fidelity round-trip" begin
    for L in 2:7
        s = siteinds("S=1/2", L)
        ps = siteinds("Pauli", L)
        vec = rand(ComplexF64, 2^L) |> normalize!
        ρ = vec * vec'
        pmps = vec2mps(pauli_list(ρ), ps)
        mpo = pmps2mpo(pmps, s)
        ψ = vec2mps(vec, s)
        @test inner(ψ', mpo, ψ) ≈ 1.0
    end
end

@testset "MPS -> PMPS round-trip" begin
    for L in 2:7
        s = siteinds("S=1/2", L)
        ps = siteinds("Pauli", L)
        vec = rand(ComplexF64, 2^L) |> normalize!
        ψ = vec2mps(vec, s)
        pmps = mps2pmps(ψ, ps)
        mpo = pmps2mpo(pmps, s)
        @test inner(ψ', mpo, ψ) ≈ 1.0
    end
end
```

- [ ] **Step 2: Run the full suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS. If these Pauli round-trips FAIL, do not "fix the test" — stop and treat it as a newly surfaced bug (the review classified this as missing coverage, not a known break); report and assess before continuing.

- [ ] **Step 3: Commit**

```bash
git add test/itensor_tests.jl
git commit -m "test: adopt orphaned Pauli MPS/PMPS round-trip tests into the suite (ext-1)"
```

---

## Phase 1 — Criticals

### Task 1.1: Fix `abelian-1` — fixed N + inversion off half-filling

**Bug:** `basis(L=6, N=2, z=1)` silently builds an incomplete basis (dim 10) with a wrong spectrum. A spin inversion maps the `N` sector to the `L−N` sector, so a fixed off-half-filling `N` sector has no inversion eigenstates — the request is ill-defined. Fix: error per the documented half-filling restriction.

**Files:**
- Create: `test/remediation_tests.jl`
- Modify: `test/runtests.jl` (add include)
- Modify: `src/Basis/AbelianBasis.jl` (insert guard before `# Dispatch between construction paths`, ~line 613)

- [ ] **Step 1: Create the regression suite with the failing test**

Create `test/remediation_tests.jl`:

```julia
# Regression tests for the 0.6.0 review-remediation roadmap.
# Self-runnable (`julia --project test/remediation_tests.jl`) and included by runtests.jl.
using EDKit, Test, LinearAlgebra

@testset "abelian-1: fixed N + inversion off half-filling errors" begin
    # Global spin-flip maps the N sector to the L-N sector. Off half-filling the
    # requested symmetry sector is empty/ill-defined and must error, not silently
    # return an incomplete basis.
    @test_throws ErrorException basis(L=6, N=2, z=1)
    @test_throws ErrorException basis(L=6, N=2, z=-1)
    # A custom inversion generator (flip all sites) off half-filling also errors.
    @test_throws ErrorException basis(L=6, N=2, symmetries=[(collect(1:6), 0, trues(6))])
    # Half-filling remains valid and builds a non-empty basis.
    B = basis(L=6, N=3, z=1)
    @test size(B, 1) > 0
end
```

- [ ] **Step 2: Wire the suite into `runtests.jl`**

In `test/runtests.jl`, add after the existing `include("analytic_physics_tests.jl")` line:

```julia
include("remediation_tests.jl")
```

- [ ] **Step 3: Run the new suite to verify it FAILS**

Run: `julia --project test/remediation_tests.jl`
Expected: FAIL — `basis(L=6, N=2, z=1)` currently returns a basis instead of throwing, so the first `@test_throws` fails.

- [ ] **Step 4: Apply the source fix**

In `src/Basis/AbelianBasis.jl`, the constructor currently reads (around lines 612-615):

```julia
    # Dispatch between construction paths
    Ndigits = isnothing(N) ? nothing : L * (base - 1) - N
    use_gosper = (base == 2 && Ndigits !== nothing && 0 <= Ndigits <= L)
```

Replace with:

```julia
    # Guard (abelian-1): a spin-inversion symmetry together with a fixed charge N
    # is only valid at half-filling. Off half-filling the inversion maps the N
    # sector to the L-N sector, so the requested symmetry sector is empty.
    if !isnothing(N) && any(any, G.inv) && 2 * N != L * (base - 1)
        error(
            "Fixed charge N together with a spin-inversion symmetry (e.g. `z`) is only " *
            "valid at half-filling (2N == L*(base-1)). Got N=$N, L=$L, base=$base. " *
            "Off half-filling the inversion maps the N sector to a different charge " *
            "sector, so the requested symmetry sector is empty."
        )
    end

    # Dispatch between construction paths
    Ndigits = isnothing(N) ? nothing : L * (base - 1) - N
    use_gosper = (base == 2 && Ndigits !== nothing && 0 <= Ndigits <= L)
```

- [ ] **Step 5: Run the new suite to verify it PASSES**

Run: `julia --project test/remediation_tests.jl`
Expected: PASS — the `abelian-1` testset is green.

- [ ] **Step 6: Run the full suite (no regressions)**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS. (The existing half-filling Abelian coverage uses `N = L÷2`, which the guard allows.)

- [ ] **Step 7: Update CHANGELOG**

In `CHANGELOG.md`, replace the `### Breaking` `_none yet_` line with:

```markdown
- `basis(...)` with a fixed charge `N` and a spin-inversion symmetry (`z` or a
  custom inversion generator) off half-filling now throws instead of silently
  returning an incomplete basis with a wrong spectrum (`abelian-1`).
```

- [ ] **Step 8: Commit**

```bash
git add test/remediation_tests.jl test/runtests.jl src/Basis/AbelianBasis.jl CHANGELOG.md
git commit -m "fix(basis): error on fixed N + inversion off half-filling (abelian-1)"
```

### Task 1.2: Fix `productstate` for reduced bases

**Bug:** `src/ToolKit.jl:107-113` mutates shared `B.dgt`, drops the orbit coefficient (`s[I]=1`), allocates `Float64` (cannot hold a complex phase), and on an out-of-sector config silently writes amplitude onto basis vector 1. Fix: local buffer, keep the coefficient, allocate with the basis element type (promoting integer eltype to `Float64` so the result is a usable state vector), and error on out-of-sector configs.

**Files:**
- Modify: `test/remediation_tests.jl` (append testset)
- Modify: `src/ToolKit.jl` (lines 107-113)

- [ ] **Step 1: Append the failing test to `test/remediation_tests.jl`**

```julia
@testset "productstate: correct on reduced bases, errors out of sector" begin
    L = 4
    v = [0, 1, 0, 1]

    # Reduced (Abelian) basis: must keep the orbit coefficient and place it at
    # the right component, with the basis element type.
    B = basis(; L, k=0)
    s = productstate(v, B)
    @test length(s) == size(B, 1)
    Texp = eltype(B) <: Integer ? Float64 : eltype(B)
    @test eltype(s) == Texp
    c, I = EDKit.index(B, collect(v))
    @test count(!iszero, s) == 1
    @test s[I] == c

    # Out-of-sector config must error, not silently write onto vector 1.
    Bp = basis(; L, N=2)                 # sum(dgt)==2 sector
    @test_throws ErrorException productstate([1, 1, 1, 1], Bp)  # sum==4, not in sector

    # Onsite basis still returns a Float64 unit vector (no behavior regression).
    Bt = TensorBasis(L=L, base=2)
    st = productstate(v, Bt)
    @test eltype(st) == Float64
    @test sum(st) == 1.0
    @test count(!iszero, st) == 1
end
```

- [ ] **Step 2: Run the suite to verify the new testset FAILS**

Run: `julia --project test/remediation_tests.jl`
Expected: FAIL — current `productstate` drops the coefficient, returns `Float64` for the reduced basis, and does not error out of sector.

- [ ] **Step 3: Apply the source fix**

In `src/ToolKit.jl`, replace the function body (lines 107-113):

```julia
function productstate(v::AbstractVector{<:Integer}, B::AbstractBasis)
    s = zeros(size(B, 1))
    B.dgt .= v 
    I = index(B)[2]
    s[I] = 1 
    s
end
```

with:

```julia
function productstate(v::AbstractVector{<:Integer}, B::AbstractBasis)
    length(v) == length(B.dgt) ||
        error("Configuration length $(length(v)) does not match basis length $(length(B.dgt)).")
    dgt = collect(v)                       # local buffer: thread-safe, no shared-state mutation
    c, I = index(B, dgt)
    iszero(c) &&
        error("Product configuration $(collect(Int, v)) is not contained in the sector spanned by the basis.")
    T = eltype(B)
    T <: Integer && (T = Float64)          # state vectors should be floating-point; keeps onsite output Float64
    s = zeros(T, size(B, 1))
    s[I] = c                               # keep the orbit phase / normalization coefficient
    s
end
```

- [ ] **Step 4: Run the suite to verify it PASSES**

Run: `julia --project test/remediation_tests.jl`
Expected: PASS — both `abelian-1` and `productstate` testsets green.

- [ ] **Step 5: Run the full suite (no regressions)**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS.

- [ ] **Step 6: Update CHANGELOG**

In `CHANGELOG.md`, replace the `### Fixed` `_none yet_` line with:

```markdown
- `productstate` now works correctly for symmetry-reduced bases: it uses a local
  digit buffer (no shared-state mutation), keeps the orbit phase/normalization
  coefficient, allocates with the basis element type, and errors on
  configurations outside the basis sector instead of silently writing amplitude
  onto basis vector 1 (`productstate`).
```

Also add under `### Breaking`:

```markdown
- `productstate` now errors on configurations outside the basis sector and
  returns a complex vector for symmetry-reduced bases (previously silent / always
  `Float64`).
```

- [ ] **Step 7: Commit**

```bash
git add test/remediation_tests.jl src/ToolKit.jl CHANGELOG.md
git commit -m "fix(toolkit): productstate keeps orbit coefficient and errors out of sector (productstate)"
```

---

## Self-Review

**Spec coverage (Phase 0 + Phase 1 of the design spec):**
- Phase 0 CI workflow → Task 0.2 ✓
- Phase 0 green baseline → Task 0.1 ✓
- Phase 0 CHANGELOG → Task 0.3 ✓
- Phase 0 `ext-1` orphaned-test adoption → Task 0.4 ✓
- Phase 1 `abelian-1` → Task 1.1 ✓
- Phase 1 `productstate` (+`productstate-2` local buffer) → Task 1.2 ✓
- Phases 2–8 → explicitly deferred to their own plans (scope note) ✓

**Placeholder scan:** No "TBD/TODO/handle edge cases" in steps. The `CHANGELOG.md` `_none yet_` markers are real initial document content, replaced by Tasks 1.1/1.2.

**Type consistency:** `index(B, dgt)` returns `(coeff, index)` for every basis (onsite via `AbstractBasis.jl:483`, Abelian via `:763`); `iszero(coeff)` is the uniform out-of-sector signal. `eltype(B)` integer→`Float64` promotion is applied identically in the test (`Texp`) and the fix (`T`). The guard uses `G.inv::Vector{BitVector}` with `any(any, G.inv)`, matching the field defined at `AbelianBasis.jl:204`.

---

## Notes for later phases (captured now, do not action here)

- After Phase 3 (eltype sweep) makes `AbelianBasis{k=0}` real, the `productstate` test's `Texp` already adapts (`eltype(B) <: Integer ? Float64 : eltype(B)`), so it will not need editing — it will assert `Float64` automatically.
- The `abelian-1` guard uses the half-filling condition `2N == L*(base-1)`, which is necessary for a global inversion. A *partial* inversion mask with fixed `N` is an untested edge case (even half-filling may not preserve `N`); revisit if such usage appears.
