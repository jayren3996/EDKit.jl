# EDKit.jl Remediation — Phase 3: eltype-correctness sweep

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Give the real-phase bases a real `Float64` element type so `index` is type-stable and real Hamiltonians on the common (k=0 / parity / flip) sectors assemble as `Float64` (half the memory and FLOPs) instead of `ComplexF64`.

**Architecture:** `ParityBasis`/`FlipBasis`/`ParityFlipBasis` carry only ±1 real phases → unconditional `eltype = Float64`. `AbelianBasis` keys `eltype` on its phase-table type parameter `Tg`: `Tg<:Real ? Float64 : ComplexF64`, mirroring `TranslationalBasis.jl:31`.

**Tech Stack:** Julia ≥1.10, `Test`, `Pkg.test()`. Branch `remediation-0.6.0`.

**Findings covered:** `parityflip-1`, `abelian-5`.

---

## Key design facts (verified against the code & by probing)

- The discrete bases' `index` returns `n * R[i]` (real: `n∈{±1}` `Int`, `R` `Float64`) in the hit branch but `zero(eltype(b))` in the miss branch. With the inherited `eltype == ComplexF64`, the return type is `Union{Tuple{Float64,Int}, Tuple{ComplexF64,Int}}` — type-unstable. Setting `eltype = Float64` makes the miss branch `0.0` and the method type-stable. (`ParityFlipBasis.index` already returns `Float64` in both branches; only its `eltype` lies.)
- For `AbelianBasis`, the phase table `c` is `ones` (`Float64`) for k=0, `[1,-1]` (`Int`) for the ±1 character (parity `p=-1`, flip `z=-1`, and `k=L/2`), and `ComplexF64` only for genuine momentum. Combined generators promote: any real-only combination yields `Tg<:Real` (`Int` or `Float64`); any genuine-momentum factor yields `Tg==ComplexF64`. Probed across k/p/z/combinations: `Tg<:Real` ⇔ all characters real ⇔ `phase(g)` real ⇔ the hit-branch coefficient is real. So `Tg<:Real ? Float64 : ComplexF64` is exact.
- The probe's apparent `ComplexF64` coefficients for `p=-1` etc. came from the *miss* branch (`zero(ComplexF64)`) on a state outside the antisymmetric sector — fixed by the same `eltype` change.
- A complex operator on a now-`Float64` basis still works: `promote_type(ComplexF64, Float64) == ComplexF64`. Only real-operator/real-sector combinations change (to `Float64`), which is the intended win.

---

## Task 3.1: `parityflip-1` — real eltype for Parity/Flip/ParityFlip

**Files:**
- Modify: `src/Basis/ParityBasis.jl` (after `order`, line 115)
- Modify: `src/Basis/FlipBasis.jl` (after `order`, line 110)
- Modify: `src/Basis/ParityFlipBasis.jl` (after `order`, line 140; and the `0.0` literal at line 122)
- Modify: `test/remediation_tests.jl` (append Phase 3 section)

- [ ] **Step 1: Append the failing test**

```julia
# ---------------------------------------------------------------------------
# Phase 3 — eltype-correctness sweep
# ---------------------------------------------------------------------------

@testset "parityflip-1: Parity/Flip/ParityFlip have real eltype and type-stable index" begin
    L = 6
    @test eltype(ParityBasis(L=L, p=1)) == Float64
    @test eltype(FlipBasis(L=L, p=1)) == Float64
    @test eltype(ParityFlipBasis(L=L, p=1, z=1)) == Float64
    # index must be type-stable: returns a concrete Float64 tuple, not a Union with ComplexF64.
    for B in (ParityBasis(L=L, p=1), FlipBasis(L=L, p=1), ParityFlipBasis(L=L, p=1, z=1))
        r = @inferred EDKit.index(B, collect(zeros(Int, L)))
        @test r isa Tuple{Float64, <:Integer}
    end
end
```

- [ ] **Step 2: Run the suite to verify it FAILS**

Run: `julia --project test/remediation_tests.jl`
Expected: FAIL — `eltype` is currently `ComplexF64`, and `@inferred` throws on the unstable `index` return type.

- [ ] **Step 3: Add `eltype = Float64` to ParityBasis**

In `src/Basis/ParityBasis.jl`, find:

```julia
order(b::ParityBasis) = 2
```

replace with:

```julia
order(b::ParityBasis) = 2
eltype(::ParityBasis) = Float64
```

- [ ] **Step 4: Add `eltype = Float64` to FlipBasis**

In `src/Basis/FlipBasis.jl`, find:

```julia
order(b::FlipBasis) = 2
```

replace with:

```julia
order(b::FlipBasis) = 2
eltype(::FlipBasis) = Float64
```

- [ ] **Step 5: Add `eltype = Float64` to ParityFlipBasis and robustify the miss-branch literal**

In `src/Basis/ParityFlipBasis.jl`, find:

```julia
order(b::ParityFlipBasis) = 4
```

replace with:

```julia
order(b::ParityFlipBasis) = 4
eltype(::ParityFlipBasis) = Float64
```

Then, in the same file, find:

```julia
    iszero(i) && return (0.0, one(b.B))
```

replace with:

```julia
    iszero(i) && return (zero(eltype(b)), one(b.B))
```

- [ ] **Step 6: Run the suite to verify it PASSES**

Run: `julia --project test/remediation_tests.jl`
Expected: PASS.

- [ ] **Step 7: Run the full suite (no regressions)**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS.

- [ ] **Step 8: Update CHANGELOG**

In `CHANGELOG.md`, under `### Performance`, append:

```markdown
- `ParityBasis`, `FlipBasis`, `ParityFlipBasis` now have a real (`Float64`)
  element type, making their `index` type-stable and letting real Hamiltonians
  assemble as real matrices instead of `ComplexF64` (`parityflip-1`).
```

- [ ] **Step 9: Commit**

```bash
git add src/Basis/ParityBasis.jl src/Basis/FlipBasis.jl src/Basis/ParityFlipBasis.jl test/remediation_tests.jl CHANGELOG.md
git commit -m "perf(basis): real Float64 eltype for parity/flip bases (parityflip-1)"
```

---

## Task 3.2: `abelian-5` — key AbelianBasis eltype on its phase type

**Files:**
- Modify: `src/Basis/AbelianBasis.jl` (after `order`, line 591)
- Modify: `test/remediation_tests.jl` (append)

- [ ] **Step 1: Append the failing test**

```julia
@testset "abelian-5: real AbelianBasis sectors are Float64 and preserve spectrum" begin
    L = 6
    # Real-character sectors -> Float64; genuine momentum stays ComplexF64.
    @test eltype(basis(L=L, k=0)) == Float64
    @test eltype(basis(L=L, p=1)) == Float64
    @test eltype(basis(L=L, p=-1)) == Float64
    @test eltype(basis(L=L, z=-1)) == Float64
    @test eltype(basis(L=L, k=0, p=-1)) == Float64
    @test eltype(basis(L=L, k=1)) == ComplexF64

    # A real, reflection-symmetric Hamiltonian assembles real on a real sector,
    # and the parity sectors still partition the full spectrum (no value corruption).
    hmat = [1.0 0 0 0; 0 -1 2 0; 0 2 -1 0; 0 0 0 1]   # swap-symmetric real 2-site term
    Bp = basis(L=L, p=1)
    Hp = trans_inv_operator(hmat, 2, Bp) |> Array
    @test eltype(Hp) == Float64

    Efull = trans_inv_operator(hmat, 2, TensorBasis(L=L)) |> Array |> Hermitian |> eigvals
    Epar = Float64[]
    for s in (1, -1)
        append!(Epar, trans_inv_operator(hmat, 2, basis(L=L, p=s)) |> Array |> Hermitian |> eigvals)
    end
    @test sort(Epar) ≈ sort(Efull)
end
```

- [ ] **Step 2: Run the suite to verify it FAILS**

Run: `julia --project test/remediation_tests.jl`
Expected: FAIL — `eltype(basis(L=6, k=0)) == Float64` fails (currently `ComplexF64`), and `eltype(Hp) == Float64` fails (assembles complex).

- [ ] **Step 3: Add the type-keyed eltype**

In `src/Basis/AbelianBasis.jl`, find:

```julia
order(b::AbelianBasis) = order(b.G)
```

replace with:

```julia
order(b::AbelianBasis) = order(b.G)
eltype(::AbelianBasis{Ti, Tg}) where {Ti, Tg} = Tg <: Real ? Float64 : ComplexF64
```

- [ ] **Step 4: Run the suite to verify it PASSES**

Run: `julia --project test/remediation_tests.jl`
Expected: PASS.

- [ ] **Step 5: Run the full suite (no regressions)**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS. Pay attention to the Abelian, entanglement, and time-evolution suites, which exercise these bases most.

- [ ] **Step 6: Update CHANGELOG**

In `CHANGELOG.md`, under `### Performance`, append:

```markdown
- `AbelianBasis` now reports a real (`Float64`) element type for real-character
  sectors (k=0, parity, spin-flip, and k=L/2), keyed on its phase-table type
  parameter, so the most common symmetry workflows assemble and diagonalize real
  matrices at half the memory/FLOPs; genuine momentum sectors remain
  `ComplexF64` (`abelian-5`).
```

- [ ] **Step 7: Commit**

```bash
git add src/Basis/AbelianBasis.jl test/remediation_tests.jl CHANGELOG.md
git commit -m "perf(basis): real eltype for real-character AbelianBasis sectors (abelian-5)"
```

---

## Self-Review

**Spec coverage (Phase 3):**
- `parityflip-1` (Parity/Flip/ParityFlip ComplexF64 eltype) → Task 3.1 ✓ (incl. the `ParityFlipBasis` `0.0`→`zero(eltype(b))` robustification)
- `abelian-5` (AbelianBasis always ComplexF64) → Task 3.2 ✓

**Placeholder scan:** None.

**Type consistency:** The `AbelianBasis` `eltype` ternary on `Tg<:Real` const-folds per concrete type (so `index`'s `zero(eltype(B))` is `Float64` for real sectors, keeping `index` type-stable). The discrete bases get unconditional `Float64`, matching their always-real `index` arithmetic. Tests assert `Tuple{Float64,<:Integer}` from `@inferred`, consistent with the fix.

---

## Risk note

If the full suite surfaces a failure where a real-sector computation actually needed complex storage (not expected, since complex operators promote), do not paper over it by reverting `eltype` — investigate which call assumed `ComplexF64` and fix that assumption. Stop and report if it is non-trivial.
