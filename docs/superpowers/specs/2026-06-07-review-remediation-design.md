# EDKit.jl Review-Remediation Roadmap — Design Spec

**Date:** 2026-06-07
**Branch:** `remediation-0.6.0`
**Source of findings:** `REVIEW_REPORT.md` (multi-agent section review, 66 findings → 50 confirmed, 7 uncertain, 9 refuted)
**Target version:** 0.6.0 (breaking allowed)

## 1. Goal & scope

Address every confirmed correctness bug from the review, do the high-leverage optimizations,
and add the highest-value extensions — as one sequenced roadmap. The review found the defects
cluster into a few recurring families rather than scattered one-offs, so the plan is organized
to fix each *root cause* once and build later work on clean foundations.

In scope: all confirmed bugs, the confirmed optimizations, the prioritized extensions, a
test/CI safety net, and triage of the 7 uncertain findings.

Out of scope: the 9 refuted findings (excluded), unrelated refactoring, and the verified-immaterial
micro-notes (parked with a one-line rationale, see §8).

## 2. Decisions & constraints

| Decision | Choice | Consequence |
|---|---|---|
| Backward-compat posture | **Breaking OK → bump to 0.6.0** | Fix correctly, no deprecation shims. Document every observable change in `CHANGELOG.md`. The behaviors being "broken" are overwhelmingly silent-wrong-answers no correct code should depend on. |
| Plan breadth | **Everything, phased roadmap** | Bugs + optimizations + extensions + test/CI infra in one sequenced spec. |
| Execution | **Plan, then start criticals** | After this spec + the implementation plan are approved, begin Phase 1 with tests, check in before continuing. |

Behavior changes that downstream users must be warned about (CHANGELOG-worthy):
- Parity/Flip/ParityFlip/Abelian real sectors now produce **real** (`Float64`) matrices instead of `ComplexF64`.
- `order()` for translational bases now returns the orbit size `L/a` (was `L`), changing symmetrizer/embedding/MPS normalization for `a>1` unit cells.
- Constructors now **error** (instead of silently overflowing) when `base^L − 1 > typemax(dtype)`.
- `productstate`, parity/flip predicate misuse, and several other previously-silent wrong-answer paths now **error** or return correct results.

## 3. Organizing approach

Chosen approach **C (hybrid)**: criticals first (isolated, highest-risk), then systemic theme
sweeps (max leverage — each root cause fixed once), then the long-tail individual bugs, then
optimizations, then extensions — all on top of a test/CI foundation laid first.

Rejected: **A (strict severity-first)** thrashes the same files repeatedly; **B (pure theme-first)**
buries the two criticals inside larger sweeps.

Dependency spine: **Phase 0 (tests)** gates everything (safety net for breaking changes).
**Phases 2–3 (integer-width + eltype)** clean the basis numeric-type core, which **Phase 7's
mixed-radix extension** then builds on. **Phase 4** unifies the operator-apply path that several
optimizations touch.

## 4. Phased roadmap

Each phase ends green (`Pkg.test()` passes) and appends to `CHANGELOG.md`. Every fix ships with
a regression test that fails before and passes after.

### Phase 0 — Test & CI foundation *(precondition)*

**Goal:** a safety net before any breaking change.

- Add `.github/workflows/CI.yml`: `Pkg.test()` on Julia 1.10 + latest, Ubuntu + macOS.
- Establish a green baseline: run the full suite, fix or quarantine anything already red (record what was quarantined and why).
- Adopt orphaned coverage: the Pauli round-trip in `test/TensorTest.jl:198` (report `ext-1`); fold relevant `MultiThreadsTest`/`DoubleBasisTest` assertions into `runtests.jl`.
- Create `CHANGELOG.md` with an `[Unreleased] → 0.6.0` section; every later phase appends.

**Exit:** CI runs the suite on PRs; baseline is green; CHANGELOG exists.

### Phase 1 — Criticals *(silent wrong answers on documented APIs)*

| Item | Location | Fix | Regression test |
|---|---|---|---|
| `abelian-1` | `src/Basis/AbelianBasis.jl:614-625` | Compute `has_inv=any(any,G.inv)`; gate `use_gosper` on `(!has_inv \|\| 2*Ndigits==L*(base-1))`. When the gate fails, **fall back to the full-scan path** (preserves the feature off half-filling rather than erroring). | `basis(L=6,N=2,z=1)` → dim 15 and spectrum matches the unsymmetrized reference. |
| `productstate` | `src/ToolKit.jl:107-113` | `c,I=index(B,collect(v))`; `iszero(c)&&error(...)`; `s=zeros(eltype(B),size(B,1))`; `s[I]=c`. Uses a local buffer (also resolves `productstate-2`). | Round-trip a product state through a reduced basis; out-of-sector input errors. |

**Exit:** both criticals fixed with tests; CHANGELOG updated.

### Phase 2 — Integer-width sweep

**Goal:** decouple integer width from the basis everywhere it silently overflows; convert overflow into a loud error.

| Item | Location | Fix |
|---|---|---|
| `bug-1` | `AbstractBasis.jl:297-310`; stored `ProjectedBasis.jl:202-203,242` | Decouple the `index` accumulator width from the *base* type. |
| `bug-2` | `ProjectedBasis.jl:36-54` (41,44,47) | Fix `_binary_fixed_weight_indices` `one(T)<<L` (returns empty basis at L=63). |
| `bug-3` | `AbstractBasis.jl:255-256` | `TensorBasis` size overflow guard. |
| `abelian-2` | `AbelianBasis.jl:648,664,681,698,705` (+ `:639`, `:724`) | Thread `dtype` through `Is`/`dgt`/`change!`, the `base` arg, and the `_abelian_select` fallback range. |
| `abelian-7` | `AbelianBasis.jl:554-569` | Decide Gosper L≤62 vs Benes L≤64 ceiling; align. |
| `parityflip-3` (triage) | `FlipBasis.jl:20`, `ParityFlipBasis.jl:21` | Make `M::Ti` for consistency (cosmetic; unreachable 2^63). |

**Cross-cutting:** one loud constructor-time guard `base^L − 1 ≤ typemax(dtype)` covers the whole family.

**Tests:** `ProjectedBasis(Int32;L=20,N=0,base=3)`; `ProjectedBasis(L=63,N=1,base=2)`; `basis(Int32;L=6,N=3,k=0)`; overflow-guard error path.

### Phase 3 — eltype-correctness sweep

**Goal:** real-phase bases produce real matrices and type-stable `index`.

| Item | Location | Fix |
|---|---|---|
| `parityflip-1` | `ParityBasis.jl:105`, `FlipBasis.jl:100`, `ParityFlipBasis.jl:121,131` | Add `eltype(...)=Float64` (mirror `TranslationParityBasis.jl:25`); change `ParityFlipBasis` literal `0.0 → zero(eltype(b))`. |
| `abelian-5` | `AbelianBasis.jl:584-590` | `eltype(::AbelianBasis{Ti,Tg}) = Tg<:Real ? Float64 : ComplexF64` (mirror `TranslationalBasis.jl:31`). |

**Tests:** `@inferred index(b,dgt)` is concrete; a real symmetric Hamiltonian on a k=0/parity sector assembles `Float64` (½ memory/FLOPs).

### Phase 4 — Cached-sparse wiring + apply-path unification

**Goal:** make `sparse!` the transparent accelerator it's documented to be, and collapse the 4-way divergence.

- **Optimal design:** introduce one `_apply!(target, opt, src, α, β)` that checks `_cached_sparse` once and dispatches SpMM vs matrix-free; make `*`, `mul`, `mul!` (3- and 5-arg) thin wrappers. Resolves `operator-2`, `operator-5`, `operator-4` (threading) together and removes the bug class.
- `operator-1` (`Operator.jl:306,354-358,370`): content-hash the cache key (or loudly document snapshot semantics).

**Tests:** `mul!` matches `*(opt,m)` with and without `sparse!`; cached path used from a block/iterative solver entry; thread-safety; stale-cache behavior under the chosen keying.

### Phase 5 — Remaining correctness bugs *(long tail)*

| Item | Location | Fix |
|---|---|---|
| `trans-1` | `TranslationalBasis.jl:352`, `TranslationalParityBasis.jl:227`, `TranslationalFlipBasis.jl:210` | `order()=ncycle(b)` (and `2*ncycle(b)` for parity/flip). Add an a>1 symmetrizer-isometry regression test. |
| `parityflip-2` | `ParityBasis.jl:46`, `FlipBasis.jl:44`, `ParityFlipBasis.jl:75` | After choosing the representative, re-check `f` on the decoded flipped/reversed orbit members; reject non-`f`-closed orbits. |
| `abelian-3` | `AbelianBasis.jl:759-761` | `index(B,dgt)=index(B,dgt,_shallow_workspace(B.G))`. Thread-safety test under threads. |
| `schmidt-1` | `Schmidt.jl:105-113,142` | Give `renyi_entropy` a `cutoff` and forward it. |
| `schmidt-2` | `Schmidt.jl:142,105-113` | Normalize over above-cutoff entries inside `renyi_entropy`; optionally widen von-Neumann dispatch to `abs(α-1)<1e-8`. |
| `expm-1` | `ToolKit.jl:58-91` | Scaling-and-squaring (`s=ceil(log2(opnorm(A)))`, Taylor of `A/2^s`, square s times); at minimum a loud `opnorm` warning above the convergence threshold. |
| ITensors `tebd4` | `TEBD.jl:163-191,185` | Clarify docstring (`exp(τ·h)`, real-time needs `τ=-i·dt`); optional EDKit-level real-time driver with truncation control. |
| `bug-4` | `AbstractBasis.jl:360` | Loosen `change!(dgt, ind; base)` to `Integer` with internal conversion. |
| `parityflip-4` | Parity/Flip/ParityFlip | Add the three missing `copy` methods. |
| `operator-6` | `Operator.jl:63-77` | `Dict{Vector{Int},Int}` O(1) dedup; guard against a silent zero-term Operator. |
| `gapratio-1` | `ToolKit.jl:15-29` | `sorted` keyword + guard; optional degeneracy `tol`; document single-sector use. |
| `qim-1` | `QIM.jl:22-31` | Document the Hermitian assumption + `ishermitian` guard (or compute the explicit anticommutator form). |

### Phase 6 — Optimizations

| Item | Location | Optimal approach |
|---|---|---|
| `optimization-1` | `ProjectedBasis.jl:195-205`; `TranslationalBasis.jl:180-182` | Bounded mixed-radix composition generator that never emits an out-of-range digit (kills 70–95% wasted iterations). **Shared machinery with Phase 7's mixed-radix extension** — build the primitive once. |
| `linearmap-1` | `LinearMap.jl:59-77,167-193,210-212` | Build `basis_embedding` via a sparse-triplet accumulator (≤1 nonzero/row) parameterized on `promote_type(eltype(B1),eltype(B2))` instead of a dense `base^L×dim`. |
| `abelian-4` | `AbelianBasis.jl:514,439-440` | Carry preallocated `ms`/`tmp` in the workspace already threaded through `colmn!`→`index`. |
| `o1` (residual) | `FermionOperator.jl:116-121` | Document/avoid the `2^span+1` Int64 `colptr` blow-up for long-range single-bond terms (endpoint + parity-coefficient construction). |

Parked micro-notes (verified immaterial — one-line rationale each, no action): `optimization-2`, `operator-7`, `o2`, `te-3`, `trans-3`, `parityflip-5`, `expm-2`, `schmidt-3`, ITensors `opt-1`.

### Phase 7 — Extensions *(ordered by value/cost)*

1. **RDM / mutual information** (`schmidt-4`, `Schmidt.jl:144-168,189-201`) — cheapest high-value: `rdm(v,Ainds,b)=S*S'` and `mutual_information(...)` reuse every existing `schmidt` dispatch. Purely additive, zero breaking surface. Ship first.
2. **Steady-state Lindblad** (`lindblad-6`, `LindbladEvolution.jl:138-151`) — `steadystate(A)`: smallest-magnitude / shift-invert eigensolve feeding the existing matrix-free `LiouvillianMap`, reshape, Hermitize, trace-normalize; dense fallback for small systems. **Dependency decision: add KrylovKit** (pure-Julia, likely already transitive via ITensorMPS) over Arpack.
3. **Spinful / multi-species fermions** (`e1`, `SpinlessFermionBasis.jl:13-17`) — 2L-mode basis with a fixed global JW ordering, `(site,spin)→mode` constructors, `(N↑,N↓)` sectors, Hubbard helper.
4. **Symmetry-resolved fermion basis** (`e2`, `FermionOperator.jl:170-179`) — momentum-resolved basis baking the boundary sign `(-1)^(N-1)` into `R`; then relax the existing JW guard.
5. **Backward / two-sided evolution** (`te-4`, `TimeEvolution.jl:556-577`) — allow signed τ (`cis(-τλ)`) with `|τ|` in the monitor and a per-cache direction flag; enables OTOC / Heisenberg-picture correlators.
6. **Mixed local dimensions** (`extension-1`, `AbstractBasis.jl:229-258,297-404`) — **most architecturally significant.** Parameterize the tensor basis on base-type (`TensorBasis{B}`, scalar `Int` keeps the allocation-free fast path; `Vector`/`NTuple` enables per-site dims), with Horner/divrem mixed-radix index math reusing the Phase 6 primitive. Builds on the clean width/eltype core from Phases 2–3.

Lower tier (demand-driven): qudit Pauli/Lindblad (`ext-3` — note `GELLMANN` already ships at `Operator.jl:830-841`), n-body fermion strings + `:right` convention (`e3`), hermitian-hopping helper (`e4`), per-jump rates/ħ in Lindblad (`lindblad-7`), SFF/number-variance/unfolding diagnostics (`diagnostics-1`).

New tests from extensions: TEBD Trotter accuracy vs exact (`ext-2`).

### Phase 8 — Uncertain triage

- `te-1` (`TimeEvolution.jl:392-393`, mirrored `LindbladEvolution.jl:384`): reproduce the defect-monitor 513-sample cap at large `ω·τ`; decide between shrinking the candidate interval / deriving the cap from `m_max` + warning, vs documenting as a known limit. Other uncertain items (`parityflip-3`, `o1` residual) are folded into Phases 2/6.

## 5. Testing strategy

- Every bug fix: a regression test that fails on the pre-fix code and passes after.
- Type-stability fixes: `@inferred` assertions.
- Behavior-changing fixes (eltype, `order()`): assert the *new* correct value and reference an unsymmetrized brute-force computation where possible.
- Thread-safety fixes (`abelian-3`, Phase 4): multi-threaded stress assertions (mirror `MultiThreadsTest`).
- All new tests wired into `runtests.jl` so CI runs them.

## 6. Versioning & changelog

- `CHANGELOG.md` created in Phase 0, appended each phase, finalized at the version bump.
- `Project.toml` version → `0.6.0` at the end of the roadmap (or at the first breaking phase to merge, if phases land as separate PRs).
- CHANGELOG groups: **Breaking**, **Fixed**, **Performance**, **Added**.

## 7. Finding → phase traceability

| Phase | Findings addressed |
|---|---|
| 0 | CI gap, `ext-1` (Pauli round-trip test), orphaned test adoption |
| 1 | `abelian-1`, `productstate-1`/`-2` |
| 2 | `bug-1`, `bug-2`, `bug-3`, `abelian-2`, `abelian-7`, `parityflip-3` |
| 3 | `parityflip-1`, `abelian-5` |
| 4 | `operator-2`, `operator-5`, `operator-4`, `operator-1` |
| 5 | `trans-1`, `parityflip-2`, `abelian-3`, `schmidt-1`, `schmidt-2`, `expm-1`, ITensors `tebd4`, `bug-4`, `parityflip-4`, `operator-6`, `gapratio-1`, `qim-1` |
| 6 | `optimization-1`, `linearmap-1`, `abelian-4`, `o1` residual; parked micro-notes |
| 7 | `schmidt-4`, `lindblad-6`, `e1`, `e2`, `te-4`, `extension-1`, `ext-3`, `e3`, `e4`, `lindblad-7`, `diagnostics-1`, `ext-2` |
| 8 | `te-1` |
| Excluded | 9 refuted findings |

## 8. Parked & excluded items

- **Parked (verified immaterial):** `optimization-2`, `operator-7`, `o2`, `te-3`, `trans-3`, `parityflip-5`, `expm-2`, `schmidt-3`, ITensors `opt-1` — recorded with a one-line rationale in the plan; revisited only if profiling motivates them.
- **Excluded (refuted):** the 9 findings the adversarial verifiers refuted (e.g. the "dense 2^span matrix" OOM framing, already built sparse). Not actioned.

## 9. Open questions for implementation-plan stage

- KrylovKit vs dense-only for `lindblad-6` (lean KrylovKit).
- Whether phases land as one PR or a stacked series (affects when `Project.toml` flips to 0.6.0).
- `expm-1`: full scaling-and-squaring vs loud-warning-only (lean full, it's small).
