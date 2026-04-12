# EDKit Agent Notes Layout Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the duplicated top-level agent guidance with a single `AGENTS.md` tree that gives short, local subsystem context for agents working in EDKit.

**Architecture:** Keep one thin root `AGENTS.md` as the package map, add short subsystem-local `AGENTS.md` files at `src/`, `src/Basis/`, `src/ITensors/`, `src/algorithms/`, `test/`, and `docs/`, and delete `CLAUDE.md` so there is one source of truth.

**Tech Stack:** Markdown, Julia repository conventions, git verification commands

---

## File Structure

| File | Responsibility | Action |
| --- | --- | --- |
| `AGENTS.md` | Thin package-level note with global invariants and local-note pointers | Rewrite |
| `CLAUDE.md` | Duplicate legacy agent note | Delete |
| `src/AGENTS.md` | Core runtime semantics for operator/map/entanglement work | Create |
| `src/Basis/AGENTS.md` | Basis invariants, orbit semantics, and 2D/3D guidance | Create |
| `src/ITensors/AGENTS.md` | ITensor bridge conventions and round-trip warnings | Create |
| `src/algorithms/AGENTS.md` | Lindblad/QIM workflow boundaries and numerical cautions | Create |
| `test/AGENTS.md` | Test map and targeted verification commands | Create |
| `docs/AGENTS.md` | Docs-layer boundaries and local docs build guidance | Create |

---

### Task 1: Replace the duplicated root guidance

**Files:**
- Modify: `AGENTS.md`
- Delete: `CLAUDE.md`

- [ ] **Step 1: Rewrite the root note as a thin package map**

Replace `AGENTS.md` with a note organized around:

```md
# AGENTS.md

## Scope
- EDKit is a Julia exact-diagonalization package built around basis selection plus operator assembly.
- This root note only covers repo-wide context.
- Read local `AGENTS.md` files when working in `src/`, `src/Basis/`, `src/ITensors/`, `src/algorithms/`, `test/`, or `docs/`.

## Package Map
- `src/Basis/`: basis types and symmetry reduction
- `src/Operator.jl`: many-body operator assembly and application
- `src/LinearMap.jl`: basis-to-basis maps
- `src/Schmidt.jl`: entanglement and Schmidt logic
- `src/ITensors/`: ITensor/MPS/MPO bridge
- `src/algorithms/`: Lindblad and QIM workflows

## Canonical Entry Points
- basis selection: `TensorBasis`, `ProjectedBasis`, `basis(...)`
- operator construction: `spin`, `operator`, `trans_inv_operator`
- maps: `DoubleBasis`, `symmetrizer`
- entanglement: `ent_S`, `schmidt`

## Repo-Wide Invariants
- Basis constructors use keyword arguments.
- `Operator` is immutable; sparse caching is external.
- Internal hot paths use explicit `dgt` buffer passing for thread safety.
- For arbitrary 2D/3D lattice symmetries, prefer `basis(...; symmetries=...)`.

## Commands
- tests: `~/.juliaup/bin/julia --project -e 'using Pkg; Pkg.test()'`
- docs: `~/.juliaup/bin/julia --project=docs docs/make.jl`
```

- [ ] **Step 2: Delete the duplicate note**

Remove `CLAUDE.md`.

Run:

```bash
cd /Users/ren/Library/CloudStorage/OneDrive-UniversityofLeeds/GitHub/EDKit.jl
rm CLAUDE.md
```

Expected: `git status --short` no longer lists `CLAUDE.md`

---

### Task 2: Add local source-layer notes

**Files:**
- Create: `src/AGENTS.md`
- Create: `src/Basis/AGENTS.md`
- Create: `src/ITensors/AGENTS.md`
- Create: `src/algorithms/AGENTS.md`

- [ ] **Step 1: Add the core source-layer note**

Create `src/AGENTS.md` with content shaped like:

```md
# AGENTS.md

## Scope
- Applies to `src/Operator.jl`, `src/LinearMap.jl`, `src/Schmidt.jl`, and `src/ToolKit.jl`.
- Read `src/Basis/AGENTS.md` before changing basis semantics.

## Core Model
- Basis objects define coordinates and symmetry sectors.
- `Operator` stores local terms plus basis embedding, not just a prebuilt matrix.
- Maps and entanglement logic assume the basis/operator split stays intact.

## Canonical APIs
- `operator`, `trans_inv_operator`, `spin`
- `DoubleBasis`, `symmetrizer`
- `ent_S`, `schmidt`

## Editing Rules
- Prefer preserving public entry points and extending internals behind them.
- When changing operator application semantics, re-check dense, sparse, and matrix-free paths.
- If you touch shared hot paths, verify basis-dependent tests as well as core operator tests.
```

- [ ] **Step 2: Add the basis-layer note**

Create `src/Basis/AGENTS.md` with content shaped like:

```md
# AGENTS.md

## Scope
- Applies to all basis constructors and basis internals.

## Basis Taxonomy
- `TensorBasis`: full tensor-product basis
- `ProjectedBasis`: constrained/fixed-charge subspace
- dedicated 1D symmetry bases: translational/parity/flip combinations
- `AbelianBasis`: generic commuting permutation symmetries
- `basis(...)`: high-level dispatcher

## Invariants
- `change!(b, i, dgt)` decodes stored coordinates into a digit buffer.
- `index(b, dgt)` interprets a digit buffer in basis coordinates and may return zero coefficient when outside the sector.
- Reduced bases store canonical representatives of symmetry orbits.
- Normalization factors and representative choice must stay consistent with `index` and operator assembly.

## Threading Rule
- Internal hot paths must use explicit external `dgt` buffers.
- Do not reintroduce shared mutable `b.dgt` use inside threaded loops.

## 2D/3D Guidance
- For arbitrary lattices, prefer `basis(...; symmetries=...)`.
- Only add a specialized 2D basis type if there is a clear ergonomics or performance gap relative to `AbelianBasis`.

## Tests
- Primary files: `test/basis_tests.jl`, `test/AbelianBasisTest.jl`, `test/abelian_overhaul_tests.jl`, `test/entanglement_tests.jl`.
```

- [ ] **Step 3: Add the ITensor and algorithm notes**

Create `src/ITensors/AGENTS.md`:

```md
# AGENTS.md

## Scope
- Applies to vector/MPS/MPO conversions, Pauli-basis helpers, and TEBD utilities.

## Working Rules
- Preserve site ordering exactly across EDKit vectors, ITensor site indices, and operator-space conversions.
- Treat EDKit coordinate conventions as the source of truth when round-tripping with ITensors.
- When changing conversions, check both plain vector round-trips and symmetry-resolved `mps2vec(..., B)` behavior.

## Main APIs
- `vec2mps`, `mps2vec`, `mat2op`, `op2mat`
- `pauli`, `pauli_list`, `commutation_mat`, `dissipation_mat`
- `mps2pmps`, `pmps2mpo`, `mpo2pmpo`, `tebd4`

## Tests
- Primary files: `test/itensor_tests.jl`, `test/TensorTest.jl`.
```

Create `src/algorithms/AGENTS.md`:

```md
# AGENTS.md

## Scope
- Applies to `Lindblad.jl` and `QIM.jl`.

## Solver Split
- `lindblad` and `DensityMatrix` cover explicit many-body density-matrix workflows.
- `quadraticlindblad`, `CovarianceMatrix`, and related helpers cover quadratic/Majorana covariance workflows.
- `qimsolve` and `covmat` are separate inverse-method utilities; do not mix their assumptions with Lindblad code casually.

## Editing Rules
- Keep input/output object conventions stable.
- Be careful with normalization, Hermiticity, and real-versus-complex assumptions.
- Verify both many-body and quadratic paths when touching shared physical semantics.

## Tests
- Primary files: `test/lindblad_tests.jl`, `test/advanced_tests.jl`.
```

---

### Task 3: Add test/docs notes and verify the note tree

**Files:**
- Create: `test/AGENTS.md`
- Create: `docs/AGENTS.md`

- [ ] **Step 1: Add the test-layer note**

Create `test/AGENTS.md` with:

```md
# AGENTS.md

## Scope
- This folder is organized by subsystem, not by one-file-per-source-file.

## Test Map
- `core_tests.jl`: operators, spin helpers, dense/sparse conversion
- `basis_tests.jl`: core basis constructors and decomposition
- `AbelianBasisTest.jl`, `abelian_overhaul_tests.jl`: Abelian and permutation-based symmetry coverage
- `entanglement_tests.jl`: entropy/Schmidt checks across sectors
- `itensor_tests.jl`, `TensorTest.jl`: ITensor and Pauli-space coverage
- `lindblad_tests.jl`, `advanced_tests.jl`: algorithms and higher-level integration

## Preferred Verification
- targeted file: `~/.juliaup/bin/julia --project test/<file>.jl`
- full suite: `~/.juliaup/bin/julia --project -e 'using Pkg; Pkg.test()'`

## Rules
- Test public behavior first, internals second.
- Add regression coverage near the subsystem most likely to break again.
```

- [ ] **Step 2: Add the docs-layer note**

Create `docs/AGENTS.md` with:

```md
# AGENTS.md

## Scope
- Applies to documentation source and build workflow.

## Documentation Layers
- `src/` docstrings: semantic reference for agents and advanced readers
- `docs/src/manual/`: concepts and decision guidance
- `docs/src/examples/`: short workflow pages
- `docs/src/reference/`: API surfacing via `@docs`

## Rules
- When adding a page, update `docs/make.jl` navigation.
- Prefer cross-links over copying the same explanation into multiple pages.
- Avoid orphan pages: discoverability matters as much as content quality.

## Build
- `~/.juliaup/bin/julia --project=docs docs/make.jl`
```

- [ ] **Step 3: Verify file placement and formatting**

Run:

```bash
cd /Users/ren/Library/CloudStorage/OneDrive-UniversityofLeeds/GitHub/EDKit.jl
find . -name AGENTS.md -o -name CLAUDE.md | sort
git diff --check
```

Expected:

- `CLAUDE.md` absent
- `AGENTS.md` present at root, `src/`, `src/Basis/`, `src/ITensors/`, `src/algorithms/`, `test/`, and `docs/`
- `git diff --check` prints no whitespace or conflict-marker errors
