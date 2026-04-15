# EDKit Full Documentation Overhaul Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rewrite EDKit's source docstrings, README, and docs pages so source documentation is AI-oriented, README is repository-orientation-oriented, and the docs site is more detailed and human-friendly.

**Architecture:** Execute the overhaul subsystem by subsystem so terminology and conventions stay consistent. Update source-level documentation first across the basis, operator, entanglement, ITensor, and algorithm layers; then rewrite README and selected docs pages to reflect the new documentation architecture; finally verify the docs build locally.

**Tech Stack:** Julia, Markdown, Documenter.jl, GitHub Pages

---

## Chunk 1: Planning And Audit

### Task 1: Record the design and target files

**Files:**
- Create: `docs/superpowers/specs/2026-03-19-full-documentation-overhaul-design.md`
- Create: `docs/superpowers/plans/2026-03-19-full-documentation-overhaul.md`

- [ ] **Step 1: Save the approved design**

Write the documentation design spec with the audience split and coverage goals.

- [ ] **Step 2: Save the implementation plan**

Write this plan file with subsystem-based tasks and verification commands.

### Task 2: Audit documentation coverage across `src/`

**Files:**
- Review: `src/Basis/*.jl`
- Review: `src/Operator.jl`
- Review: `src/LinearMap.jl`
- Review: `src/Schmidt.jl`
- Review: `src/ToolKit.jl`
- Review: `src/ITensors/*.jl`
- Review: `src/algorithms/*.jl`

- [ ] **Step 1: Read each subsystem source file**

Run file reads on each source file and identify thin docstrings, missing argument/return explanations, and internal helpers lacking semantic guidance.

- [ ] **Step 2: Group gaps by subsystem**

Organize the rewrite into basis, operator/map, entanglement, ITensor, and algorithm layers.

## Chunk 2: Basis Layer

### Task 3: Rewrite basis docstrings comprehensively

**Files:**
- Modify: `src/Basis/AbstractBasis.jl`
- Modify: `src/Basis/ProjectedBasis.jl`
- Modify: `src/Basis/TranslationalBasis.jl`
- Modify: `src/Basis/TranslationalParityBasis.jl`
- Modify: `src/Basis/TranslationalFlipBasis.jl`
- Modify: `src/Basis/ParityBasis.jl`
- Modify: `src/Basis/FlipBasis.jl`
- Modify: `src/Basis/ParityFlipBasis.jl`
- Modify: `src/Basis/AbelianBasis.jl`

- [ ] **Step 1: Expand base abstractions**

Document `AbstractBasis`, index/change conventions, digit conventions, and `TensorBasis`.

- [ ] **Step 2: Expand projected and translational basis docs**

Clarify predicates, charge sectors, momentum conventions, normalization, and helper roles.

- [ ] **Step 3: Expand discrete-symmetry and Abelian basis docs**

Document parity/flip combinations, compatibility assumptions, and the high-level `basis(...)` constructor.

## Chunk 3: Operator, Maps, And Entanglement

### Task 4: Rewrite operator and mapping docstrings

**Files:**
- Modify: `src/Operator.jl`
- Modify: `src/LinearMap.jl`
- Modify: `src/ToolKit.jl`

- [ ] **Step 1: Rewrite `Operator` construction and conversion docs**

Document matrices, indices, basis semantics, merging behavior, matrix-free usage, and helper routines.

- [ ] **Step 2: Rewrite map and symmetrizer docs**

Document `DoubleBasis`, `symmetrizer`, embedding helpers, and compatibility constraints.

- [ ] **Step 3: Rewrite utilities docs**

Clarify return semantics, numerical caveats, and intended usage for statistics and small linear algebra helpers.

### Task 5: Rewrite entanglement docstrings

**Files:**
- Modify: `src/Schmidt.jl`

- [ ] **Step 1: Rewrite Schmidt matrix and entropy helper docs**

Document bipartition conventions, basis dependence, return meanings, and internal helper roles.

- [ ] **Step 2: Rewrite basis-specific Schmidt docs**

Explain how translational, parity, flip, and Abelian bases alter the decomposition logic.

## Chunk 4: ITensor And Algorithm Layers

### Task 6: Rewrite ITensor and Pauli-layer docstrings

**Files:**
- Modify: `src/ITensors/ITensorsKit.jl`
- Modify: `src/ITensors/PauliBasis.jl`
- Modify: `src/ITensors/TEBD.jl`

- [ ] **Step 1: Expand conversion helper docs**

Document row-major conventions, shape assumptions, and basis restrictions.

- [ ] **Step 2: Expand Pauli-space docs**

Document operator-space conventions, coefficient meaning, and tensor-network conversion semantics.

- [ ] **Step 3: Expand TEBD docs**

Document sweep construction, expected inputs, and output interpretation.

### Task 7: Rewrite algorithm-layer docstrings

**Files:**
- Modify: `src/algorithms/Lindblad.jl`
- Modify: `src/algorithms/QIM.jl`

- [ ] **Step 1: Expand Lindblad docs**

Document both many-body and quadratic workflows, matrix conventions, and return objects.

- [ ] **Step 2: Expand inverse-method docs**

Document covariance interpretation, assumptions on operator lists, and solution semantics.

## Chunk 5: README And Human-Facing Docs

### Task 8: Rewrite README for agent onboarding

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Add architecture and navigation guidance**

Explain subsystem ownership, canonical entry points, and common task-to-file mapping.

- [ ] **Step 2: Keep installation and quick-start concise**

Preserve basic user-facing setup while making the README more useful to AI agents.

### Task 9: Upgrade human-facing docs pages

**Files:**
- Modify: `docs/src/index.md`
- Modify: `docs/src/getting-started.md`
- Modify: `docs/src/manual/*.md`
- Modify: `docs/src/reference/*.md`
- Modify: `docs/src/developer.md`

- [ ] **Step 1: Update landing and manual pages**

Make explanations more detailed and explicitly human-oriented.

- [ ] **Step 2: Update developer page**

Explain the documentation split between source docstrings, README, and docs pages.

## Chunk 6: Verification

### Task 10: Verify that the documentation site still builds

**Files:**
- Verify: `docs/Project.toml`
- Verify: `docs/make.jl`
- Verify: `docs/src/`
- Verify: `README.md`
- Verify: `src/`

- [ ] **Step 1: Instantiate docs environment**

Run: `julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'`
Expected: docs dependencies resolve successfully.

- [ ] **Step 2: Build the docs site**

Run: `julia --project=docs docs/make.jl`
Expected: successful Documenter build with no fatal errors.

- [ ] **Step 3: Review git state**

Run: `git status --short` and `git diff --stat`
Expected: changes are limited to documentation source files and deliberate docstring edits.

- [ ] **Step 4: Summarize remaining gaps**

If any lower-priority dev files or unusually repetitive helpers remain less detailed than the rest, call that out explicitly.
