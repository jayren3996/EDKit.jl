# EDKit Documentation Site Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a manual-first `Documenter.jl` site for `EDKit.jl`, deploy it with GitHub Actions to GitHub Pages, and add a documentation link to the README.

**Architecture:** Add a dedicated docs environment and build script, create a structured set of manual and reference pages under `docs/src`, wire deployment through a standard Documenter GitHub Actions workflow, and verify the site builds locally before closing the task.

**Tech Stack:** Julia 1.9+, Documenter.jl, GitHub Actions, Markdown

---

## Chunk 1: Docs Infrastructure

### Task 1: Add the Documenter environment and build entrypoint

**Files:**
- Create: `docs/Project.toml`
- Create: `docs/make.jl`

- [ ] **Step 1: Review the repository package metadata**

Run: `sed -n '1,220p' Project.toml`
Expected: confirm package name, compat, and base Julia version for docs setup.

- [ ] **Step 2: Write the docs environment**

Create `docs/Project.toml` with:
- `Documenter`
- local package dependency via `EDKit`
- any stdlib or package dependencies needed by documentation examples

- [ ] **Step 3: Write the docs build script**

Create `docs/make.jl` that:
- activates the docs environment,
- loads `Documenter` and `EDKit`,
- configures `makedocs` with site name, module list, source format, and navigation,
- configures `deploydocs` for the GitHub repository.

- [ ] **Step 4: Run the docs build once to observe the initial failure or missing-pages state**

Run: `julia --project=docs docs/make.jl`
Expected: fail because content pages are not all present yet, or reveal missing docs structure that must be added next.

## Chunk 2: Manual Pages

### Task 2: Create the top-level docs pages

**Files:**
- Create: `docs/src/index.md`
- Create: `docs/src/getting-started.md`
- Create: `docs/src/developer.md`

- [ ] **Step 1: Draft the landing page**

Write `docs/src/index.md` using README material, but reshape it into a docs homepage with:
- package overview,
- installation,
- quick navigation,
- a short “what EDKit covers” section.

- [ ] **Step 2: Draft the getting-started page**

Write `docs/src/getting-started.md` around the package’s core workflow and include one compact example that constructs and applies an `Operator`.

- [ ] **Step 3: Draft the developer page**

Write `docs/src/developer.md` describing docs layout, build commands, and how to extend docstrings and manual pages.

### Task 3: Create the manual section pages

**Files:**
- Create: `docs/src/manual/architecture.md`
- Create: `docs/src/manual/bases.md`
- Create: `docs/src/manual/operators.md`
- Create: `docs/src/manual/maps.md`
- Create: `docs/src/manual/entanglement.md`
- Create: `docs/src/manual/itensors.md`
- Create: `docs/src/manual/lindblad.md`
- Create: `docs/src/manual/utilities.md`

- [ ] **Step 1: Review the relevant source files and examples**

Run:
- `sed -n '1,260p' src/Operator.jl`
- `sed -n '1,260p' src/LinearMap.jl`
- `sed -n '1,260p' src/ToolKit.jl`
- `sed -n '1,260p' src/Schmidt.jl`
- inspect the basis and ITensor source files as needed

Expected: gather enough material to explain the public behavior of each subsystem accurately.

- [ ] **Step 2: Write the architecture and basis pages**

Focus on mental model, subsystem boundaries, and basis selection guidance.

- [ ] **Step 3: Write the operator, map, and entanglement pages**

Use short examples that connect directly to central workflows.

- [ ] **Step 4: Write the ITensor, Lindblad, and utilities pages**

Explain advanced workflows clearly and distinguish the different solver/tool regimes.

## Chunk 3: Examples And Reference

### Task 4: Add curated example pages

**Files:**
- Create: `docs/src/examples/basic-workflows.md`
- Create: `docs/src/examples/symmetry-workflows.md`
- Create: `docs/src/examples/tensor-network-workflows.md`
- Create: `docs/src/examples/open-system-workflows.md`

- [ ] **Step 1: Review example README files and existing example topics**

Run:
- `sed -n '1,220p' examples/README.md`
- `sed -n '1,220p' examples/Basic/README.md`
- `sed -n '1,220p' examples/Maps/README.md`
- `sed -n '1,220p' examples/TensorNetworks/README.md`
- `sed -n '1,220p' examples/Lindblad/README.md`

- [ ] **Step 2: Convert the highest-value example ideas into docs pages**

Keep each page focused, shorter than the raw notebook set, and aligned with the manual.

### Task 5: Add grouped API reference pages

**Files:**
- Create: `docs/src/reference/api-overview.md`
- Create: `docs/src/reference/bases.md`
- Create: `docs/src/reference/operators.md`
- Create: `docs/src/reference/maps.md`
- Create: `docs/src/reference/entanglement.md`
- Create: `docs/src/reference/itensors.md`
- Create: `docs/src/reference/lindblad.md`
- Create: `docs/src/reference/utilities.md`

- [ ] **Step 1: Identify the most important public names for each subsystem**

Use `rg` on `src/` and docstrings to group the exported or user-facing APIs.

- [ ] **Step 2: Write the reference pages with `@docs` blocks**

Each page should include a short orientation paragraph and grouped docs blocks.

- [ ] **Step 3: Build the docs to catch missing names or undocumented references**

Run: `julia --project=docs docs/make.jl`
Expected: either success or actionable warnings/errors about bad references.

## Chunk 4: Deployment And README

### Task 6: Add GitHub Actions docs deployment

**Files:**
- Create: `.github/workflows/documentation.yml`

- [ ] **Step 1: Review current GitHub workflows**

Run: `find .github/workflows -maxdepth 1 -type f | sort`
Expected: understand current CI footprint and avoid conflicting workflow assumptions.

- [ ] **Step 2: Write the documentation workflow**

Create a workflow that:
- runs on pushes to the default branch and pull requests,
- installs Julia,
- instantiates the docs environment,
- builds the docs,
- deploys on pushes to the default branch.

### Task 7: Add the README documentation link

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Choose a visible location near the top of the README**

Prefer the title/badge region or immediately after the opening summary.

- [ ] **Step 2: Add the docs link and a brief mention of the manual**

Keep the change small and clear.

## Chunk 5: Verification

### Task 8: Verify the documentation build and repository state

**Files:**
- Verify: `docs/`
- Verify: `.github/workflows/documentation.yml`
- Verify: `README.md`

- [ ] **Step 1: Run the complete docs build**

Run: `julia --project=docs docs/make.jl`
Expected: successful build with generated site output under `docs/build`.

- [ ] **Step 2: Inspect the resulting pages list or build output**

Confirm the navigation includes manual, examples, and reference sections.

- [ ] **Step 3: Review git diff**

Run: `git status --short` and `git diff --stat`
Expected: only the intended docs, workflow, and README changes.

- [ ] **Step 4: Summarize any residual risks**

If there are remaining doc gaps due to sparse docstrings or untested example snippets, note them explicitly in the handoff.
