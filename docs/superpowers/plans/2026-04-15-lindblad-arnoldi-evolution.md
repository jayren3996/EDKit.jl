# Lindblad Arnoldi Evolution Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an adaptive many-body Lindblad time-evolution solver that reuses Arnoldi subspaces and exposes approximation-quality diagnostics.

**Architecture:** Keep the legacy dense `Lindblad` Taylor stepper unchanged and add a new `LiouvillianMap` plus `LindbladArnoldiCache` in `src/algorithms/LindbladEvolution.jl`. Internally the solver uses vectorized density matrices and a cached complex-Schur reduced model, while the public API stays `DensityMatrix`-based.

**Tech Stack:** Julia, `LinearAlgebra`, `SparseArrays`, EDKit `DensityMatrix`/`Lindblad`, adaptive Arnoldi with cached Schur reduced-space propagation.

---

### Task 1: Add failing public-API tests

**Files:**
- Create: `test/lindblad_timeevolve_tests.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Write the failing tests**

Add tests that call:

```julia
lb = lindblad(H, jumps)
rho0 = densitymatrix(ComplexF64[0.0, 1.0])
rho1 = lindblad_timeevolve(lb, rho0, 0.2; tol = 1e-11, m_init = 6, m_max = 16)
```

and:

```julia
A = LiouvillianMap(sparse(H), sparse.(jumps))
cache = LindbladArnoldiCache(A, rho0; tol = 1e-11, m_init = 6, m_max = 16)
rho2 = lindblad_timeevolve!(cache, 0.2)
```

- [ ] **Step 2: Run tests to verify they fail**

Run:

```bash
~/.juliaup/bin/julia --project -e 'include("test/TestHelpers.jl"); include("test/lindblad_timeevolve_tests.jl")'
```

Expected: failures for undefined `LiouvillianMap`, `LindbladArnoldiCache`, and `lindblad_timeevolve`.

### Task 2: Implement the new Arnoldi solver

**Files:**
- Create: `src/algorithms/LindbladEvolution.jl`

- [ ] **Step 1: Define the new types**

Add `LiouvillianMap`, `LindbladArnoldiCache`, `LindbladArnoldiDiagnostics`, and internal workspace helpers.

- [ ] **Step 2: Implement matrix-free Liouvillian action**

Support dense/sparse matrix backends and `Lindblad` input through:

```julia
LiouvillianMap(lb::Lindblad)
LiouvillianMap(H::AbstractMatrix, jumps::AbstractVector{<:AbstractMatrix})
```

- [ ] **Step 3: Implement adaptive Arnoldi + cached Schur reduced model**

Add Arnoldi build/extend/restart logic, Schur cache rebuilds, reduced-space step caching, and state reconstruction back to `DensityMatrix`.

- [ ] **Step 4: Implement diagnostics**

Track raw per-output `trace_drifts`, `hermiticity_drifts`, and optional `min_hermitian_eigs`, plus summary extrema and cache activity counters.

### Task 3: Make tests green incrementally

**Files:**
- Modify: `src/algorithms/LindbladEvolution.jl`
- Test: `test/lindblad_timeevolve_tests.jl`

- [ ] **Step 1: Run the new targeted tests**

Run:

```bash
~/.juliaup/bin/julia --project -e 'include("test/TestHelpers.jl"); include("test/lindblad_timeevolve_tests.jl")'
```

Expected: PASS.

- [ ] **Step 2: Run neighboring integration tests**

Run:

```bash
~/.juliaup/bin/julia --project -e 'include("test/TestHelpers.jl"); include("test/lindblad_tests.jl"); include("test/timeevolve_tests.jl")'
```

Expected: PASS.

### Task 4: Document the new solver

**Files:**
- Modify: `docs/src/manual/lindblad.md`
- Modify: `docs/src/reference/lindblad.md`
- Modify: `docs/src/examples/open-system-workflows.md`

- [ ] **Step 1: Add user-facing docs**

Document the new `lindblad_timeevolve` API, the Schur/Arnoldi reuse model, and the caveats around trace, Hermiticity, and positivity.

- [ ] **Step 2: Build docs**

Run:

```bash
~/.juliaup/bin/julia --project=docs docs/make.jl
```

Expected: docs build succeeds.

### Task 5: Final verification

**Files:**
- Modify: `src/algorithms/LindbladEvolution.jl`
- Modify: `test/lindblad_timeevolve_tests.jl`
- Modify: docs files above

- [ ] **Step 1: Run the relevant full verification**

Run:

```bash
~/.juliaup/bin/julia --project test/runtests.jl
```

Expected: PASS.

- [ ] **Step 2: Review the diff and summarize remaining risks**

Review:

```bash
git diff -- src/algorithms/LindbladEvolution.jl test/lindblad_timeevolve_tests.jl docs/src/manual/lindblad.md docs/src/reference/lindblad.md docs/src/examples/open-system-workflows.md test/runtests.jl
```
