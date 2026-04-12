# EDKit Agent Notes Layout Design

**Date:** 2026-04-12
**Status:** Design
**Audience:** AI agents working in the repository
**Goal:** Make it easy for an agent to pick up local context quickly without loading a large root note or reading human-facing documentation first.

## Motivation

EDKit already has useful top-level orientation in `README.md` and `AGENTS.md`, but
the current guidance is too centralized:

- an agent entering `src/Basis/` still has to infer basis-specific invariants
  from code and scattered docs,
- an agent entering `test/` or `docs/` has no local operating guidance,
- the repository currently duplicates the same guidance in `AGENTS.md` and
  `CLAUDE.md`, which increases maintenance cost and drift risk.

The new design should optimize for agent context efficiency, not human prose.

## Scope

This design covers repository-local agent notes only.

In scope:

- slimming and refocusing the root `AGENTS.md`,
- adding a small set of folder-local `AGENTS.md` files at subsystem boundaries,
- deleting `CLAUDE.md` and standardizing on `AGENTS.md` as the only agent-note
  format,
- documenting maintenance rules so the note tree stays short and useful.

Out of scope:

- changing package runtime behavior,
- changing the public API,
- rewriting the user manual for human readers,
- adding note files to every subfolder,
- treating `examples/` or `dev/` as note-bearing ownership boundaries.

## Design Principles

### 1. Locality over completeness

An agent should learn the minimum needed for the folder it is currently in.
Folder-local notes are preferred over one large global note.

### 2. Stable invariants only

Agent notes should capture stable semantics, ownership boundaries, hot-path
warnings, and test expectations. They should not duplicate long examples or
copy large chunks of the user manual.

### 3. One source of truth

The repository should use `AGENTS.md` only. Parallel `CLAUDE.md` files create
unnecessary duplication and will be removed.

### 4. Subsystem boundaries, not directory saturation

Only add local notes where the working model changes materially. A note tree
that is too dense wastes context and becomes hard to maintain.

## Proposed Note Layout

### Root

`AGENTS.md`

Purpose:

- package-wide orientation,
- top-level module map,
- canonical build/test/docs commands,
- global invariants that apply everywhere,
- pointers to local `AGENTS.md` files.

This file should become thinner than it is today. Detailed subsystem guidance
should move down to local notes.

### Core source layer

`src/AGENTS.md`

Purpose:

- explain the core EDKit abstraction split between basis, operator, maps, and
  entanglement,
- identify canonical entry points in `src/`,
- describe how to reason about public API additions versus internal helpers,
- flag performance-sensitive verification expectations after touching core files.

This note covers `src/Operator.jl`, `src/LinearMap.jl`, `src/Schmidt.jl`, and
`src/ToolKit.jl`.

### Basis layer

`src/Basis/AGENTS.md`

Purpose:

- basis taxonomy and when each basis type is the correct target,
- invariants for `index`, `change!`, `content`, `order`, representative states,
  and normalization,
- explicit-buffer-passing and thread-safety expectations,
- guidance for 2D/3D lattice support through `basis(...; symmetries=...)`,
- basis-specific test placement guidance.

This is the most important local note because basis semantics drive operator
application, maps, and entanglement logic.

### ITensor bridge

`src/ITensors/AGENTS.md`

Purpose:

- explain what is source-of-truth between EDKit vectors/operators and ITensor
  objects,
- record round-trip expectations,
- warn about site ordering and basis mismatch failures,
- list the tests to run after changing conversion helpers.

### Algorithm layer

`src/algorithms/AGENTS.md`

Purpose:

- distinguish many-body Lindblad workflows from quadratic/covariance workflows,
- identify numerically sensitive helper families,
- state expected input/output conventions,
- warn against mixing conceptual models between algorithm families.

### Test layer

`test/AGENTS.md`

Purpose:

- map test files to subsystems,
- give minimal targeted test commands before `Pkg.test()`,
- state the rule that new user-visible behavior should be tested at the public
  API level first,
- identify where slower regression/performance coverage belongs.

### Docs layer

`docs/AGENTS.md`

Purpose:

- explain the split between source docstrings, manual pages, examples, and
  reference pages,
- provide local docs build commands,
- remind agents to update docs navigation when adding pages,
- emphasize discoverability so pages do not become orphaned.

## Explicit Non-Layout Decisions

The design intentionally does **not** add:

- `examples/AGENTS.md`
- `dev/AGENTS.md`
- one-note-per-subfolder inside `src/Basis/`
- duplicate `AGENTS.md` and `CLAUDE.md` trees

These areas are useful reference material, but they are not the right
boundaries for persistent agent operating guidance.

## Content Rules Per Note

Each `AGENTS.md` should be:

- short,
- bullet-oriented,
- operational rather than tutorial-style,
- scoped to decisions an agent actually needs in that folder.

Target size:

- root note: roughly 40-70 lines,
- local notes: roughly 25-60 lines.

Preferred content:

- ownership boundaries,
- invariants,
- hot-path and performance warnings,
- subsystem-specific conventions,
- test/build commands,
- “if you are doing X, read Y first” pointers.

Avoid:

- long package introductions,
- repeated examples that already exist in docs,
- detailed physics exposition,
- repeating the entire root note in local notes.

## Proposed Content by File

### `AGENTS.md`

Keep:

- package overview,
- subsystem map,
- global build/test commands,
- global performance notes that truly apply repo-wide.

Add:

- a short note tree map telling agents which local note to read in `src/`,
  `src/Basis/`, `src/ITensors/`, `src/algorithms/`, `test/`, and `docs/`.

Trim:

- detailed subsystem content that belongs in local notes.

### `src/AGENTS.md`

Should include:

- the core abstraction graph,
- public entry points in `src/`,
- expectations when editing exported APIs,
- verification hints for operator/map/entanglement changes.

### `src/Basis/AGENTS.md`

Should include:

- basis family overview,
- representative/orbit semantics,
- `index` and `change!` invariants,
- normalization and compatibility assumptions,
- explicit-buffer-passing rule,
- specialization guidance for new 2D basis work.

### `src/ITensors/AGENTS.md`

Should include:

- round-trip expectations,
- site index ordering rules,
- basis consistency requirements,
- which tests are the canonical checks.

### `src/algorithms/AGENTS.md`

Should include:

- split between solver families,
- numerical caution points,
- expected verification routes.

### `test/AGENTS.md`

Should include:

- subsystem-to-test-file map,
- fast targeted commands,
- when to add regression tests versus smoke tests.

### `docs/AGENTS.md`

Should include:

- documentation layer boundaries,
- docs build commands,
- navigation/update expectations.

## Migration Plan

1. Rewrite the root `AGENTS.md` to be thinner and pointer-oriented.
2. Add local `AGENTS.md` files in `src/`, `src/Basis/`, `src/ITensors/`,
   `src/algorithms/`, `test/`, and `docs/`.
3. Delete `CLAUDE.md`.
4. Ensure there is no duplicated large guidance block between root and local
   notes.
5. Verify that each local note can be read independently and still make sense.

## Risks

### Risk: Drift between root and local notes

Mitigation:

- keep the root note short,
- move detailed rules down,
- avoid duplicating the same bullets in multiple places.

### Risk: Notes become mini-manuals

Mitigation:

- keep them operational,
- cap their size,
- link outward instead of copying tutorial content.

### Risk: Missing a useful boundary

Mitigation:

- start with the six-note tree above,
- only add finer-grained notes later if repeated agent failures justify it.

## Success Criteria

This design is successful when:

- an agent entering a subsystem can get enough local context from one short
  nearby note,
- the root note no longer has to explain every subsystem in detail,
- there is only one agent-note format (`AGENTS.md`),
- maintenance burden is lower than the current duplicated setup,
- the note tree improves agent onboarding without noticeably increasing
  repository noise.
