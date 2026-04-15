# EDKit Full Documentation Overhaul Design

## Goal

Rework the package documentation at three layers so that each layer serves a distinct audience and does so well:

- source-level docstrings in `src/` should primarily help AI agents and advanced developers understand code semantics precisely,
- `README.md` should primarily help AI agents and repository newcomers orient themselves quickly inside the codebase,
- the Documenter site under `docs/src/` should primarily help human readers learn the package and its workflows in a more guided and explanatory way.

The work should be broad rather than selective. The target is not merely public API polish, but near-complete documentation coverage across types, functions, and important internal helpers.

## Audience Split

### Source Docstrings

Source docstrings are the machine-facing and implementation-facing documentation layer.

They should help an AI agent answer questions such as:

- What problem does this function solve?
- What does each argument mean semantically, not just syntactically?
- What assumptions or invariants must already hold?
- What does the return value represent?
- Is the function mutating, allocating, normalizing, projecting, or reindexing?
- Where does this function sit in the larger algorithmic pipeline?

Therefore, source docstrings should be explicit, structured, and semantically dense.

### README

The README should function as a compact agent onboarding guide for the repository.

It should answer:

- What is this package for?
- What are the main subsystems?
- Which files are the canonical entry points?
- Which APIs are stable and central?
- Which conventions does the package use for bases, digits, symmetries, and operators?
- Where should an agent look first for examples or workflow guidance?

This README can still serve human readers, but it should lean toward orienting an intelligent code-reading collaborator rather than trying to be the complete user manual.

### Documentation Site

The Documenter site remains the human-facing manual.

It should emphasize:

- conceptual explanations,
- subsystem overviews,
- careful examples,
- decision guidance such as when to use one basis or solver instead of another,
- and navigation across the package.

Compared with the README, the docs site should be more pedagogical and more willing to spend words on intuition and interpretation.

## Source Docstring Standard

The codebase should adopt a consistent docstring style across public APIs and internal helpers.

For user-facing functions and types, docstrings should generally include:

- a short summary line,
- a fuller explanation of purpose,
- an arguments section,
- a returns section,
- notes on conventions, invariants, and edge cases,
- related functions when useful,
- and a small example when it genuinely clarifies usage.

For internal helpers, docstrings should still be substantive, but examples are optional. Internal docstrings should prioritize:

- the helper's role in a larger algorithm,
- required preconditions on mutable buffers or basis state,
- normalization or indexing conventions,
- and the reason it exists as a separate helper.

The guiding principle is that an AI agent should be able to infer intent from the docstring without having to reconstruct all semantics from implementation details.

## Coverage Targets

The overhaul should cover:

- exported types and functions,
- non-exported but important algorithmic helpers,
- internal abstractions that encode package conventions,
- matrix/basis conversion utilities,
- threading and helper routines whose behavior is non-obvious,
- and important structs/constants involved in Lindblad, Pauli, Schmidt, and basis workflows.

The work should focus on files in `src/` that define the package's actual behavior. Development scratch files in `src/dev/` are lower priority and may be treated separately if time permits.

## Documentation Architecture

The codebase is already organized roughly by subsystem. The rewrite should follow that structure:

1. basis layer
2. operator and mapping layer
3. entanglement layer
4. ITensor/Pauli layer
5. algorithm layer
6. top-level orientation materials such as README and docs pages

This subsystem-first organization makes it easier to keep terminology consistent within a cluster of related functions before moving to the next cluster.

## Terminology Alignment

The rewrite should consistently explain recurring concepts using stable vocabulary:

- digit representation
- local base / on-site dimension
- basis representative
- orbit / normalization factor
- symmetry sector
- tensor-product embedding
- matrix-free application
- explicit dense or sparse representation
- covariance-matrix vs many-body Lindblad evolution

This matters because the package has many interrelated abstractions, and inconsistent wording makes both human learning and agent reasoning harder.

## README Upgrade Direction

The README should be extended to include:

- a short architecture map,
- a section explicitly titled for repository navigation or agent orientation,
- file and subsystem pointers,
- recommended entry points for common tasks,
- and stronger explanation of which abstractions matter most.

It should still contain installation and quick-start material, but the upgraded emphasis should be orientation and guidance.

## Documentation Site Upgrade Direction

The docs site should be expanded and sharpened rather than radically restructured.

Likely upgrades include:

- more detailed explanations on manual pages,
- clearer cross-links between conceptual pages and API reference pages,
- explicit explanations of conventions used by source docstrings,
- and a developer page that explains the intended split between README, docstrings, and human-facing docs.

## Risks

### Risk: Inconsistency across many files

Mitigation:

- rewrite by subsystem,
- keep a clear docstring template,
- and use recurring wording intentionally.

### Risk: Overly verbose or repetitive docs

Mitigation:

- let source docstrings carry detailed semantics,
- let README carry orientation,
- let docs pages carry pedagogy,
- and avoid copying identical paragraphs between layers.

### Risk: Internal helper documentation becoming noisy

Mitigation:

- keep helper docstrings focused on role, assumptions, and return semantics,
- not full tutorial prose.

## Success Criteria

The overhaul is successful when:

- most important types, functions, and helpers in `src/` have explicit high-quality docstrings,
- docstrings document arguments and return semantics in detail,
- README serves as a strong orientation guide for AI agents and repository newcomers,
- the Documenter manual is clearer and more detailed for human readers,
- and the docs site still builds successfully after the rewrite.
