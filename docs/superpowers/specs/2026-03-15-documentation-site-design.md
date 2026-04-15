# EDKit Documentation Site Design

## Goal

Build a full `Documenter.jl` documentation site for `EDKit.jl`, deploy it with the standard GitHub Pages workflow, and link the published manual from the repository README.

The site should teach users how to think about the package, explain the purpose of each major subsystem, document important functions with examples, and provide a maintainable API reference generated from the code.

## Why This Matters

`EDKit.jl` already has substantial functionality spread across basis construction, operator assembly, symmetry-aware mappings, entanglement helpers, ITensor tooling, and Lindblad workflows. The current repository has a strong README and example folders, but it does not yet provide a coherent manual that helps new users discover:

- what the package is for,
- which abstractions are central,
- when to use each basis or algorithmic layer,
- and how to translate common research workflows into package calls.

The new docs should close that gap while staying close to the source through docstrings and generated API pages.

## Documentation Strategy

The site will follow a manual-first approach.

This means:

- narrative pages come before raw API listings,
- conceptual explanations are grouped by subsystem,
- examples are curated and shortened into documentation-friendly walkthroughs,
- and a generated API reference is provided as a secondary resource rather than the primary entry point.

This is the best fit for the package because EDKit’s value lies not only in individual functions, but in how basis objects, operators, mappings, and algorithmic helpers work together.

## Site Structure

### Home

The landing page should present:

- a concise description of EDKit,
- installation and compatibility details,
- the main user pathways,
- and quick links into the most useful manual sections.

It should not duplicate the full README verbatim. Instead, it should adapt the strongest parts of the README into a docs-native landing page.

### Getting Started

This page should teach the core workflow:

1. define local operators,
2. specify where they act,
3. choose a basis,
4. build an `Operator`,
5. apply it as a linear map or convert it into explicit matrix forms.

It should include one compact example, ideally a small spin-chain Hamiltonian, and give users a first successful workflow without requiring deep prior knowledge of the package.

### Manual

The manual should be organized by subsystem:

- `architecture`: the package mental model and how the major parts relate
- `bases`: full-space and symmetry-reduced bases, plus when to use each
- `operators`: operator construction, translation-invariant helpers, and matrix conversions
- `maps`: `DoubleBasis`, inter-basis maps, and symmetrization workflows
- `entanglement`: Schmidt and entropy tools
- `itensors`: vector/MPS/MPO conversion, Pauli-space helpers, and TEBD tooling
- `lindblad`: many-body and quadratic open-system workflows
- `utilities`: smaller helper functions worth surfacing to users

Each page should explain both purpose and practice:

- what this subsystem is for,
- which objects or functions matter most,
- typical inputs and outputs,
- a small example,
- and links to the relevant API reference.

### Worked Examples

Examples should be curated into a few pages rather than exposing the raw notebook tree directly. The worked example pages should adapt ideas from `examples/` into shorter walkthroughs that are easier to read in a browser.

Planned example pages:

- `basic-workflows`
- `symmetry-workflows`
- `tensor-network-workflows`
- `open-system-workflows`

These pages should complement the manual, not replace it.

### API Reference

The API reference should be generated with `@docs` and grouped by subsystem. It should provide discoverability and exact signatures while leaving long explanations to the manual pages.

Reference pages should include short orientation text at the top so users understand what each group of functions is for.

### Developer Notes

A small developer page should explain:

- how the docs site is built,
- where source pages live,
- how the navigation is organized,
- and how to improve docstrings so the API reference remains useful.

## Content Sources

The docs should be assembled from four sources:

1. `README.md` for package overview and initial workflow framing
2. existing docstrings in `src/`
3. current examples in `examples/`
4. new explanatory prose written specifically for the manual

The new prose is necessary because the existing materials are informative but distributed across code, examples, and the README. The docs site should unify those sources into a clear teaching path.

## Information Architecture Principles

The site should follow these principles:

- Concept first, API second
- Prefer short runnable examples over long notebook dumps
- Explain why a user would choose one basis or workflow over another
- Keep pages scoped so each has one clear teaching goal
- Cross-link manual pages and reference pages heavily
- Preserve accuracy by reusing docstrings where possible

## Deployment Design

The documentation will use the standard `Documenter.jl` GitHub deployment model:

- a dedicated `docs/Project.toml`
- a `docs/make.jl` build script
- a GitHub Actions workflow for build and deploy
- publishing to the repository’s `gh-pages` branch

This should run on pushes to the default branch and on pull requests for build validation.

## README Integration

The README should be updated to include a prominent documentation link, ideally near the top and in the examples or quick-start area, so repository visitors can immediately find the full manual.

## Non-Goals For This First Version

To keep the first release reliable and maintainable, the following should be excluded:

- executing every notebook as part of the docs build,
- building a notebook-to-docs pipeline,
- exhaustive documentation for every internal helper,
- broad refactors unrelated to documentation quality.

The first version should focus on producing a polished, coherent, and deployable manual for the current package version.

## Risks And Mitigations

### Risk: Missing or uneven docstrings

Mitigation:

- write new narrative manual pages where docstrings are thin,
- improve targeted docstrings only where it materially helps generated reference pages.

### Risk: Docs build instability from heavy examples

Mitigation:

- use compact examples in Markdown pages,
- avoid live execution of large or dependency-heavy notebooks in the initial version.

### Risk: Confusing package organization

Mitigation:

- explicitly explain that EDKit is a single Julia module with several functional subsystems,
- describe those subsystems in the architecture page.

## Success Criteria

The work is successful when:

- the repository contains a complete `Documenter.jl` site structure,
- the manual explains the main EDKit subsystems with examples,
- the API reference is organized and generated from code,
- GitHub Actions can build and deploy the docs,
- and the README links to the published documentation.
