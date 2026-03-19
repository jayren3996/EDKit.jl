# Developer Notes

This page explains how the documentation is organized and how to extend it.

EDKit deliberately uses different documentation layers for different audiences.

## Build The Docs Locally

From the repository root:

```julia
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The built site is written to `docs/build/`.

## Docs Layout

The documentation source lives under `docs/src/`.

Main sections:

- `index.md` and `getting-started.md` for entry points,
- `manual/` for narrative explanations,
- `examples/` for short workflow pages,
- `reference/` for grouped `@docs` pages.

## Documentation Roles

The repository now uses three complementary documentation layers:

- Source docstrings in `src/`
  This is the most detailed semantic layer. It is written to help advanced
  code readers and AI agents understand argument meaning, return semantics,
  invariants, internal helper roles, and package conventions.
- `README.md`
  This is the repository-orientation layer. It should answer where to start,
  which files own which concepts, and which entry points are canonical.
- The Documenter site in `docs/src/`
  This is the human-learning layer. It should prioritize explanation, examples,
  and conceptual flow.

## Where Content Comes From

The docs combine:

- package overview material from the README,
- existing docstrings from `src/`,
- curated ideas from `examples/`,
- and new manual prose written specifically for the documentation site.

## How To Improve The Docs

When adding or changing functionality:

- update the relevant manual page if the workflow changes,
- add or improve docstrings for both user-facing functions and important
  internal helpers,
- keep examples short enough that they explain one idea clearly,
- prefer cross-links between manual pages and reference pages over repeating the same explanation in many places,
- use README additions sparingly and only when they genuinely improve
  repository-level orientation.

## Deployment

GitHub Actions builds and deploys the site through `Documenter.jl`. The deployment workflow lives in `.github/workflows/documentation.yml`.
