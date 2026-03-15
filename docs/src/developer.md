# Developer Notes

This page explains how the documentation is organized and how to extend it.

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

## Where Content Comes From

The docs combine:

- package overview material from the README,
- existing docstrings from `src/`,
- curated ideas from `examples/`,
- and new manual prose written specifically for the documentation site.

## How To Improve The Docs

When adding or changing functionality:

- update the relevant manual page if the workflow changes,
- add or improve docstrings for user-facing functions,
- keep examples short enough that they explain one idea clearly,
- prefer cross-links between manual pages and reference pages over repeating the same explanation in many places.

## Deployment

GitHub Actions builds and deploys the site through `Documenter.jl`. The deployment workflow lives in `.github/workflows/documentation.yml`.
