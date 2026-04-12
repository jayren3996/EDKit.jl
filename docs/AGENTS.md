# AGENTS.md

## Scope

- Applies to documentation source and the local docs build workflow.

## Documentation Layers

- `src/` docstrings: semantic reference for agents and advanced readers
- `docs/src/manual/`: concepts and decision guidance
- `docs/src/examples/`: short workflow pages
- `docs/src/reference/`: API surfacing via `@docs`

## Rules

- When adding or moving a page, update `docs/make.jl` navigation.
- Prefer cross-links over copying the same explanation into multiple pages.
- Avoid orphan pages: discoverability matters as much as content quality.
- Agent notes are local workflow guidance; do not copy them into the user docs.

## Build

- docs build:
  `~/.juliaup/bin/julia --project=docs docs/make.jl`
- built site output:
  `docs/build/`
