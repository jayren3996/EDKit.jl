# Tests

This folder contains the current baseline test suite for `EDKit`.

## How to run

From the repository root, run:

```bash
julia --project=. test/runtests.jl
```

This uses the local package environment from the repository `Project.toml`.

## Structure

- `runtests.jl`: main entrypoint
- `TestHelpers.jl`: shared helpers used by multiple test files
- `core_tests.jl`: operators, spin helpers, dense/sparse conversion, multiplication, entropy, and toolkit helpers
- `basis_tests.jl`: projected bases, symmetry sectors, `basis(...)`, and translational sectors including `a > 1`
- `advanced_tests.jl`: linear maps and algorithm-level helpers such as Lindblad and QIM routines
- `itensor_tests.jl`: ITensor conversion utilities and `mps2vec(psi, B)` on symmetry sectors

## Scope

The new suite is intended to be:

- compact enough to run regularly during development
- broad enough to catch regressions in the main user-facing API

It does not yet exhaustively cover every edge case or every large-system behavior. Those can be added later as slower, more targeted tests if needed.
