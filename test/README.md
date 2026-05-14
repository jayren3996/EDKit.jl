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
- `core_tests.jl`: operators, spin helpers, dense/sparse conversion, multiplication variants, sparse cache behavior, entropy, and toolkit helpers
- `basis_tests.jl`: projected bases, symmetry sectors, `basis(...)`, threaded constructor parity, and translational sectors including `a > 1`
- `entanglement_tests.jl`: entropy and Schmidt checks across tensor, translational, parity, flip, and Abelian sectors
- `advanced_tests.jl`: linear maps and algorithm-level helpers such as Lindblad and QIM routines
- `lindblad_tests.jl`: Lindblad and quadratic Lindblad regression coverage
- `itensor_tests.jl`: ITensor conversion utilities and `mps2vec(psi, B)` on symmetry sectors
- `abelian_overhaul_tests.jl`: Abelian operator internals, custom symmetries, two-dimensional lattice sectors, and small performance sanity checks
- `timeevolve_tests.jl`: adaptive Krylov time evolution
- `lindblad_timeevolve_tests.jl`: adaptive Lindblad Arnoldi time evolution

Legacy one-off tests with CamelCase names are not included by `runtests.jl`.
`MultiThreadsTest.jl` is a standalone threaded-basis check; by default it uses a
small CI-safe PXP system, and the original `L = 28` case runs only with:

```bash
EDKIT_SLOW_TESTS=1 julia --project=. test/MultiThreadsTest.jl
```

## Scope

The new suite is intended to be:

- compact enough to run regularly during development
- broad enough to catch regressions in the main user-facing API

It does not yet exhaustively cover every edge case or every large-system behavior. Those can be added later as slower, more targeted tests if needed.
