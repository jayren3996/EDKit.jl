# AGENTS.md

## Scope

- This folder is organized by subsystem, not by one-file-per-source-file.

## Test Map

- `core_tests.jl`: operators, spin helpers, dense/sparse conversion, toolkit
- `basis_tests.jl`: core basis constructors, threaded constructor parity,
  empty-sector behavior, and symmetry decomposition
- `abelian_overhaul_tests.jl`: Abelian internals, permutation-based symmetry
  coverage, custom symmetries, and two-dimensional lattice sectors
- `entanglement_tests.jl`: entropy and Schmidt checks across sectors
- `advanced_tests.jl`: linear maps, symmetrizers, and higher-level integration
- `lindblad_tests.jl`: Lindblad and quadratic Lindblad dynamics
- `itensor_tests.jl`: ITensor, Pauli-space, and MPS coverage
- `timeevolve_tests.jl`: adaptive Krylov time evolution
- `lindblad_timeevolve_tests.jl`: adaptive Lindblad Arnoldi time evolution
- `MultiThreadsTest.jl`: standalone threaded PXP check; the default case is
  small, and the original heavy `L = 28` case is gated by
  `EDKIT_SLOW_TESTS=1`

## Preferred Verification

- targeted pattern:
  `~/.juliaup/bin/julia --project -e 'include("test/TestHelpers.jl"); include("test/<file>.jl")'`
- full suite:
  `~/.juliaup/bin/julia --project test/runtests.jl`
- package-level:
  `~/.juliaup/bin/julia --project -e 'using Pkg; Pkg.test()'`

## Rules

- Test public behavior first, internals second.
- Add regression coverage near the subsystem most likely to break again.
- If a change spans basis and operator logic, run both the basis tests and at
  least one operator/integration file.
