# AGENTS.md

## Scope

- This folder is organized by subsystem, not by one-file-per-source-file.

## Test Map

- `core_tests.jl`: operators, spin helpers, dense/sparse conversion, toolkit
- `basis_tests.jl`: core basis constructors and symmetry decomposition
- `AbelianBasisTest.jl`, `abelian_overhaul_tests.jl`: Abelian and
  permutation-based symmetry coverage
- `entanglement_tests.jl`: entropy and Schmidt checks across sectors
- `itensor_tests.jl`, `TensorTest.jl`: ITensor, Pauli-space, and MPS coverage
- `lindblad_tests.jl`, `advanced_tests.jl`: algorithms and higher-level
  integration

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
