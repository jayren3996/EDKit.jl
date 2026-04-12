# AGENTS.md

## Scope

- EDKit is a Julia exact-diagonalization package built around basis selection
  plus operator assembly.
- This root note covers repo-wide context only.
- Read local `AGENTS.md` files when working in `src/`, `src/Basis/`,
  `src/ITensors/`, `src/algorithms/`, `test/`, or `docs/`.

## Package Map

- `src/Basis/`: basis types, symmetry reduction, and `basis(...)`
- `src/Operator.jl`: many-body operator assembly and application
- `src/LinearMap.jl`: `DoubleBasis`, `symmetrizer`, and overlap maps
- `src/Schmidt.jl`: entanglement and Schmidt decomposition
- `src/ITensors/`: ITensor/MPS/MPO bridge and Pauli-space tools
- `src/algorithms/`: Lindblad and QIM workflows

## Canonical Entry Points

- basis selection: `TensorBasis`, `ProjectedBasis`, `basis(...)`
- operator construction: `spin`, `operator`, `trans_inv_operator`
- maps: `DoubleBasis`, `symmetrizer`
- entanglement: `ent_S`, `schmidt`
- docs for arbitrary lattice symmetries: `docs/src/abelian_basis.md`

## Repo-Wide Invariants

- Basis constructors use keyword arguments.
- `Operator` is immutable; sparse caching lives outside the struct.
- Internal hot paths use explicit `dgt` buffer passing for thread safety.
- If you touch `index`, `change!`, or operator application internals, preserve
  structural validity of existing `@inbounds` loops.
- For arbitrary 2D/3D lattice symmetries, prefer
  `basis(...; symmetries=...)` over adding a new specialized basis by default.

## Commands

- full tests:
  `~/.juliaup/bin/julia --project -e 'using Pkg; Pkg.test()'`
- direct test entrypoint:
  `~/.juliaup/bin/julia --project test/runtests.jl`
- docs build:
  `~/.juliaup/bin/julia --project=docs docs/make.jl`

## Local Note Tree

- `src/AGENTS.md`: core runtime semantics outside basis internals
- `src/Basis/AGENTS.md`: basis invariants and symmetry-sector rules
- `src/ITensors/AGENTS.md`: ITensor bridge conventions
- `src/algorithms/AGENTS.md`: solver-family boundaries
- `test/AGENTS.md`: test map and targeted verification patterns
- `docs/AGENTS.md`: docs-layer boundaries and navigation rules
