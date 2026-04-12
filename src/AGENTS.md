# AGENTS.md

## Scope

- Applies to `src/Operator.jl`, `src/LinearMap.jl`, `src/Schmidt.jl`, and
  `src/ToolKit.jl`.
- Read `src/Basis/AGENTS.md` before changing basis semantics.

## Core Model

- Basis objects define coordinates and symmetry sectors.
- `Operator` stores local terms plus basis embedding logic, not just a prebuilt
  matrix.
- Maps and entanglement helpers assume the basis/operator split stays intact.

## Canonical APIs

- `spin`, `operator`, `trans_inv_operator`
- `DoubleBasis`, `symmetrizer`
- `ent_S`, `schmidt`

## Editing Rules

- Prefer preserving public entry points and extending internals behind them.
- When changing operator application semantics, re-check dense, sparse, and
  matrix-free behavior.
- Repeated `Operator * matrix` workflows should still benefit from
  `sparse!(opt)` and `clear_sparse_cache!()`.
- If you touch shared hot paths, verify basis-dependent tests as well as core
  operator tests.

## Verification Targets

- primary tests: `test/core_tests.jl`, `test/advanced_tests.jl`
- also relevant when semantics change: `test/basis_tests.jl`,
  `test/entanglement_tests.jl`
