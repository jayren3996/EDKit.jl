# AGENTS.md

## Scope

- Applies to all basis constructors and basis internals in this folder.

## Basis Taxonomy

- `TensorBasis`: full tensor-product basis
- `ProjectedBasis`: constrained or fixed-charge subspace
- dedicated 1D symmetry bases: translational, parity, flip, and combinations
- `AbelianBasis`: generic commuting permutation symmetries
- `basis(...)`: high-level dispatcher across the above

## Invariants

- `change!(b, i, dgt)` decodes stored coordinates into a digit buffer.
- `index(b, dgt)` interprets a digit buffer in basis coordinates and may return
  a zero coefficient when the state lies outside the sector.
- `content(B, i)` returns the stored representative index.
- Reduced bases store canonical representatives of symmetry orbits.
- Normalization factors and representative choice must stay consistent with
  `index`, operator assembly, and Schmidt logic.

## Threading Rule

- Internal hot paths must use explicit external `dgt` buffers.
- Do not reintroduce shared mutable `b.dgt` use inside threaded loops.
- Be careful when editing `@inbounds` sections in `AbstractBasis.jl` and
  orbit-search helpers.

## 2D/3D Guidance

- For arbitrary finite lattices, prefer `basis(...; symmetries=...)`.
- Symmetry generators are permutation-based and must be mutually commuting.
- Only add a specialized 2D basis type if there is a clear ergonomics or
  performance gap relative to `AbelianBasis`.

## Verification Targets

- primary tests: `test/basis_tests.jl`, `test/AbelianBasisTest.jl`,
  `test/abelian_overhaul_tests.jl`
- also relevant: `test/entanglement_tests.jl`, `test/MultiThreadsTest.jl`
