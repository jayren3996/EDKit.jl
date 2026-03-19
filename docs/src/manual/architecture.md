# Architecture

`EDKit` is a single Julia module with several internal subsystems that work together.

The important thing for users is not the file layout, but the workflow boundary between these subsystems.

For contributors and AI agents, the file layout is still important as a
navigation tool, so this page also points back to the implementation layers
where each concept lives.

## Mental Model

A typical EDKit workflow looks like this:

1. choose a basis,
2. build a many-body operator from local terms,
3. act on vectors or convert to explicit matrices,
4. optionally move between bases or compute diagnostics,
5. optionally switch into tensor-network or Lindblad-specific tooling.

The package is designed so that the basis and the operator representation stay separate. You do not hard-code a Hamiltonian as one huge matrix first and then try to add symmetry later. Instead, the basis decides how local actions are embedded into the state space.

## Main Subsystems

### Bases

Basis objects describe the state space in which operators act.

EDKit includes:

- `TensorBasis` for the full tensor-product Hilbert space,
- `ProjectedBasis` for constrained or fixed-charge subspaces,
- translational, parity, flip, and combined symmetry bases,
- `basis(...)` as the high-level constructor for common symmetry combinations.

These are covered in [Bases and Sectors](bases.md).

Implementation files:

- `src/Basis/AbstractBasis.jl`
- `src/Basis/ProjectedBasis.jl`
- `src/Basis/TranslationalBasis.jl`
- `src/Basis/ParityBasis.jl`
- `src/Basis/FlipBasis.jl`
- `src/Basis/ParityFlipBasis.jl`
- `src/Basis/TranslationalParityBasis.jl`
- `src/Basis/TranslationalFlipBasis.jl`
- `src/Basis/AbelianBasis.jl`

### Operators

The `Operator` abstraction stores a many-body operator as local matrices plus the sites they act on, together with a basis object. This lets one operator description serve several numerical roles.

These are covered in [Operators](operators.md).

Implementation files:

- `src/Operator.jl`

### Maps

`DoubleBasis` and `symmetrizer` let you build linear maps between related bases, especially between full-space and symmetry-reduced descriptions.

These are covered in [Maps and Symmetrizers](maps.md).

Implementation files:

- `src/LinearMap.jl`

### Entanglement

The Schmidt and entropy tools work on both full and symmetry-aware bases. They provide a practical bridge between state vectors and bipartite diagnostics.

These are covered in [Entanglement](entanglement.md).

Implementation files:

- `src/Schmidt.jl`

### ITensor Tools

EDKit also ships utilities for moving between vectors, operators, MPS, MPO, and Pauli-space representations. These are not the entry point for most users, but they become valuable once exact-diagonalization and tensor-network workflows need to meet.

These are covered in [ITensor Workflows](itensors.md).

Implementation files:

- `src/ITensors/ITensorsKit.jl`
- `src/ITensors/PauliBasis.jl`
- `src/ITensors/TEBD.jl`

### Lindblad And Quadratic Solvers

The package includes both many-body Lindblad evolution and quadratic-fermion covariance-matrix workflows. They solve related physical problems with different scaling and assumptions.

These are covered in [Lindblad Workflows](lindblad.md).

Implementation files:

- `src/algorithms/Lindblad.jl`
- `src/algorithms/QIM.jl`

## How The Pieces Fit

The architecture is easiest to understand through examples:

- use `basis(...)` or a concrete basis type to define the state space,
- use `operator` or `trans_inv_operator` to construct the model,
- use `H * psi`, `mul(H, psi)`, `Array(H)`, or `sparse(H)` depending on the calculation,
- use `DoubleBasis` when you need overlap or projection maps between bases,
- use `ent_S` or `schmidt` when you need bipartite diagnostics,
- move into ITensor or Lindblad workflows only when the problem demands them.

This separation is what makes EDKit flexible: model specification stays stable while representation choice changes with the numerical task.
