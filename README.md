# EDKit.jl

[![Documentation](https://img.shields.io/badge/docs-live-blue.svg)](https://jayren3996.github.io/EDKit.jl/)

`EDKit.jl` is a Julia package for exact diagonalization and symmetry-resolved many-body calculations, with additional ITensor-based tools for MPS, Pauli-space operator representations, and Lindblad or quadratic-fermion workflows.

This README is intentionally written as a repository-orientation guide as much
as a package overview. It should help both human contributors and AI agents
answer:

- what the package is fundamentally for,
- which abstractions are central,
- where the implementation of each subsystem lives,
- and which files or APIs are the right entry points for a given task.

The package is built around one core idea:

- describe local terms once,
- choose a basis or symmetry sector,
- construct an `Operator`,
- use it as a linear map, a dense matrix, or a sparse matrix.

## What EDKit covers

- Generic many-body operator construction with `operator` and `trans_inv_operator`
- Tensor-product and symmetry-reduced bases:
  `TensorBasis`, `TranslationalBasis`, `TranslationParityBasis`, `ParityBasis`, `FlipBasis`, `ProjectedBasis`, and the high-level `basis(...)` helper
- Entanglement utilities such as `ent_S`
- Pauli-basis and operator-space utilities:
  `pauli`, `pauli_list`, `commutation_mat`, `dissipation_mat`
- ITensor/MPS helpers:
  `vec2mps`, `mps2vec`, `mat2op`, `op2mat`, `mps2pmps`, `pmps2mpo`, `mpo2pmpo`, `tebd4`
- Lindblad and quadratic-fermion solvers:
  `lindblad`, `qimsolve`, `quadraticlindblad`, and related helpers

## Repository Orientation

If you are an AI agent or a new contributor, the most important structural fact
about EDKit is that it is a single Julia module with several clearly separated
subsystems.

### Main source layout

- `src/EDKit.jl`
  Top-level module that includes and re-exports everything.
- `src/Basis/`
  Basis representations and symmetry reduction logic. Start here for questions
  about digits, representatives, momentum sectors, parity, spin flip, or the
  high-level `basis(...)` helper.
- `src/Operator.jl`
  The core many-body `Operator` abstraction, local term assembly, matrix-free
  application, and symbolic spin-string helpers.
- `src/LinearMap.jl`
  Inter-basis maps such as `DoubleBasis` and `symmetrizer`.
- `src/Schmidt.jl`
  Entanglement and Schmidt-decomposition logic across full and reduced bases.
- `src/ToolKit.jl`
  Small numerical helpers such as gap-ratio statistics and convenience
  exponentials.
- `src/ITensors/`
  Vector/MPS/MPO conversion, Pauli-basis utilities, and TEBD support.
- `src/algorithms/`
  Lindblad and quantum-inverse-method workflows.

### Canonical entry points

For most tasks, prefer these APIs first:

- Basis selection:
  `TensorBasis`, `ProjectedBasis`, `TranslationalBasis`, `ParityBasis`,
  `FlipBasis`, `ParityFlipBasis`, `TranslationParityBasis`,
  `TranslationFlipBasis`, `basis(...)`
- Operator construction:
  `spin`, `operator`, `trans_inv_operator`
- Basis-to-basis maps:
  `DoubleBasis`, `symmetrizer`
- Entanglement:
  `ent_S`, `ent_spec`, `schmidt`
- ITensor bridge:
  `vec2mps`, `mps2vec`, `mat2op`, `op2mat`, `pauli`, `pauli_list`, `tebd4`
- Open-system and inverse-method tools:
  `lindblad`, `densitymatrix`, `quadraticlindblad`, `covariancematrix`,
  `majoranaform`, `fermioncorrelation`, `covmat`, `qimsolve`

### Package conventions that matter

- Product states are represented internally by digit vectors such as
  `[0, 1, 0, 1]`.
- Most reduced bases store one canonical representative per symmetry orbit.
- `index(...)` usually means “interpret the current digit buffer in basis
  coordinates.”
- `change!(...)` usually means “load basis state `i` into the mutable digit
  buffer.”
- `Operator` stores local terms plus the basis embedding logic, not just one
  prebuilt matrix.
- Many helper functions intentionally return zero coefficients instead of
  throwing when a state lies outside a symmetry sector, because this keeps
  matrix assembly and projections composable.

### How to navigate the docs

- Source docstrings in `src/` are the most precise semantic reference and are
  written to be highly useful for code-reading agents.
- The documentation site at [https://jayren3996.github.io/EDKit.jl/](https://jayren3996.github.io/EDKit.jl/) is more human-oriented and more tutorial/manual in style.
- The examples under `examples/` are useful for end-to-end workflows, but they
  are not the best place to learn package architecture.

## Installation

From the Julia package manager:

```julia
pkg> add EDKit
```

Or install the GitHub version directly:

```julia
pkg> add https://github.com/jayren3996/EDKit.jl
```

Current compat:

- Julia `1.9+`
- `ITensors.jl` `0.7` to `0.9`
- `ITensorMPS.jl` `0.3`

## Quick Start

The simplest workflow is:

1. define local operators,
2. specify where they act,
3. choose a basis,
4. construct the many-body operator.

Example: open XXZ chain on `L = 8` sites.

```julia
using EDKit, LinearAlgebra

L = 8
Delta = 1.5

mats = [
    fill(spin("XX"), L - 1);
    fill(spin("YY"), L - 1);
    fill(Delta * spin("ZZ"), L - 1);
]

inds = vcat(
    [[i, i + 1] for i in 1:L-1],
    [[i, i + 1] for i in 1:L-1],
    [[i, i + 1] for i in 1:L-1],
)

H = operator(mats, inds, L)
```

`H` is an `Operator`, not just a matrix. You can:

```julia
psi = normalize(randn(ComplexF64, 2^L))
Hpsi = H * psi

Hdense = Array(H)
Hsparse = sparse(H)

vals = eigvals(Hermitian(Hdense))
```

The key design idea is that operator definition and representation choice are
separate. The same `Operator` can remain matrix-free, be turned into a dense
matrix, or be turned into a sparse matrix depending on what the next algorithm
needs.

## Symmetry Sectors

EDKit becomes most useful when you stop working in the full Hilbert space.

Example: momentum-zero, half-filling sector of a spin-1/2 chain:

```julia
using EDKit

L = 12
B = basis(L = L, N = L ÷ 2, k = 0)

h2 = spin((1.0, "xx"), (1.0, "yy"), (1.0, "zz"))
H = trans_inv_operator(h2, 1:2, B)
```

You can then diagonalize the reduced matrix:

```julia
vals, vecs = eigen(Hermitian(Array(H)))
```

Available symmetry-oriented basis types include:

- `TranslationalBasis`
- `TranslationParityBasis`
- `ParityBasis`
- `FlipBasis`
- `ProjectedBasis`
- `AbelianBasis`
- `basis(...)` for common combinations of quantum numbers

## ITensor and Operator-Space Tools

EDKit also includes a second layer of tooling for tensor-network workflows.

### MPS conversion

```julia
using EDKit, ITensors, ITensorMPS, LinearAlgebra

L = 6
s = siteinds(2, L)
psi = randn(ComplexF64, 2^L) |> normalize

psi_mps = vec2mps(psi, s)
psi_back = mps2vec(psi_mps)
```

### Pauli-space evolution

For operator growth or Lindbladian calculations, EDKit can work in a Pauli basis:

```julia
using EDKit, ITensors, ITensorMPS

L = 5
ps = siteinds("Pauli", L)
h2 = spin((1.0, "xx"), (1.0, "yy"), (1.0, "zz"))
superop = commutation_mat(h2)
gates = tebd4(fill(superop, L - 1), ps, 0.05)

O = productMPS(ps, ["X", fill("I", L - 1)...])
O = apply(gates, O)
normalize!(O)
```

## Examples

The main documentation style in this repository is example-driven.

Start here:

- [`examples/README.md`](examples/README.md)
- [`examples/Basic/README.md`](examples/Basic/README.md)
- [`examples/Maps/README.md`](examples/Maps/README.md)
- [`examples/Symmetries/README.md`](examples/Symmetries/README.md)
- [`examples/TensorNetworks/README.md`](examples/TensorNetworks/README.md)
- [`examples/Lindblad/README.md`](examples/Lindblad/README.md)
- [`examples/Models/README.md`](examples/Models/README.md)

Basic notebooks:

- [`examples/Basic/OperatorConstruction.ipynb`](examples/Basic/OperatorConstruction.ipynb): construct a many-body Hamiltonian, apply it as a linear map, and compare dense and sparse forms.
- [`examples/Basic/SymmetryReduction.ipynb`](examples/Basic/SymmetryReduction.ipynb): build the same model in full and symmetry-reduced bases, then verify sector recombination.
- [`examples/Basic/MPSAndPauli.ipynb`](examples/Basic/MPSAndPauli.ipynb): convert vectors to MPS, move into Pauli-space MPS/MPO form, and inspect bond dimensions.

Maps notebooks:

- [`examples/Maps/DoubleBasisBasics.ipynb`](examples/Maps/DoubleBasisBasics.ipynb): introduce `DoubleBasis(Btarget, Bsource)`, map direction, and the agreement between `T(v)` and `symmetrizer(T) * v`.
- [`examples/Maps/Symmetrizers.ipynb`](examples/Maps/Symmetrizers.ipynb): project full-space states into parity or Abelian symmetry sectors and lift them back to the full Hilbert space.
- [`examples/Maps/InterbasisOperators.ipynb`](examples/Maps/InterbasisOperators.ipynb): build non-square operators with `DoubleBasis` and compare them with explicit symmetrization.

Lindblad notebooks:

- [`examples/Lindblad/DissipativeXXChain.ipynb`](examples/Lindblad/DissipativeXXChain.ipynb): small many-body open-system walkthrough for a two-site XX chain with local loss.
- [`examples/Lindblad/QuadraticXXLoss.ipynb`](examples/Lindblad/QuadraticXXLoss.ipynb): free-fermion loss benchmark comparing `quadraticlindblad` with full many-body evolution.
- [`examples/Lindblad/PauliSuperoperators.ipynb`](examples/Lindblad/PauliSuperoperators.ipynb): Pauli-basis walkthrough for `pauli`, `pauli_list`, `commutation_mat`, and `dissipation_mat`.

Other topic folders:

- [`examples/Symmetries/MPSProjection.ipynb`](examples/Symmetries/MPSProjection.ipynb)
- [`examples/Symmetries/SectorCatalogue.ipynb`](examples/Symmetries/SectorCatalogue.ipynb)
- [`examples/TensorNetworks/Tensors.ipynb`](examples/TensorNetworks/Tensors.ipynb)
- [`examples/TensorNetworks/MPSBasics.ipynb`](examples/TensorNetworks/MPSBasics.ipynb)
- [`examples/Models/PXPDrive.ipynb`](examples/Models/PXPDrive.ipynb)
- [`examples/Models/ConstrainedPXP.ipynb`](examples/Models/ConstrainedPXP.ipynb)
- [`examples/Models/GapRatioXXZ.ipynb`](examples/Models/GapRatioXXZ.ipynb)

## Design Notes

`Operator` is the main abstraction. It lets you keep model construction separate from representation choice:

- stay matrix-free for large sparse calculations,
- convert to `Array` when exact diagonalization is affordable,
- convert to `sparse` when you want explicit sparse storage,
- reuse the same local terms across different bases or symmetry sectors.

That separation is what makes EDKit flexible enough to handle both textbook ED workflows and the symmetry or MPS-based utilities added later.

## Documentation Layers

EDKit now uses three documentation layers on purpose:

- Source docstrings in `src/`
  These are the most detailed semantic descriptions of functions and types and
  are intended to help code readers and AI agents reason correctly about
  arguments, return values, conventions, and invariants.
- `README.md`
  This file is the repository-orientation layer: architecture map, entry points,
  conventions, and where to look next.
- Documenter manual
  This is the human-learning layer with more narrative explanations and guided
  examples.

## Development

The repository examples usually default to loading the local source tree with:

```julia
const DEV = true
```

Set `DEV = false` inside those example files if you want them to run against an installed package instead.
