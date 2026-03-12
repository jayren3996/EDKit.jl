# EDKit.jl

`EDKit.jl` is a Julia package for exact diagonalization and symmetry-resolved many-body calculations, with additional ITensor-based tools for MPS, operator-space evolution, and Lindblad or quadratic-fermion workflows.

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
- Operator-space truncation MPOs:
  `daoe` and `fdaoe`
- Lindblad and quadratic-fermion solvers:
  `lindblad`, `qimsolve`, `quadraticlindblad`, and related helpers

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
- `ITensors.jl` `0.7`
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

### DAOE and fDAOE

The package includes two MPO filters for operator-space truncation:

```julia
D = daoe(ps, 3, 0.4)
FD = fdaoe(ps, 2, 0.3)
```

See the dedicated notebooks in [`examples/DAOE`](examples/DAOE) for small physical benchmarks.

## Examples

The main documentation style in this repository is example-driven.

Start here:

- [`examples/README.md`](examples/README.md)
- [`examples/DAOE/README.md`](examples/DAOE/README.md)

Basic notebooks:

- [`examples/Basic/OperatorConstruction.ipynb`](examples/Basic/OperatorConstruction.ipynb): construct a many-body Hamiltonian, apply it as a linear map, and compare dense and sparse forms.
- [`examples/Basic/SymmetryReduction.ipynb`](examples/Basic/SymmetryReduction.ipynb): build the same model in full and symmetry-reduced bases, then verify sector recombination.
- [`examples/Basic/MPSAndPauli.ipynb`](examples/Basic/MPSAndPauli.ipynb): convert vectors to MPS, move into Pauli-space MPS/MPO form, and inspect bond dimensions.

Notable notebooks and scripts:

- [`examples/Basics.ipynb`](examples/Basics.ipynb)
- [`examples/Tensors.ipynb`](examples/Tensors.ipynb)
- [`examples/Symmetries/MPSProjection.ipynb`](examples/Symmetries/MPSProjection.ipynb)
- [`examples/DAOE/XXZOperatorGrowth.ipynb`](examples/DAOE/XXZOperatorGrowth.ipynb)
- [`examples/DAOE/XXMajoranaGrowth.ipynb`](examples/DAOE/XXMajoranaGrowth.ipynb)
- [`examples/ConstrainedPXP.jl`](examples/ConstrainedPXP.jl)
- [`examples/GapRatioXXZ.jl`](examples/GapRatioXXZ.jl)

## Design Notes

`Operator` is the main abstraction. It lets you keep model construction separate from representation choice:

- stay matrix-free for large sparse calculations,
- convert to `Array` when exact diagonalization is affordable,
- convert to `sparse` when you want explicit sparse storage,
- reuse the same local terms across different bases or symmetry sectors.

That separation is what makes EDKit flexible enough to handle both textbook ED workflows and the operator-space or MPS-based utilities added later.

## Development

The repository examples usually default to loading the local source tree with:

```julia
const DEV = true
```

Set `DEV = false` inside those example files if you want them to run against an installed package instead.
