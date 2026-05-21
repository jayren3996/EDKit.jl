<p align="center">
  <img src="docs/src/assets/logo.svg" width="200" alt="EDKit.jl logo">
</p>

<h1 align="center">EDKit.jl</h1>

<p align="center">
  <em>Exact diagonalization for symmetry-resolved many-body quantum systems in Julia.</em>
</p>

<p align="center">
  <a href="https://jayren3996.github.io/EDKit.jl/"><img src="https://img.shields.io/badge/docs-live-blue.svg" alt="Documentation"></a>
  <a href="https://github.com/jayren3996/EDKit.jl/actions/workflows/documentation.yml"><img src="https://github.com/jayren3996/EDKit.jl/actions/workflows/documentation.yml/badge.svg" alt="Docs build"></a>
  <a href="https://julialang.org"><img src="https://img.shields.io/badge/julia-1.10%2B-blueviolet.svg" alt="Julia 1.10+"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-yellow.svg" alt="License: MIT"></a>
</p>

---

## 60-Second Demo

Build a Heisenberg chain Hamiltonian once; reuse the same description in the full Hilbert space and in a symmetry sector:

```julia
using EDKit, LinearAlgebra

L  = 12
h2 = spin((1.0, "xx"), (1.0, "yy"), (1.0, "zz"))

H_full = trans_inv_operator(h2, 1:2, L)              # 4096 states
H_red  = trans_inv_operator(h2, 1:2,
            basis(L = L, N = L ÷ 2, k = 0, p = 1))   # 50 states, same GS

E_full = eigvals(Hermitian(Array(H_full)))[1]
E_red  = eigvals(Hermitian(Array(H_red )))[1]
@show E_full, E_red, abs(E_full - E_red)
# (E_full, E_red, abs(E_full - E_red)) = (-5.387391, -5.387391, 1.15e-14)
```

The same `Operator` survives the change of basis: pick the sector that fits the problem, the matrix shrinks (here 82×), the model description stays unchanged.

## What You Can Do With EDKit

- **Many-body operators** from local terms via `operator` and `trans_inv_operator`, evaluated matrix-free or as `Array` / `SparseMatrixCSC` on demand
- **Symmetry sectors** in 1D, 2D, and 3D: translation, parity, spin-flip, custom abelian permutations, and combinations
- **Spinless fermions** with automatic left-string Jordan-Wigner via `fermion_operator`
- **Entanglement diagnostics** with `ent_S`, `ent_spec`, and `schmidt`
- **ITensor / MPS** bridge: `vec2mps`, `mps2vec`, `op2mat`, `pauli`, `tebd4`
- **Open-system dynamics**: adaptive Krylov `lindblad`, quadratic-fermion `quadraticlindblad`, dense step loops, and quantum-inverse-method solvers
- **Closed-system real-time** evolution via `timeevolve` with adaptive Lanczos

## Choose Your Path

| If you want to … | Start with |
|------------------|------------|
| Build and diagonalize a many-body Hamiltonian | [Getting Started](https://jayren3996.github.io/EDKit.jl/getting-started/) → [Operators](https://jayren3996.github.io/EDKit.jl/manual/operators/) |
| Reduce into symmetry sectors (1D, 2D, 3D) | [Bases and Sectors](https://jayren3996.github.io/EDKit.jl/manual/bases/) → [General Abelian Symmetries](https://jayren3996.github.io/EDKit.jl/abelian_basis/) |
| Run closed-system real-time dynamics | [Time Evolution](https://jayren3996.github.io/EDKit.jl/manual/time-evolution/) |
| Work with MPS or Pauli-space representations | [ITensor Workflows](https://jayren3996.github.io/EDKit.jl/manual/itensors/) |
| Solve open-system or quadratic-fermion problems | [Lindblad Workflows](https://jayren3996.github.io/EDKit.jl/manual/lindblad/) |
| Use spinless fermion operators with auto JW | [Spinless Fermions](https://jayren3996.github.io/EDKit.jl/manual/fermions/) |

## Installation

```julia
pkg> add EDKit
```

Or directly from GitHub:

```julia
pkg> add https://github.com/jayren3996/EDKit.jl
```

Compat:

- Julia `1.10+`
- ITensors.jl `0.7`–`0.9`
- ITensorMPS.jl `0.3`

## Repository Map

| Path | What lives there |
|------|-------------------|
| [`src/EDKit.jl`](src/EDKit.jl) | Top-level module, includes and re-exports everything |
| [`src/Basis/`](src/Basis) | Basis types: momentum, parity, flip, projected, abelian, fermion |
| [`src/Operator.jl`](src/Operator.jl) | `Operator` type, local-term assembly, sparse caching |
| [`src/FermionOperator.jl`](src/FermionOperator.jl) | Spinless fermion helpers with automatic Jordan-Wigner |
| [`src/LinearMap.jl`](src/LinearMap.jl) | `DoubleBasis`, `symmetrizer` — inter-basis maps |
| [`src/Schmidt.jl`](src/Schmidt.jl) | Entanglement, Schmidt decomposition |
| [`src/ITensors/`](src/ITensors) | MPS/MPO bridge, Pauli-space utilities, TEBD |
| [`src/algorithms/`](src/algorithms) | Lindblad, quantum-inverse-method solvers |
| [`docs/`](docs) | Documenter manual, worked examples, API reference |
| [`examples/`](examples) | End-to-end Jupyter notebooks |

## Conventions Worth Knowing

- Product states are digit vectors: `[0, 1, 0, 1]`. Most reduced bases store one canonical representative per orbit.
- `index(b, dgt)` reads the digit buffer into a basis coordinate; `change!(b, i, dgt)` writes basis state `i` back into the buffer.
- `Operator` stores **local terms plus a basis**, not a prebuilt matrix. The same object can be applied matrix-free, converted with `Array(opt)`, or with `sparse(opt)`.
- For repeated multiplications, call `sparse!(opt)` once to switch onto a cached `SparseMatrixCSC` (10–1000× faster). Vector products `opt * v` stay matrix-free regardless.
- Many helpers return zero coefficients instead of throwing when a state lies outside a sector — keeps assembly composable.

For deeper conventions and the full entry-point catalog, see [AGENTS.md](AGENTS.md).

## Examples

End-to-end Jupyter notebooks live under [`examples/`](examples). Three good starting points:

- [`Basic/OperatorConstruction.ipynb`](examples/Basic/OperatorConstruction.ipynb) — build a Hamiltonian, apply it matrix-free, compare dense and sparse forms.
- [`Basic/SymmetryReduction.ipynb`](examples/Basic/SymmetryReduction.ipynb) — full vs. reduced bases, sector recombination.
- [`Lindblad/DissipativeXXChain.ipynb`](examples/Lindblad/DissipativeXXChain.ipynb) — small open-system walkthrough.

See [`examples/README.md`](examples/README.md) for the full index.

## Documentation

- **Live docs:** [https://jayren3996.github.io/EDKit.jl/](https://jayren3996.github.io/EDKit.jl/) — manual, worked examples, API reference.
- **Source docstrings** under `src/` are the most precise semantic reference.
- **[AGENTS.md](AGENTS.md)** carries repository-level orientation for AI agents and code-reading contributors.

## Development

The notebooks under `examples/` default to loading the local source tree via:

```julia
const DEV = true
```

Set `DEV = false` inside an example to run it against an installed package instead. To run the test suite:

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

## License

MIT — see [LICENSE](LICENSE).
