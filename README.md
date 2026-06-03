<div align="center">

<img src="docs/src/assets/logo.svg" alt="EDKit.jl logo" width="170"/>

# EDKit.jl

**Exact diagonalization for symmetry-resolved quantum many-body systems.**

From the full Hilbert space to translation, parity, spin-flip, and custom abelian sectors in 1D, 2D, and 3D.

[![Docs](https://img.shields.io/badge/docs-latest-9558B2.svg)](https://jayren3996.github.io/EDKit.jl/) [![Docs Build](https://github.com/jayren3996/EDKit.jl/actions/workflows/documentation.yml/badge.svg)](https://github.com/jayren3996/EDKit.jl/actions/workflows/documentation.yml) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE) [![Julia](https://img.shields.io/badge/Julia-1.10%2B-389826.svg)](https://julialang.org)

</div>

---

EDKit.jl builds many-body operators from local terms, reduces them into symmetry sectors, and
diagonalizes them for spectra, entanglement, and real- or imaginary-time dynamics. The same
`Operator` works whether you keep the full Hilbert space or restrict to a momentum, parity, or
spin-flip sector, so the model you write down does not change when the basis does. It is written
for researchers who already know the physics and want a precise Julia interface with explicit
conventions.

## ✨ Features

|  |  |
| --- | --- |
| 🧩 **Operators from local terms** | `operator` and `trans_inv_operator` assemble many-body Hamiltonians from small local matrices, applied matrix-free or materialized as `Array` / `SparseMatrixCSC` on demand. |
| 🔁 **Symmetry-resolved bases** | Translation, reflection parity, spin-flip, and custom abelian permutation sectors in 1D, 2D, and 3D, all from one `basis` constructor. |
| ⚛️ **Spinless fermions** | `fermion_operator` inserts the left-string Jordan-Wigner prefix automatically, so `c†ᵢcⱼ` terms read as plainly as spin terms. |
| ✂️ **Entanglement** | `ent_S`, `ent_spec`, and `schmidt` give entropies, entanglement spectra, and Schmidt matrices, with the orbit-weight bookkeeping for reduced bases handled for you. |
| 🔗 **ITensor / MPS bridge** | `vec2mps`, `mps2vec`, `op2mat`, and `tebd4` move states and operators between EDKit vectors and ITensor MPS/MPO objects. |
| 🌊 **Closed & open dynamics** | Adaptive-Lanczos `timeevolve` for unitary evolution, plus the `lindblad` stepper and the quadratic-fermion `quadraticlindblad` solver for open systems. |

## 📦 Installation

```julia
pkg> add EDKit
```

## 🚀 Quick Start

Build a translation-invariant Heisenberg chain and read off its ground-state energy:

```julia
using EDKit, LinearAlgebra

L  = 12
h2 = spin((1.0, "xx"), (1.0, "yy"), (1.0, "zz"))   # one Heisenberg bond

H  = trans_inv_operator(h2, 1:2, L)                # Σᵢ h2 on bond (i, i+1)
E0 = eigvals(Hermitian(Array(H)))[1]
```

## 🧭 Choosing a Basis

| Modeling need | Construct with |
| --- | --- |
| Full Hilbert space, no symmetry | `basis(L=…)` → `TensorBasis` |
| Fixed particle number or magnetization | `basis(L=…, N=…)` → `ProjectedBasis` |
| A local constraint (e.g. Rydberg blockade) | `basis(L=…, f=pred)` |
| A momentum (translation) sector | `basis(L=…, k=…)` |
| Momentum together with reflection parity | `basis(L=…, k=…, p=±1)` |
| A spin-flip (Z₂) sector | `basis(L=…, z=±1)` |
| 2D / 3D lattices or custom abelian symmetries | `basis(L=…, symmetries=[(perm, q), …])` |
| Spinless fermions, optionally fixed-N | `SpinlessFermionBasis(L=…, N=…)` |

Most rows call the same `basis` constructor, which returns whichever concrete type the requested
symmetries need. The operator does not care which one you pick: an `Operator` stores local terms
and a basis, not a prebuilt matrix, so one Hamiltonian description carries across every row of
this table. Change the basis and the matrix shrinks to that sector while the model stays the same.

## 🛠 Usage

<details>
<summary><b>Symmetry-sector reduction</b></summary>

<br>

The same Heisenberg description, diagonalized in the full space and in one momentum/parity sector.
The ground-state energies agree to machine precision while the sector matrix is far smaller (here
4096 → 50 states).

```julia
using EDKit, LinearAlgebra

L  = 12
h2 = spin((1.0, "xx"), (1.0, "yy"), (1.0, "zz"))

H_full = trans_inv_operator(h2, 1:2, L)                # 4096 states
H_red  = trans_inv_operator(h2, 1:2,
            basis(L = L, N = L ÷ 2, k = 0, p = 1))     # 50 states

E_full = eigvals(Hermitian(Array(H_full)))[1]
E_red  = eigvals(Hermitian(Array(H_red )))[1]
abs(E_full - E_red)   # ≈ 1e-14
```

</details>

<details>
<summary><b>Spinless fermions with automatic Jordan-Wigner</b></summary>

<br>

`trans_inv_fermion_operator` builds a translation-invariant fermion Hamiltonian and inserts the
long-way Jordan-Wigner chain on the wrap-around bond. This is the half-filled t–V chain
`H = -t Σᵢ (c†ᵢc_{i+1} + h.c.) + V Σᵢ nᵢn_{i+1}`.

```julia
using EDKit, LinearAlgebra

L, N, t, V = 8, 4, 1.0, 2.0
B = SpinlessFermionBasis(L = L, N = N)

H_hop = trans_inv_fermion_operator("+-", [1, 2], B)
H_int = trans_inv_fermion_operator("nn", [1, 2], B)
H = -t * (H_hop + adjoint(H_hop)) + V * H_int

E0 = eigvals(Hermitian(Array(H)))[1]
```

</details>

<details>
<summary><b>Entanglement entropy across a cut</b></summary>

<br>

`ent_S` takes a state vector, the sites in subsystem `A`, and the basis the vector lives in. It
behaves the same way in a symmetry-reduced basis: EDKit accounts for the orbit weights and
symmetry phases when it assembles the Schmidt matrix.

```julia
using EDKit, LinearAlgebra

B   = basis(L = 10, N = 5, k = 0)
psi = normalize(randn(ComplexF64, size(B, 1)))

S = ent_S(psi, 1:5, B)        # von Neumann entropy of the half-chain cut
```

</details>

<details>
<summary><b>Open-system Lindblad dynamics</b></summary>

<br>

The dense stepper materializes the Hamiltonian and jump operators and advances a density matrix by
one truncated-Taylor step of size `dt`. For many output times use `lindblad_timeevolve` (adaptive
Arnoldi); for quadratic models use `quadraticlindblad` on the covariance matrix.

```julia
using EDKit

H     = zeros(ComplexF64, 2, 2)
jumps = [sqrt(0.4) * ComplexF64[0 1; 0 0]]      # decay at rate 0.4

lb = lindblad(H, jumps)
ρ0 = densitymatrix(ComplexF64[0.0, 1.0])        # start in the excited state
ρ1 = lb(ρ0, 0.01; order = 8)                    # one step of dt = 0.01
```

</details>

## 📚 Documentation

Full documentation lives at **[jayren3996.github.io/EDKit.jl](https://jayren3996.github.io/EDKit.jl/)**.

- [Getting Started](https://jayren3996.github.io/EDKit.jl/getting-started/) — install, build a first operator, diagonalize.
- [Operators](https://jayren3996.github.io/EDKit.jl/manual/operators/) and [Bases and Sectors](https://jayren3996.github.io/EDKit.jl/manual/bases/) — the core assembly and symmetry-reduction APIs.
- [General Abelian Symmetries](https://jayren3996.github.io/EDKit.jl/abelian_basis/) — 2D/3D lattices and custom permutation sectors.
- [Time Evolution](https://jayren3996.github.io/EDKit.jl/manual/time-evolution/) and [Lindblad Workflows](https://jayren3996.github.io/EDKit.jl/manual/lindblad/) — closed- and open-system dynamics.
- [Spinless Fermions](https://jayren3996.github.io/EDKit.jl/manual/fermions/) and [ITensor Workflows](https://jayren3996.github.io/EDKit.jl/manual/itensors/) — Jordan-Wigner operators and MPS conversions.

## 🗂 Examples

Runnable Jupyter notebooks live in [`examples/`](examples/):

- [`Basic/OperatorConstruction.ipynb`](examples/Basic/OperatorConstruction.ipynb) — build a Hamiltonian, apply it matrix-free, compare dense and sparse forms.
- [`Basic/SymmetryReduction.ipynb`](examples/Basic/SymmetryReduction.ipynb) — full versus reduced bases and sector recombination.
- [`Lindblad/DissipativeXXChain.ipynb`](examples/Lindblad/DissipativeXXChain.ipynb) — a small open-system walkthrough.

See [`examples/README.md`](examples/README.md) for the full index, which also covers symmetry catalogues, MPS workflows, and worked physics models. Repository-level orientation for contributors and AI agents lives in [AGENTS.md](AGENTS.md).

## 📄 License

[MIT](LICENSE) © 2021 JieRen and contributors.
