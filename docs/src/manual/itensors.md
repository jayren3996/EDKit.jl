# ITensor Workflows

EDKit includes a practical bridge between exact-diagonalization style data structures and ITensor-based tensor-network workflows.

This layer is useful when you want to:

- convert state vectors into MPS form,
- convert dense operators into ITensor operators,
- move between ordinary MPS/MPO objects and Pauli-space representations,
- build TEBD gates from local matrices.

## Vector And MPS Conversion

The simplest entry points are `vec2mps` and `mps2vec`.

```julia
using EDKit, ITensors, ITensorMPS, LinearAlgebra

L = 6
s = siteinds(2, L)
psi = normalize(randn(ComplexF64, 2^L))

psi_mps = vec2mps(psi, s)
psi_back = mps2vec(psi_mps)
```

`mps2vec(psi, B)` also supports extracting amplitudes in a symmetry-aware basis `B` when the MPS already lies in that sector.

That last point matters: `mps2vec(psi, B)` is exact when the MPS already lies
in the symmetry sector described by `B`. If you need a true projection of a
general MPS into a sector, use an explicit basis map instead of treating this as
an orthogonal projector.

## Operator Conversion

`mat2op` and `op2mat` convert between dense matrices and ITensor operator objects:

```julia
using EDKit, ITensors

L = 4
s = siteinds(2, L)
h = spin((1.0, "xx"), (1.0, "yy"), (1.0, "zz"))

O = mat2op(h, s[1], s[2])
h_back = op2mat(O, s[1], s[2])
```

This is useful when you already have a local matrix from the EDKit side and want to feed it into an ITensor workflow.

## MPS State Constructors And Entanglement

EDKit also provides two small MPS-side conveniences:

- `productstate(sites, states)` builds a bond-dimension-1 product MPS from local state vectors,
- `pbcmps(sites, tensors)` builds an MPS from a periodic-boundary-style tensor list.

For entanglement diagnostics on MPS objects:

- `ent_S(psi, b)` computes the entropy across bond `b`,
- `ent_S!(psi, b)` is the in-place variant,
- `ent_specs!(psi, b)` returns the singular values at that cut,
- `ent_S(psi)` returns the entropy profile along the chain.

These are the tensor-network-side companions to the vector-based entanglement
tools described in [Entanglement](entanglement.md).

## Pauli-Space Tools

EDKit also supports operator-space calculations in a Pauli basis.

The main pieces are:

- `pauli` and `pauli_list` for converting between matrices and Pauli coefficients,
- `commutation_mat` for the superoperator `-i[H, \cdot]`,
- `dissipation_mat` for the dissipator `D[L]`,
- `mps2pmps`, `pmps2mpo`, and `mpo2pmpo` for Pauli-space tensor-network conversions.

This is particularly useful for operator growth and Lindbladian tensor-network workflows.

## TEBD Gates

`tebd4` builds fourth-order Trotter gate sequences from a list of local Hamiltonian pieces:

```julia
using EDKit, ITensors

L = 5
ps = siteinds("Pauli", L)
h2 = commutation_mat(spin((1.0, "xx"), (1.0, "yy"), (1.0, "zz")))
gates = tebd4(fill(h2, L - 1), ps, 0.05)
```

The lower-level `tebd_n!` routine applies an n-site TEBD update directly to an MPS.

## When To Use This Layer

Use the ITensor layer when exact diagonalization and tensor networks need to meet. If you only need small-system dense or sparse exact diagonalization, you can stay entirely inside the basis and operator layers.
