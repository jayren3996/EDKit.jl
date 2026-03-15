# Tensor-Network Workflows

These examples show how EDKit can connect exact-diagonalization-style objects with ITensor workflows.

## Convert A Vector To MPS And Back

```julia
using EDKit, ITensors, ITensorMPS, LinearAlgebra

L = 6
s = siteinds(2, L)
psi = normalize(randn(ComplexF64, 2^L))

psi_mps = vec2mps(psi, s)
psi_back = mps2vec(psi_mps)
```

## Build Pauli-Space Gates

```julia
using EDKit, ITensors

L = 5
ps = siteinds("Pauli", L)
h2 = commutation_mat(spin((1.0, "xx"), (1.0, "yy"), (1.0, "zz")))
gates = tebd4(fill(h2, L - 1), ps, 0.05)
```

This is a good starting point for operator-space MPS evolution.

## Convert Between MPO And Pauli MPO

```julia
using EDKit, ITensors, ITensorMPS

L = 4
s = siteinds(2, L)
S = siteinds("Pauli", L)

h2 = spin((1.0, "xx"), (1.0, "yy"), (1.0, "zz"))
local_op = mat2op(h2, s[1], s[2])
```

In full workflows, `mps2pmps`, `pmps2mpo`, and `mpo2pmpo` let you move between standard tensor-network operators and their Pauli-space representations.
