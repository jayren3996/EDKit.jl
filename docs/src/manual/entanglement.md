# Entanglement

EDKit includes tools for constructing Schmidt decompositions and entanglement entropies directly from state vectors in either full or symmetry-aware bases.

## Core Functions

The main user-facing entry point is `ent_S`.

For a state vector `v`, subsystem `Ainds`, and basis `B`, it computes the entanglement entropy across the chosen bipartition:

```julia
using EDKit, LinearAlgebra

L = 8
B = TensorBasis(L = L, base = 2)
psi = normalize(randn(ComplexF64, 2^L))

S = ent_S(psi, 1:4, B)
```

This works with symmetry-reduced bases too:

```julia
B = basis(L = 10, N = 5, k = 0)
psi = normalize(randn(ComplexF64, size(B, 1)))
S = ent_S(psi, 1:5, B)
```

## Schmidt Matrices

When you need more than the entropy, `schmidt` constructs the Schmidt matrix associated with a bipartition.

That is useful when you want:

- singular values directly,
- Schmidt vectors,
- or a custom post-processing step beyond the built-in entropy functions.

The helper `ent_spec` returns the singular values of the Schmidt matrix, and `entropy` converts squared singular values into von Neumann or Renyi entropies.

## Why Basis-Aware Entanglement Is Useful

In a full tensor-product basis, Schmidt decomposition is straightforward. In symmetry-reduced bases, however, the state coordinates no longer correspond directly to plain product states.

EDKit handles that bookkeeping for you. The basis object determines how each representative contributes, and the Schmidt matrix is assembled accordingly.

That makes the entanglement tools especially useful when you diagonalize in momentum, parity, or combined symmetry sectors but still want real-space bipartite diagnostics.

## Practical Notes

- `Ainds` lists the sites in subsystem `A`; all other sites belong to subsystem `B`.
- `α = 1` gives the von Neumann entropy.
- other `α` values give Renyi entropies.
- `cutoff` controls how tiny Schmidt values are handled in entropy calculations.
