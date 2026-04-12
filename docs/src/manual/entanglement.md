# Entanglement

EDKit includes two related entanglement layers:

- vector-and-basis tools built around `schmidt`, `ent_spec`, and `ent_S`,
- ITensor/MPS cut entanglement helpers described in [ITensor Workflows](itensors.md).

For state vectors in EDKit bases, the workflow is:

1. use `schmidt` when you want the actual bipartite matrix,
2. use `ent_spec` when you want the Schmidt singular values,
3. use `ent_S` when you want the final entropy directly.

## Entropy From A Vector In A Basis

The main user-facing entry point is `ent_S`.

For a state vector `v`, subsystem `Ainds`, and basis `B`, it computes the entanglement entropy across the chosen bipartition:

```julia
using EDKit, LinearAlgebra

L = 8
B = TensorBasis(L = L, base = 2)
psi = normalize(randn(ComplexF64, 2^L))

S = ent_S(psi, 1:4, B)
```

There is also a convenience overload for full-space vectors when you only know
the system size:

```julia
using EDKit, LinearAlgebra

L = 8
psi = normalize(randn(ComplexF64, 2^L))

S = ent_S(psi, 1:4, L)
```

The same interface works with symmetry-reduced bases:

```julia
B = basis(L = 10, N = 5, k = 0)
psi = normalize(randn(ComplexF64, size(B, 1)))
S = ent_S(psi, 1:5, B)
```

## What Each Function Returns

Use the lower-level helpers when you want more control:

- `schmidt(v, Ainds, B)` returns the Schmidt matrix itself,
- `ent_spec(v, Ainds, B)` returns the singular values of that matrix,
- `ent_S(v, Ainds, B; α=..., cutoff=...)` squares those singular values into
  Schmidt probabilities and applies `EDKit.entropy`.

`EDKit.entropy` also supports direct Renyi calculations:

- `α = 1` gives the von Neumann entropy,
- `α = 0` gives the support-size variant used by the package,
- other `α` values give Renyi entropies.

## Custom Subsystem Bases

By default, `schmidt` interprets the `A` and `B` subsystems in full
tensor-product bases. When you need something else, you can provide subsystem
bases explicitly through `B1` and `B2`.

This is useful when:

- you want to resolve a conserved quantity inside a subsystem,
- you want a constrained Schmidt matrix,
- or you want to post-process the bipartition in a custom reduced basis.

Example:

```julia
using EDKit, LinearAlgebra

L = 8
B = basis(L = L, N = 4, k = 0)
psi = normalize(randn(ComplexF64, size(B, 1)))

BA = ProjectedBasis(L = 4, N = 2)
BB = ProjectedBasis(L = 4, N = 2)
M = EDKit.schmidt(psi, 1:4, B; B1=BA, B2=BB)
```

## Why Symmetry-Aware Entanglement Matters

In a full tensor-product basis, Schmidt decomposition is straightforward. In
symmetry-reduced bases, however, the state coordinates no longer correspond
directly to plain product states.

EDKit handles that bookkeeping for you. The basis object determines how each
representative contributes, including orbit weights and symmetry phases, and
the Schmidt matrix is assembled accordingly.

That makes the entanglement tools especially useful when you diagonalize in momentum, parity, or combined symmetry sectors but still want real-space bipartite diagnostics.

## MPS Entanglement

For MPS objects, EDKit provides bond-cut entanglement helpers in the ITensor
layer:

- `ent_S(psi, b)` computes the entropy across bond `b`,
- `ent_S!(psi, b)` is the in-place version that may move the orthogonality
  center,
- `ent_specs!(psi, b)` returns the singular values at that cut,
- `ent_S(psi)` returns the full entropy profile along the chain.

See [ITensor Workflows](itensors.md) for that side of the API.

## Practical Notes

- `Ainds` lists the sites in subsystem `A`; all other sites belong to subsystem `B`.
- `cutoff` controls how tiny Schmidt values are treated when converting spectra to entropies.
- For symmetry-reduced bases, the state vector should already live in the basis you pass to `ent_S` or `schmidt`.
