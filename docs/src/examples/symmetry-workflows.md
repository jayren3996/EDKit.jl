# Symmetry Workflows

EDKit is especially strong when a problem has good quantum numbers you can exploit.

## Work In A Fixed Sector

```julia
using EDKit, LinearAlgebra

L = 12
B = basis(L = L, N = L ÷ 2, k = 0, p = 1)

h2 = spin((1.0, "xx"), (1.0, "yy"), (1.0, "zz"))
H = trans_inv_operator(h2, 1:2, B)
vals = eigvals(Hermitian(Array(H)))
```

This combines charge, momentum, and parity into one symmetry-reduced basis.

## Compare A Reduced Basis With Full Space

```julia
using EDKit

L = 8
Bfull = TensorBasis(L = L, base = 2)
Bred = basis(L = L, N = 4, p = 1)
T = DoubleBasis(Bred, Bfull)
P = symmetrizer(T)
```

This is a useful pattern when you want to validate a reduced-basis calculation against a full-space reference.

## Build Constrained Hilbert Spaces

```julia
using EDKit

Bpxp = ProjectedBasis(L = 10, f = x -> all(x[i] + x[i + 1] <= 1 for i in 1:9))
```

This kind of projected basis is a natural fit for kinetically constrained or Rydberg-blockade models.
