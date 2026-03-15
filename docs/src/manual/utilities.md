# Utilities

A few small helpers in EDKit are worth knowing even though they are not the package's headline features.

## Spectral Statistics

`gapratio` and `meangapratio` compute the adjacent-gap ratio statistics often used in level-statistics studies.

```julia
using EDKit

E = sort(randn(20))
r = gapratio(E)
rbar = meangapratio(E)
```

The input spectrum should already be sorted.

## Product States In A Basis

`productstate(v, B)` returns the basis vector associated with a product configuration `v` inside the basis `B`.

```julia
using EDKit

B = basis(L = 6, N = 3)
psi = productstate([0, 1, 0, 1, 0, 1], B)
```

This is useful when you want a simple initial state inside a constrained or symmetry-aware basis without manually working out the representative index.

There is also an ITensor-side `productstate(sites, states)` method for building product MPS objects from local state vectors.

## Matrix Exponentials For Quick Experiments

`expm` and `expv` provide lightweight truncated-Taylor approximations for matrix exponentials and exponential actions on vectors or matrices.

They are convenient for quick experiments and small problems, but they are not intended to replace more robust specialized algorithms for demanding numerical work.
