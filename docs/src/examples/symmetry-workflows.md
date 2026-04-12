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

## Work In A 2D Momentum Sector

```julia
using EDKit, LinearAlgebra

Lx, Ly = 4, 3
L = Lx * Ly
sites = [(x, y) for y in 0:Ly-1 for x in 0:Lx-1]

T_x = [mod(x + 1, Lx) + Lx * y + 1 for (x, y) in sites]
T_y = [x + Lx * mod(y + 1, Ly) + 1 for (x, y) in sites]

J = spin((1.0, "xx"), (1.0, "yy"), (1.0, "zz"))
bonds = Tuple{Int,Int}[]
for (x, y) in sites
    i = x + Lx * y + 1
    push!(bonds, (i, mod(x + 1, Lx) + Lx * y + 1))
    push!(bonds, (i, x + Lx * mod(y + 1, Ly) + 1))
end

B2d = basis(L = L, N = L ÷ 2, symmetries = [(T_x, 0), (T_y, 0)])
H2d = operator(fill(J, length(bonds)), [[i, j] for (i, j) in bonds], B2d)
vals2d = eigvals(Hermitian(Array(H2d)))
```

For higher-dimensional lattices, `basis(...; symmetries=...)` is the main
symmetry entry point. In that setting you usually build the bond list
explicitly, as above, instead of relying on the 1D ring helper
`trans_inv_operator`.

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

For more on permutation-defined symmetries, including 2D reflections, 3D
translations, and spin-inverting generators, see
[General Abelian Symmetries](../abelian_basis.md).
