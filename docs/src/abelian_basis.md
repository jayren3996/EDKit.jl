# AbelianBasis: General Symmetry-Reduced Bases

`AbelianBasis` is EDKit's most general symmetry-reduction mechanism. It constructs
a basis of states that are simultaneous eigenstates of an arbitrary collection of
commuting discrete symmetries, each specified as a permutation of lattice sites
(optionally with spin inversion).

The main user-facing entry point is [`basis`](@ref). For 1D chains, the
keywords `k`, `p`, and `z` provide shorthand for translation, reflection, and
spin-flip sectors. For arbitrary finite lattices, including 2D and 3D systems,
the same constructor accepts explicit symmetry generators through
`symmetries=...`.

## What AbelianBasis Is For

Use `AbelianBasis` when you want to work in a sector defined by several
commuting discrete symmetries without having to build a dedicated basis type for
each geometry.

This is the right tool when:

- you want to combine several symmetry quantum numbers in one construction,
- your lattice is not covered by the dedicated 1D basis types,
- or you want to describe symmetries directly as site permutations.

If you only need the full Hilbert space or a simple constrained subspace,
start with [Bases and Sectors](manual/bases.md) and return here once you need
explicit symmetry reduction.

## 2D Bases In EDKit

There is no separate "2D basis" type in EDKit. A 2D symmetry-reduced basis is
normally an `AbelianBasis` built with `basis(...; symmetries=...)`.

The key idea is that you first choose a linear ordering of the lattice sites,
then express each lattice symmetry as a permutation of `1:L` in that ordering.
Translations in `x` and `y`, reflections, and sublattice inversions all fit
into the same representation. Once that ordering is fixed, use it consistently
for:

- the permutation arrays in `symmetries`,
- the site indices used when you build operators or bond lists,
- and any custom geometry-dependent post-processing.

## Quick Start

### 1D Chain (via convenience keywords)

For 1D periodic chains, the `basis()` constructor provides shorthand keywords:

```julia
# Heisenberg chain, L=12, half-filling, zero momentum, even parity
B = basis(; L=12, N=6, k=0, p=1)

# With spin-flip symmetry
B = basis(; L=12, N=6, k=0, p=1, z=1)
```

### 2D Lattice (via `symmetries` keyword)

For arbitrary lattice geometries, pass symmetry generators as permutation arrays.
This is the recommended route for 2D and 3D systems:

```julia
Lx, Ly = 4, 3
L = Lx * Ly
sites = [(x, y) for y in 0:Ly-1 for x in 0:Lx-1]

# Define symmetry generators as permutation arrays
T_x = [mod(x+1, Lx) + Lx*y + 1 for (x,y) in sites]  # translation in x
T_y = [x + Lx*mod(y+1, Ly) + 1 for (x,y) in sites]   # translation in y
P_x = [Lx-1-x + Lx*y + 1 for (x,y) in sites]          # reflection in x

# Build basis with (kx=0, ky=0, px=+1) symmetry sector
B = basis(; L, N=L÷2, base=2,
    symmetries=[(T_x, 0), (T_y, 0), (P_x, 0)])
```

Each symmetry tuple is `(perm, quantum_number)` where:
- `perm`: 1-indexed permutation array of length L (`perm[i] = j` means site i maps to site j)
- `quantum_number`: integer labeling the irreducible representation (0 to period-1)

The group order (period) is auto-computed from the permutation.

### 3D Lattice

The same approach works for any dimension:

```julia
Lx, Ly, Lz = 2, 2, 3
L = Lx * Ly * Lz
sites = [(x,y,z) for z in 0:Lz-1 for y in 0:Ly-1 for x in 0:Lx-1]

T_x = [mod(x+1,Lx) + Lx*y + Lx*Ly*z + 1 for (x,y,z) in sites]
T_y = [x + Lx*mod(y+1,Ly) + Lx*Ly*z + 1 for (x,y,z) in sites]
T_z = [x + Lx*y + Lx*Ly*mod(z+1,Lz) + 1 for (x,y,z) in sites]

B = basis(; L, symmetries=[(T_x, 0), (T_y, 0), (T_z, 0)])
```

## Spin Inversion

Some symmetries also flip spins. Specify with a third element in the tuple:

```julia
# Spin-flip: identity permutation + full inversion
z_flip = (collect(1:L), 0, trues(L))

# Sublattice inversion: flip spins on even sites only
even_inv = BitVector([iseven(i) for i in 1:L])
z_sub = (collect(1:L), 0, even_inv)

B = basis(; L, symmetries=[z_flip])
```

## Performance

For `base=2` systems, several optimizations accelerate basis construction:

- **Benes networks**: Permutations are compiled into bit-manipulation circuits for fast integer-state operations, avoiding digit-buffer overhead.
- **Integer orbit search**: Canonical representatives are found by operating directly on integer states.
- **Gosper's hack**: When particle number N is fixed, only states with exactly N particles are enumerated, reducing the search space by a factor of `2^L / binomial(L, N)`.
- **Multi-threading**: Basis construction is parallelized across available threads.

## Full Example: 2D Heisenberg Model

```julia
using EDKit, LinearAlgebra

Lx, Ly = 4, 3
L = Lx * Ly
sites = [(x, y) for y in 0:Ly-1 for x in 0:Lx-1]

# Symmetry generators
T_x = [mod(x+1, Lx) + Lx*y + 1 for (x,y) in sites]
T_y = [x + Lx*mod(y+1, Ly) + 1 for (x,y) in sites]

# Nearest-neighbor bonds
J = spin((1.0, "xx"), (1.0, "yy"), (1.0, "zz"))
bonds = Tuple{Int,Int}[]
for (x, y) in sites
    i = x + Lx * y + 1
    push!(bonds, (i, mod(x+1, Lx) + Lx*y + 1))
    push!(bonds, (i, x + Lx*mod(y+1, Ly) + 1))
end

# Ground state in (kx=0, ky=0, N=L/2) sector
B = basis(; L, N=L÷2, base=2, symmetries=[(T_x, 0), (T_y, 0)])
H = operator([J for _ in bonds], [[b[1], b[2]] for b in bonds], B)
E, V = eigen(Hermitian(Array(H)))
println("Ground state energy per site: ", E[1] / L)
println("Reduced Hilbert space dimension: ", size(B, 1))
```
