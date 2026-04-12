# AbelianBasis: General Symmetry-Reduced Bases

`AbelianBasis` is EDKit's low-level basis for commuting discrete symmetry
actions on digit strings. The object that actually generates those actions is
`EDKit.AbelianOperator`, and `EDKit.AbelianBasis` builds a reduced basis from
the orbits of that action.

This page focuses on the direct construction route:

- define one or more cyclic generators as actions on a digit buffer,
- combine them into one commuting Abelian action,
- build an `EDKit.AbelianBasis` from that action,
- then use `index`, `change!`, and operator construction exactly as with other
  basis types.

If you only want the convenience wrapper, [Bases and Sectors](manual/bases.md)
covers `basis(...; symmetries=...)`.

## The Low-Level Objects

The two important objects are:

- `EDKit.AbelianOperator(order, k, perm; inv=falses(length(perm)))`
- `EDKit.AbelianBasis(; L, G, base=2, f=x->true, N=nothing, threaded=true)`

`EDKit.AbelianOperator` is one cyclic generator together with its phase
representation. `EDKit.AbelianBasis` stores the canonical orbit
representatives, normalization factors, and the combined group action.

## Designing A Generator On Digits

Each generator acts on a digit buffer in two steps:

1. permute the sites,
2. optionally complement selected local digits.

The permutation convention is:

- `perm[i] = j` means the digit currently at site `i` moves to site `j`.

For example, on a 4-site buffer:

```julia
using EDKit

perm = [2, 3, 4, 1]
g = EDKit.AbelianOperator(4, 0, perm)

dgt = [0, 1, 2, 3]
dgt2 = copy(dgt)
EDKit.init!(g)
g(dgt2, 4)

# dgt2 is now [3, 0, 1, 2]
```

The `order` argument must match the actual period of the action. If applying
the generator `order` times does not return the digits to themselves, the basis
construction is not describing the symmetry you think it is.

The `k` argument selects the character of that cyclic generator:

- `k = 0` gives the trivial phase sector,
- `k = order ÷ 2` gives a sign sector when `order` is even,
- more generally the phase is `exp(2πimk / order)` after `m` applications.

## Optional Digit Inversion

After the permutation, `EDKit.AbelianOperator` can also complement selected
digits through the `inv` bit mask:

```julia
using EDKit

perm = collect(1:4)
inv = BitVector([true, false, true, false])
g = EDKit.AbelianOperator(2, 0, perm; inv=inv)
```

At a site marked by `inv[i] = true`, the local digit is replaced by
`base - 1 - dgt[i]`. For spin-1/2 systems this is the usual `0 <-> 1` flip.

## Combining Generators

Multiple commuting cyclic generators are combined with `+`:

```julia
using EDKit

Lx, Ly = 2, 2
L = Lx * Ly
sites = [(x, y) for y in 0:Ly-1 for x in 0:Lx-1]

T_x = [mod(x + 1, Lx) + Lx * y + 1 for (x, y) in sites]
T_y = [x + Lx * mod(y + 1, Ly) + 1 for (x, y) in sites]

Gx = EDKit.AbelianOperator(Lx, 0, T_x)
Gy = EDKit.AbelianOperator(Ly, 0, T_y)
G = Gx + Gy
```

This direct-product construction assumes the generators commute. In practice,
that means applying them in different orders should produce the same action on
every digit string.

## Building AbelianBasis Directly

Once the action is defined, build the reduced basis directly:

```julia
using EDKit

Lx, Ly = 2, 2
L = Lx * Ly
sites = [(x, y) for y in 0:Ly-1 for x in 0:Lx-1]

T_x = [mod(x + 1, Lx) + Lx * y + 1 for (x, y) in sites]
T_y = [x + Lx * mod(y + 1, Ly) + 1 for (x, y) in sites]

G = EDKit.AbelianOperator(Lx, 0, T_x) +
    EDKit.AbelianOperator(Ly, 0, T_y)

B = EDKit.AbelianBasis(; L, G, base=2, N=L÷2)
```

This is the low-level route behind the higher-level `basis(...; symmetries=...)`
wrapper. It is the better choice when you want direct control over the group
action itself.

## Checking A Custom Action

Two helpers are especially useful when debugging a user-defined generator:

- `EDKit.check_min(dgt, G; base=...)` checks whether `dgt` is already the
  canonical orbit representative,
- `EDKit.shift_canonical!(dgt, G; base=...)` finds the representative reached by
  the symmetry action.

Example:

```julia
using EDKit

perm = [2, 3, 4, 1]
G = EDKit.AbelianOperator(4, 0, perm)
dgt = [0, 1, 1, 0]

EDKit.check_min(copy(dgt), deepcopy(G); base=2)
EDKit.shift_canonical!(copy(dgt), deepcopy(G); base=2)
```

When the basis is built, `index(B, dgt)` uses the same canonicalization logic.
It returns the coefficient induced by the symmetry phase and normalization,
together with the basis position of the representative.

## 2D And 3D Lattices

There is no separate "2D basis" type in EDKit. A higher-dimensional
symmetry-reduced basis is still an `EDKit.AbelianBasis`; the only extra work is
choosing a site ordering and translating lattice symmetries into permutations of
`1:L`.

Once that ordering is fixed, use it consistently for:

- the permutation arrays that define the generators,
- the operator site indices and bond lists,
- and any geometry-dependent post-processing.

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

The same direct construction can also be written with explicit generators:

```julia
using EDKit

G = EDKit.AbelianOperator(Lx, 0, T_x) +
    EDKit.AbelianOperator(Ly, 0, T_y) +
    EDKit.AbelianOperator(Lz, 0, T_z)

B = EDKit.AbelianBasis(; L, G, base=2)
```

## Convenience Wrapper

If you do not need direct control over `EDKit.AbelianOperator`, the same ideas
can be passed through `basis(...; symmetries=...)`. Each tuple is
`(perm, quantum_number)` or `(perm, quantum_number, inv)`:

```julia
# Spin-flip: identity permutation + full inversion
z_flip = (collect(1:L), 0, trues(L))

# Sublattice inversion: flip spins on even sites only
even_inv = BitVector([iseven(i) for i in 1:L])
z_sub = (collect(1:L), 0, even_inv)

B = basis(; L, symmetries=[z_flip])
```

Here the period is inferred from the action automatically, which is convenient
for quick lattice constructions. The direct `EDKit.AbelianOperator` route is
still the better documentation target when you are designing the symmetry action
itself.

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
