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
the generator `order` times does not return the digits to themselves,
`EDKit.AbelianOperator` throws an error instead of constructing a misleading
generator.

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

Two helpers are useful when debugging a user-defined generator:

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

Here the period is inferred from the action automatically. The wrapper
validates that every `perm` is a permutation of `1:L`, that each `inv` mask has
length `L`, and that the custom generators commute pairwise. To deliberately use a known compatible sector of a
non-commuting generator set, pass `allow_noncommuting_symmetries=true` and check
the resulting sector against a full-space projection.

For `base=2`, the fixed-`N` convention is `sum(dgt) == L - N`; that is, `N`
counts digit-`0` sites. Keep this explicit when comparing against other
libraries whose sector labels count the opposite digit or use magnetization.

## Performance

For `base=2` systems, several optimizations accelerate basis construction:

- **Benes networks**: Permutations are compiled into bit-manipulation circuits for fast integer-state operations, avoiding digit-buffer overhead.
- **Integer orbit search**: Canonical representatives are found by operating directly on integer states.
- **Gosper's hack**: When `N` is fixed and `L <= 62`, only states in the requested fixed-`N` digit sector are enumerated, reducing the search space by a factor of `2^L / binomial(L, N)`.
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

# Unique undirected nearest-neighbor bonds
J = spin((1.0, "xx"), (1.0, "yy"), (1.0, "zz"))
bonds = Set{Tuple{Int,Int}}()
for (x, y) in sites
    i = x + Lx * y + 1
    for j in (
        mod(x+1, Lx) + Lx*y + 1,
        x + Lx*mod(y+1, Ly) + 1,
    )
        push!(bonds, minmax(i, j))
    end
end
bonds = collect(bonds)

# Ground state in (kx=0, ky=0, N=L/2) sector
B = basis(; L, N=L÷2, base=2, symmetries=[(T_x, 0), (T_y, 0)])
H = operator([J for _ in bonds], [[b[1], b[2]] for b in bonds], B)
E, V = eigen(Hermitian(Array(H)))
println("Ground state energy per site: ", E[1] / L)
println("Reduced Hilbert space dimension: ", size(B, 1))
```

## Non-Abelian Point Groups Via Abelian Subgroups

The dihedral group `D_4` of a square lattice has eight elements: four rotations `e, C_4, C_4^2, C_4^3` and four reflections `σ_h, σ_v, σ_d, σ_{d'}`. The dihedral relation `σ_h C_4 σ_h = C_4^{-1}` says rotations and reflections do not commute, so `D_4` is non-abelian and `AbelianBasis` cannot resolve it as a single combined sector. What it can resolve is any abelian subgroup, and a subgroup chain in which a non-abelian generator happens to commute with the projector onto an earlier sector.

The two largest abelian subgroups of `D_4` are `C_4` (rotations only, cyclic of order 4) and `Z_2 × Z_2 = {e, C_4^2, σ_h, σ_v}` (the inversion together with two perpendicular reflections). `C_4` is the more common choice because the rotation eigenvalue `e^{i π k / 2}` is a direct quantum number.

### Why `σ_h` is compatible only with `k = 0` and `k = 2`

The relation `σ_h C_4 = C_4^{-1} σ_h` means `σ_h` takes a `C_4`-eigenstate with momentum `k` to one with momentum `-k mod 4`. The sector closes under `σ_h` exactly when `k ≡ -k mod 4`, which selects `k = 0` and `k = 2`. The remaining sectors `k = 1` and `k = 3` are paired by `σ_h` into a 2D `E` representation of `D_4`; `σ_h` cannot split them.

### Construction on a 3×3 lattice

```julia
using EDKit, LinearAlgebra

Lx = Ly = 3
L = Lx * Ly
sites = [(x, y) for y in 0:Ly-1 for x in 0:Lx-1]
idx(x, y) = mod(x, Lx) + Lx * mod(y, Ly) + 1

# C_4 rotates around the center (1, 1): (x, y) -> (2 - y, x).
C4 = [idx(2 - y, x) for (x, y) in sites]
# σ_h reflects across the middle row: (x, y) -> (x, 2 - y).
sh = [idx(x, 2 - y) for (x, y) in sites]
```

By default, `basis(...; symmetries = [(C4, 0), (sh, 0)])` errors with the message "Custom symmetry generators must commute pairwise". Pass `allow_noncommuting_symmetries = true` to opt into a compatible sector:

```julia
B_A1 = basis(L = L, base = 2,
             symmetries = [(C4, 0), (sh, 0)],
             allow_noncommuting_symmetries = true)
```

For the transverse-field Ising model `H = -J Σ_⟨ij⟩ Z_i Z_j - h Σ_i X_i` on the 3×3 torus, the four compatible `D_4` sectors have:

| representation | `(k, p)` | dimension |
|----------------|----------|-----------|
| `A_1`          | `(0, 0)` | 102       |
| `A_2`          | `(0, 1)` | 38        |
| `B_1`          | `(2, 0)` | 66        |
| `B_2`          | `(2, 1)` | 66        |

These four sectors sum to 272 states. The remaining 240 sit in the 2D `E` representation `{k = 1} ⊕ {k = 3}`. Calling `basis(...; symmetries = [(C4, 1), (sh, p)], allow_noncommuting_symmetries = true)` does not raise an error but returns a basis that is not an irreducible-representation sector and yields a meaningless spectrum. Keep the second-step reflection only in `k = 0` and `k = 2`.

### Verifying against the full space

```julia
ZZ, X = spin((1.0, "zz")), spin((1.0, "x"))
bonds = Set{Tuple{Int,Int}}()
for (x, y) in sites
    i = idx(x, y)
    push!(bonds, minmax(i, idx(x + 1, y)))
    push!(bonds, minmax(i, idx(x, y + 1)))
end
bonds = collect(bonds)
bond_pairs   = [[b[1], b[2]] for b in bonds]
field_sites  = [[i] for i in 1:L]

J, h = 1.0, 0.7
H_full = -J * operator([ZZ for _ in bonds], bond_pairs, L) -
         h  * operator([X  for _ in 1:L],   field_sites, L)
H_A1   = -J * operator([ZZ for _ in bonds], bond_pairs, B_A1) -
         h  * operator([X  for _ in 1:L],   field_sites, B_A1)

E0_full = eigvals(Hermitian(Array(H_full)))[1]
E0_A1   = eigvals(Hermitian(Array(H_A1  )))[1]
@assert abs(E0_full - E0_A1) < 1e-10   # ground state sits in A_1
```

Excited states with other `D_4` quantum numbers are reached by changing `(k, p)` to `(0, 1)`, `(2, 0)`, or `(2, 1)`.
