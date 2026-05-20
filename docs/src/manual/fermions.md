# Spinless fermions

EDKit provides a `SpinlessFermionBasis` type and a small helper layer that
inserts the Jordan-Wigner string automatically when you build local fermion
operators.

## Basis

```julia
B  = SpinlessFermionBasis(L = 8)                # full Hilbert space
Bn = SpinlessFermionBasis(L = 8, N = 4)         # one fixed-N sector
Bf = SpinlessFermionBasis(L = 8, nf = 1/2)      # density shorthand (= N=4)
Bm = SpinlessFermionBasis(L = 8, N = [3, 4, 5]) # union of several N sectors
```

Each digit `dgt[i]` is `0` for an empty site and `1` for an occupied site.
The basis is internally a `base=2` projected occupation basis; the distinct
type lets `fermion_operator` dispatch on it.

The `N` keyword counts **occupations** (`dgt=1` entries), unlike
[`ProjectedBasis`](@ref) where `N` follows the spin convention (`dgt=0`
entries). The constructor compensates internally so users can write `N`
directly as the particle number. Pass `N = [n1, n2, …]` to construct the
union of several fixed-N sectors (representatives are sorted), or `nf` as a
filling fraction shorthand — the constructor errors if `nf * L` is not an
integer.

## Local operators

The `fermion(op, span)` helper returns the local matrix for an operator on a
contiguous block of `span` sites, with the Jordan-Wigner `σ_z` chain
inserted on the intermediate sites:

| `op`   | meaning                                       | span |
|--------|-----------------------------------------------|------|
| `"n"`  | number operator `n_i = c†_i c_i`              | 1    |
| `"z"`  | `n_i - 1/2` (eigenvalues `±1/2`)              | 1    |
| `"I"`  | identity                                      | 1    |
| `"+"`  | bare creation `c†_i` (no JW)                  | 1    |
| `"-"`  | bare annihilation `c_i` (no JW)               | 1    |
| `"+-"` | hop `c†_1 c_span`                             | ≥ 2  |
| `"-+"` | `c_1 c†_span`                                 | ≥ 2  |
| `"++"` | pair creation `c†_1 c†_span`                  | ≥ 2  |
| `"--"` | pair annihilation `c_1 c_span`                | ≥ 2  |
| `"nn"` | density-density `n_1 n_span` (no JW)          | ≥ 2  |

The default Jordan-Wigner convention is left-string:
`c_j = (σ_z^1 σ_z^2 ⋯ σ_z^{j-1}) σ⁺_j`. EDKit identifies `dgt=0` with empty
and `dgt=1` with occupied, so `c†` maps to `σ⁻` and `c` maps to `σ⁺`.

**Sign convention reminder.** The left Jordan-Wigner string puts `σ_z`
operators to the left of each fermion. For two-fermion products the
on-site `σ⁺ σ_z = -σ⁺` versus `σ⁻ σ_z = +σ⁻` asymmetry causes the
overall sign of `"--"` and `"-+"` to carry a `-1` even at nearest-neighbor
sites, while `"++"` and `"+-"` do not. This is required for the canonical
fermion anticommutation `{c_i, c_j} = 0`.

## Full operators

`fermion_operator(op, sites, B)` embeds the local matrix into the many-body
basis:

```julia
B = SpinlessFermionBasis(L = 6)

# Number on site 3
N3 = fermion_operator("n", [3], B)

# Nearest-neighbor hop c†_1 c_2
hop = fermion_operator("+-", [1, 2], B)

# Long-range hop c†_1 c_5 with σ_z chain on sites 2, 3, 4
long_hop = fermion_operator("+-", [1, 5], B)

# Swapped sites: c†_3 c_1 (= (c†_1 c_3)†, no overall sign)
swapped = fermion_operator("+-", [3, 1], B)

# Bare creation c†_3 — the σ_z prefix on sites 1, 2 is added automatically
cdag3 = fermion_operator("+", [3], B)

# Density-density n_2 n_5 (n's commute; order of sites doesn't matter)
nn25 = fermion_operator("nn", [2, 5], B)
```

A simple tight-binding chain with PBC:

```julia
L = 8
t = 1.0
B = SpinlessFermionBasis(L = L)

H = sum(
    -t * (
        fermion_operator("+-", [i, mod1(i + 1, L)], B) +
        fermion_operator("+-", [mod1(i + 1, L), i], B)
    )
    for i in 1:L
)
```

## Limitations

Supported:

- Single-site `"n"`, `"z"`, `"I"` (diagonal operators, no JW)
- Single-site `"+"` and `"-"` embedded with the full left-string JW prefix
- Two-fermion products `"+-"`, `"-+"`, `"++"`, `"--"` (with JW chain)
- Density-density `"nn"` (no JW; diagonal at both endpoints)
- Fixed-N sectors, multi-sector unions (`N = [n1, n2, …]`), and density
  shorthand `nf = n / L`
- Custom predicate filter `f(dgt) -> Bool` on basis construction

Not yet supported:

- Spinful fermion basis with separate `(N↑, N↓)` sectors
- Symmetry-reduced fermion bases (translation, parity, particle-hole)
- Operator-string parsing for ≥ 4-fermion products (e.g. `"++--"` in one
  call) — build these by composing two-fermion or `"nn"` operators
- Majorana operators `"x"`, `"y"`
- Right-string JW convention
