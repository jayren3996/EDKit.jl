# Spinless fermions

EDKit provides a `SpinlessFermionBasis` type and a small helper layer that
inserts the Jordan-Wigner string automatically when you build local fermion
operators.

## Basis

```julia
B = SpinlessFermionBasis(L = 8)               # full Hilbert space
Bn = SpinlessFermionBasis(L = 8, N = 4)       # fixed-N sector
```

Each digit `dgt[i]` is `0` for an empty site and `1` for an occupied site.
The basis is internally a `base=2` projected occupation basis; the distinct
type lets `fermion_operator` dispatch on it.

The `N` keyword counts **occupations** (`dgt=1` entries), unlike
[`ProjectedBasis`](@ref) where `N` follows the spin convention (`dgt=0`
entries). The constructor compensates internally so users can write `N`
directly as the particle number.

## Local operators

The `fermion(op, span)` helper returns the local matrix for an operator on a
contiguous block of `span` sites, with the Jordan-Wigner `σ_z` chain
inserted on the intermediate sites:

| `op`     | meaning                                  | span |
|----------|------------------------------------------|------|
| `"n"`    | number operator `n_i = c†_i c_i`         | 1    |
| `"+"`    | bare creation `c†_i` (no JW)             | 1    |
| `"-"`    | bare annihilation `c_i` (no JW)          | 1    |
| `"+-"`   | hop `c†_1 c_span`                        | ≥ 2  |
| `"-+"`   | `c_1 c†_span`                            | ≥ 2  |
| `"++"`   | pair creation `c†_1 c†_span`             | ≥ 2  |
| `"--"`   | pair annihilation `c_1 c_span`           | ≥ 2  |

The default Jordan-Wigner convention is left-string:
`c_j = (σ_z^1 σ_z^2 ⋯ σ_z^{j-1}) σ⁺_j`. EDKit identifies `dgt=0` with empty
and `dgt=1` with occupied, so `c†` maps to `σ⁻` and `c` maps to `σ⁺`.

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

## Limitations of the MVP

The current implementation supports:

- Single-site `"n"` (number operator embedded into the basis)
- Single-site `"+"` and `"-"` via `fermion()` only (raw local matrix; no
  JW string embedding through `fermion_operator`)
- Two-fermion products: `"+-"`, `"-+"`, `"++"`, `"--"` (with JW string)

Not yet supported:

- Spinful fermion basis (separate ↑/↓ tracking)
- Four-fermion interaction terms in a single call (build them as products of
  `"n"` operators or as sums of `"+-"` calls)
- Symmetry-reduced fermion bases (translation, parity, particle-hole)
- Embedding bare `c†_i` / `c_i` through `fermion_operator` — this would
  require an explicit Jordan-Wigner string and is not exposed in this MVP.
  Use `fermion("+")` / `fermion("-")` to get the raw 2×2 matrix and embed it
  yourself with [`operator`](@ref) if you need it.
