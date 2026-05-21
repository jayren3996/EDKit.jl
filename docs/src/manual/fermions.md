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
filling fraction shorthand; the constructor errors if `nf * L` is not an
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
| `"+"`  | bare creation `c†_i` (no JW — see warning)    | 1    |
| `"-"`  | bare annihilation `c_i` (no JW — see warning) | 1    |
| `"+-"` | hop `c†_1 c_span`                             | ≥ 2  |
| `"-+"` | `c_1 c†_span`                                 | ≥ 2  |
| `"++"` | pair creation `c†_1 c†_span`                  | ≥ 2  |
| `"--"` | pair annihilation `c_1 c_span`                | ≥ 2  |
| `"nn"` | density-density `n_1 n_span` (no JW)          | ≥ 2  |

!!! warning
    The single-character entries `"+"` and `"-"` return the **bare** on-site
    matrix (`σ⁻` / `σ⁺`) without any Jordan-Wigner prefix. Passing the result
    of `fermion("+")` directly to [`operator`](@ref) at a site `i ≥ 2` will
    produce a wrong operator. Always build `c†_i` / `c_i` via
    [`fermion_operator`](@ref), which inserts the full `σ_z^1 ⋯ σ_z^{i-1}`
    prefix automatically.

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

The same Hamiltonian using the [`trans_inv_fermion_operator`](@ref) helper,
which inserts the long-way Jordan-Wigner chain on the wrap-around bond
automatically:

```julia
H_hop = trans_inv_fermion_operator("+-", [1, 2], B)
H = -t * (H_hop + adjoint(H_hop))
```

!!! warning
    Do **not** use the spin-flavored [`trans_inv_operator`](@ref) for fermion
    `c†c` hops on a ring. It duplicates the same 2-site matrix on every
    translation, including the boundary bond `[L, 1]`, with no Jordan-Wigner
    string on the interior sites. The result happens to match the correct
    fermion Hamiltonian in the single-particle sector but produces wrong
    eigenvalues at `N ≥ 2`. Use `trans_inv_fermion_operator` instead.

## Diagonalization and observables

`Operator` objects materialize to dense or sparse matrices for direct
diagonalization. For the canonical spinless-fermion chain
`H = -t Σ_i (c†_i c_{i+1} + h.c.) + V Σ_i n_i n_{i+1}`:

```julia
using LinearAlgebra
L, N, t, V = 8, 4, 1.0, 2.0
B = SpinlessFermionBasis(L = L, N = N)

H_hop = trans_inv_fermion_operator("+-", [1, 2], B)
H_int = trans_inv_fermion_operator("nn", [1, 2], B)
H = -t * (H_hop + adjoint(H_hop)) + V * H_int

vals, vecs = eigen(Hermitian(Array(H)))
E0, ψ = vals[1], vecs[:, 1]
```

For larger systems, call [`sparse!`](@ref) on the operator and use a Krylov
solver such as `KrylovKit.eigsolve` to find low-lying eigenstates without
forming the dense matrix.

Single-site expectation values and density-density correlators follow the
standard `dot(ψ, opt * ψ)` pattern:

```julia
ni  = real(dot(ψ, fermion_operator("n",  3,      B) * ψ))   # ⟨n_3⟩
nij = real(dot(ψ, fermion_operator("nn", [3, 5], B) * ψ))   # ⟨n_3 n_5⟩
```

(`fermion_operator("n", 3, B)` is shorthand for `fermion_operator("n", [3], B)`.)

## Entanglement entropy

[`SpinlessFermionBasis`](@ref) is an `AbstractOnsiteBasis`, so the standard
[`ent_S`](@ref) helper works without any fermion-specific glue:

```julia
S = ent_S(ψ, 1:L÷2, B)              # von Neumann entropy of the half-chain cut
```

If you want the underlying singular values or the Schmidt matrix itself,
use [`ent_spec`](@ref) or [`EDKit.schmidt`](@ref). Note that `schmidt(...)` returns
the bipartite Schmidt matrix (not a vector of singular values), so to compute
entropy by hand use `svdvals` on its result.

!!! warning "Schmidt convention for fermions"
    [`EDKit.schmidt`](@ref) and [`ent_S`](@ref) compute the decomposition of `ψ`
    viewed as a Jordan-Wigner *spin* state. For a **contiguous** real-space
    bipartition (such as `1:k`), this coincides with the fermionic Schmidt
    decomposition and the entanglement entropies agree. For a
    **non-contiguous** cut (e.g. `Ainds = [1, 3, 5]`), the fermionic
    anticommutation signs are *not* inserted automatically; interpret the
    result with care.

## Symmetry caveats

The Jordan-Wigner transform is *not* a local map: `c†_j` and `c_j` each carry a
string `σ_z^1 ⋯ σ_z^{j-1}` whose length and sign depend on the site index. As
a consequence, spatial spin symmetries (translation, parity/reflection,
spin-flip) do **not** commute with `c†` and `c`. Equivalently, "spin-translation
on the JW representation" and "fermion-translation on the physical operators"
are different group actions:

- A site-permutation `σ : i ↦ π(i)` applied to a spin basis sends
  `σ⁻_i ↦ σ⁻_{π(i)}` (unchanged otherwise).
- The same permutation applied to the fermion operators picks up additional
  fermionic signs from re-ordering the JW string. These signs depend on
  the total particle number `N` and the parity of the permutation.

For this reason EDKit currently does **not** support symmetry-reduced fermion
bases. Calling [`fermion_operator`](@ref) with a [`TranslationalBasis`](@ref),
[`ParityBasis`](@ref), [`FlipBasis`](@ref), [`ParityFlipBasis`](@ref), or
`AbelianBasis` on a JW-bearing operator (`"+"`, `"-"`, `"+-"`, `"-+"`, `"++"`,
`"--"`) raises an error rather than silently returning the wrong matrix.

Diagonal operators (`"n"`, `"z"`, `"I"`, `"nn"`) do **not** carry a JW string
and commute with any onsite permutation, so they remain valid on every basis.

**Recommended workflow:** for fermion problems use
[`SpinlessFermionBasis`](@ref) with `N=…` or `nf=…` (and any extra predicate
filter via `f=…`). Build translation-invariant Hamiltonians with
[`trans_inv_fermion_operator`](@ref) so the wrap-around bond carries the
correct Jordan-Wigner chain.

## ITensor / MPS interface

EDKit's MPS conversion routines (`mps2vec`, `vec2mps`) treat a
[`SpinlessFermionBasis`](@ref) as an ordinary `S = 1/2` chain in the
Jordan-Wigner spin representation, i.e. they expect MPS site indices
built with `siteinds("S=1/2", L)`. ITensor's native
`siteinds("Fermion", L)` site type inserts additional anticommuting
fermionic-parity phases that EDKit does **not** account for, so the two
are incompatible without an explicit conversion. When round-tripping
fermion states between EDKit and ITensor, stick with `"S=1/2"` site
indices.

## Limitations

Supported:

- Single-site `"n"`, `"z"`, `"I"` (diagonal operators, no JW)
- Single-site `"+"` and `"-"` embedded with the full left-string JW prefix
- Two-fermion products `"+-"`, `"-+"`, `"++"`, `"--"` (with JW chain)
- Density-density `"nn"` (no JW; diagonal at both endpoints)
- Translation-invariant assembly via [`trans_inv_fermion_operator`](@ref),
  with the correct long-way JW chain on the PBC wrap-around bond
- Fixed-N sectors, multi-sector unions (`N = [n1, n2, …]`), and density
  shorthand `nf = n / L`
- Custom predicate filter `f(dgt) -> Bool` on basis construction

Not yet supported:

- Spinful fermion basis with separate `(N↑, N↓)` sectors
- Symmetry-reduced fermion bases (translation, parity, particle-hole). Using
  [`fermion_operator`](@ref) with a JW-bearing operator on these bases now
  raises a loud error rather than returning a silently wrong matrix.
- 2D / 3D fermion lattices via [`AbelianBasis`](@ref). The Jordan-Wigner
  string does not commute with arbitrary lattice permutations, so
  [`fermion_operator`](@ref) on an `AbelianBasis` is rejected. A
  fermion-aware lattice basis is future work.
- Operator-string parsing for ≥ 4-fermion products (e.g. `"++--"` in one
  call); build these by composing two-fermion or `"nn"` operators
- Majorana operators `"x"`, `"y"`
- Right-string JW convention
