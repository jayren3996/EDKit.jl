# Spinful / Multi-Species Fermions (e1) — Design Spec

**Date:** 2026-06-08
**Branch:** `remediation-0.6.0`
**Phase:** 7 (Extensions), item `e1`
**Status:** approved (mode ordering = blocked; scope = general multi-species with spin-½ ergonomics)

## 1. Goal

Add a spinful / multi-species fermion basis on top of the existing spinless
fermion machinery, so users can build Hubbard-type models with `(site, spin)`
labels and per-species number sectors, without hand-managing Jordan-Wigner mode
indices. Purely additive — no breaking surface.

## 2. Decisions

| Decision | Choice | Consequence |
|---|---|---|
| Mode ordering | **Blocked**: `mode(i,σ) = (σ−1)·L + i` | Same-spin hopping is nearest-neighbour (span 2); Hubbard `U` is JW-free (`nn`). Best for the canonical Hubbard kinetic term. |
| Species scope | **General `S`-species**, spin-½ sugar | Core supports `S` species with per-species sectors; `S=2` gets `:↑`/`:↓` labels and a `hubbard` helper. |
| Representation | Typed `AbstractOnsiteBasis` over `M=S·L` modes | Inherits `index`/`change!`/`size`/`length` from the generic onsite methods (fields `dgt`,`I`,`B`); only `copy` is added. Onsite (not permute) ⇒ JW operators allowed. |
| Enumeration | Reuse `SpinlessFermionBasis(M)` with a per-species Hamming-weight predicate | No new enumeration code. |

## 3. Public API

```julia
struct SpinfulFermionBasis{T<:Integer} <: AbstractOnsiteBasis
    dgt::Vector{T}   # length M = S·L
    I::Vector{T}     # sorted representatives
    B::T             # == 2
    L::Int           # sites
    S::Int           # species
end

SpinfulFermionBasis(dtype=Int64; L, S=2, N=nothing, f=nothing, kwargs...)
```

- `N`:
  - `nothing` — all particle-number sectors.
  - `Integer` — total particle number across all species.
  - `NTuple{S,Int}` / `Vector` of length `S` — per-species numbers `(N₁,…,N_S)`,
    e.g. `(N↑, N↓)`. Each `0 ≤ N_σ ≤ L`.
- `f` — optional extra predicate on the `M`-mode occupation digit string.
- `kwargs` — forwarded to `SpinlessFermionBasis` (`alloc`, `threaded`, `small_N`).

Helpers / operators:
- `fermionmode(B, site, spin) → Int` — `(σ−1)·L + site`; `spin` is an `Integer`
  species in `1:S`, or for `S=2` the symbols `:↑`/`:↓` (`:up`/`:down`).
- `fermion_operator(op, sitespins, B::SpinfulFermionBasis)` — `sitespins` is a
  vector of `(site, spin)` tuples (or a single tuple); maps to modes and
  delegates to the existing `fermion_operator(op, modes, B)`. Raw-mode calls
  still work.
- `hubbard(B; t=1.0, U=0.0, μ=0.0, boundary=:periodic)` — **spin-½ only**:
  `H = −t Σ_{⟨ij⟩,σ}(c†_{iσ}c_{jσ}+h.c.) + U Σ_i n_{i↑}n_{i↓} − μ Σ_{iσ} n_{iσ}`.
  `boundary ∈ (:periodic, :open)`; `L=2` periodic collapses to a single bond.

## 4. Internals

- `_spinful_sector(N, S, L) → (totalN, Nvec_or_nothing)`.
- `_species_weights_ok(dgt, Nvec, L, S)` — each block `[(σ−1)L+1 : σL]` has
  Hamming weight `Nvec[σ]`.
- Constructor builds `SpinlessFermionBasis(dtype; L=M, N=totalN, f=pred)` and
  reuses its `.I`/`.B`.

## 5. Tests (TDD, brute-force / analytic references)

1. **Sector dimension** `= C(L,N↑)·C(L,N↓)` (and total-N, all-sectors counts).
2. **Mode map** `fermionmode(B,i,σ) == (σ−1)L+i`; spin-symbol equivalence.
3. **Canonical anticommutation** `{c_p, c†_q} = δ_pq I` on the full basis.
4. **2-site half-filled Hubbard** ground state `E₀ = (U−√(U²+16t²))/2` (exact),
   `U=0 ⇒ −2t` limit.
5. **Cross-check**: spin-aware construction `==` raw-mode construction on a
   `SpinlessFermionBasis(M)`.
6. **Number conservation** `[N̂, H] = 0`; Hubbard preserves `(N↑,N↓)` sectors.
7. **Errors**: out-of-range `N_σ`, `:↓` with `S<2`, `hubbard` with `S≠2`.

## 6. Files

- `src/Basis/SpinfulFermionBasis.jl` (new) — basis type + constructor + helpers.
- `src/FermionOperator.jl` — `fermion_operator` tuple overload + `hubbard`.
- `src/EDKit.jl` — `include` the new basis file.
- `test/extension_tests.jl` — `e1` testset.
