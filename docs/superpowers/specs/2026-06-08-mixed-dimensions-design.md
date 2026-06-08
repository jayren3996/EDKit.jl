# Mixed Local Dimensions (extension-1) — Design Spec

**Date:** 2026-06-08
**Branch:** `remediation-0.6.0`
**Phase:** 7 (Extensions), item `extension-1`
**Status:** approved (separate `MixedTensorBasis` type)

## 1. Goal

A full tensor-product basis with **per-site local dimensions**, so EDKit can
diagonalize systems mixing different local Hilbert spaces — alternating
spin-½/spin-1 chains, spin-boson models, bosons with different per-site cutoffs.

## 2. Decision

| Decision | Choice | Consequence |
|---|---|---|
| Representation | **Separate `MixedTensorBasis`** (not `TensorBasis{B}`) | The scalar `TensorBasis` and its base-2 fast kernel are bit-for-bit untouched (zero hot-path risk). |
| Mixed-radix encoding | Extend `index`/`change!` with a `base isa AbstractVector` branch | Keyword-arg specialization folds the branch to a compile-time constant per call site, so the scalar path is unaffected. |
| Operators | Reuse the existing matrix-free path unchanged | `colmn!` already calls `index(dgt, sites; base=b.B)`; with `b.B` a vector it takes the mixed-radix branch. No `colmn!` changes. |

## 3. Implementation

- `MixedTensorBasis <: AbstractOnsiteBasis` with fields `dgt::Vector{Int64}`,
  `B::Vector{Int64}` (per-site dims). Constructor `MixedTensorBasis(; dims)`
  validates `dims ≥ 1` and `∏ dims < typemax(Int64)`.
- Mixed-radix helpers in `AbstractBasis.jl` (`_mixed_index`, `_mixed_change!`,
  full and `sites`-subset), selected by the early `base isa AbstractVector`
  return inside `index`/`change!`. Big-endian convention (site 1 most
  significant), matching `TensorBasis`.
- `index(b, dgt)`/`size`/`content`/`norm`/`eltype`/`int_type`/`copy` mirror
  `TensorBasis` but use the per-site vector.
- Entanglement: `schmidtmatrix`'s default subsystem basis is made base-aware
  (`_schmidt_subbasis`) so `schmidt`/`rdm`/`ent_S`/`mutual_information` build a
  `MixedTensorBasis` over the selected sites' dimensions.

## 4. Verification (TDD)

1. `index`/`change!` round-trip over the whole space vs a brute-force big-endian
   mixed-radix reference; digits stay within per-site ranges.
2. A spin-½ ⊗ spin-1 ⊗ spin-½ Hamiltonian built with `operator` equals the
   explicit `kron` embedding (matrix and spectrum).
3. `dims = [2,2,…]` reproduces the uniform `TensorBasis` operator.
4. `rdm`/`schmidt` across a mixed bipartition: site-2 (dim 3) RDM matches a
   brute-force partial trace; Hermitian, unit trace, PSD.
5. Errors: empty `dims`, a zero dimension.

## 5. Deferred (future work)

- Mixed dimensions under symmetry sectors (`ProjectedBasis`/momentum) — fixed-`N`
  and translation are ill-defined or subtle with heterogeneous local spaces.
- A `find_base`-style auto-inference for mixed local-operator dimensions.

## 6. Files

- `src/Basis/MixedTensorBasis.jl` (new).
- `src/Basis/AbstractBasis.jl` — mixed-radix branch in the four `index`/`change!`
  methods.
- `src/Schmidt.jl` — base-aware `_schmidt_subbasis`.
- `src/EDKit.jl` — `include` after `AbstractBasis.jl`.
- `test/extension_tests.jl` — `extension-1` testset.
