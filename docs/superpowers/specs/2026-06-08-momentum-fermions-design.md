# Momentum-Resolved Spinless Fermions (e2) — Design Spec

**Date:** 2026-06-08
**Branch:** `remediation-0.6.0`
**Phase:** 7 (Extensions), item `e2`
**Status:** approved (scope = spinless translation only)

## 1. Goal

A momentum-resolved basis for **spinless fermions** at fixed particle number `N`,
the fermionic analogue of `TranslationalBasis`, so Hamiltonians can be
block-diagonalized by lattice momentum.

## 2. Physics — the translation sign

The many-body translation operator on the Jordan-Wigner occupation
representation differs from the bosonic site shift by a fermion sign. Under one
single-site translation (site `i → i+1`, site `L → 1` in EDKit's big-endian
`dgt`), the occupied boundary site `dgt[L]` wraps past the other `N−1` particles:

    sign = (−1)^((N−1)·dgt[L])

- **odd `N`** ⇒ sign `≡ +1` (periodic boundary for the single-particle problem),
- **even `N`** ⇒ sign `= −1` when the boundary site is occupied (antiperiodic),
  which **shifts the compatible momenta**.

Over a full period of `L` translations every particle wraps once, total sign
`(−1)^(N(N−1)) = +1`, consistent with `Tᴸ = 1`.

*Validated empirically before implementation:* the explicit `T_f` built with this
rule satisfies `[T_f, H] = 0` and is unitary for a translation-invariant fermion
Hamiltonian, for both even and odd `N`.

## 3. Design

`TranslationalFermionBasis{Ti,T} <: AbstractPermuteBasis` mirrors
`TranslationalBasis` (fields `dgt, I, R, C, A, B`) plus `evenN::Bool`.

- **Judge** (`TranslationFermionJudge`): cycles the orbit, accumulates the
  fermion sign `σ`, and at period `p` accepts iff the combined phase
  `φ_p·σ = 1`, where `φ_p = exp(−i2πkp/L)`. Concretely
  `_fermion_check_momentum(p,k,L,σ)`: with `r = mod(2kp, 2L)`, accept when
  `r==0 ∧ σ==+1` or `r==L ∧ σ==−1`.
- **`index(b, dgt)`**: integer-translates `dgt` to its representative,
  accumulating `σ`, and returns `C[M+1]·σ·R[i]` (momentum phase × fermion sign ×
  normalization).
- Restrictions (MVP): `base = 2`, single-site translation (`a = 1`), fixed `N`.

Operators: the Jordan-Wigner guard in `fermion_operator` is relaxed for
`TranslationalFermionBasis`; `trans_inv_fermion_operator` accepts it and builds
translation-invariant fermion operators per momentum sector.

## 4. Verification (decisive test)

For a translation-invariant fermion Hamiltonian (hopping + `nn` interaction):
the **union of all `k`-sector spectra equals the full `SpinlessFermionBasis(L,N)`
spectrum**, and the sector dimensions partition the `N`-sector. Checked for
several `(L, N, V)` with even and odd `N` — matches to machine precision.

## 5. Deferred (future work)

- Fermion **parity/reflection** momentum sectors (reflection reverses the JW
  order → extra global sign).
- **Spinful** momentum (translation permutes modes across spin blocks under the
  blocked ordering; sign depends on total `N`).
- Unit cells `a > 1`.

## 6. Files

- `src/Basis/TranslationalFermionBasis.jl` (new).
- `src/FermionOperator.jl` — relaxed JW guard; `trans_inv_fermion_operator`
  accepts the new basis.
- `src/EDKit.jl` — `include` after `TranslationalBasis.jl`.
- `test/extension_tests.jl` — `e2` testset.
