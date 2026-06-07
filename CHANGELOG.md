# Changelog

All notable changes to EDKit.jl are recorded here, following
[Keep a Changelog](https://keepachangelog.com/) and (loosely) SemVer.

## [Unreleased] — targeting 0.6.0

This release is the outcome of a full section-by-section review. It includes
breaking changes; see each entry below.

### Breaking
- `basis(...)` with a fixed charge `N` and a spin-inversion symmetry (`z` or a
  custom inversion generator) off half-filling now throws instead of silently
  returning an incomplete basis with a wrong spectrum (`abelian-1`).
- `productstate` now errors on configurations outside the basis sector and
  returns a complex vector for symmetry-reduced bases (previously silent / always
  `Float64`).
- Basis constructors (`TensorBasis`, `ProjectedBasis`, `AbelianBasis`, and the
  symmetry bases) now error when `base^L` does not fit the chosen index `dtype`,
  instead of silently overflowing to wrong/negative/empty representatives
  (`bug-1`, `bug-2`, `bug-3`).
- `order(b)` for the translational bases now returns the orbit size `L/a` (and
  `2L/a` for the parity/flip variants) instead of `L`/`2L`. This corrects
  `symmetrizer`/`basis_embedding`/`mps2vec` normalization for unit cells `a>1`
  (results were off by a factor of `a`); `a=1` is unaffected (`trans-1`).

### Fixed
- `productstate` now works correctly for symmetry-reduced bases: it uses a local
  digit buffer (no shared-state mutation), keeps the orbit phase/normalization
  coefficient, allocates with the basis element type, and errors on
  configurations outside the basis sector instead of silently writing amplitude
  onto basis vector 1 (`productstate`).
- `AbelianBasis` / `basis(...)` with a non-default index `dtype` (e.g. `Int32`)
  no longer throws a `MethodError`; representatives and base are converted to the
  requested type at construction (`abelian-2`). Note: the base-2 integer fast
  path remains bounded by 64-bit orbit machinery, so `L ≤ 62`/`64` still applies
  regardless of `dtype` (`abelian-7`; `_gosper_enumerate` already errors above
  this, and the new capacity guard reports it with a clearer message).
- `mul!(target, opt::Operator, m::AbstractMatrix, …)` now uses the `sparse!`
  cache (SpMM) like `*` and `mul` already did, so block and iterative-solver code
  that calls the standard `mul!` gets the documented acceleration instead of
  silently falling back to the matrix-free path (`operator-2`, `operator-5`).
- Documented `sparse!`'s snapshot semantics: the cache is keyed by object
  identity, so mutating an operator's stored matrices in place after `sparse!`
  requires a fresh `sparse!` call to refresh the cache (`operator-1`).
- The public 2-arg `index(B::AbelianBasis, dgt)` no longer mutates the basis's
  shared group odometer `B.G`, making it (and `index_nocheck`) thread-safe like
  the internal hot paths (`abelian-3`).

### Performance
- `ParityBasis`, `FlipBasis`, `ParityFlipBasis` now have a real (`Float64`)
  element type, making their `index` type-stable and letting real Hamiltonians
  assemble as real matrices instead of `ComplexF64` (`parityflip-1`).
- `AbelianBasis` now reports a real (`Float64`) element type for real-character
  sectors (k=0, parity, spin-flip, and k=L/2), keyed on its phase-table type
  parameter, so the most common symmetry workflows assemble and diagonalize real
  matrices at half the memory/FLOPs; genuine momentum sectors remain
  `ComplexF64` (`abelian-5`).

### Added
- Continuous-integration workflow running the test suite on Julia 1.10 and
  latest across Ubuntu and macOS.
