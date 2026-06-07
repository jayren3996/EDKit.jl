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

### Fixed
- `productstate` now works correctly for symmetry-reduced bases: it uses a local
  digit buffer (no shared-state mutation), keeps the orbit phase/normalization
  coefficient, allocates with the basis element type, and errors on
  configurations outside the basis sector instead of silently writing amplitude
  onto basis vector 1 (`productstate`).

### Performance
- _none yet_

### Added
- Continuous-integration workflow running the test suite on Julia 1.10 and
  latest across Ubuntu and macOS.
