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

### Fixed
- _none yet_

### Performance
- _none yet_

### Added
- Continuous-integration workflow running the test suite on Julia 1.10 and
  latest across Ubuntu and macOS.
