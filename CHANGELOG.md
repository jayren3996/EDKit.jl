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
- `ParityBasis`/`FlipBasis`/`ParityFlipBasis` now re-check a custom predicate `f`
  on each symmetry-orbit partner (reflection/flip), rejecting orbits that are not
  `f`-closed. A non-symmetry-invariant `f` previously produced an over-counted,
  invalid basis; it now yields the correct (smaller, possibly empty) basis
  (`parityflip-2`).

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
- Added `copy` methods for `ParityBasis`, `FlipBasis`, and `ParityFlipBasis`
  (previously a `MethodError`), matching the other bases (`parityflip-4`).
- Documented that `covmat` computes the anticommutator covariance only for
  **Hermitian** operators (`qim-1`).
- Documented `tebd4`'s `exp(τ·h)` convention: real-time evolution by `dt` needs
  `τ = -im*dt`; a real `dt` performs imaginary-time evolution (`tebd4`).
- `change!(dgt, ind; base)` accepts `ind`/`base` integer types that differ from
  the digit buffer's element type (e.g. an `Int64` index into an `Int32` buffer),
  converting internally instead of throwing a `MethodError` (`bug-4`).
- `entropy(s; α, cutoff)` now forwards `cutoff` to the Rényi branch (α∉{0,1}) and
  `renyi_entropy` normalizes its input over above-cutoff entries, so noise
  Schmidt values are dropped and unnormalized inputs give correct entropies
  (`schmidt-1`, `schmidt-2`).
- `expm`/`expv` now use scaling-and-squaring, so they stay accurate for matrices
  with large norm (the previous fixed-order Taylor truncation diverged for
  `‖A‖ ≳ 5`) (`expm-1`).
- `gapratio`/`meangapratio` gained a `sorted` keyword (sorts the spectrum when
  `false`) and now return `NaN` for fully-degenerate `0/0` gaps instead of a
  spurious `1.0`; `meangapratio` filters those out (`gapratio-1`).
- The adaptive Krylov time-evolution defect monitor now **shrinks the candidate
  interval** instead of under-sampling it when full Nyquist resolution would
  exceed the sample cap. Previously a very large `ω·τ` could leave a long
  interval sampled below the Nyquist rate, letting the integrated defect
  `∫η` (which bounds the error) exceed `tol` between samples; the monitor is now
  always fully resolved on every accepted interval (`te-1`).

### Performance
- `ParityBasis`, `FlipBasis`, `ParityFlipBasis` now have a real (`Float64`)
  element type, making their `index` type-stable and letting real Hamiltonians
  assemble as real matrices instead of `ComplexF64` (`parityflip-1`).
- `AbelianBasis` now reports a real (`Float64`) element type for real-character
  sectors (k=0, parity, spin-flip, and k=L/2), keyed on its phase-table type
  parameter, so the most common symmetry workflows assemble and diagonalize real
  matrices at half the memory/FLOPs; genuine momentum sectors remain
  `ComplexF64` (`abelian-5`).
- `operator(...)` deduplicates repeated local-term supports in `O(num)` via a
  `Dict` instead of an `O(num²)` linear scan, speeding up construction of
  operators with many terms (`operator-6`).
- Fixed-charge construction for `base>2` (`ProjectedBasis`/`TranslationalBasis`
  with `N`) now enumerates only in-range digit strings via a bounded mixed-radix
  generator, instead of generating all compositions and discarding 70–95% of
  them — markedly faster for large qudit bases (`optimization-1`).
- `basis_embedding` (and therefore `symmetrizer`) now builds a **sparse** matrix
  with the basis's natural element type instead of a dense `base^L × dim`
  `ComplexF64` matrix, so these maps scale to larger systems and real sectors
  stay real (`linearmap-1`).

### Added
- `MixedTensorBasis(dims=[…])` — full tensor-product basis with **per-site local
  dimensions** (e.g. alternating spin-½/spin-1 chains, spin-boson models, bosons
  with per-site cutoffs). Operators (`operator(mat, sites, B)`) and entanglement
  (`schmidt`/`rdm`/`ent_S`/`mutual_information`) work via mixed-radix
  `index`/`change!`. The scalar `TensorBasis` and its base-2 fast kernel are
  untouched; the mixed-radix branch is folded away at compile time on the scalar
  path (`extension-1`).
- `TranslationalFermionBasis(L=…, N=…, k=…)` — momentum-resolved basis for
  spinless fermions at fixed `N`. The many-body translation sign
  `(−1)^((N−1)·n_wrap)` (periodic for odd `N`, antiperiodic for even `N`) is
  baked into the orbit phase, so `trans_inv_fermion_operator` builds correct
  Hamiltonians per momentum sector. Verified: the union of all `k`-sector spectra
  equals the full `N`-sector spectrum. Real (`Float64`) for `k=0` and `k=L/2`
  sectors. The Jordan-Wigner guard is relaxed for this one basis type (`e2`).
- `SpinfulFermionBasis(L=…, S=2, N=…)` — spinful / multi-species fermion
  occupation basis over `M = S·L` modes (blocked ordering `mode(i,σ)=(σ−1)L+i`),
  with per-species `(N₁,…,N_S)` or total-`N` sectors. Spin-aware operators via
  `fermion_operator(op, [(site, spin)…], B)` (spin `:↑`/`:↓` or integer species)
  and a `fermionmode(B, site, spin)` helper. A `hubbard(B; t, U, μ, boundary)`
  helper builds the spin-½ Fermi–Hubbard Hamiltonian (`e1`).
- `timeevolve`/`timeevolve!` now support **backward and two-sided** evolution: a
  negative time evolves by `exp(+i|t|H)`, and a forward pass followed by a
  backward pass enables OTOC / Heisenberg-picture workflows. Each cache evolves
  in one direction (fixed by its first motion); mixing positive and negative
  times in one stateless call is rejected (`te-4`).
- `steadystate(A)` computes a Lindblad steady state `ρ_ss` (`𝓛[ρ_ss]=0`) from a
  `LiouvillianMap`, `Lindblad`, or `(H, jumps)`. Uses dense diagonalization for
  small systems (`d ≤ 16`) and matrix-free Arnoldi (`KrylovKit.eigsolve`) for
  large ones, returning a trace-normalized Hermitian `DensityMatrix`
  (`lindblad-6`). This adds **KrylovKit** as a direct dependency (previously
  transitive via ITensorMPS).
- `rdm(v, Ainds, b)` returns the reduced density matrix `ρ_A = S S†` of a
  subsystem from a state vector, reusing every existing `schmidt` dispatch
  (works for tensor, projected, translational, parity/flip, and Abelian bases);
  real-phase bases give a real `ρ_A` (`schmidt-4`).
- `mutual_information(v, Ainds, Cinds, b)` computes `I(A:C) = S(A)+S(C)−S(A∪C)`
  between two disjoint subsystems, with `α`/`cutoff` forwarded to `ent_S`
  (`schmidt-4`).
- Continuous-integration workflow running the test suite on Julia 1.10 and
  latest across Ubuntu and macOS.
