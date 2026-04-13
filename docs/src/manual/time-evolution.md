# Time Evolution

This page describes EDKit's closed-system real-time evolution layer for state
vectors. It is the recommended route for unitary dynamics of demanding
many-body systems, especially when you want to stay matrix-free and keep using
the same `Operator` object you already built for diagonalization.

## What It Solves

Given a Hermitian Hamiltonian `H` and an initial state `ψ₀`, the
time-evolution layer computes

```
ψ(t) = exp(-i t H) ψ₀
```

without ever forming `exp(-i t H)` explicitly, and without requiring `H` to be
a dense or explicit sparse matrix.

The typical use cases are:

- quench dynamics of a product state under an EDKit `Operator`,
- evolution of a wavepacket or random state in a symmetry sector,
- sweeping a time grid to compute expectation values or overlaps at many
  output times,
- any situation where you already have a matrix-free `Operator` and want to
  reuse it for dynamics.

## Why It Belongs In EDKit

The rest of EDKit is built around a simple workflow:

1. describe local terms,
2. choose a basis,
3. build an `Operator`,
4. reuse that same operator for whichever numerical task you need.

Until now, that last step has covered dense/sparse conversion, matrix-free
application, diagonalization, entanglement diagnostics, and Lindblad-style
evolution. Closed-system real-time dynamics is the natural missing piece.
[`timeevolve`](@ref) consumes exactly the `Operator` objects you already build
with [`operator`](../manual/operators.md) or
[`trans_inv_operator`](../manual/operators.md), in whichever basis you chose,
and it never requires you to materialize `Array(H)` or `sparse(H)` first.

## Public API

The main user-facing entry points are:

```julia
timeevolve(H, ψ0, t::Real; kwargs...)             -> Vector
timeevolve(H, ψ0, ts::AbstractVector; kwargs...)  -> Matrix
timeevolve!(out, H, ψ0, t; kwargs...)             -> out
```

and, when you want to drive the Lanczos cache explicitly rather than rebuild
it on each call:

```julia
cache = KrylovEvolutionCache(H, ψ0; kwargs...)
ψa = timeevolve!(cache, 0.5)
ψb = timeevolve!(cache, 1.5)    # reuses or advances the same basis
timeevolve!(out_matrix, cache, ts)
```

`H` can be any object that supports `LinearAlgebra.mul!(y, H, x)`, including an
EDKit `Operator`, an `AbstractMatrix`, a `SparseMatrixCSC`, or a
`Hermitian` wrapper. The matrix-free `Operator` path is the reason this layer
exists; the plain matrix paths are convenient for tests and small-system
validation.

The single-time form returns `exp(-i t H) ψ₀`. The multi-time form takes a
**sorted** vector of times and returns a matrix whose `k`-th column is the
state at `ts[k]`. A single Krylov basis serves as many requested times as
possible before the adaptive controller decides to restart or extend.

Pass `return_diagnostics = true` to the multi-time form (or inspect
`cache.diagnostics` directly) to see how many basis builds, restarts, and
Hamiltonian applications were used. The diagnostics structure is
[`KrylovEvolutionDiagnostics`](@ref).

## The Mental Model

The central idea is simple. A Lanczos process anchored at a normalized state
`v₁ = ψ₀ / ‖ψ₀‖` produces an orthonormal basis `V_m = [v₁, …, v_m]` and a
small symmetric tridiagonal matrix `T_m` satisfying

```
H V_m = V_m T_m + β_m v_{m+1} eₘᵀ.
```

Inside the Krylov subspace the propagator is cheap: once `T_m` is
diagonalized, every reduced evaluation

```
c(τ) = exp(-i τ T_m) e₁
ψ(τ) ≈ ‖ψ₀‖ V_m c(τ)
```

is just a tiny dense operation, no matter how large the full Hilbert space
is. The expensive part is building the Lanczos basis in the first place,
because that needs `m` matrix-free Hamiltonian applications.

The optimization that matters is therefore basis *reuse*. Instead of rebuilding
a fresh Krylov basis for every small propagation step, `timeevolve` treats
the current basis as a reusable reduced dynamical model anchored at a state
and a time. As long as the reduced model stays accurate, every requested
output time is served by cheap small-matrix work. Only when the basis is no
longer good enough does the propagator either extend it in place (a few more
Lanczos steps) or restart it from the current physical state.

Validity is decided by a cheap defect monitor tied to the Lanczos boundary
term,

```
η_m(τ) = | β_m · eₘᵀ exp(-i τ T_m) e₁ |,
```

which can be evaluated in `O(m)` work per time once the reduced propagator is
cached. The monitor is sampled at several interior points of the candidate
interval and a conservative safe prefix is accepted.

The practical consequence is that one Lanczos build can often serve an entire
output time window with a handful of matrix-vector multiplies, which is why
this API is noticeably faster than a naive call-`expv`-per-step loop for
demanding problems.

## Controlling Accuracy

The most important controls are:

- `tol`            — the defect-monitor threshold that decides whether the
  current basis remains valid over a candidate interval. Smaller `tol` means
  earlier restarts and more accurate propagation. Default `1e-10`.
- `m_init`         — the starting Lanczos dimension. Default `30`. A larger
  value builds a fatter reduced model up front, which amortizes better across
  long output windows but costs more matrix-free applications at build time.
- `m_max`          — the largest Lanczos dimension the cache will ever reach
  before restarting. Default `60`. For very stiff problems a larger `m_max`
  can help reuse survive longer.
- `extend_basis`   — if `true` (default), the cache tries to add more Lanczos
  vectors in place before giving up on the current anchor.
- `extend_step`    — how many vectors to add per in-place extension. Default
  `10`.
- `normalize_output` — opt-in post-processing that renormalizes every
  reconstructed state. Default `false`. Setting this to `true` can **mask
  accumulated propagation error** because it hides norm drift: accuracy is
  delivered by `tol`, not by renormalization. Use this flag only when you
  explicitly want a normalized output (e.g. for plotting) and never as a
  substitute for tightening `tol`.

Operationally, `tol` is the main knob you should touch. It roughly controls
the worst-case local truncation error of the reduced model on each validated
sub-interval. If your problem fails to converge at the default tolerance,
tighten `tol` before you start pushing `m_init` and `m_max`.

## Diagnostics

Pass `return_diagnostics = true` to recover a
[`KrylovEvolutionDiagnostics`](@ref) struct alongside the states. It tracks:

- `basis_builds`       — how many Lanczos bases were built (initial build
  plus any restarts),
- `basis_extensions`   — how many times a basis was extended in place,
- `restarts`           — how many times the anchor moved forward,
- `matvecs`            — total Hamiltonian applications,
- `total_times_served` — how many requested output times were returned,
- `max_dim_used`       — the largest Lanczos dimension ever reached,
- `accepted_intervals` — the interval lengths consumed per restart.

These counters are particularly useful for confirming that adaptive reuse is
actually happening. For a short closely-spaced time window you typically
expect `basis_builds == 1`, `restarts == 0`, and `matvecs ≈ max_dim_used`.

## When To Use This Method

Use [`timeevolve`](@ref) when:

- you want real-time dynamics of a large Hermitian `Operator`,
- you want to stay matrix-free and avoid materializing `Array(H)`,
- you need the state at many output times and want one Krylov build to serve
  several of them,
- or you want accurate long-time propagation without manually managing the
  time step of a low-order integrator.

The older truncated-Taylor helpers `EDKit.expm` and `EDKit.expv` are still
documented in [Utilities](utilities.md), but they are only appropriate for
quick experiments on small matrices. They are not a replacement for
`timeevolve` on demanding many-body problems.

## What It Does Not Do

- **Non-Hermitian generators** are not supported. The current implementation
  assumes `H` is Hermitian and uses a Lanczos (rather than Arnoldi)
  tridiagonalization. Hermiticity is **not checked at runtime**: passing a
  non-Hermitian operator will silently produce meaningless output. There is
  no `hermitian` kwarg — the propagator is Hermitian-only by construction.
- **Backward time evolution is not supported.** All requested times must be
  non-negative, and a persistent `KrylovEvolutionCache` rejects any request
  earlier than its last served time. The stateless `timeevolve(H, ψ0, ts)`
  form accepts an unsorted `ts`, sorts it internally, and returns the
  columns in the caller's original order; the stateful
  `timeevolve!(out, cache, ts)` form requires `ts` to be sorted ascending so
  the anchor cannot overshoot. Start a fresh cache if you need to evolve a
  new state.
- **Density-matrix / open-system evolution** is a different problem entirely;
  use [Lindblad Workflows](lindblad.md) for those.
- **MPS/TEBD-based propagation** is available through the ITensor bridge; see
  [ITensor Workflows](itensors.md). `timeevolve` operates on state vectors in
  an EDKit basis, not on MPS.
- `timeevolve` is a state-vector propagator. It does not compute full spectra
  or Floquet quasi-energies — use explicit diagonalization for those.

## A Note On The Defect Monitor

The defect monitor is a trigonometric sum over the eigenvalues of the
projected tridiagonal matrix, so its oscillation frequency is bounded by the
spectral spread of that small matrix. The adaptive interval selector uses
this bound: it picks the sample count so that the spacing resolves oscillations
up to the Nyquist rate, with a safety factor of 4. A free early-exit also
accepts the full candidate interval whenever the Lanczos boundary coefficient
|β_m| is already below `tol`, since in that case the monitor cannot exceed
`tol` anywhere. This densification keeps the sampled check robust on long
intervals and wide-spectrum operators, where a naive fixed sample grid could
miss narrow peaks between samples.

If you suspect a specific problem is slipping through the monitor, the safest
response is to tighten `tol`. You can also bump `nsample` as a fallback — the
adaptive densifier only raises the count above the user-supplied floor, so a
larger `nsample` is always respected.

## Example 1: Full Basis Quench

A half-filled Néel state quenched under a short XXZ chain:

```julia
using EDKit, LinearAlgebra

L = 10
B = TensorBasis(L = L, base = 2)
H = trans_inv_operator(spin((1.0, "xx"), (1.0, "yy"), (0.7, "zz")), 1:2, B)

ψ0 = productstate([iseven(i) ? 0 : 1 for i in 1:L], B) .+ 0.0im

# Single-time propagation.
ψt = timeevolve(H, ψ0, 1.5; tol = 1e-10)

# Multi-time propagation with diagnostics.
ts = collect(range(0.0, 3.0; length = 61))
ψs, diag = timeevolve(H, ψ0, ts;
                      tol = 1e-10, m_init = 25, m_max = 50,
                      return_diagnostics = true)
@show diag.basis_builds, diag.restarts, diag.matvecs
```

The `ψs` array holds one column per requested time. The diagnostics object
tells you how many Lanczos bases were needed to cover the full window.

## Example 2: Symmetry Sector Dynamics

`timeevolve` works directly on a symmetry-reduced basis. There is nothing
special to set up — the cache asks the operator how to apply itself, and the
operator already knows how to act inside the sector.

```julia
using EDKit, LinearAlgebra

L = 12
B = basis(L = L, N = L ÷ 2, k = 0)      # half-filling, k = 0 sector
bond = spin((1.0, "+-"), (1.0, "-+"), (0.5, "zz"))
H = trans_inv_operator(bond, 1:2, B)

ψ0 = randn(ComplexF64, size(B, 1))
normalize!(ψ0)

ts = [0.2, 0.5, 1.0, 2.0]
ψs = timeevolve(H, ψ0, ts; tol = 1e-11, m_init = 25, m_max = 50)
```

Because the operator lives in the reduced basis, every Lanczos application
runs in the symmetry sector without ever touching the full Hilbert space.
This is the main reason `timeevolve` belongs in EDKit rather than in a
general-purpose linear algebra package.

## Example 3: Reusing A Cache Across Calls

For interactive workflows where you decide the next time to query based on
earlier results, the explicit cache form is convenient:

```julia
using EDKit, LinearAlgebra

L = 10
B = TensorBasis(L = L, base = 2)
H = trans_inv_operator(spin((1.0, "xx"), (1.0, "yy"), (0.7, "zz")), 1:2, B)
ψ0 = productstate([iseven(i) ? 0 : 1 for i in 1:L], B) .+ 0.0im

cache = KrylovEvolutionCache(H, ψ0; tol = 1e-10, m_init = 25, m_max = 50)

ψa = timeevolve!(cache, 0.3)
ψb = timeevolve!(cache, 0.6)        # often served by the same Lanczos basis
ψc = timeevolve!(cache, 2.0)        # may restart the basis internally

@show cache.diagnostics
```

Each call to `timeevolve!(cache, t)` advances forward only; asking for an
earlier time errors by design.

## See Also

- [Operators](operators.md) — building the `Operator` objects that
  [`timeevolve`](@ref) consumes.
- [Bases and Sectors](bases.md) — choosing the basis in which the dynamics
  happens.
- [Time Evolution Workflows](../examples/time-evolution-workflows.md) — a
  longer worked example page with full-basis, symmetry-sector, and
  method-choice comparisons.
- [Utilities](utilities.md) — the older `EDKit.expm` / `EDKit.expv` helpers,
  and when you would still use them.
- [Lindblad Workflows](lindblad.md) — the open-system counterpart for
  density-matrix evolution.
- [Time Evolution API Reference](../reference/time-evolution.md) — the full
  docstring-driven API.
