# Lindblad Workflows

This page is the main routing guide for open-system dynamics in EDKit. It
covers the two many-body density-matrix paths, the quadratic
covariance-matrix path, and how they relate to the closed-system Krylov page
[Time Evolution](time-evolution.md).

## Which Dynamics Path Should I Use?

EDKit now has four distinct dynamics routes. Choosing the right one depends
first on whether the system is closed or open, and then on whether you need a
general many-body density matrix or a quadratic special-case representation.

| Physical setting | Evolved object | Recommended API | Best fit |
| --- | --- | --- | --- |
| Closed-system unitary dynamics with Hermitian `H` | state vector `ψ` | [`timeevolve`](@ref) | Matrix-free `Operator` workflows, symmetry sectors, large Hermitian many-body problems |
| General many-body open-system dynamics | density matrix `ρ` | [`lindblad_timeevolve`](@ref) | Adaptive Arnoldi reuse across many output times, explicit dense or sparse Hamiltonian/jump matrices |
| Small explicit many-body Lindblad stepping | density matrix `ρ` | `lb = lindblad(H, jumps)` followed by `lb(dm, dt; order=...)` | Small dense systems, custom step loops, quick checks against a simple reference integrator |
| Quadratic or free-fermion Lindblad dynamics | covariance matrix `Γ` | [`quadraticlindblad`](@ref) | Scalable special case when the model stays quadratic in Majorana/fermionic variables |

Two practical consequences are worth stating explicitly:

- The open-system many-body routes evolve density matrices, not state vectors.
- The closed-system [`timeevolve`](@ref) path is the only one that directly
  accepts an EDKit `Operator` as a matrix-free Hamiltonian. The open-system
  many-body routes work with explicit Hamiltonian and jump matrices through
  the dense `Lindblad` container or [`LiouvillianMap`](@ref).

If the jump set is empty and the problem is genuinely closed, prefer
[`timeevolve`](@ref) on the state vector rather than propagating a density
matrix.

## Adaptive Krylov / Arnoldi Many-Body Lindblad Evolution

This is EDKit's main general-purpose solver for demanding many-body Lindblad
dynamics. It is the open-system analogue of the closed-system
[`timeevolve`](@ref) layer: instead of a Lanczos basis for a Hermitian
Hamiltonian acting on a state vector, it builds an Arnoldi basis for the
Liouvillian acting on a vectorized density matrix.

### What It Solves

Given a Hamiltonian `H`, jump operators `L_j`, and an initial density matrix
`ρ₀`, the solver propagates

```
∂ₜρ = -i[H, ρ] + Σ_j L_j ρ L_j† - 1/2 Σ_j {L_j† L_j, ρ}.
```

The public inputs are either:

- a legacy dense `Lindblad` object created by [`lindblad`](@ref), or
- a low-level [`LiouvillianMap`](@ref) built from explicit dense or sparse
  Hamiltonian and jump matrices.

The outputs are always `DensityMatrix` objects. The solver is
forward-only: requested times must be non-negative, and a persistent cache can
only be advanced to later times.

Unlike [`timeevolve`](@ref), this layer does not consume an EDKit `Operator`
directly. You first assemble explicit matrices for `H` and the jumps, then
pass those into [`lindblad`](@ref) or [`LiouvillianMap`](@ref).

### Public API

The user-facing entry points are:

```julia
lindblad_timeevolve(A, dm, t::Real; kwargs...)            -> DensityMatrix
lindblad_timeevolve(A, dm, ts::AbstractVector; kwargs...) -> Vector{DensityMatrix}

cache = LindbladArnoldiCache(A, dm; kwargs...)
ρa = lindblad_timeevolve!(cache, 0.5)
ρb = lindblad_timeevolve!(cache, 1.5)
```

with `A` equal to either a `Lindblad` or a [`LiouvillianMap`](@ref), and `dm`
equal to either a `DensityMatrix` or an explicit matrix.

The main supporting types are:

- [`LiouvillianMap`](@ref): low-level linear map representing the Lindblad
  generator on `vec(ρ)`.
- [`LindbladArnoldiCache`](@ref): persistent adaptive Arnoldi workspace,
  anchored at one density matrix and one time.
- [`LindbladArnoldiDiagnostics`](@ref): counters and raw diagnostics collected
  during propagation.

The stateless multi-time form sorts the requested times internally and returns
the densities in the caller's original order. The cache-based
[`lindblad_timeevolve!`](@ref) form is stateful and must be called in
non-decreasing time order.

Initial conditions are usually prepared with [`densitymatrix`](@ref), either
from an explicit matrix, a pure-state vector, or a basis-state index. The
returned outputs are `DensityMatrix` wrappers, so
[`expectation`](@ref) remains the standard way to evaluate observables along
the trajectory.

### Mental Model

At user level, the solver works as follows:

1. Take the current anchor density matrix `ρ_anchor` and vectorize it.
2. Build an Arnoldi basis around that vector, producing a small reduced
   Hessenberg model of the Liouvillian.
3. Cache a reduced propagator by taking a complex Schur decomposition of that
   small matrix.
4. Use the reduced model to serve as many requested output times as possible
   without touching the full many-body generator again.
5. Monitor a cheap Arnoldi defect estimate. When the current basis stops being
   trustworthy over the requested interval, either extend the basis in place
   or restart from a newly reconstructed anchor density matrix.

The important optimization is basis reuse. A successful Arnoldi build is not
thrown away after one time step. It becomes a reusable reduced dynamical model
anchored at the current density matrix and time. Closely spaced output times
can often be served by reduced-space work alone, while long intervals trigger
either basis extension or restart.

This is the same architectural idea as [`timeevolve`](@ref), but with two key
differences:

- the generator is a generally non-Hermitian Liouvillian rather than a
  Hermitian Hamiltonian,
- and the evolved object is a density matrix rather than a state vector.

### Diagnostics

Ask for `return_diagnostics = true`, or inspect `cache.diagnostics` directly,
to recover a [`LindbladArnoldiDiagnostics`](@ref) object.

It tracks two kinds of information.

Operational counters:

- `basis_builds`: fresh Arnoldi builds, including the initial build and every
  restart.
- `basis_extensions`: in-place growth of an existing Arnoldi basis.
- `restarts`: how many times the anchor moved forward and the basis was rebuilt.
- `liouvillian_applies`: total Liouvillian applications.
- `total_times_served`: number of requested output times returned.
- `max_dim_used`: largest Arnoldi dimension reached.
- `accepted_intervals`: validated interval lengths consumed before restarting.

Raw output diagnostics:

- `trace_drifts`: `abs(tr(ρ_raw) - trace_target)` at each returned time.
- `hermiticity_drifts`: `norm(ρ_raw - ρ_raw')` at each returned time.
- `min_hermitian_eigs`: optional `eigmin(Hermitian((ρ_raw + ρ_raw') / 2))`,
  recorded only when `record_min_eig = true`.
- `max_trace_drift`, `max_hermiticity_drift`, and `min_min_hermitian_eig`:
  aggregate extrema over the full run.

These are diagnostics of approximation quality, not guarantees of physicality.
Near-zero trace and Hermiticity drift are reassuring, but they do not prove
that the reduced model is accurate in every other respect. A negative value in
`min_hermitian_eigs` is a useful warning sign, but it is still only a
diagnostic on the Hermitian part of the returned matrix, not a proof of exact
positivity or complete positivity.

One implementation detail matters in practice: the diagnostics are recorded on
the raw reconstructed density matrix before any optional `normalize_trace`
post-processing is applied. That keeps the reported drifts honest even when
the returned density is rescaled for presentation.

### Controlling Accuracy

The main controls are:

- `tol`: defect-monitor tolerance. This is the primary accuracy knob.
- `m_init`: starting Arnoldi dimension.
- `m_max`: largest Arnoldi dimension allowed before a restart becomes
  necessary.
- `extend_basis`: whether the cache may grow the basis in place before
  restarting.
- `extend_step`: number of Arnoldi vectors added per extension.
- `nsample`: minimum number of defect-monitor sample points used when checking
  a candidate interval.
- `bisect_iters`: number of bisection refinements after the first failed
  sample.
- `normalize_trace`: opt-in post-processing that rescales returned density
  matrices back to the initial trace.
- `record_min_eig`: opt-in diagnostic that also tracks the minimum eigenvalue
  of the Hermitian part.

The practical rule is the same as on the closed-system page: touch `tol`
first. Tightening `tol` forces earlier extension or restart and is the right
response when you want more reliable output. Adjust `m_init` and `m_max` when
the problem needs a larger reduced model to amortize well across the time
window.

`nsample` and `bisect_iters` are more advanced controls. The solver already
adapts the effective sample count upward when the reduced model looks harder
to certify, so `nsample` is best read as the minimum sampling density rather
than the exact number of checks.

`normalize_trace = true` should be treated as presentation-oriented
post-processing only. It rescales the returned density matrix; it does not
repair a poor approximation, and it is not a substitute for tightening `tol`.

### When To Use This Method

Prefer [`lindblad_timeevolve`](@ref) when:

- you need general many-body open-system dynamics without assuming quadratic
  structure,
- you want the density matrix at many output times,
- you want adaptive basis reuse rather than manually choosing a small time step,
- you already have explicit dense or sparse matrices for `H` and the jumps,
- or you want raw diagnostics for trace and Hermiticity drift during the run.

It is also the right many-body path when the pure-Hamiltonian limit is only
part of a larger open-system workflow, for example when you want to compare
zero-jump and finite-jump runs with the same density-matrix interface.

### What It Does Not Do

- It is not exactly positivity-preserving at finite tolerance. Generic Arnoldi
  projection does not enforce complete positivity on the returned density
  matrix.
- It does not eliminate the basic cost of many-body density-matrix evolution.
  The reduced model helps, but the underlying object still scales like the full
  Hilbert-space density matrix.
- It is not the same as the closed-system Krylov solver. The generator is
  non-Hermitian and often non-normal, so Arnoldi reuse can be harder to
  maintain than Lanczos reuse for Hermitian Hamiltonians.
- It is not a state-vector propagator. If the system is closed and you only
  need `ψ(t)`, use [`timeevolve`](@ref) instead.
- It is not EDKit's matrix-free `Operator` workflow. The many-body open-system
  layer needs explicit Hamiltonian and jump matrices.
- It is forward-only. A persistent cache rejects requests earlier than its
  last served time.

### Relationship To The Closed-System Krylov Page

The unitary and Lindbladian Krylov layers share the same package-level design,
but they solve different equations:

| Closed system | Open system |
| --- | --- |
| [`timeevolve`](@ref) | [`lindblad_timeevolve`](@ref) |
| state vector `ψ` | density matrix `ρ` |
| Hermitian Hamiltonian | generally non-Hermitian Liouvillian |
| Lanczos basis | Arnoldi basis |
| reduced tridiagonal model | reduced Hessenberg model with cached Schur form |
| norm-focused interpretation | trace / Hermiticity / minimum-Hermitian-eigenvalue diagnostics |

If you already understand the closed-system page, the open-system mental model
is the same broad story with the extra caveats that come from density-matrix
propagation and non-normal generators. For the unitary counterpart, see
[Time Evolution](time-evolution.md).

## Legacy Dense Lindblad Stepper

EDKit's original many-body Lindblad path is still available and still useful:

```julia
using EDKit

H = zeros(ComplexF64, 2, 2)
jumps = [sqrt(0.4) * ComplexF64[0 1; 0 0]]

lb = lindblad(H, jumps)
ρ0 = densitymatrix(ComplexF64[0.0, 1.0])

ρ1 = lb(ρ0, 0.01; order = 8)
dρ = lb * ρ0.ρ
```

This route is intentionally simple:

- [`lindblad`](@ref) materializes the Hamiltonian and jumps as dense matrices,
- `lb * ρ` applies the Lindblad right-hand side directly,
- `lb(dm, dt; order=...)` performs one truncated-Taylor step of size `dt`.

Conceptually, this is very different from the Arnoldi solver. The dense
stepper is a local explicit integrator that knows nothing about future output
times. You control accuracy by choosing `dt` and Taylor `order`, then driving
the step loop yourself. There is no adaptive reduced model, no basis reuse,
and no built-in restart/extension logic.

This is still the right tool when:

- the system is small and fully dense anyway,
- you want the simplest possible explicit stepping workflow,
- you are embedding the update inside your own custom loop,
- or you want a transparent baseline for validating a more sophisticated setup.

Prefer [`lindblad_timeevolve`](@ref) instead when the run asks for many output
times, long trajectories, or enough manual stepping that a reusable reduced
model becomes worthwhile.

## Quadratic Lindblad / Covariance-Matrix Workflows

[`quadraticlindblad`](@ref) is EDKit's scalable special case for quadratic or
free-fermion open systems. Instead of evolving a full many-body density
matrix, it evolves a Majorana covariance matrix:

```julia
using EDKit

A = [0.0 1.0; 1.0 0.0]
B = zeros(2, 2)
Hmaj = majoranaform(A, B)
Ljump = zeros(4, 1)

ql = quadraticlindblad(Hmaj, Ljump)
cm0 = covariancematrix([1, 0])
cm1 = ql(cm0, 0.05; order = 8)
Ablock, Bblock = fermioncorrelation(cm1)
```

The related helpers are:

- [`covariancematrix`](@ref) for constructing covariance matrices,
- [`majoranaform`](@ref) for converting quadratic fermionic data,
- [`fermioncorrelation`](@ref) for recovering correlation blocks from the
  covariance matrix.

This path is not a general replacement for many-body Lindblad evolution. It is
appropriate only when the model really stays within the quadratic/Majorana
formalism. When that structure is available, it is usually much more scalable
than the many-body density-matrix routes. When it is not, return to either the
adaptive Arnoldi many-body solver or the dense many-body stepper.

## Practical Decision Summary

- Use [`timeevolve`](@ref) when the problem is closed, Hermitian, and naturally
  expressed in terms of state vectors.
- Use [`lindblad_timeevolve`](@ref) for general many-body Lindblad dynamics
  when you need adaptive reuse across many output times or longer intervals.
- Use [`lindblad`](@ref) with `lb(dm, dt; order=...)` when the system is small,
  explicit, and you want a simple dense step loop.
- Use [`quadraticlindblad`](@ref) when the physics is quadratic and a
  covariance-matrix description is valid.

## See Also

- [Time Evolution](time-evolution.md) for the closed-system Krylov/Lanczos
  counterpart.
- [Open-System Workflows](../examples/open-system-workflows.md) for worked
  examples of all three open-system routes.
- [Lindblad Reference](../reference/lindblad.md) for the full API surface.
