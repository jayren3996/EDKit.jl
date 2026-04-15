# Time Evolution Reference

EDKit ships an adaptive Krylov/Lanczos real-time propagator for Hermitian
Hamiltonians. The conceptual overview, accuracy controls, diagnostics, and
worked walkthroughs live in the manual and example pages:

- [Time Evolution](../manual/time-evolution.md) — the manual chapter for
  users.
- [Time Evolution Workflows](../examples/time-evolution-workflows.md) —
  end-to-end examples.
- [Lindblad Workflows](../manual/lindblad.md) — the open-system analogue for
  density-matrix evolution.
- [Lindblad Reference](lindblad.md) — the many-body and quadratic
  open-system APIs.
- [Operators](../manual/operators.md) — how the matrix-free `Operator`
  objects consumed here are built.
- [Utilities](../manual/utilities.md) — the older lightweight `EDKit.expv` /
  `EDKit.expm` helpers and when you would still use them instead.

The central idea is to treat the Krylov subspace as a *reusable reduced
dynamical model*: a Lanczos basis is built once from an anchor state, the
small projected matrix is diagonalized to cache a cheap reduced propagator,
and a defect-monitor tied to the Lanczos boundary term decides when the
current basis is still good enough to serve more output times. Only when the
monitor exceeds tolerance is the basis extended in place or restarted from a
new anchor.

The public entry points live in `src/algorithms/TimeEvolution.jl` and work
with any object that supports `LinearAlgebra.mul!(y, H, x)`, including EDKit
`Operator` instances, `AbstractMatrix`, and `SparseMatrixCSC`.

## Example

```julia
using EDKit, LinearAlgebra

L = 10
B = TensorBasis(L = L, base = 2)
H = trans_inv_operator(spin((1.0, "xx"), (1.0, "yy"), (0.7, "zz")), 1:2, B)

ψ0 = productstate([iseven(i) ? 0 : 1 for i in 1:L], B) .+ 0.0im

# single-time evolution
ψt = timeevolve(H, ψ0, 1.5; tol = 1e-10)

# multi-time evolution: one Krylov basis serves many closely-spaced times
ts = collect(range(0.0, 2.0; length = 41))
ψs, diag = timeevolve(H, ψ0, ts;
                      tol = 1e-10, m_init = 25, m_max = 50,
                      return_diagnostics = true)
@show diag.basis_builds, diag.restarts, diag.matvecs
```

For workflows that want explicit control over the Lanczos cache, use
[`KrylovEvolutionCache`](@ref) directly and advance it in place with
`timeevolve!`.

```@docs
timeevolve
timeevolve!
KrylovEvolutionCache
KrylovEvolutionDiagnostics
```
