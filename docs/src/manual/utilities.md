# Utilities

A few smaller helpers in EDKit are worth knowing even though they are not the
package's headline features. They cover four main jobs:

- spectral statistics,
- basis-space product-state preparation,
- lightweight exponential actions for small problems,
- and the quantum inverse method utilities.

For serious closed-system real-time dynamics, use the dedicated
[Time Evolution](time-evolution.md) layer rather than the lightweight
exponential helpers on this page.

## Spectral Statistics

`gapratio` and `meangapratio` compute the adjacent-gap ratio statistics often used in level-statistics studies.

```julia
using EDKit

E = sort(randn(20))
r = gapratio(E)
rbar = meangapratio(E)
```

The input spectrum should already be sorted.

## Product States In A Basis

`productstate(v, B)` returns the basis vector associated with a product configuration `v` inside the basis `B`.

```julia
using EDKit

B = basis(L = 6, N = 3)
psi = productstate([0, 1, 0, 1, 0, 1], B)
```

This is useful when you want a simple initial state inside a constrained or
symmetry-aware basis without manually working out the representative index.

The configuration should belong to the basis sector you chose. In other words,
if `B` fixes particle number, parity, or other constraints, the digits you pass
to `productstate` should respect those constraints too.

There is also an ITensor-side `productstate(sites, states)` method for building
product MPS objects from local state vectors.

## Matrix Exponentials For Quick Experiments

`EDKit.expm` and `EDKit.expv` provide lightweight truncated-Taylor
approximations for matrix exponentials and exponential actions on vectors or
matrices.

```julia
using EDKit

A = randn(4, 4)
v = randn(4)

U = EDKit.expm(A; order=12)
w = EDKit.expv(A, v; order=12, λ=0.1)
```

These helpers are deliberately simple: a fixed-order Taylor expansion with no
adaptive error control. They are convenient for quick experiments, for
computing a small dense exponential in a test, or for one-off evolution of a
low-dimensional auxiliary operator.

They are **not** the right tool for real-time dynamics of a many-body
Hamiltonian:

- there is no adaptive error control,
- there is no Krylov subspace that would let you reuse work across output times,
- and the fixed-order Taylor series loses accuracy quickly as `‖λA‖` grows.

For demanding closed-system many-body dynamics, prefer the adaptive
Krylov/Lanczos propagator documented in [Time Evolution](time-evolution.md).
It works directly on matrix-free EDKit `Operator`s, reuses one Lanczos basis
across many output times, and exposes a principled tolerance control.

Use `EDKit.expm` / `EDKit.expv` when:

- the operator is small and already explicit,
- you want a quick one-liner,
- or you are writing a reference against which to check a heavier method.

Use [`timeevolve`](@ref) for everything else.

## Quantum Inverse Method Helpers

`covmat` and `qimsolve` support inverse problems where you want to recover or
constrain a Hamiltonian from target state data.

`covmat(ops, v)` builds the covariance matrix of a list of operators on one
target state or an ensemble of target states. `qimsolve(ops, v; tol=...)` then
extracts the small-variance directions from that covariance matrix.

```julia
using EDKit

ops = [randn(4, 4) for _ in 1:3]
psi = randn(4)

C = covmat(ops, psi)
sol = qimsolve(ops, psi; tol=1e-7)
```

The operator list may contain dense matrices or EDKit `Operator` objects, and
the target data may be one vector or a matrix whose columns represent several
target states.
