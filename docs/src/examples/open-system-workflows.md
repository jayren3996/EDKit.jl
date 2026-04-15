# Open-System Workflows

This page collects runnable open-system workflows for EDKit's three distinct
Lindblad layers:

- adaptive many-body Arnoldi propagation with [`lindblad_timeevolve`](@ref),
- the legacy dense many-body Taylor stepper built around [`lindblad`](@ref),
- the quadratic covariance-matrix path built around [`quadraticlindblad`](@ref).

For the conceptual routing page, see [Lindblad Workflows](../manual/lindblad.md).
For the closed-system state-vector counterpart, see
[Time Evolution Workflows](time-evolution-workflows.md).

## Example 1: Many-Body Lindblad With Adaptive Arnoldi

This is the main general-purpose open-system workflow. We build an explicit
many-body Hamiltonian and local jump matrices, propagate one initial density
matrix to many output times, and inspect the diagnostics.

```julia
using EDKit, LinearAlgebra

L = 4
B = TensorBasis(L = L, base = 2)

bond = spin((1.0, "xx"), (1.0, "yy"), (0.4, "zz"))
H = Array(trans_inv_operator(bond, 1:2, B))
jumps = [sqrt(0.15) * Array(operator(spin("+"), [i], L)) for i in 1:L]

lb = lindblad(H, jumps)

ψ0 = productstate([iseven(i) ? 0 : 1 for i in 1:L], B) .+ 0.0im
ρ0 = densitymatrix(ψ0)

ts = collect(range(0.0, 1.5; length = 31))
ρs, diag = lindblad_timeevolve(lb, ρ0, ts;
                               tol = 1e-10,
                               m_init = 8,
                               m_max = 24,
                               return_diagnostics = true,
                               record_min_eig = true)

n1 = Hermitian(Array(operator([1.0 0.0; 0.0 0.0], [1], L)))
occ1 = [real(expectation(n1, ρ)) for ρ in ρs]
traces = [real(tr(ρ.ρ)) for ρ in ρs]

@show size(ρs)
@show diag.basis_builds, diag.restarts, diag.liouvillian_applies
@show maximum(diag.trace_drifts), maximum(diag.hermiticity_drifts)
@show minimum(diag.min_hermitian_eigs)
```

Useful observations:

- `ρs` is a vector of `DensityMatrix` objects, one per requested time.
- `diag.basis_builds`, `diag.restarts`, and `diag.liouvillian_applies` tell
  you whether the reduced Arnoldi model is being reused effectively.
- `diag.trace_drifts` and `diag.hermiticity_drifts` are raw approximation
  diagnostics, not guarantees.
- `record_min_eig = true` adds a simple positivity-related warning signal by
  tracking the minimum eigenvalue of the Hermitian part of each returned matrix.

If your Hamiltonian and jumps are already stored as sparse matrices, build a
[`LiouvillianMap`](@ref) directly and use the same API:

```julia
using SparseArrays

A = LiouvillianMap(sparse(H), sparse.(jumps))
ρs = lindblad_timeevolve(A, ρ0, ts; tol = 1e-10, m_init = 8, m_max = 24)
```

## Example 2: Legacy Dense Taylor-Step Lindblad

The legacy many-body workflow is still convenient when the system is small and
you want an explicit step loop that you control yourself.

```julia
using EDKit, LinearAlgebra

γ = 0.4
jump = sqrt(γ) * ComplexF64[0 1; 0 0]
lb = lindblad(zeros(ComplexF64, 2, 2), [jump])

ρ = densitymatrix(ComplexF64[0.0, 1.0])
projector_excited = ComplexF64[0 0; 0 1]

dt = 0.02
steps = 50
pop_excited = Float64[]

for _ in 1:steps
    push!(pop_excited, real(expectation(projector_excited, ρ)))
    ρ = lb(ρ, dt; order = 8)
end

@show real(tr(ρ.ρ))
@show minimum(eigvals(Hermitian(ρ.ρ)))
```

This route is useful when:

- the density matrix is small enough that dense stepping is cheap,
- you want to integrate inside your own loop with your own stopping logic,
- or you want a transparent baseline for checking a more elaborate setup.

Remember that this is a fixed-step truncated-Taylor update. You choose `dt`
and `order`; EDKit does not adapt them for you.

## Example 3: Quadratic Lindblad

When the model stays quadratic, the covariance-matrix path is usually the more
natural and more scalable choice.

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

Use the three-argument `quadraticlindblad(H, L, M)` form when the dissipator
also contains the additional quadratic terms represented by `M`.

## Example 4: Choosing Between Methods

EDKit's dynamics stack is easiest to navigate if you start from the object you
actually need to evolve.

| Situation | Recommended call | Why |
| --- | --- | --- |
| Closed-system unitary dynamics of a state vector | [`timeevolve`](@ref) | Matrix-free Lanczos on a Hermitian Hamiltonian; the right path for `Operator` workflows and symmetry sectors |
| General many-body Lindblad dynamics at many output times | [`lindblad_timeevolve`](@ref) | Adaptive Arnoldi reuse, diagnostics, and restart/extension control |
| Small explicit many-body Lindblad stepping | `lb(dm, dt; order=...)` after [`lindblad`](@ref) | Simple dense update inside a user-managed loop |
| Quadratic/free-fermion Lindblad dynamics | [`quadraticlindblad`](@ref) | Covariance-matrix scaling rather than full density-matrix scaling |

Rule of thumb:

- if the system is actually closed, move back to [`timeevolve`](@ref),
- if it is open and general, start with [`lindblad_timeevolve`](@ref),
- if it is open but tiny and you want a simple dense loop, use [`lindblad`](@ref),
- if it is quadratic, use [`quadraticlindblad`](@ref).

## See Also

- [Lindblad Workflows](../manual/lindblad.md) for method selection, mental
  model, diagnostics, and caveats.
- [Lindblad Reference](../reference/lindblad.md) for the full API list.
- [Time Evolution Workflows](time-evolution-workflows.md) for the closed-system
  Krylov/Lanczos analogue.
