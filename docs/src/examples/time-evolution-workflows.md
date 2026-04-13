# Time Evolution Workflows

This page walks through end-to-end closed-system real-time dynamics with
EDKit. The conceptual overview lives in [Time Evolution](../manual/time-evolution.md);
the snippets below are meant to be runnable.

All three examples share the same pattern: build an `Operator` in the basis
you care about, prepare an initial state vector, and call
[`timeevolve`](@ref) on that operator. The propagator is matrix-free, adaptive,
and reuses a single Krylov basis across as many output times as its defect
monitor allows.

## Example 1: Full-Basis XXZ Dynamics

A periodic-boundary XXZ chain with a random normalized initial state, evolved
over a time grid. We also read back the diagnostics to confirm that the
Lanczos basis is actually being reused.

```julia
using EDKit, LinearAlgebra

L = 10
B = TensorBasis(L = L, base = 2)

bond = spin((1.0, "xx"), (1.0, "yy"), (0.7, "zz"))
H = trans_inv_operator(bond, 1:2, B)

ψ0 = randn(ComplexF64, size(B, 1))
normalize!(ψ0)

ts = collect(range(0.0, 2.0; length = 41))
ψs, diag = timeevolve(H, ψ0, ts;
                      tol = 1e-10,
                      m_init = 25, m_max = 50,
                      return_diagnostics = true)

# ψs is a (dim × 41) matrix with one column per requested time.
@show size(ψs)
@show diag.basis_builds, diag.restarts
@show diag.matvecs, diag.max_dim_used
```

A useful sanity check is to track a physical observable, for example the
magnetization of the first site:

```julia
Sz1 = operator(spin("z"), [1], B)
mz  = [real(dot(ψs[:, k], Sz1 * ψs[:, k])) for k in eachindex(ts)]
```

For a closely-spaced time grid like this you typically expect
`basis_builds == 1` and `restarts == 0`: one Lanczos basis serves the whole
window.

## Example 2: Dynamics In A Symmetry Sector

`timeevolve` does not need a new code path for symmetry-reduced bases. Build
the operator in the sector and feed it to the propagator.

```julia
using EDKit, LinearAlgebra

L = 12
B = basis(L = L, N = L ÷ 2, k = 0)     # half-filling, momentum 0

bond = spin((1.0, "+-"), (1.0, "-+"), (0.5, "zz"))
H = trans_inv_operator(bond, 1:2, B)

# A representative state inside the sector.
ψ0 = randn(ComplexF64, size(B, 1))
normalize!(ψ0)

ts = [0.2, 0.5, 1.0, 2.0]
ψs = timeevolve(H, ψ0, ts; tol = 1e-11, m_init = 25, m_max = 50)

# Check that the propagator preserves norm inside the sector.
@show [norm(ψs[:, k]) for k in eachindex(ts)]
```

Because `H` lives in the reduced basis, every Lanczos application runs
entirely inside the momentum-zero half-filling sector. The full `2^L`
Hilbert space is never touched.

A product state inside a sector can be prepared with `productstate`:

```julia
ψN = productstate([i <= L ÷ 2 ? 1 : 0 for i in 1:L], B)
ψt = timeevolve(H, ComplexF64.(ψN), 0.5; tol = 1e-11, m_init = 20, m_max = 40)
```

## Example 3: Driving The Cache Explicitly

If you want to make scheduling decisions between calls, keep the cache alive
across calls. The anchor moves forward adaptively inside
`timeevolve!(cache, t)`; your code just asks for one time at a time.

```julia
using EDKit, LinearAlgebra

L = 10
B = TensorBasis(L = L, base = 2)
H = trans_inv_operator(spin((1.0, "xx"), (1.0, "yy"), (0.7, "zz")), 1:2, B)
ψ0 = productstate([iseven(i) ? 0 : 1 for i in 1:L], B) .+ 0.0im

cache = KrylovEvolutionCache(H, ψ0; tol = 1e-10, m_init = 25, m_max = 50)

ψ1 = timeevolve!(cache, 0.3)
ψ2 = timeevolve!(cache, 0.8)   # typically same Lanczos basis as ψ1
ψ3 = timeevolve!(cache, 3.0)   # may restart the anchor internally

@show cache.diagnostics
```

The diagnostics record how many basis builds and restarts happened over the
whole session. You can use them to set `m_init` / `m_max` to a size where
reuse actually pays off for your specific problem.

## Example 4: Choosing Between Methods

EDKit exposes several ways to move a state under a Hamiltonian. Picking the
right one matters:

| Situation | Recommended call |
|-----------|------------------|
| Closed-system real-time dynamics on a many-body `Operator` | `timeevolve(H, ψ0, t)` |
| Closed-system dynamics on a small dense or sparse matrix | `timeevolve(H, ψ0, t)` — works uniformly |
| One-off toy evolution on a small matrix | `EDKit.expv(A, v; λ = -1im * t)` |
| Full spectrum of a small system | `eigen(Hermitian(Array(H)))` |
| Explicit sparse linear algebra | `sparse(H)` or `sparse!(H)` |
| Open-system density-matrix dynamics | `lindblad(...)` — see [Lindblad Workflows](../manual/lindblad.md) |
| Quadratic free-fermion dynamics | `quadraticlindblad(...)` |
| MPS-based time evolution | `tebd4` / `tebd_n!` — see [ITensor Workflows](../manual/itensors.md) |

Concrete rule of thumb:

- if you find yourself writing a `for t in ts` loop that rebuilds a
  truncated-Taylor approximation on each iteration, replace the loop body
  with a single `timeevolve(H, ψ0, ts)` call,
- if you need the state vector at hundreds of times along a trajectory,
  passing the whole grid to `timeevolve` lets one Lanczos basis serve many
  of them,
- if you only need the state at one time and the system is tiny,
  `exp(-1im * t * Array(H)) * ψ0` is already fine.

The matrix-free path becomes essential when `Array(H)` is too large to hold
in memory, which is the regime EDKit is designed for in the first place.

## See Also

- [Time Evolution](../manual/time-evolution.md) — concepts, controls, and
  limitations of the adaptive Krylov/Lanczos layer.
- [Time Evolution API Reference](../reference/time-evolution.md) — full
  public API.
- [Operators](../manual/operators.md) — constructing the `Operator` objects
  used here.
- [Lindblad Workflows](../manual/lindblad.md) — open-system counterpart for
  density matrices.
