# Getting Started

This page introduces the core EDKit workflow:

1. define local terms,
2. specify where they act,
3. choose a basis,
4. build an operator,
5. apply it or convert it to an explicit matrix.

## The Main Abstraction

The central user-facing object in EDKit is `Operator`.

An `Operator` stores a many-body operator as a sum of local terms on a chosen basis. That means you can describe your Hamiltonian once and then decide how you want to use it:

- apply it directly to a state vector,
- convert it to `Array(H)` for dense exact diagonalization,
- convert it to `sparse(H)` for sparse workflows.

## A First Example

Here is a compact open-chain XXZ Hamiltonian on `L = 8` sites:

```julia
using EDKit, LinearAlgebra, SparseArrays

L = 8
Delta = 1.5

mats = [
    fill(spin("XX"), L - 1);
    fill(spin("YY"), L - 1);
    fill(Delta * spin("ZZ"), L - 1);
]

inds = vcat(
    [[i, i + 1] for i in 1:L-1],
    [[i, i + 1] for i in 1:L-1],
    [[i, i + 1] for i in 1:L-1],
)

H = operator(mats, inds, L)
```

`H` is now an `Operator` acting on the full `TensorBasis` of a spin-1/2 chain.

## Using The Operator

Apply it to a state:

```julia
psi = normalize(randn(ComplexF64, 2^L))
Hpsi = H * psi
```

Convert it to explicit matrix representations:

```julia
Hdense = Array(H)
Hsparse = sparse(H)
```

Diagonalize the dense version:

```julia
vals = eigvals(Hermitian(Hdense))
```

## A First Time-Evolution Example

The same `Operator` you just built can be fed directly into the closed-system
real-time propagator, with no intermediate dense conversion:

```julia
using EDKit, LinearAlgebra

L = 8
B = TensorBasis(L = L, base = 2)
H = trans_inv_operator(spin((1.0, "xx"), (1.0, "yy"), (0.7, "zz")), 1:2, B)

ψ0 = productstate([iseven(i) ? 0 : 1 for i in 1:L], B) .+ 0.0im

ψt = timeevolve(H, ψ0, 1.0)                 # state at t = 1
ts = collect(range(0.0, 2.0; length = 21))
ψs = timeevolve(H, ψ0, ts)                  # one column per requested time
```

[`timeevolve`](@ref) uses an adaptive Krylov/Lanczos propagator that reuses
one Lanczos basis across as many requested times as the error tolerance
allows. It is the recommended route for demanding real-time dynamics; the
older truncated-Taylor helpers `EDKit.expm` and `EDKit.expv` are only meant
for quick experiments on small matrices.

See [Time Evolution](manual/time-evolution.md) for the full manual page and
[Time Evolution Workflows](examples/time-evolution-workflows.md) for longer
worked examples.

For open systems, the routing is slightly different:

- use [Lindblad Workflows](manual/lindblad.md) when you need density-matrix
  evolution,
- use [`lindblad_timeevolve`](@ref) for the adaptive many-body Arnoldi path,
- use `lindblad(...)(ρ, dt; order=...)` for small dense explicit stepping,
- use [`quadraticlindblad`](@ref) when the problem is quadratic and a
  covariance-matrix description applies.

## Working In A Symmetry Sector

EDKit becomes especially useful when you do not want the full Hilbert space.

For example, the half-filling momentum-zero sector can be built with the high-level `basis(...)` helper:

```julia
using EDKit, LinearAlgebra

L = 12
B = basis(L = L, N = L ÷ 2, k = 0)

h2 = spin((1.0, "xx"), (1.0, "yy"), (1.0, "zz"))
H = trans_inv_operator(h2, 1:2, B)

vals, vecs = eigen(Hermitian(Array(H)))
```

The construction logic is the same as before, but the basis is now symmetry reduced.

## What To Read Next

- [Bases and Sectors](manual/bases.md) explains which basis to use.
- [Operators](manual/operators.md) explains `operator`, `trans_inv_operator`, and `spin`.
- [Maps and Symmetrizers](manual/maps.md) explains how to move between bases.
- [Time Evolution](manual/time-evolution.md) explains how to propagate state
  vectors with the adaptive Krylov/Lanczos API.
- [Lindblad Workflows](manual/lindblad.md) explains how to choose between the
  many-body Arnoldi, legacy dense Lindblad, and quadratic open-system routes.
