# Operators

The operator layer is the heart of EDKit.

Instead of asking you to build a full matrix immediately, EDKit lets you specify local terms and where they act, then stores the result as an `Operator`.

## Why `Operator` Matters

An `Operator` keeps three things together:

- the local matrices,
- the sites each local matrix acts on,
- the basis in which the operator should act.

That means the same object can support:

- matrix-free application with `H * psi`,
- explicit dense conversion with `Array(H)`,
- explicit sparse conversion with `sparse(H)`.

## Constructing Operators

The main constructor is `operator(mats, inds, basis_or_length)`.

Example:

```julia
using EDKit

L = 6
h2 = spin((1.0, "xx"), (1.0, "yy"), (0.5, "zz"))

mats = fill(h2, L - 1)
inds = [[i, i + 1] for i in 1:L-1]

H = operator(mats, inds, L)
```

Repeated site patterns are merged automatically, which is convenient when you build a Hamiltonian by concatenating several contributions.

## Translation-Invariant Construction

For periodic translation-invariant models, `trans_inv_operator` is usually the cleanest interface:

```julia
using EDKit

L = 10
h2 = spin((1.0, "xx"), (1.0, "yy"), (1.0, "zz"))
H = trans_inv_operator(h2, 1:2, L)
```

This takes one seed local term and translates it around the ring with periodic boundary conditions.

## The `spin` Helper

`spin` is the standard way to build local matrices.

Examples:

```julia
spin("X")
spin("zz")
spin((1.0, "xx"), (1.0, "yy"), (0.3, "zz"))
```

For spin-1/2 systems, uppercase Pauli-style labels like `"X"` and `"Z"` are supported, along with lowercase operator strings such as `"x"`, `"y"`, `"z"`, `"+"`, `"-"`, and `"1"`.

Multi-site strings are interpreted as Kronecker products in the written order.

## Acting With Operators

For many workflows you do not need the full matrix:

```julia
using EDKit, LinearAlgebra

L = 8
H = trans_inv_operator(spin("zz"), 1:2, L)
psi = normalize(randn(ComplexF64, 2^L))

Hpsi = H * psi
```

If you want threaded multiplication, use `mul(H, psi)`.

## Converting To Matrices

Use dense matrices when exact diagonalization is affordable:

```julia
Hdense = Array(H)
vals = eigvals(Hermitian(Hdense))
```

Use sparse matrices when you need explicit sparse storage or want to call sparse linear algebra routines:

```julia
using SparseArrays

Hsparse = sparse(H)
```

If you need control over the output buffer, `addto!` writes into an existing matrix.

## Working In Symmetry Sectors

The same construction functions work with reduced bases:

```julia
using EDKit

B = basis(L = 12, N = 6, k = 0)
h2 = spin((1.0, "xx"), (1.0, "yy"), (1.0, "zz"))
H = trans_inv_operator(h2, 1:2, B)
```

This is one of the main advantages of the EDKit design: the model description stays the same while the basis changes.

## When To Use Which Interface

- Use `spin` to define local building blocks.
- Use `operator` when terms act on specific site lists.
- Use `trans_inv_operator` when one local term is repeated by translation symmetry.
- Use `Array(H)` or `sparse(H)` only when you really want an explicit matrix.
- Use `H * psi` or `mul(H, psi)` when you want to stay matrix free.
