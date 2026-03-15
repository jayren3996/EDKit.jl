# Basic Workflows

This page collects a few compact workflows that cover the most common first uses of EDKit.

## Build A Hamiltonian From Local Terms

```julia
using EDKit

L = 8

mats = [
    fill(spin("XX"), L - 1);
    fill(spin("YY"), L - 1);
    fill(1.2 * spin("ZZ"), L - 1);
]

inds = vcat(
    [[i, i + 1] for i in 1:L-1],
    [[i, i + 1] for i in 1:L-1],
    [[i, i + 1] for i in 1:L-1],
)

H = operator(mats, inds, L)
```

This is the basic “define terms once, then assemble” pattern used throughout the package.

## Use The Same Model In Different Representations

```julia
using EDKit, LinearAlgebra, SparseArrays

Hdense = Array(H)
Hsparse = sparse(H)
vals = eigvals(Hermitian(Hdense))
```

This is one of the package's main conveniences: the same `Operator` can serve as a matrix-free object or as an explicit matrix.

## Build A Translation-Invariant Model

```julia
using EDKit

L = 10
h2 = spin((1.0, "xx"), (1.0, "yy"), (0.5, "zz"))
H = trans_inv_operator(h2, 1:2, L)
```

Use this when the same local interaction is repeated around a periodic chain.
