# Maps and Symmetrizers

Not every workflow stays inside one basis. Sometimes you want to move amplitudes between a full tensor-product basis and a symmetry-reduced basis, or between two different symmetry sectors that share the same local Hilbert space.

That is what `DoubleBasis` is for.

## `DoubleBasis`

`DoubleBasis(B1, B2)` creates a basis-like object representing maps from coordinates in `B2` to coordinates in `B1`.

Typical uses:

- project a full-space vector into a symmetry sector,
- lift sector amplitudes back into a less reduced basis,
- compare two sector descriptions through their common embedding.

The two bases must have:

- the same chain length,
- the same local on-site base.

The convention is:

- `B1` is the target basis,
- `B2` is the source basis.

So `size(DoubleBasis(B1, B2)) == (dim(B1), dim(B2))`.

## Symmetrizers

The most direct basis-only map is `symmetrizer(B)`, where `B` is a `DoubleBasis`.

It returns the overlap map between the embeddings of the two bases into the full tensor-product basis.

In the common case where `B2` is less symmetric than `B1`, it behaves like a projection or symmetrization matrix into the target sector.

Use `symmetrizer(T)` when you want the explicit matrix. Use `T(v)` when you
only need the action of that map on a vector.

## Example: Full Space To A Symmetry Sector

```julia
using EDKit

L = 8
Bfull = TensorBasis(L = L, base = 2)
Bsector = basis(L = L, N = 4, p = 1)
T = DoubleBasis(Bsector, Bfull)

P = symmetrizer(T)
v_full = randn(size(Bfull, 1))
v_sector = T(v_full)
```

`P` now maps coordinates from the full basis into the chosen symmetry sector,
and `T(v_full)` applies the same overlap map without materializing the matrix.

The reverse construction also makes sense:

```julia
Tback = DoubleBasis(Bfull, Bsector)
v_sector = randn(size(Bsector, 1))
v_full = Tback(v_sector)
```

That gives the overlap-based embedding of sector amplitudes back into the full
space.

## Applying A `DoubleBasis`

`DoubleBasis` itself can also act on a vector, as in the example above:
`v_sector = T(v_full)`.

This is convenient when you want the basis map as an operation, not just as an
explicit matrix. It also avoids building a potentially large dense matrix when
you only need one application.

## How The Map Is Defined

The map is defined through the common embedding of both bases into the full
tensor-product basis. That matters because two reduced bases do not need to
share the same representative set or orbit normalization conventions.

`DoubleBasis` therefore computes overlaps through the underlying physical
product states rather than by matching stored representative indices directly.

For `ProjectedBasis`, off-sector states are treated as zero contributions in the
map construction rather than as hard errors. That behavior is what makes the
projection and embedding workflows practical.

## When To Reach For This Layer

Use `DoubleBasis` and `symmetrizer` when:

- you want to compare states across basis choices,
- you need an explicit projection map,
- you want to verify that reduced-basis calculations agree with full-space calculations,
- or you want to build inter-basis operators.

If you only need a Hamiltonian inside one symmetry sector, you usually do not need this layer at all.
