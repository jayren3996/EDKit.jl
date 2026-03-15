# Bases and Sectors

Basis objects are the foundation of EDKit. They tell the package how to interpret digit strings, how large the Hilbert space is, and how local operators should be embedded into that space.

## The Full-Space Starting Point

`TensorBasis` is the simplest basis:

```julia
using EDKit

B = TensorBasis(L = 6, base = 2)
size(B)
```

This is the full tensor-product basis with no symmetry reduction. It is the right starting point when:

- system sizes are still modest,
- you want the most transparent interpretation,
- or you want to compare a reduced basis against the full space.

## Projected Subspaces

`ProjectedBasis` keeps only basis states that satisfy a predicate or a fixed charge sector.

Typical use cases:

- fixed particle number or magnetization,
- constrained Hilbert spaces such as PXP/Rydberg blockade,
- custom local constraints.

Examples:

```julia
using EDKit

B_fixed = ProjectedBasis(L = 8, N = 4)
B_pxp = ProjectedBasis(L = 8, f = x -> all(x[i] + x[i + 1] <= 1 for i in 1:7))
```

Use `ProjectedBasis` when you want a physically meaningful subspace but do not need momentum or parity quantum numbers.

## Translation Sectors

`TranslationalBasis` organizes states by lattice momentum.

```julia
using EDKit

B = TranslationalBasis(L = 12, k = 0, N = 6)
```

This is useful when:

- the Hamiltonian is translation invariant,
- you want block-diagonalization by momentum,
- or you want representative-orbit states rather than full product states.

The keyword `a` controls the unit-cell length. When `a > 1`, translations are taken by one unit cell instead of one site.

## Reflection And Flip Sectors

EDKit also provides discrete symmetry bases:

- `ParityBasis` for reflection parity,
- `FlipBasis` for spin-flip symmetry,
- `ParityFlipBasis` for simultaneous reflection and flip sectors,
- `TranslationParityBasis` and `TranslationFlipBasis` for translation plus a discrete symmetry.

These are useful when the model preserves the corresponding symmetry and you want a smaller Hilbert space or a cleaner quantum-number classification.

## High-Level Constructor

For most users, the easiest entry point is `basis(...)`:

```julia
using EDKit

B = basis(L = 12, N = 6, k = 0, p = 1)
```

This helper chooses the concrete basis type for you based on the requested symmetries and constraints.

It is a good default choice when:

- you know the quantum numbers you want,
- you do not need to manually control the concrete basis type,
- you want one consistent entry point across many workflows.

## Choosing A Basis

Use this rule of thumb:

| Goal | Recommended basis |
| --- | --- |
| No symmetry reduction, simplest interpretation | `TensorBasis` |
| Fixed charge or custom local constraint | `ProjectedBasis` |
| Momentum sector | `TranslationalBasis` |
| Reflection parity only | `ParityBasis` |
| Spin-flip only | `FlipBasis` |
| Multiple commuting symmetries | `basis(...)` |
| Explicit basis-to-basis map construction | a pair of bases plus `DoubleBasis` |

## Basis Internals You Will See Often

Even if you do not implement new basis types, a few operations matter to users:

- `size(B)` gives the Hilbert-space dimension,
- `change!(B, i)` loads the `i`th basis representative into the working digit buffer,
- `index(B)` returns the coefficient and position associated with the current digit configuration,
- `content(B, i)` returns the stored representative index.

Those low-level operations are what allow operator application, basis conversion, and Schmidt decomposition to work uniformly across symmetry sectors.

## Next Step

Once the basis is chosen, the next task is usually to build an `Operator`. See [Operators](operators.md).
