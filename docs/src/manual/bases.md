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

## The Digit Buffer and Thread Safety

Every basis object carries a mutable digit buffer `b.dgt`, a `Vector{Int}` of length `L` that stores the local state label for each site. The two core operations — `change!` (decode a basis index into digits) and `index` (interpret digits as basis coordinates) — read from and write to this buffer.

Because `b.dgt` is shared mutable state, naive concurrent use of the same basis across threads would cause data races. EDKit solves this through **explicit buffer passing**: all hot-path functions accept an optional `dgt::AbstractVector` argument that replaces `b.dgt`:

```julia
# Thread-safe variants (used internally by mul, *, schmidt, etc.)
change!(b, i, dgt)        # decode into external buffer
index(b, dgt)             # interpret external buffer
colmn!(target, M, I, b, dgt, coeff)  # apply local term using external buffer
```

Inside EDKit, every function that iterates over basis states allocates a thread-local buffer at the top of its loop:

```julia
dgt = similar(b.dgt)    # one allocation, reused for all iterations
for j in 1:length(v)
    colmn!(target, opt, j, dgt, v[j])
end
```

In the multi-threaded `mul(H, psi)`, each thread allocates its own `dgt`:

```julia
Threads.@threads for i in 1:nthreads
    dgt = similar(opt.B.dgt)   # thread-local, no sharing
    for j in range_for_thread_i
        colmn!(Ms[i], opt, j, dgt, v[j])
    end
end
```

This design means:

- **Thread-safe by construction** — no locks, no atomics, no synchronization overhead.
- **Zero performance cost** — benchmarks show the explicit-buffer path is slightly faster than the old struct-based path (~6%), because the compiler can optimize local variables better than struct field accesses.
- **Backward compatible** — `b.dgt` still exists and the old no-argument forms of `index(b)` and `change!(b, i)` still work by delegating to the buffer-passing variants with `b.dgt`.

!!! note "For custom code using basis objects"
    If you write your own loops over basis states and want thread safety, allocate
    a local buffer with `dgt = similar(b.dgt)` and use the `dgt`-accepting
    overloads. The built-in `mul`, `*`, `addto!`, and `schmidt` functions already
    do this automatically.

## Next Step

Once the basis is chosen, the next task is usually to build an `Operator`. See [Operators](operators.md).
