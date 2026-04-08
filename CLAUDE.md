# CLAUDE.md

## Project Overview

EDKit.jl is a Julia package for exact diagonalization of quantum many-body systems. It provides symmetry-resolved basis construction, operator assembly from local terms, entanglement diagnostics, ITensor conversions, and open-system workflows.

## Key Architecture

- `src/Basis/` -- Basis types (TensorBasis, ProjectedBasis, TranslationalBasis, ParityBasis, FlipBasis, etc.)
- `src/Operator.jl` -- Central `Operator` type: stores local sparse matrices + site indices + basis; supports matrix-free and cached-sparse multiplication
- `src/LinearMap.jl` -- Inter-basis maps (DoubleBasis, symmetrizer)
- `src/Schmidt.jl` -- Entanglement and Schmidt decomposition
- `src/ToolKit.jl` -- Utilities (gap ratios, exponentials)
- `src/ITensors/` -- ITensor/MPS/MPO integration
- `src/algorithms/` -- Lindblad and QIM solvers

## Build and Test

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Julia binary may be at `~/.juliaup/bin/julia`. The project requires Julia >= 1.9.

## Performance Notes

### Operator multiplication

`Operator * matrix` has two paths:

1. **Matrix-free** (default): Iterates over basis states, applying local terms via digit manipulation. Cost is dominated by `index()` and `change!()` calls in `src/Basis/`. TranslationalBasis is the most expensive (O(L^2) orbit search per state).

2. **Cached sparse** (after `sparse!(opt)`): Stores the explicit `SparseMatrixCSC` in an LRU cache. Subsequent `opt * matrix` and `mul(opt, matrix)` calls use SparseArrays SpMM, which is 10-1000x faster. Vector multiplication (`opt * v`) always uses the matrix-free path regardless.

When writing code that multiplies the same operator by matrices repeatedly, always call `sparse!(opt)` first. Call `clear_sparse_cache!()` when done to release memory. For very large systems (dim > 100k), the sparse matrix may be too large to cache -- in that case use the matrix-free path.

### Hot path functions

The innermost loops in `src/Basis/AbstractBasis.jl` (`index`, `change!`) and `src/Operator.jl` (`colmn!`) use `@inbounds` for performance. When modifying these functions, ensure array accesses remain structurally valid.

## Conventions

- All basis constructors use keyword arguments: `TensorBasis(L=10, base=2)`, `ProjectedBasis(L=12, base=2, f=...)`, `TranslationalBasis(L=12, base=2, k=0, f=...)`
- `Operator` is an immutable struct; the sparse cache is stored externally in a module-level LRU
- Digit buffers (`b.dgt`) are mutable working state inside basis objects -- each thread needs its own copy
- Tests are in `test/` and run via `Pkg.test()`; no separate test runner needed

## 2D/3D Lattice Support

`AbelianBasis` supports arbitrary lattice geometries via the `symmetries` keyword in `basis()`. Each symmetry generator is a `(perm, quantum_number)` tuple where `perm` is a 1-indexed permutation array:

```julia
# 2D square lattice example
Lx, Ly = 4, 3; L = Lx * Ly
sites = [(x, y) for y in 0:Ly-1 for x in 0:Lx-1]
T_x = [mod(x+1, Lx) + Lx*y + 1 for (x,y) in sites]
T_y = [x + Lx*mod(y+1, Ly) + 1 for (x,y) in sites]
B = basis(; L, N=L÷2, symmetries=[(T_x, 0), (T_y, 0)])
```

For base=2 systems, permutation networks provide fast integer-state manipulation and Gosper's hack accelerates particle-conserving enumeration. See `docs/src/abelian_basis.md` for full documentation.
