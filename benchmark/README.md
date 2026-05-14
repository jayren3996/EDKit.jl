# Benchmarks

This folder is for local performance checks, not mandatory CI tests.

The package does not add `BenchmarkTools` as a dependency for these notes. For
repeatable local runs, start Julia with the project environment and add
temporary benchmark tooling in a throwaway environment if needed.

Suggested first-pass measurements:

- `index(dgt; base=2)` and `index(dgt; base>2)`
- `change!(dgt, ind; base=2)` and `change!(dgt, ind; base>2)`
- `operator(mats, inds, B)` with many repeated site patterns
- `Operator * vector`
- `Operator * matrix`
- `EDKit.mul(H, vector)` and `EDKit.mul(H, matrix)`
- `sparse!(H)` followed by repeated `H * matrix`
- `basis(; L, N, symmetries=...)` for small custom two-dimensional symmetries
- `schmidt` and `ent_S` on small translational and Abelian sectors

The TensorBasis base-2 operator application benchmark is:

```sh
julia -e 'using Pkg; Pkg.activate(; temp=true); Pkg.add("BenchmarkTools"); Pkg.develop(path=pwd()); include("benchmark/tensor_base2_operator_apply.jl")'
```

Record Julia version, thread count, exact command, and representative output in
any benchmark note or PR description. Do not claim a speedup without before and
after measurements from the same machine and Julia session shape.
