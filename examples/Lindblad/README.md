# Lindblad examples

This folder contains notebook-based walkthroughs for the open-system solvers in `EDKit.jl`.

- `DissipativeXXChain.ipynb`: many-body Lindblad evolution for a two-site XX chain with local loss, tracking occupations, entropy, and trace preservation.
- `QuadraticXXLoss.ipynb`: four-site free-fermion XX chain with local loss, comparing `quadraticlindblad` against the full many-body `lindblad` evolution.
- `PauliSuperoperators.ipynb`: Pauli-basis introduction to `pauli`, `pauli_list`, `commutation_mat`, and `dissipation_mat`.

Both notebooks default to loading the local source tree with:

```julia
const DEV = true
```

Set `DEV = false` if you want them to run against an installed `EDKit` package instead.
