# Basic examples

This folder contains small notebook-based examples for the core EDKit workflows.

- `OperatorConstruction.ipynb`: construct a many-body operator from local terms and compare linear-map, dense, and sparse representations.
- `SymmetryReduction.ipynb`: compare full-space and symmetry-reduced constructions and check parity-sector recombination.
- `MPSAndPauli.ipynb`: basic `vec2mps`/`mps2vec` usage plus Pauli-space MPS and MPO conversion.
- `Basics.ipynb`: older broad introductory notebook that mixes several first-use workflows.

All notebooks default to loading the local source tree with:

```julia
const DEV = true
```

Set `DEV = false` if you want them to run against an installed `EDKit` package instead.
