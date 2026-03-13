# Examples

The `examples/` directory is now organized by topic. Start with the folder that matches the workflow you want to learn.

- `Basic/`: first-use walkthroughs for operator construction, symmetry reduction, and simple MPS/Pauli conversions.
- `Maps/`: basis-to-basis maps, `DoubleBasis`, `symmetrizer`, and rectangular operators.
- `Symmetries/`: sector-construction examples, Abelian basis combinations, and MPS projection workflows.
- `TensorNetworks/`: tensor-network and ITensor-oriented examples beyond the basic introduction.
- `Lindblad/`: open-system examples for both many-body and quadratic Lindblad evolution.
- `Models/`: model-specific examples where multiple package features come together in one physics problem.

Suggested starting points:

- `Basic/OperatorConstruction.ipynb`
- `Basic/SymmetryReduction.ipynb`
- `Basic/MPSAndPauli.ipynb`
- `Maps/DoubleBasisBasics.ipynb`
- `Maps/Symmetrizers.ipynb`
- `Lindblad/PauliSuperoperators.ipynb`
- `Symmetries/SectorCatalogue.ipynb`
- `Lindblad/DissipativeXXChain.ipynb`

Each script or notebook defaults to loading the local source tree:

```julia
const DEV = true
```

Set `DEV = false` if you want to run them against an installed `EDKit` package instead.
