# Examples

The notebooks in this folder are the main walkthroughs for the package. The scripts below are smaller terminal-friendly examples that exercise the same API without requiring Jupyter.

Additional ITensor notebooks:

- `DAOE/XXZOperatorGrowth.ipynb`: short-chain XXZ/Heisenberg benchmark comparing exact operator growth, raw Pauli-MPS TEBD, and DAOE-filtered TEBD.
- `DAOE/XXMajoranaGrowth.ipynb`: short-chain XX benchmark comparing exact free-fermion string growth against `daoe` and `fdaoe`.

- `ConstrainedPXP.jl`: projected Hilbert spaces, `productstate`, diagonalization, and entanglement entropy.
- `GapRatioXXZ.jl`: random-field spin chain in a symmetry sector and level-statistics analysis with `meangapratio`.
- `AbelianBasisSectors.jl`: the high-level `basis(...)` helper for combining `N`, `k`, and `p` quantum numbers.

Each script defaults to loading the local source tree:

```julia
const DEV = true
```

Set `DEV = false` if you want to run them against an installed `EDKit` package instead.
