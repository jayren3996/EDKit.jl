# Examples

The notebooks in this folder are the main walkthroughs for the package.

Basic notebook examples:

- `Basic/OperatorConstruction.ipynb`: construct a nearest-neighbor spin Hamiltonian, apply it as a linear map, and compare dense and sparse representations.
- `Basic/SymmetryReduction.ipynb`: build the same model in a full basis and reduced symmetry sectors, then check parity-sector recombination.
- `Basic/MPSAndPauli.ipynb`: convert vectors to MPS, move to Pauli-space MPS/MPO objects, and inspect bond dimensions.

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
