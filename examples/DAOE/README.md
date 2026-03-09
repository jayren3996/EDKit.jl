# DAOE examples

This folder contains the MPS-based operator-growth notebooks for `daoe` and `fdaoe`.

- `XXZOperatorGrowth.ipynb`: short-chain XXZ/Heisenberg benchmark comparing exact operator growth, raw Pauli-MPS TEBD, and DAOE-filtered TEBD.
- `XXMajoranaGrowth.ipynb`: short-chain XX benchmark comparing exact free-fermion string growth against `daoe` and `fdaoe`.

Both notebooks default to loading the local source tree:

```julia
const DEV = true
```

If `Plots` is installed in the active Julia environment, the notebooks render figures. Otherwise they still run and return the underlying arrays and summaries.
