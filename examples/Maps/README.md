# Maps

This folder is about basis-to-basis linear maps in `EDKit.jl`.

- `DoubleBasisBasics.ipynb`: what `DoubleBasis(Btarget, Bsource)` stores, how to read the map direction, and how the fast `T(v)` action matches `symmetrizer(T) * v`.
- `Symmetrizers.ipynb`: project full-space states into symmetry sectors and lift them back to the full Hilbert space.
- `InterbasisOperators.ipynb`: build non-square operators with `DoubleBasis` and interpret their matrix shapes.

Recommended order:

1. `DoubleBasisBasics.ipynb`
2. `Symmetrizers.ipynb`
3. `InterbasisOperators.ipynb`
