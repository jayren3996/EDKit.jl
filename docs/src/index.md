# EDKit.jl

Julia package for general many-body exact diagonalization calculation for spin systems. 

All information we need to specify a *local operator* in the many-body Hilbert space are:

1. Matrix form of the local operators;
2. Indices the local operator act on;
3. Basis for the many-body operator.

For example, the AKLT model (with size ``L=8``)
```math
H = \sum_{i=1}^{8}\left[\vec S_i \cdot \vec S_{i+1} + \frac{1}{3}\left(\vec S_i \cdot \vec S_{i+1}\right)^2\right],
```
can be easyly constructed using the command:
```julia
L = 8
SS = spin((1, "xx"), (1, "yy"), (1, "zz"), D=3)
mat = SS + 1/3 * SS^2
H = trans_inv_operator(mat, 1:2, L)
```

The `EDKit.jl` provide a general function `operator(mats, inds, basis)` which helps to create local operator. Especially when we are doing exact diagonalization calculation on a specific *symmetry sector* (e.g., sector with total ``S^z=0`` and total momentum ``k=0``). The functionalities can be extended with user-defined bases.

## Installation

Run the following script in the ```Pkg REPL``` environment:

```julia
pkg> add EDKit
```

