# Symmetry Examples

Each `.jl` file in this folder constructs the same symmetry-preserving XXZ Hamiltonian,

```math
H = \sum_i \left(S_i^xS_{i+1}^x + S_i^yS_{i+1}^y + 0.6 S_i^zS_{i+1}^z\right),
```

but in a different symmetry sector selected with `basis(...)`.

- `01_N.jl`: `N`
- `02_k.jl`: `k`
- `03_p.jl`: `p`
- `04_z.jl`: `z`
- `05_Nk.jl`: `N + k`
- `06_Np.jl`: `N + p`
- `07_Nz.jl`: `N + z`
- `08_kp.jl`: `k + p`
- `09_kz.jl`: `k + z`
- `10_pz.jl`: `p + z`
- `11_Nkp.jl`: `N + k + p`
- `12_Nkz.jl`: `N + k + z`
- `13_Npz.jl`: `N + p + z`
- `14_kpz.jl`: `k + p + z`
- `15_Nkpz.jl`: `N + k + p + z`

All files share `_common.jl`, which defines the model and a small helper for printing sector size and ground-state energy.

Compatibility rules:

- `N + z` requires half filling for spin-1/2 systems.
- `k + p` requires `k = 0` or `k = L/2`.
- Any combination containing both `k` and `p` uses `k = 0` here for that reason.
