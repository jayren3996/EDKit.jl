# Symmetry Examples

This folder collects examples centered on symmetry-resolved bases and sector projections.

- `Basics.ipynb`: introductory symmetry notebook.
- `MPSProjection.ipynb`: symmetry projection workflow with MPS objects.
- `AbelianBasisSectors.ipynb`: compact notebook using the high-level `basis(...)` helper and checking sector recombination.
- `SectorCatalogue.ipynb`: notebook replacement for the old numbered sector scripts, showing representative symmetry combinations and compatibility rules.

Compatibility rules:

- `N + z` requires half filling for spin-1/2 systems.
- `k + p` requires `k = 0` or `k = L/2`.
- Any combination containing both `k` and `p` uses `k = 0` here for that reason.
