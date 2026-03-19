# API Overview

The reference pages group the main EDKit APIs by subsystem.

Use the manual when you want conceptual guidance or workflow advice. Use this
section when you need exact names, signatures, or a quick reminder of the
relevant entry points.

For the most detailed semantics of arguments, return values, invariants, and
internal helper roles, read the source-level docstrings in `src/`. The
Documenter reference pages expose those docstrings, but the source files remain
the authoritative implementation-facing layer.

## Reference Map

- [Bases](bases.md): basis types, indexing, and high-level sector constructors
- [Operators](operators.md): operator construction, local spin helpers, matrix conversion
- [Maps](maps.md): `DoubleBasis` and `symmetrizer`
- [Entanglement](entanglement.md): entropy and Schmidt-related tools
- [ITensors](itensors.md): MPS, MPO, Pauli, and TEBD helpers
- [Lindblad](lindblad.md): open-system and quadratic workflows
- [Utilities](utilities.md): spectral statistics and small helper routines
