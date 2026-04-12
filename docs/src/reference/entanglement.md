# Entanglement Reference

These APIs cover the vector-and-basis side of EDKit entanglement:

- `schmidt` builds the bipartite matrix,
- `ent_spec` returns Schmidt singular values,
- `ent_S` returns entropies for EDKit basis states and full-space vectors,
- `EDKit.entropy` is the shared entropy helper used underneath.

For MPS bond-cut helpers such as `ent_S!(psi, b)` and `ent_specs!(psi, b)`, see
[ITensors](itensors.md).

```@docs
ent_S
EDKit.ent_spec
EDKit.entropy
EDKit.schmidt
```
