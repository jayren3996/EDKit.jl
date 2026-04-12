# Utilities Reference

These APIs collect smaller but frequently useful helper routines:

- spectral statistics (`gapratio`, `meangapratio`),
- lightweight exponential helpers (`EDKit.expm`, `EDKit.expv`),
- and the quantum inverse method tools (`covmat`, `qimsolve`).

The basis-side `productstate(v, B)` helper is described in the manual
[Utilities](../manual/utilities.md), while the ITensor-dispatched `productstate`
constructor appears in [ITensors](itensors.md).

```@docs
gapratio
meangapratio
EDKit.expm
EDKit.expv
EDKit.covmat
EDKit.qimsolve
```
