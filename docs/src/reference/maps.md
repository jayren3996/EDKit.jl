# Maps Reference

These APIs handle linear maps between bases and symmetry embeddings.

`DoubleBasis(B1, B2)` always means "target basis `B1`, source basis `B2`".
Use the callable object when you want the action on a vector, and use
`symmetrizer` when you want the explicit overlap matrix.

```@docs
DoubleBasis
symmetrizer
```
