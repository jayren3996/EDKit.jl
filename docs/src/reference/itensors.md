# ITensor Reference

These APIs connect EDKit models and states with ITensor and Pauli-space
workflows.

They include:

- vector/MPS and matrix/operator conversion,
- MPS constructors and MPS entanglement helpers,
- Pauli-space superoperator tools,
- TEBD gate construction.

```@docs
EDKit.vec2mps
mps2vec
mat2op
op2mat
pbcmps
productstate
EDKit.ent_specs
EDKit.ent_specs!
EDKit.ent_S!
pauli
pauli_list
commutation_mat
dissipation_mat
mps2pmps
pmps2mpo
mpo2pmpo
tebd_n!
EDKit.tebd4
```
