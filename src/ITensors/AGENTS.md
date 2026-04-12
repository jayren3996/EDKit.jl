# AGENTS.md

## Scope

- Applies to vector/MPS/MPO conversions, Pauli-basis helpers, and TEBD
  utilities in this folder.

## Working Rules

- Preserve site ordering exactly across EDKit vectors, ITensor site indices,
  and operator-space conversions.
- Treat EDKit coordinate conventions as the source of truth when round-tripping
  with ITensors.
- When changing conversions, check both plain vector round-trips and
  symmetry-resolved `mps2vec(..., B)` behavior.

## Main APIs

- `vec2mps`, `mps2vec`, `mat2op`, `op2mat`
- `pauli`, `pauli_list`, `commutation_mat`, `dissipation_mat`
- `mps2pmps`, `pmps2mpo`, `mpo2pmpo`, `tebd4`

## Failure Modes

- site ordering drift between EDKit and ITensor indices
- basis mismatch during symmetry-resolved conversion
- normalization drift across density-operator and Pauli-space conversions

## Verification Targets

- primary tests: `test/itensor_tests.jl`, `test/TensorTest.jl`
