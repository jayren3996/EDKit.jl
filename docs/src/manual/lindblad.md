# Lindblad Workflows

EDKit includes two different open-system layers:

- a many-body Lindblad representation acting on density matrices,
- a quadratic-fermion covariance-matrix representation for free or quadratic problems.

They solve related physical questions, but they are meant for different regimes.

## Many-Body Lindblad Evolution

Use `lindblad` when you already have a Hamiltonian and jump operators as explicit matrices and want to evolve a density matrix directly.

```julia
using EDKit

H = Array(trans_inv_operator(spin((1.0, "xx"), (1.0, "yy")), 1:2, 4))
Lops = [Array(operator(spin("+" ), [i], 4)) for i in 1:4]

lb = lindblad(H, Lops)
rho0 = densitymatrix(1, 4)
rho1 = lb(rho0, 0.01)
```

This layer is the natural open-system companion to exact diagonalization, but it scales with the full density-matrix dimension.

## Density Matrices And Observables

The helper `densitymatrix` constructs density matrices from:

- an explicit matrix,
- a pure-state vector,
- or a basis-state index.

Use `expectation(O, dm)` to evaluate observables and `entropy(dm)` for density-matrix entropy.

## Quadratic Lindblad Evolution

Use `quadraticlindblad` when the model can be expressed in quadratic Majorana or fermionic form and you want covariance-matrix evolution rather than full density-matrix evolution.

The related helpers are:

- `covariancematrix` to construct covariance matrices,
- `majoranaform` to build Majorana quadratic forms from fermionic data,
- `fermioncorrelation` to recover correlation blocks.

This layer is usually far more scalable than the many-body Lindblad path, but only applies when the model structure permits it.

## Which One Should You Use

Use the many-body `lindblad` path when:

- your system is small enough for explicit density matrices,
- you want a direct, general representation,
- or your problem is not quadratic.

Use `quadraticlindblad` when:

- the model is quadratic,
- covariance-matrix evolution is the natural language,
- you need much larger system sizes than many-body density matrices allow.

## Practical Advice

If you are unsure, start from the physical representation of your problem:

- explicit finite Hilbert-space Hamiltonian and jump matrices suggests `lindblad`,
- quadratic fermion or Majorana data suggests `quadraticlindblad`.
