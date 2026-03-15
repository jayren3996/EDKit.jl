# Open-System Workflows

EDKit offers both direct many-body Lindblad evolution and quadratic covariance-matrix evolution.

## Many-Body Lindblad Example

```julia
using EDKit

L = 4
H = Array(trans_inv_operator(spin((1.0, "xx"), (1.0, "yy")), 1:2, L))
J = [Array(operator(spin("+"), [i], L)) for i in 1:L]

lb = lindblad(H, J)
rho0 = densitymatrix(1, L)
rho1 = lb(rho0, 0.01)
```

This is the most direct route when you already have explicit Hamiltonian and jump matrices.

## Quadratic Covariance-Matrix Example

```julia
using EDKit

A = [0.0 1.0; 1.0 0.0]
B = zeros(2, 2)
Hmaj = majoranaform(A, B)

L = zeros(4, 4)
ql = quadraticlindblad(Hmaj, L)
cm = covariancematrix([1, 0])
cm_next = ql(cm, 0.05)
```

This route is appropriate when the problem is naturally quadratic and you want the scaling advantages of covariance-matrix evolution.
