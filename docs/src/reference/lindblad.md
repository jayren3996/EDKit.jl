# Lindblad Reference

These APIs cover both open-system layers in EDKit:

- the many-body density-matrix path built around `lindblad` and `DensityMatrix`,
- the quadratic covariance-matrix path built around `quadraticlindblad` and
  `CovarianceMatrix`.

The many-body layer is the general exact-diagonalization workflow. The
quadratic layer is the scalable special case for free or quadratic fermionic
systems.

```@docs
lindblad
EDKit.densitymatrix
EDKit.expectation
EDKit.quadraticlindblad
EDKit.covariancematrix
majoranaform
EDKit.fermioncorrelation
EDKit.Lindblad
EDKit.DensityMatrix
EDKit.CovarianceMatrix
```
