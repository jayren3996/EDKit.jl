# Lindblad Reference

These APIs cover EDKit's open-system dynamics layers:

- the adaptive many-body Arnoldi solver for general density-matrix evolution,
- the legacy dense many-body Lindblad stepper and its density-matrix helpers,
- the quadratic covariance-matrix workflow for free or quadratic systems.

For the conceptual routing page, see [Lindblad Workflows](../manual/lindblad.md).
For worked open-system examples, see
[Open-System Workflows](../examples/open-system-workflows.md). For the
closed-system analogue based on Lanczos and state vectors, see
[Time Evolution](../manual/time-evolution.md) and
[Time Evolution Reference](time-evolution.md).

## Many-Body Adaptive Arnoldi APIs

These are the primary APIs for general many-body open-system propagation when
you want adaptive basis reuse, restart/extension logic, and propagation
diagnostics.

```@docs
EDKit.lindblad_timeevolve
EDKit.lindblad_timeevolve!
EDKit.LiouvillianMap
EDKit.LindbladArnoldiCache
EDKit.LindbladArnoldiDiagnostics
```

## Legacy Dense Lindblad And Density-Matrix APIs

These APIs define the original many-body density-matrix workflow. The same
`DensityMatrix` wrapper is also used by the adaptive Arnoldi solver's returned
outputs.

```@docs
EDKit.lindblad
EDKit.Lindblad
EDKit.DensityMatrix
EDKit.densitymatrix
EDKit.expectation
```

## Quadratic Lindblad APIs

These APIs cover the covariance-matrix route for quadratic/Majorana Lindblad
systems. This layer is structurally different from the many-body density-matrix
solvers above.

```@docs
EDKit.quadraticlindblad
EDKit.CovarianceMatrix
EDKit.covariancematrix
EDKit.majoranaform
EDKit.fermioncorrelation
```
