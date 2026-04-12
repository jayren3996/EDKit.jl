# AGENTS.md

## Scope

- Applies to `Lindblad.jl` and `QIM.jl`.

## Solver Split

- `lindblad`, `DensityMatrix`, and `expectation` cover explicit many-body
  density-matrix workflows.
- `quadraticlindblad`, `CovarianceMatrix`, `majoranaform`, and
  `fermioncorrelation` cover quadratic/Majorana covariance workflows.
- `qimsolve` and `covmat` are inverse-method utilities; do not casually mix
  their assumptions with Lindblad code.

## Editing Rules

- Keep input/output object conventions stable.
- Be careful with normalization, Hermiticity, and real-versus-complex
  assumptions.
- Verify both many-body and quadratic paths when touching shared physical
  semantics.

## Verification Targets

- primary tests: `test/lindblad_tests.jl`, `test/advanced_tests.jl`
- also relevant: examples and docs when changing solver-facing APIs
