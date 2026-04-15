# EDKit.jl

`EDKit.jl` is a Julia package for exact diagonalization, symmetry-resolved many-body calculations, and selected tensor-network and open-system workflows.

This documentation site is the human-facing manual for the package. It focuses
on concepts, workflows, and examples. If you are reading the code directly, the
source-level docstrings in `src/` are intentionally more exhaustive and more
implementation-oriented.

It is built around a simple pattern:

1. describe local terms once,
2. choose the basis or symmetry sector you care about,
3. construct an operator,
4. use that same object as a linear map, a dense matrix, or a sparse matrix.

## What EDKit Covers

EDKit has a few major subsystems:

- basis construction for full Hilbert spaces, projected subspaces, and symmetry sectors,
- arbitrary commuting lattice symmetries, including 2D and 3D permutation-defined sectors,
- many-body operator construction from local terms,
- inter-basis maps and symmetrizers,
- entanglement and Schmidt-decomposition helpers,
- closed-system real-time evolution of state vectors via an adaptive Krylov/Lanczos propagator,
- ITensor/MPS utilities, including Pauli-space workflows,
- Lindblad and quadratic-fermion open-system tools.

If you are new to the package, start with [Getting Started](getting-started.md), then move to the manual pages for the subsystem you need.

## Installation

Install the registered package with Julia's package manager:

```julia
pkg> add EDKit
```

Or use the GitHub repository directly:

```julia
pkg> add https://github.com/jayren3996/EDKit.jl
```

Current package compat:

- Julia `1.9+`
- `ITensors.jl` `0.7` to `0.9`
- `ITensorMPS.jl` `0.3`

## Choose A Path

### I want to build and diagonalize many-body Hamiltonians

Start with:

- [Getting Started](getting-started.md)
- [Operators](manual/operators.md)
- [Bases and Sectors](manual/bases.md)

### I want to reduce into symmetry sectors

Start with:

- [Bases and Sectors](manual/bases.md)
- [General Abelian Symmetries](abelian_basis.md)
- [Maps and Symmetrizers](manual/maps.md)
- [Symmetry Workflows](examples/symmetry-workflows.md)

### I want closed-system real-time dynamics

Start with:

- [Getting Started](getting-started.md)
- [Time Evolution](manual/time-evolution.md)
- [Operators](manual/operators.md)
- [Time Evolution Workflows](examples/time-evolution-workflows.md)

### I want to work with ITensors, MPS, or Pauli-space representations

Start with:

- [ITensor Workflows](manual/itensors.md)
- [Tensor-Network Workflows](examples/tensor-network-workflows.md)

### I want open-system or quadratic-fermion evolution

Start with:

- [Lindblad Workflows](manual/lindblad.md)
- [Open-System Workflows](examples/open-system-workflows.md)

That path now covers four cases:

- closed-system Krylov/Lanczos dynamics belongs on [Time Evolution](manual/time-evolution.md),
- general many-body open-system dynamics at many times belongs on
  [the adaptive many-body Lindblad path](manual/lindblad.md),
- small dense many-body step loops belong on
  [the legacy dense Lindblad path](manual/lindblad.md),
- quadratic/free-fermion open systems belong on
  [the quadratic covariance-matrix path](manual/lindblad.md).

## Design Philosophy

EDKit keeps model construction separate from representation choice.

That separation matters in practice:

- you can stay matrix-free when systems are too large for explicit diagonalization,
- convert to dense arrays when exact eigensolvers are appropriate,
- convert to sparse matrices for sparse linear algebra or interoperability,
- feed the same matrix-free `Operator` into [`timeevolve`](manual/time-evolution.md) for closed-system real-time dynamics,
- and reuse the same local terms across different bases and symmetry sectors.

## Manual Layout

The docs are organized in three layers:

- the manual explains concepts and workflows for human readers,
- the worked examples show compact end-to-end use cases,
- the API reference lists the core functions and types grouped by subsystem,
- the source docstrings in `src/` carry the most detailed function-by-function
  semantic documentation.

The package is a single Julia module, `EDKit`, but it contains several functional layers. The [Architecture](manual/architecture.md) page explains how those layers fit together.
