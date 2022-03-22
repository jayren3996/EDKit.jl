# Basis

In `EDKit.jl`, the fundamental objects are basis and operator. The `AbstractBasis` is the abstract type of basis. Currently there are 4 concrete basis:

1. `TensorBasis`: Ordinary basis without any symmetry.
2. `ProjectedBasis`: Basis for subspace that is spanned only by product states.
3. `TranslationalBasis`: Basis for translational symmetric Hamiltonian.
4. `TranslationalParityBasis` : Basis with translational as well as reflection symmetry. The momensum should be 0 or π.

The basis object can be extended. To construct linear operation, we need to define 3 functions for a new basis type:

1. `size(b::AbstractBasis)` : Size of matrix representations of the operators in this subspace.
2. `change!(b::AbstractBasis, i::Integer)` : Return the normalization of ith state and change the digits to ith states in this subspace.
3. `index(b::AbstractBasis)` : Return the coefficient and index of the digits.

Optionally, we can define `eltype` for a basis object (default is `ComplexF64`).

If the calculation is done on the entire Hilbert space, the basis object need not be explicitly constructed. The `Operator` will use `TensorBasis` by default. The construction of other basis with symmetry concern are discussed below.

In addition, if the entaglement entropy is needed, the user-defined basis should implement a function `schmidt!(target, v, Ainds, b::AbstractBasis)`.

Here we introduce 3 concrete implementation of `AbstractBasis`.

## TensorBasis

The type `TensorBasis` has the fields:

```julia
struct TensorBasis <: AbstractBasis
    dgt::Vector{Int}
    B::Int
end
```

We can use the function

```julia
tensorbasis(L::Integer; base::Integer=2)
```

to construct a basis, though in most cases it is not necessary.

## ProjectedBasis

The type `ProjectedBasis` has the fields:

```julia
struct ProjectedBasis <: AbstractBasis
    dgt::Vector{Int}
    I::Vector{Int}
    B::Int
end
```

We can use the function

```julia
projectedbasis(f, L::Integer; base::Integer=2)
```

In the definition, `f` is a function acting on digits that tells whether a given digits is valid in this basis. For eaxample, consider a S=1/2 chain with L=6 (conserve megnetic quantum number). If we consider the subspace spaned by those states with 3 up-spins, we can create the basis for the subspace by

```julia
projectedbasis(x -> sum(x)==3, 6; base=2)
```

## TranslationalBasis

The type `TranslationalBasis` has the fields:

```julia
struct TranslationalBasis <: AbstractBasis
    dgt::Vector{Int}
    I::Vector{Int}
    R::Vector{Float64}
    C::ComplexF64
    B::Int
end
```

We can use the function

```julia
translationalbasis(k::Integer, L::Integer; base::Integer=2)
translationalbasis(f, k::Integer, L::Integer; base::Integer=2)
```

In the definition, `k` is the momentum sector we are interested, and `f` is a function acting on digits that tells whether a given digits is valid in this basis. Note that `k` is an integer, representing the momentum at `2πk/L`. For example, consider a `S=1/2` chain with `L=6` (with translational symmetry and conserve magnetic quantum number). If we consider the subspace spaned by those states with 3 up-spins with zero momentum, we can create the basis for the subspace by

```julia
translationalbasis(x -> sum(x)==3, 0, 6; base=2)
```

