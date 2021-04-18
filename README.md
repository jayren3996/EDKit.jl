# EDKit.jl

 Julia package for general many-body exact diagonalization calculation. The package provide a general Hamiltonian constructing routine for specific symmetry sectors. The functionalities can be extended with user-defined bases.

## Installation & Update

Run the following script in the ```Pkg REPL``` :

```julia
pkg> add EDKit
```

To update the package, run:

```julia
pkg> update EDKit
```

## Basis Object

In `EDKit.jl`, the fundamental objects are basis and operator. The `AbstractBasis` is the abstract type of basis. Currently there are 3 concrete basis:

1. `TensorBasis`: This is ordinary basis without any symmetry.
2. `ProjectedBasis`: This is a basis for subspace that is spanned only by product states.
3. `TranslationalBasis`: This is a basis for translational symmetric Hamiltonian.

The basis object can be expanded, by defining 3 functions:

- `size(b::AbstractBasis)`;
- `change!(b::AbstractBasis, i::Integer)`;
- `index(b::AbstractBasis)`.

Optionally, we can define `eltype` for a basis object (default is `ComplexF64`).

If the calculation is done on the entire Hilbert space, the basis object need not be explicitly constructed. The `Operator` will use `TensorBasis` by default. The construction of other basis with symmetry concern are discussed below.

## Operator Object

In `EDKit.jl` , a many-body operator is represented by the type `Operator`:

```julia
struct Operator{Tv<:Number, Tb<:AbstractBasis}
    M::Vector{SparseMatrixCSC{Tv, Int}}
    I::Vector{Vector{Int}}
    B::Tb
end
```

In this definition, `M` is the list of matrix representations of local operators, `I` is the list of indices of sites it acts on.

### Construction

To construct an `Operator` object, we need 3 inputs:

```julia
M <: AbstractVector{<:AbstractMatrix}
I <: AbstractVector{<:AbstractVector{<:Integer}}
B <: AbstractBasis
```

where `M` is the list of matrices representing the local operators, `I` is the list of vectors representing the sites it acts on. `B` is a basis object. If use `TensorBasis`, we can directly using the constructing method

```julia
operator(mats::AbstractVector{<:AbstractMatrix}, inds::AbstractVector{<:AbstractVector}, L::Integer)
```

### Convert to matrix

An `Operator` object is basically like a matrix, and it can be convert to dense matrix using the function

```julia
Array(opt::Operation)
```

It can also write to a given matrix with correct dimension using function

```julia
addto!(M::AbstractMatrix, opt::Operator)
```

Note that to get correct answer, `M` should de initialized as a zero matrix.

### Multiply to vector or matrix

We can directly using

```julia
O::Operator * M::AbstractVecOrMat
```

to do the multiplycation. Or, use the function

```julia
mul!(target::AbstractVecOrMat, opt::Operator, v::AbstractVecOrMat)
```

to modify `target` (similarly, `target` should be initialized as a zero vector/matrix).

## Concrete Implementations of Basis

Here we introduce 3 concrete implementation of `AbstractBasis`.

### TensorBasis

The type `TensorBasis` has the fields:

```julia
struct TensorBasis <: AbstractBasis
    dgt::Vector{Int}
    B::Int
end
```

We can use the function

```julia
tensorbasis(L::Integer; base::Integer=2) = TensorBasis(zeros(Int, L), base)
```

to construct a basis, though in most cases it is not necessary.

### ProjectedBasis

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

### TranslationalBasis

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

In the definition, `k` is the momentum sector we are interested, and `f` is a function acting on digits that tells whether a given digits is valid in this basis. Note that `k` is an integer, representing the momentum at `2πk/L`. For example, consider a S=1/2 chain with L=6 (with translational symmetry and conserve magnetic quantum number). If we consider the subspace spaned by those states with 3 up-spins with zero momentum, we can create the basis for the subspace by

```julia
translationalbasis(x -> sum(x)==3, 0, 6; base=2)
```

## Spin tools

We also provide a helper function to create spin-s operators (represented by matrices):

```julia
function spin(spins...; D::Integer=2)
```

In the definition, `spins` are arbituary number of tuples such as `(1.0, "xzx")`. The supported characters are:

```julia
"x", "y", "z", "1", "+", "-"
```

The other input `D` is the dimension of the matrix (`D = 2s+1`).
