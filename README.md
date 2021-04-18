# EDKit.jl

 Julia package for general many-body exact diagonalization calculation. The package provide a general Hamiltonian constructing routine for specific symmetry sectors. The functionalities can be extended providing user-defined bases.

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

The basis object can be expanded, by defining 3 functions 

- `size(b::AbstractBasis)`;
- `change!(b::AbstractBasis, i::Integer)`;
- `index(b::AbstractBasis)`.

Optionally, we can define `eltype` for a basis object (default is `ComplexF64`).

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

## Spin tools

We also provide a helper function to create spin-s operators (represented by matrices):

```julia
function spin(spins...; D::Integer=2)
```

In the definition, `spins` are arbituary number of tuples such as `(1.0, "xzx")`. The supported characters are

```julia
"x", "y", "z", "1", "+", "-", "Y",
```

where `Y=iSʸ`. The other input `D` is the dimension of the matrix (`D = 2s+1`).

## Concrete Implementations of Basis

Here we introduce 3 concrete implementation of `AbstractBasis`.

### TensorBasis

The type `TensorBasis`

```julia
struct TensorBasis <: AbstractBasis
    dgt::Vector{Int}
    B::Int
end
```

### ProjectedBasis

The type `ProjectedBasis`

```julia
struct ProjectedBasis <: AbstractBasis
    dgt::Vector{Int}
    I::Vector{Int}
    B::Int
end
```

### TranslationalBasis

The type `TranslationalBasis`

```julia
struct TranslationalBasis <: AbstractBasis
    dgt::Vector{Int}
    I::Vector{Int}
    R::Vector{Float64}
    C::ComplexF64
    B::Int
end
```



