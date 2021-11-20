# EDKit.jl

Julia package for general many-body exact diagonalization calculation. The package provide a general Hamiltonian constructing routine for specific symmetry sectors. The functionalities can be extended with user-defined bases.

## Installation

Run the following script in the ```Pkg REPL``` environment:

```julia
pkg> add EDKit
```

## Examples

### Transverse Field Ising Model

Consider the transverse insing model with the parameters:

```julia
J, h, L = 1.0, 0.5, 10
```

The Hamiltonian operator can be constructed as:

```julia
julia> H = J * trans_inv_operator("XX", L) + h * trans_inv_operator("Z", L)
Operator of size (1024, 1024) with 20 terms.
```

Use the function `Array` to create the dense matrix:

```julia
julia> Array(H)
1024×1024 Matrix{Float64}:
 5.0  0.0  0.0  1.0  …   0.0   0.0   0.0
 0.0  4.0  1.0  0.0      0.0   0.0   0.0
 0.0  1.0  4.0  0.0      0.0   0.0   0.0
 1.0  0.0  0.0  3.0      0.0   0.0   0.0
 0.0  0.0  1.0  0.0      0.0   0.0   0.0
 0.0  0.0  0.0  1.0  …   0.0   0.0   0.0
 ⋮                   ⋱              
 0.0  0.0  0.0  0.0      0.0   0.0   0.0
 0.0  0.0  0.0  0.0      1.0   0.0   0.0
 0.0  0.0  0.0  0.0  …   0.0   0.0   1.0
 0.0  0.0  0.0  0.0     -4.0   1.0   0.0
 0.0  0.0  0.0  0.0      1.0  -4.0   0.0
 0.0  0.0  0.0  0.0      0.0   0.0  -5.0
```

Or use the function `sparse` to create the sparse matrix (requires the module `SparseArrays` being imported):

```julia
julia> sparse(H)
1024×1024 SparseMatrixCSC{Float64, Int64} with 11264 stored entries:
⣿⣿⣷⡘⠦⠀⠀⠳⣄⠘⢦⡀⠀⠀⠳⣄⠀⠀⠀
⣙⠻⣿⣿⣶⢦⡀⠀⠈⠃⠀⠙⢦⡀⠀⠈⠳⣄⠀
⠈⠃⠸⣟⢻⣶⣽⠆⢦⡰⣄⠀⠀⠙⢦⡀⠀⠈⠓
⢤⡀⠀⠈⠳⠟⣿⣿⣦⡅⠈⠳⣄⠀⠀⠙⢦⡀⠀
⣀⠙⠦⠀⢈⡳⠌⠿⠿⢇⣀⣀⢈⡳⠄⠀⣀⠙⠦
⠈⠳⣄⠀⠀⠙⢦⡀⠀⢸⣿⣿⣦⢙⡂⠀⠈⠳⣄
⠀⠀⠈⠳⣄⠀⠀⠙⢦⡰⣌⢛⣿⣿⡝⢦⡀⠀⠀
⠙⢦⡀⠀⠈⠳⣄⠀⠀⠁⠈⠈⠳⣍⣿⣿⣧⡘⠦
⠀⠀⠙⢦⡀⠀⠈⠳⣄⠘⢦⡀⠀⠈⣉⠻⣿⣿⣷
⠀⠀⠀⠀⠙⠀⠀⠀⠈⠃⠀⠙⠀⠀⠈⠃⠙⠛⠛
```

## Basis Object

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

### Compute entaglement entropy

After obtaining Hamiltonian in a symmetry sector. We can calculate the entaglement entropy of an eigenvector `v` (assume the system size is `2L`, and the entropy cut is at the middel of the chain) by

```julia
ent_S(v::AbstractVector, 1:L, b::AbstractBasis)
```

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
tensorbasis(L::Integer; base::Integer=2)
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

In the definition, `k` is the momentum sector we are interested, and `f` is a function acting on digits that tells whether a given digits is valid in this basis. Note that `k` is an integer, representing the momentum at `2πk/L`. For example, consider a `S=1/2` chain with `L=6` (with translational symmetry and conserve magnetic quantum number). If we consider the subspace spaned by those states with 3 up-spins with zero momentum, we can create the basis for the subspace by

```julia
translationalbasis(x -> sum(x)==3, 0, 6; base=2)
```
