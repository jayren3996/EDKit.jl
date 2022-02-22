# Operators

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

