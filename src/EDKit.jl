module EDKit

using LinearAlgebra, SparseArrays, JLD2

import Base: +, -, *, /, Array, Vector, Matrix, size, length, eltype, digits
import LinearAlgebra: norm, mul!
import SparseArrays: sparse

export content, base, index, change!

"""
    AbstractBasis

Abstract type for all concrete `Basis` types. 
A typical concrete `Basis` type usually has the field:

- `dgt`: Digits representation of a state. For a symmetric basis like translational-invariant basis, `dgt` is just a representation of the state.
- `I`  : List of indices in this basis.
- `R`  : List of normalization of each state in the basis.
- `B`  : Base of the basis.

Correspondingly, there are property functions for a `Basis` type:

- `digits` : Return digits representation.
- `content`: Return list of indices in this basis.
- `norm`   : Return list of norm.
- `base`   : Return the base of the basis.

In addition, we have the following property functions:

- `eltype` : Used for type promotion when constructing matrix representation of an `Operator`.
- `length` : Return the length of the digits of the basis.
- `size`   : Return the size of the matrix that can be constructed for an `Operator`.

The most inportant function for a `Basis` type is the index of a given digits and the digits representation of its ith content:

- `index`  : Return the norm and the index of the digits of the given basis.
- `change!`: Change the digits to the target index, and return the new normalization.

If we want to calculate the schmidt decomposition (in real space), we can use the `schmidt` function.
"""
abstract type AbstractBasis end

include("ToolKit.jl")
include("Search.jl")
include("Schmidt.jl")
include("Symmetries.jl")
include("Basis.jl")
include("Operator.jl")
include("LinearOperation.jl")
include("Miscellaneous.jl")
include("algorithms/BlockDiagonal.jl")

# Default function definitions:
@inline digits(b::AbstractBasis) = b.dgt
@inline content(b::AbstractBasis) = b.I
@inline content(b::AbstractBasis, i::Integer) = b.I[i]
@inline norm(b::AbstractBasis) = b.R
@inline norm(b::AbstractBasis, i::Integer) = b.R[i]
@inline base(b::AbstractBasis) = b.B
@inline eltype(::AbstractBasis) = ComplexF64
@inline length(b::AbstractBasis) = length(b.dgt)
@inline size(b::AbstractBasis, i::Integer) = isone(i) || isequal(i, 2) ? length(b.I) : 1
@inline size(b::AbstractBasis) = (size(b, 1), size(b, 2))
@inline change!(b::AbstractBasis, i::Integer) = (change!(b.dgt, content(b, i), base=b.B); norm(b, i))

end # module
