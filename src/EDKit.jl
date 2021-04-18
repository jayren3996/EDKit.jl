module EDKit

using LinearAlgebra
using SparseArrays
#-----------------------------------------------------------------------------------------------------
import Base: +, -, *, /
import Base: eltype, view, length, size
import Base: Array
import LinearAlgebra: norm, mul!
import SparseArrays: sparse
#-----------------------------------------------------------------------------------------------------
export change!, index, operator
#-----------------------------------------------------------------------------------------------------
# Abstract Basis Type
#-----------------------------------------------------------------------------------------------------
"""
    AbstractBasis

Abstract type for basis object. Any concrete implementation should have the following functions:

Functions:
----------
- size(basis)       : Size of the block.
- change!(basis, i) : Change the bit and return normalization.
- index(basis)      : Get normalization and index of new state.

Optional:
---------
- eltype(basis)     : Default is ComplexF64
"""
abstract type AbstractBasis end
eltype(::AbstractBasis) = ComplexF64
length(b::AbstractBasis) = size(b, 2)
#-----------------------------------------------------------------------------------------------------
# INCLUDE
#-----------------------------------------------------------------------------------------------------
include("basis/TensorBasis.jl")
include("basis/ProjectedBasis.jl")
include("basis/TranslationalBasis.jl")
include("Operator.jl")
include("LinearOperation.jl")
include("Spin.jl")

end # module
