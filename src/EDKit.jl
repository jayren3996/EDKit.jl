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

In `EDKit.jl`, the fundamental objects are basis and operator. The `AbstractBasis` is the abstract type of basis. Currently there are 3 concrete basis:

1. `TensorBasis`: Ordinary basis without any symmetry.
2. `ProjectedBasis`: Basis for subspace that is spanned only by product states.
3. `TranslationalBasis`: Basis for translational symmetric Hamiltonian.
    
The basis object can be extended. To construct linear operation, we need to define 4 functions for a new basis type:
    
1. `size(b::AbstractBasis)` : Size of matrix representations of the operators in this subspace.
2. `change!(b::AbstractBasis, i::Integer)` : Change the digits to ith states in this subspace.
3. `index(b::AbstractBasis)` : Return the coefficient and index of the digits.
4. `norm(b::AbstractBasis, i)` : Normalization of the given basis.
    
Optionally, we can define `eltype` for a basis object (default is `ComplexF64`).
    
If the calculation is done on the entire Hilbert space, the basis object need not be explicitly constructed. The `Operator` will use `TensorBasis` by default. The construction of other basis with symmetry concern are discussed below.
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
include("basis/DoubleBasis.jl")
include("Operator.jl")
include("LinearOperation.jl")
include("Spin.jl")

end # module
