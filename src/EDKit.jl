module EDKit

using LinearAlgebra, SparseArrays, JLD2
#-----------------------------------------------------------------------------------------------------
import Base: +, -, *, /
import Base: eltype, view, length, size, append!
import Base: Array, Vector
import LinearAlgebra: norm, mul!
import SparseArrays: sparse
#-----------------------------------------------------------------------------------------------------
export change!, index, operator
const BITTYPE = UInt64
#-----------------------------------------------------------------------------------------------------
# INCLUDE
#-----------------------------------------------------------------------------------------------------
include("basis/TensorBasis.jl")
include("Operator.jl")
include("LinearOperation.jl")
include("Entanglement.jl")
include("Miscellaneous.jl")
include("algorithms/BlockDiagonal.jl")
end # module
