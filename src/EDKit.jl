module EDKit

using LinearAlgebra
using SparseArrays
#-----------------------------------------------------------------------------------------------------
import Base: +, -, *, /
import Base: eltype, view, length, size, append!
import Base: Array, Vector
import LinearAlgebra: norm, mul!
import SparseArrays: sparse
#-----------------------------------------------------------------------------------------------------
export change!, index, operator

#-----------------------------------------------------------------------------------------------------
# INCLUDE
#-----------------------------------------------------------------------------------------------------
include("basis/BasisCore.jl")
include("basis/ToolKit.jl")
include("basis/TensorBasis.jl")
include("basis/ProjectedBasis.jl")
include("basis/TranslationalBasis.jl")
include("basis/TranslationParityBasis.jl")

include("basis/DoubleBasis.jl")
include("Operator.jl")
include("LinearOperation.jl")
include("Entanglement.jl")
include("Miscellaneous.jl")

end # module
