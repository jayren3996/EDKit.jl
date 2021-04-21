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
include("basis/TensorBasis.jl")
include("basis/ProjectedBasis.jl")
include("basis/TranslationalBasis.jl")
include("basis/ParityBasis.jl")
include("basis/DoubleBasis.jl")
include("Operator.jl")
include("LinearOperation.jl")
include("Spin.jl")

end # module
