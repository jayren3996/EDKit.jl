module EDKit

using LinearAlgebra, SparseArrays, Random

import Base: +, -, *, /, Array, Vector, Matrix, size, length, eltype, digits
import LinearAlgebra: norm, mul!
import SparseArrays: sparse

include("Basis.jl")
include("LinearMap.jl")
include("Schmidt.jl")
include("Symmetries.jl")
include("Operator.jl")
include("ToolKit.jl")

for file in readdir("$(@__DIR__)/algorithms/")
    if file[end-2:end] == ".jl"
        include("$(@__DIR__)/algorithms/$file")
    end
end

end # module
