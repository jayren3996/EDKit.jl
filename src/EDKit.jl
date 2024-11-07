module EDKit

using LinearAlgebra, SparseArrays, Random, Combinatorics
using ITensors, ITensorMPS, LRUCache

import Base: +, -, *, /, Array, Vector, Matrix, size, length, eltype, digits, copy
import LinearAlgebra: norm, mul!

include("Basis/AbstractBasis.jl")
include("Basis/ProjectedBasis.jl")
include("Basis/TranslationalBasis.jl")
include("Basis/TranslationalParityBasis.jl")
include("Basis/TranslationalFlipBasis.jl")
include("Basis/ParityBasis.jl")
include("Basis/FlipBasis.jl")
include("Basis/ParityFlipBasis.jl")
include("Basis/AbelianBasis.jl")


include("LinearMap.jl")
include("Schmidt.jl")
include("Operator.jl")
include("ToolKit.jl")

for file in readdir("$(@__DIR__)/algorithms/")
    if file[end-2:end] == ".jl"
        include("$(@__DIR__)/algorithms/$file")
    end
end

for file in readdir("$(@__DIR__)/ITensors/")
    if file[end-2:end] == ".jl"
        include("$(@__DIR__)/ITensors/$file")
    end
end

end # module
