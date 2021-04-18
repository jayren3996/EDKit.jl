module OldED

using LinearAlgebra
using SparseArrays
import Base: view, copy, eltype, size, Array
import Base: +, -, *, /, sum
import LinearAlgebra: mul!

include("Basis.jl")
include("Operation.jl")
include("Operation_construct.jl")
include("Measure.jl")
include("Spin.jl")
#-----------------------------------------------------------------------------------------------------
# Krylov space
#-----------------------------------------------------------------------------------------------------
function reducespace(
    vector_space::AbstractMatrix; 
    rtol::Real=1e-3
)
    mat = vector_space' * vector_space
    vals, vecs = eigen(Hermitian(mat))
    pos = vals .> rtol
    vals_pos = vals[pos]
    vecs_pos = vecs[:, pos]
    null_vector_space = vector_space * vecs_pos
    for i = 1:length(vals_pos)
        null_vector_space[:, i] ./= sqrt(vals_pos[i])
    end
    return null_vector_space
end
# Special case for single vector
reducespace(v::AbstractVector) = normalize(v)
#-----------------------------------------------------------------------------------------------------
export krylovspace
function krylovspace(
    H, 
    vector_space::AbstractMatrix; 
    rtol::Real=1e-3
)
    vs_rank = rank(vector_space)
    cat_space = hcat(vector_space, H * vector_space)
    new_space = reducespace(cat_space, rtol=rtol)
    new_rank = size(new_space, 2)
    while new_rank > vs_rank
        vs_rank = new_rank
        cat_space = hcat(new_space, H * new_space)
        new_space = reducespace(cat_space, rtol=rtol)
        new_rank = size(new_space, 2)
    end
    return new_space
end
# Special case for vector
krylovspace(H, v::AbstractVector; rtol::Real=1e-3) = krylovspace(H, reshape(v, :, 1), rtol=rtol)


end # module ExactDiagonalization
