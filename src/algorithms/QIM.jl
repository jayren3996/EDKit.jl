"""
Quantum Inverse Method

Calculate Hamiltonian from eigen state(s).
"""


@inline inner_product(v1::AbstractVector{<:Number}, v2::AbstractVector{<:Number}) = dot(v1, v2)
@inline inner_product(v1::AbstractMatrix{<:Number}, v2::AbstractMatrix{<:Number}) = dot(v1, v2) / size(v1, 2)

export covmat
"""
covmat(ol::AbstractVector, v::AbstractVecOrMat{<:Number})

Return the covariant matrix for a given list of operators.

Inputs:
-------
- `ol`: List of operators, represented by matrices or `Operator`s.
- `V` : Target state represented by a vector, or target states represented by a matrix.

Outputs:
--------
- cm: (Hermitian) matrix.
"""
function covmat(ol::AbstractVector, v::AbstractVecOrMat{<:Number})
    n = length(ol)
    vs = [oi * v for oi in ol]
    am = [real(inner_product(v, vsi)) for vsi in vs]
    cm = Matrix{Float64}(undef, n, n)
    for i=1:n
        for j=i:n
            cm[i, j] = real(inner_product(vs[i], vs[j])) - am[i] * am[j] 
        end
    end
    Hermitian(cm)
end
