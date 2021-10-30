"""
Quantum Inverse Method

Calculate Hamiltonian from eigen state(s).
"""


inner_product(v1::AbstractVector{<:Number}, v2::AbstractVector{<:Number}) = dot(v1, v2)
inner_product(v1::AbstractMatrix{<:Number}, v2::AbstractMatrix{<:Number}) = sum(dot(v1[:, i], v2[:, i]) for i in 1:size(v1, 2)) / size(v1, 2)

export covmat
function covmat(ol::AbstractVector{T}, v::AbstractVecOrMat{<:Number}) where T <: Union{<:AbstractMatrix, <:Operator}
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
