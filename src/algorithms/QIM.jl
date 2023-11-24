#---------------------------------------------------------------------------------------------------
# Quantum Inverse Method
# Calculate Hamiltonian from eigen state(s).
#---------------------------------------------------------------------------------------------------
export covmat
"""
covmat(ol::AbstractVector, v::AbstractVecOrMat{<:Number})

Return the covariant matrix for a given list of operators:
    Cᵢⱼ = 1/2⟨v|{hᵢ,hⱼ}|v⟩ - ⟨v|hᵢ|v⟩⟨v|hⱼ|v⟩,
or for a set of vector,
    Cᵢⱼ = 1/2Tr[ρ{hᵢ,hⱼ}] - Tr[ρhᵢ]Tr[ρhⱼ],
where ρ = 1/N∑ₙ|vₙ⟩⟨vₙ|.

Inputs:
-------
- `ol`: List of operators, represented by matrices or `Operator`s.
- `V` : Target state represented by a vector, or target states represented by a matrix.

Outputs:
--------
- cm: (Symmetric real) matrix.
"""
function covmat(ol::AbstractVector, v::AbstractVecOrMat{<:Number})
    n = length(ol)
    vs = [oi * v for oi in ol]
    am = [real(inner_product(v, vsi)) for vsi in vs]
    cm = Matrix{Float64}(undef, n, n)
    for i=1:n, j=i:n
        cm[i, j] = real(inner_product(vs[i], vs[j])) - am[i] * am[j] 
    end
    Symmetric(cm)
end
#---------------------------------------------------------------------------------------------------
export qimsolve
function qimsolve(ol::AbstractVector, v::AbstractVecOrMat{<:Number}; tol::Real=1e-7)
    cm = covmat(ol, v)
    e, v = eigen(cm)
    i = findfirst(x -> x > tol, e) - 1
    v[:, 1:i] |> simplify
end

#---------------------------------------------------------------------------------------------------
# Helper
#---------------------------------------------------------------------------------------------------
"""
Inner product ⟨v₁|v₂⟩. For two set of vector.
"""
function inner_product(v1::AbstractVecOrMat, v2::AbstractVecOrMat)
    dot(v1, v2) / size(v1, 2)
end
#---------------------------------------------------------------------------------------------------
"""
Simplify the solutions.
"""
function simplify(h)
    n = size(h, 2)
    for i in 1:n
        j = findfirst(x -> abs(x) > 1e-7, h[:,i])
        for k in 1:n
            isequal(k, i) && continue
            h[:,k] .-= h[j,k] / h[j,i] * h[:,i]
        end
    end
    for i in 1:n
        j = findfirst(x -> abs(x) > 1e-7, h[:,i])
        h[:,i] ./= h[j,i]
    end
    for i in eachindex(h)
        abs(h[i]) < 1e-10 && (h[i] = 0.0)
    end
    h
end
