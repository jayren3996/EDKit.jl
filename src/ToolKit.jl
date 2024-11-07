#-------------------------------------------------------------------------------------------------------------------------
# Level statistics
#-------------------------------------------------------------------------------------------------------------------------
export gapratio, meangapratio
function gapratio(E::AbstractVector{<:Real})
    dE = diff(E)
    r = zeros(length(dE)-1)
    for i = 1:length(r)
        if dE[i] < dE[i+1]
            r[i] = dE[i]/dE[i+1]
        elseif dE[i] > dE[i+1]
            r[i] = dE[i+1]/dE[i]
        else
            r[i] = 1.0
        end
    end
    r
end

meangapratio(E::AbstractVector{<:Real}) = sum(gapratio(E)) / (length(E) - 2)

#-------------------------------------------------------------------------------------------------------------------------
# LinearAlgebra
#-------------------------------------------------------------------------------------------------------------------------
"""
    expm(A, order::Integer=10)

Matrix exponential using Taylor expansion.
"""
function expm(A; order::Integer=10)
    mat = I + A / order
    order -= 1
    while order > 0
        mat = A * mat
        mat ./= order
        mat += I
        order -= 1
    end
    mat
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    expv(A, v::AbstractVecOrMat; order=10, λ=1)

Compute exp(λA)*v using Taylor expansion
"""
function expv(A, v::AbstractVecOrMat; order::Integer=10, λ::Number=1)
    vec = v + λ * A * v / order
    order -= 1
    while order > 0
        vec = λ * A * vec
        vec ./= order
        vec += v
        order -= 1
    end
    vec
end

#-----------------------------------------------------------------------------------------------------
# State
#-----------------------------------------------------------------------------------------------------
export productstate
"""
    productstate(v::AbstractVector{<:Integer}, B::AbstractBasis)

Construction for product state.

Inputs:
-------
- `v`: Vector representing the product state.
- `B`: Basis.

Outputs:
- `s`: Vector representing the many-body state.
"""
function productstate(v::AbstractVector{<:Integer}, B::AbstractBasis)
    s = zeros(size(B, 1))
    B.dgt .= v 
    I = index(B)[2]
    s[I] = 1 
    s
end