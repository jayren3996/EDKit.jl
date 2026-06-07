#-------------------------------------------------------------------------------------------------------------------------
# Level statistics
#-------------------------------------------------------------------------------------------------------------------------
export gapratio, meangapratio
"""
    gapratio(E::AbstractVector{<:Real})

Compute adjacent-gap ratios from a sorted energy spectrum `E`.

For consecutive level spacings `δ_n = E[n+1] - E[n]`, this returns the vector
with entries `min(δ_n, δ_{n+1}) / max(δ_n, δ_{n+1})`.

The input should already be sorted in ascending order.
"""
function gapratio(E::AbstractVector{<:Real})
    dE = diff(E)
    length(dE) < 2 && return Float64[]
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

"""
    meangapratio(E::AbstractVector{<:Real})

Return the mean adjacent-gap ratio of a sorted spectrum `E`.

This is a convenience wrapper around [`gapratio`](@ref).
"""
function meangapratio(E::AbstractVector{<:Real})
    r = gapratio(E)
    isempty(r) ? NaN : sum(r) / length(r)
end

#-------------------------------------------------------------------------------------------------------------------------
# LinearAlgebra
#-------------------------------------------------------------------------------------------------------------------------
"""
    expm(A, order::Integer=10)

Approximate `exp(A)` using a truncated Taylor expansion of the given order.

This routine is lightweight and convenient for small matrices, but it is not
intended to replace the more numerically robust algorithms in dedicated linear
algebra packages.

Returns:
- An approximation to the matrix exponential of `A`.
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

Approximate `exp(λA) * v` using a truncated Taylor expansion.

This avoids explicitly constructing `exp(λA)` and is useful for quick tests or
small problems.

Returns:
- An approximation to the action of `exp(λA)` on `v`.
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

Construct the basis vector corresponding to the product configuration `v`.

`v` is a digit representation of local states in base `B.B`, for example
`[0, 1, 0, 1]` for a spin-1/2 chain. The returned vector lives in the Hilbert
space defined by `B`, so this works for tensor-product, projected, and
symmetry-reduced bases alike.
"""
function productstate(v::AbstractVector{<:Integer}, B::AbstractBasis)
    length(v) == length(B.dgt) ||
        error("Configuration length $(length(v)) does not match basis length $(length(B.dgt)).")
    dgt = collect(v)                       # local buffer: thread-safe, no shared-state mutation
    c, I = index(B, dgt)
    iszero(c) &&
        error("Product configuration $(collect(Int, v)) is not contained in the sector spanned by the basis.")
    T = eltype(B)
    T <: Integer && (T = Float64)          # state vectors should be floating-point; keeps onsite output Float64
    s = zeros(T, size(B, 1))
    s[I] = c                               # keep the orbit phase / normalization coefficient
    s
end
