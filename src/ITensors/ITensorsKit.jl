#=---------------------------------------------------------------------------------------------------
Conversion
Due to julia's column major convention, the conversion between tensor and vector is not straight-
forward. It is basically a reshape and permute indices.
---------------------------------------------------------------------------------------------------=#
"""
vec2tensor!(dest, v)

Convert vector `v` to tensor, with row-major convention.

Inputs:
-------
dest: Tensor to write
v   : Vector
"""
function vec2tensor!(dest::Array, v::AbstractVector)
    T = reshape(v, size(dest))
    permutedims!(dest, T, ndims(dest):-1:1)
end
#----------------------------------------------------------------------------------------------------
"""
Convert vector `v` to tensor, with row-major convention.
"""
function vec2tensor(v::AbstractVector; base::Integer=2)
    L = round(Int, log(base, length(v)))
    @assert base^L == length(v) "Dimension not right."
    T = Array{eltype(v), L}(undef, fill(base, L)...)
    vec2tensor!(T, v)
end
#----------------------------------------------------------------------------------------------------
export mps2vec
"""
mps2vec(psi::MPS)

Convert MPS to a vector
"""
function mps2vec(psi::MPS)
    v = ITensor(1.0)
    for t in psi 
        v *= t 
    end
    T = Array(v, reverse(siteinds(psi)))
    reshape(T, :)
end
#----------------------------------------------------------------------------------------------------
"""
mps2vec(psi::MPS, B::AbstractBasis)

Convert MPS to a vector, with a given basis.
"""
function mps2vec(psi::MPS, B::AbstractBasis)
    L = length(psi)
    s = siteinds(psi)
    v = Vector{eltype(psi[1])}(undef, size(B, 1))
    for i in eachindex(v)
        R = change!(B, i)
        val = ITensor(1.0)
        for j in 1:L
            val *= psi[j] * state(s[j], B.dgt[j]+1)
        end
        v[i] = scalar(val) * L / R
    end
    v
end
#----------------------------------------------------------------------------------------------------
export vec2mps
function vec2mps(v::AbstractVector, s::AbstractVector{<:Index})
    T = vec2tensor(v; base=space(s[1]))
    MPS(T, s)
end
#----------------------------------------------------------------------------------------------------
export mat2op
"""
Convert a matrix to an ITensor operator
"""
function mat2op(mat::AbstractMatrix, s::Index...)
    L = length(s)
    d = space(s[1])
    T = reshape(mat, fill(d, 2L)...)
    Tp = permutedims(T, [L:-1:1; 2L:-1:L+1])
    op(Tp, s...)
end
#----------------------------------------------------------------------------------------------------
export op2mat
"""
Convert a ITensor operator to a matrix
"""
function op2mat(o::ITensor, s::Index...)
    rs = reverse(s)
    T = Array(o, adjoint.(rs), rs) 
    N = space(s[1]) ^ length(s)
    reshape(T, N, N)
end


#=---------------------------------------------------------------------------------------------------
MPS constructions
---------------------------------------------------------------------------------------------------=#
export pbcmps
"""
    pbcmps(sites, tensors)

Construct an MPS from a list of tensor. 
"""
function pbcmps(
    sites::AbstractVector,
    tensors::AbstractVector{<:AbstractArray}
)
    L = length(sites)
    link = map(1:L-1) do i 
        χ = size(tensors[i], 2)
        Index(χ, "Link,l=$i")
    end
    ltensor, rtensor = tensors[1], tensors[end]
    χ = size(ltensor, 1)

    ψ = MPS(L)
    for i in 2:L-1 
        ψ[i] = ITensor(tensors[i], link[i-1], link[i], sites[i])
    end

    # Fix boundary |1⟩⟨1|
    ψ[1] = ITensor(ltensor[1, :, :], link[1], sites[1])
    ψ[L] = ITensor(rtensor[:, 1, :], link[L-1], sites[L])

    res = deepcopy(ψ)
    for n in 2:χ    # Fix boundary |n⟩⟨n|
        ψ[1] = ITensor(ltensor[n, :, :], link[1], sites[1])
        ψ[L] = ITensor(rtensor[:, n, :], link[L-1], sites[L])
        res += ψ
    end
    
    orthogonalize!(res, 1)
    normalize!(res)
end
#----------------------------------------------------------------------------------------------------
export productstate
ITensors.state(s::Index, v::AbstractVector) = ITensor(v, s)
"""
product_state(s, states)

Return a product MPS 

Inputs:
-------
s     : Vector of indices
states: Vector of vector representing local states 
"""
function productstate(s::AbstractVector{<:Index}, states::AbstractVector{<:AbstractVector})
    L, d = length(s), length(states[1])
    link = [Index(1, "Link,l=$i") for i in 1:L-1]
    ψ = MPS(L)
    ψ[1] = ITensor(reshape(states[1], d, 1), s[1], link[1])
    ψ[L] = ITensor(reshape(states[L], d, 1), s[L], link[L-1])
    for i in 2:L-1 
        ψ[i] = ITensor(reshape(states[i], d, 1, 1), s[i], link[i-1], link[i])
    end
    ψ
end

#=---------------------------------------------------------------------------------------------------
MPS properties
---------------------------------------------------------------------------------------------------=#
export ent_spec, ent_specs!, ent_S!
function ent_specs(ψ::MPS, b::Integer)
    ψ = orthogonalize(ψ, b)
    svd(ψ[b], (linkind(ψ, b-1), siteind(ψ, b))).spec.eigs
end
#----------------------------------------------------------------------------------------------------
function ent_specs!(ψ::MPS, b::Integer)
    orthogonalize!(ψ, b)
    svd(ψ[b], (linkind(ψ, b-1), siteind(ψ, b))).spec.eigs
end
#----------------------------------------------------------------------------------------------------
function ent_S(ψ::MPS, b::Integer)
    spec = ent_specs(ψ, b)
    ITensors.entropy(spec)
end
#----------------------------------------------------------------------------------------------------
function ent_S!(ψ::MPS, b::Integer)
    spec = ent_specs!(ψ, b)
    ITensors.entropy(spec)
end
#----------------------------------------------------------------------------------------------------
function ent_S(ψ::MPS)
    L = length(ψ)
    S = Vector{Float64}(undef, L)
    ψ = orthogonalize(ψ, 1)
    for i in 1:L 
        S[i] = ent_S!(ψ, i)
    end
    S
end
