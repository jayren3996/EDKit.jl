#=---------------------------------------------------------------------------------------------------
Conversion
Due to julia's column major convention, the conversion between tensor and vector is not straight-
forward. It is basically a reshape and permute indices.
---------------------------------------------------------------------------------------------------=#
"""
    vec2tensor!(dest, v)

Convert vector `v` to tensor, with row-major convention.

Arguments:
- `dest`: destination tensor to overwrite.
- `v`: flat input vector.

Returns:
- The mutated tensor `dest`.

Notes:
- EDKit interprets vectors in row-major physical-site order before converting to
  Julia's column-major array layout.
"""
function vec2tensor!(dest::Array, v::AbstractVector)
    T = reshape(v, size(dest))
    permutedims!(dest, T, ndims(dest):-1:1)
end
#----------------------------------------------------------------------------------------------------
"""
    vec2tensor(v; base=2)

Convert vector `v` to tensor, with row-major convention.

Arguments:
- `v`: flat vector of length `base^L`.
- `base`: local dimension used to infer `L`.

Returns:
- An `L`-index dense tensor whose entries correspond to `v` in EDKit's ordering
  convention.
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

Convert an MPS to a dense state vector.

Returns:
- A vector in EDKit's row-major site-order convention.
"""
function mps2vec(psi::MPS)
    v = ITensor(1.0)
    for t in psi 
        v *= t 
    end
    T = Array(v, reverse(siteinds(psi))...)
    reshape(T, :)
end
#----------------------------------------------------------------------------------------------------
"""
    mps2vec(psi::MPS, B)

Convert an MPS to coordinates in the symmetry-reduced basis `B`.

This routine is exact when `psi` already lies in the symmetry sector described
by `B`. In that case it reconstructs the sector amplitudes from the amplitude on
each representative product state together with the normalization/orbit data of
`B`.

It is not a generic orthogonal projector for arbitrary MPS outside the sector:
for a general state, use an explicit projection constructed from the basis if
you need the true projected coefficients.
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
        v[i] = scalar(val) * order(B) / R
    end
    v
end
#----------------------------------------------------------------------------------------------------
export vec2mps
"""
    vec2mps(v::AbstractVector, s)

Convert a state vector `v` into an MPS on the ITensor site indices `s`.

The vector is interpreted in EDKit's row-major site ordering convention, then
reshaped into a tensor and wrapped as an `MPS`.
"""
function vec2mps(v::AbstractVector, s::AbstractVector{<:Index})
    T = vec2tensor(v; base=space(s[1]))
    MPS(T, s)
end
#----------------------------------------------------------------------------------------------------
export mat2op
"""
Convert a matrix to an ITensor operator

Arguments:
- `mat`: dense local operator matrix.
- `s...`: ITensor site indices describing the target Hilbert spaces.

Returns:
- An `ITensor` operator with bra/ket site structure compatible with `op`.
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

Returns:
- The dense matrix representation of `o` in the site ordering supplied by
  `s...`.
"""
function op2mat(o::ITensor, s::Index...)
    rs = reverse(s)
    T = Array(o, adjoint.(rs)..., rs...)
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

Arguments:
- `sites`: physical site indices.
- `tensors`: local array data for a periodic-boundary-style construction.

Returns:
- An MPS normalized after summing over the boundary-sector contributions.
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
productstate(s, states)

Return a product MPS 

Arguments:
- s     : Vector of indices
- states: Vector of vector representing local states 

Returns:
- A product-state `MPS` with bond dimension `1`.
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
"""
    ent_specs(ψ::MPS, b::Integer)

Return entanglement spectrum between site `b` and `b+1`.

Arguments:
- ψ: MPS
- b: Link index

Returns:
- The singular values across the cut at bond `b`.
"""
function ent_specs(ψ::MPS, b::Integer)
    ψ = orthogonalize(ψ, b)
    isone(b) && return svd(ψ[b], siteind(ψ, b)).spec.eigs
    svd(ψ[b], (linkind(ψ, b-1), siteind(ψ, b))).spec.eigs
end
#----------------------------------------------------------------------------------------------------
"""
    ent_specs!(ψ::MPS, b::Integer)

In-place version of [`ent_specs`](@ref) that may change the orthogonality center
of `ψ`.
"""
function ent_specs!(ψ::MPS, b::Integer)
    orthogonalize!(ψ, b)
    isone(b) && return svd(ψ[b], siteind(ψ, b)).spec.eigs
    svd(ψ[b], (linkind(ψ, b-1), siteind(ψ, b))).spec.eigs
end
#----------------------------------------------------------------------------------------------------
"""
    ent_S(ψ::MPS, b::Integer)

Return entanglement entropy between site `b` and `b+1`.

Arguments:
- ψ: MPS
- b: Link index

Returns:
- The entropy associated with the cut at bond `b`.
"""
function ent_S(ψ::MPS, b::Integer)
    spec = ent_specs(ψ, b)
    entropy(spec)
end
#----------------------------------------------------------------------------------------------------
"""
    ent_S!(ψ::MPS, b::Integer)

In-place version of [`ent_S(ψ::MPS, b)`](@ref) that may move the orthogonality
center of `ψ`.
"""
function ent_S!(ψ::MPS, b::Integer)
    spec = ent_specs!(ψ, b)
    entropy(spec)
end
#----------------------------------------------------------------------------------------------------
"""
    ent_specs(ψ::MPS)

Return entanglement entropies along the full chain.

Arguments:
- ψ: MPS

Returns:
- S: Vector of entropy
"""
function ent_S(ψ::MPS)
    L = length(ψ)
    S = Vector{Float64}(undef, L)
    ψ = orthogonalize(ψ, 1)
    for i in 1:L 
        S[i] = ent_S!(ψ, i)
    end
    S
end
