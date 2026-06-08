#-------------------------------------------------------------------------------------------------------------------------
"""
    SchmidtMatrix{Tm, Ta<:SubArray, Tb<:SubArray, Ti<:Integer}

Struct that is used to construct the Schmidt matrix:

`|Ψ⟩ = ∑ᵢⱼ Mᵢᵢ|ψᵢ⟩ ⊗ |ψⱼ⟩`

Fields:
- `M`: Schmidt matrix being accumulated.
- `A`: view onto subsystem-`A` digits inside the working basis buffer.
- `B`: view onto subsystem-`B` digits inside the working basis buffer.
- `B1`: basis used for subsystem `A`.
- `B2`: basis used for subsystem `B`.
"""
struct SchmidtMatrix{Tm <: Number, Ta <: SubArray, Tb <: SubArray, TB1 <: AbstractBasis, TB2 <: AbstractBasis, Td1 <: AbstractVector, Td2 <: AbstractVector}
    M::Matrix{Tm}
    A::Ta
    B::Tb
    B1::TB1
    B2::TB2
    dgt1::Td1
    dgt2::Td2
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    schmidtmatrix(T, b::AbstractBasis, Ainds::AbstractVector)

Construct the mutable bookkeeping object used to assemble a Schmidt matrix.

Arguments:
- `T`: element type of the matrix entries.
- `b`: full-system basis whose digit buffer will be partitioned.
- `Ainds`: site indices belonging to subsystem `A`.
- `B1`, `B2`: optional subsystem bases; when omitted, tensor-product bases are
  used.

Returns:
- A [`SchmidtMatrix`](@ref) whose matrix `M` is initialized to zeros.
"""
# Default subsystem basis for the Schmidt bipartition: a uniform `TensorBasis`
# for scalar-base bases, a `MixedTensorBasis` carrying the per-site dimensions of
# the selected sites for a mixed-dimension basis.
_schmidt_subbasis(b::AbstractBasis, sites) = TensorBasis(L=length(sites), base=b.B)
_schmidt_subbasis(b::MixedTensorBasis, sites) = MixedTensorBasis(dims=b.B[sites])

function schmidtmatrix(
    T::DataType, b::AbstractBasis, Ainds::AbstractVector{Ta},
    B1=nothing, B2=nothing;
    dgt::AbstractVector=b.dgt
) where Ta <: Integer
    L = length(b)
    Binds = Vector{Ta}(undef, L-length(Ainds))
    P = 1
    for i in range(one(Ta), stop=convert(Ta, L))
        if !in(i, Ainds)
            Binds[P] = i
            P += 1
        end
    end
    B1 = isnothing(B1) ? _schmidt_subbasis(b, Ainds) : B1
    B2 = isnothing(B2) ? _schmidt_subbasis(b, Binds) : B2
    M = zeros(T, size(B1, 1), size(B2, 1))
    dgt1 = similar(B1.dgt)
    dgt2 = similar(B2.dgt)
    SchmidtMatrix(M, view(dgt, Ainds), view(dgt, Binds), B1, B2, dgt1, dgt2)
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    addto!(S::SchmidtMatrix, val)

Accumulate a contribution `val` into the Schmidt matrix entry selected by the
current subsystem digits stored in `S.A` and `S.B`.
"""
function addto!(S::SchmidtMatrix, val::Number)
    S.dgt1 .= S.A
    S.dgt2 .= S.B
    _, ia = index(S.B1, S.dgt1)
    _, ib = index(S.B2, S.dgt2)
    S.M[ia, ib] += val
end

#-------------------------------------------------------------------------------------------------------------------------
# Entanglement Entropy
#-------------------------------------------------------------------------------------------------------------------------
"""
ent_spec(v::AbstractVector, Aind::AbstractVector{<:Integer}, b::AbstractBasis)

Compute the Schmidt singular values of a state across a bipartition.

Returns:
- The singular values of the Schmidt matrix produced by [`schmidt`](@ref).
"""
ent_spec(v::AbstractVector, Aind::AbstractVector{<:Integer}, b::AbstractBasis) = svdvals(schmidt(v, Aind, b))
#-------------------------------------------------------------------------------------------------------------------------
"""
entropy(s::AbstractVector{<:Real}; α::Real=1, cutoff::Real=1e-20)

Compute the entropy of Schmidt values.  

Inputs:
-------
- `s`     : Schmidt values.
- `α`     : Renyi index.
- `cutoff`: Cutoff of the Schmidt values. 

Outputs:
--------
- `S`: Entanglement entropy.
"""
function entropy(s::AbstractVector{<:Real}; α::Real=1, cutoff::Real=1e-20)
    if isone(α)
        shannon_entropy(s, cutoff=cutoff)
    elseif iszero(α)
        renyi_zero_entropy(s, cutoff=cutoff)
    else
        renyi_entropy(s, α, cutoff=cutoff)
    end
end
#-------------------------------------------------------------------------------------------------------------------------
"""
Compute the Shannon/von Neumann entropy of a probability vector.
"""
function shannon_entropy(s::AbstractVector{<:Real}; cutoff::Real=1e-20)
    ent = 0.0
    for si in s
        si > cutoff || continue
        ent -= si * log(si)
    end
    ent
end
#-------------------------------------------------------------------------------------------------------------------------
"""
Compute Renyi-0 entropy, that is, the log-support size.
"""
function renyi_zero_entropy(s::AbstractVector{<:Real}; cutoff::Real=1e-20)
    N = 0
    for si in s
        si > cutoff || continue
        N += 1
    end
    iszero(N) ? 0.0 : log(N)
end
#-------------------------------------------------------------------------------------------------------------------------
"""
Compute the Renyi entropy of order `α` for a probability vector `s`.
"""
function renyi_entropy(s::AbstractVector{<:Real}, α::Real; cutoff::Real=1e-20)
    Z = 0.0
    for si in s
        si > cutoff && (Z += si)
    end
    iszero(Z) && return 0.0
    acc = 0.0
    for si in s
        si > cutoff && (acc += (si / Z)^α)
    end
    log(acc) / (1 - α)
end
#-------------------------------------------------------------------------------------------------------------------------
export ent_S
"""
ent_S(v::AbstractVector, Aind::AbstractVector{<:Integer}, b::AbstractBasis; α::Real=1, cutoff::Real=1e-20)

Compute the bipartite entanglement entropy of a state represented in basis `b`.

The singular values returned by [`ent_spec`](@ref) are squared into Schmidt
probabilities before being passed to [`entropy`](@ref).
"""
function ent_S(v::AbstractVector, Aind::AbstractVector{<:Integer}, b::AbstractBasis; α::Real=1, cutoff::Real=1e-20)
    s = ent_spec(v, Aind, b) .^ 2
    entropy(s, α=α, cutoff=cutoff)
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    ent_S(v, Aind, L; α=1, cutoff=1e-20)

Convenience overload that infers a `TensorBasis` from the vector length and
system size `L`.
"""
function ent_S(v::AbstractVector, Aind::AbstractVector{<:Integer}, L::Integer; α::Real=1, cutoff::Real=1e-20)
    b = TensorBasis(L=L, base=round(Int, length(v)^(1/L)))
    s = ent_spec(v, Aind, b) .^ 2
    entropy(s, α=α, cutoff=cutoff)
end
#-------------------------------------------------------------------------------------------------------------------------
export rdm
"""
    rdm(v, Ainds, b::AbstractBasis; B1=nothing, B2=nothing)
    rdm(v, Ainds, L::Integer; B1=nothing, B2=nothing)

Reduced density matrix of subsystem `A` for the state `v` represented in basis
`b`.

Writing the Schmidt decomposition `|v⟩ = Sᵢⱼ |Aᵢ⟩|Bⱼ⟩` (see [`schmidt`](@ref)),
the subsystem-`A` reduced density matrix is `ρ_A = S S†`. It is Hermitian,
positive semidefinite, and has unit trace when `v` is normalized.

Arguments:
- `v`     : state vector in basis `b`.
- `Ainds` : site indices of subsystem `A`; the rest form subsystem `B`.
- `b`/`L` : the basis, or a system size `L` from which a [`TensorBasis`](@ref) is inferred.
- `B1`,`B2` : optional subsystem bases forwarded to [`schmidt`](@ref).

The matrix is expressed in subsystem `A`'s tensor-product basis (`Float64` for
real-phase bases, `ComplexF64` for momentum sectors), so for a symmetry-reduced
`b` the result is the physical RDM in the full local Hilbert space of `A`.

See also [`ent_S`](@ref), [`mutual_information`](@ref).
"""
function rdm(v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::AbstractBasis; B1=nothing, B2=nothing)
    S = schmidt(v, Ainds, b; B1, B2)
    S * S'
end
function rdm(v::AbstractVector, Ainds::AbstractVector{<:Integer}, L::Integer; B1=nothing, B2=nothing)
    b = TensorBasis(L=L, base=round(Int, length(v)^(1/L)))
    rdm(v, Ainds, b; B1, B2)
end
#-------------------------------------------------------------------------------------------------------------------------
export mutual_information
"""
    mutual_information(v, Ainds, Cinds, b::AbstractBasis; α=1, cutoff=1e-20)
    mutual_information(v, Ainds, Cinds, L::Integer; α=1, cutoff=1e-20)

Quantum mutual information

    I(A:C) = S(A) + S(C) − S(A∪C)

between two **disjoint** subsystems `A` and `C` (which need not be
complementary), where `S(·)` is the entanglement entropy [`ent_S`](@ref) of the
indicated region. `α` selects the Rényi index and `cutoff` the Schmidt-value
floor, both forwarded to `ent_S`.

For a pure state and `C = Ā` (the complement of `A`), this reduces to
`I(A:Ā) = 2 S(A)`. An `ArgumentError` is thrown if `A` and `C` overlap.

See also [`rdm`](@ref), [`ent_S`](@ref).
"""
function mutual_information(
    v::AbstractVector, Ainds::AbstractVector{<:Integer}, Cinds::AbstractVector{<:Integer},
    b::AbstractBasis; α::Real=1, cutoff::Real=1e-20
)
    A = sort(collect(Ainds))
    C = sort(collect(Cinds))
    overlap = intersect(A, C)
    isempty(overlap) || throw(ArgumentError("mutual_information requires disjoint subsystems; sites $overlap appear in both."))
    AC = sort(vcat(A, C))
    ent_S(v, A, b; α, cutoff) + ent_S(v, C, b; α, cutoff) - ent_S(v, AC, b; α, cutoff)
end
function mutual_information(
    v::AbstractVector, Ainds::AbstractVector{<:Integer}, Cinds::AbstractVector{<:Integer},
    L::Integer; α::Real=1, cutoff::Real=1e-20
)
    b = TensorBasis(L=L, base=round(Int, length(v)^(1/L)))
    mutual_information(v, Ainds, Cinds, b; α, cutoff)
end

#-------------------------------------------------------------------------------------------------------------------------
# Specific Bases
#-------------------------------------------------------------------------------------------------------------------------
"""
schmidt(v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::AbstractOnsiteBasis)

Schmidt decomposition of state `v`, with respect to given lattice bipartition.

Inputs:
-------
- `v`    : State represented by a (abstract) vector. 
- `Ainds`: List of indices in subsystem `A`, the remaining indices are regarded as subsystem `B`.
- `b`    : Basis.

Outputs:
--------
- `S`: Matrix S in the decomposition: |v⟩ = Sᵢⱼ |Aᵢ⟩|Bⱼ⟩.
"""
function _generic_onsite_schmidt(v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::AbstractOnsiteBasis, B1, B2)
    dgt = similar(b.dgt)
    S = schmidtmatrix(eltype(v), b, Ainds, B1, B2; dgt)
    for i = 1:length(v)
        change!(b, i, dgt)
        addto!(S, v[i])
    end
    S.M
end

function schmidt(v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::AbstractOnsiteBasis; B1=nothing, B2=nothing)
    _generic_onsite_schmidt(v, Ainds, b, B1, B2)
end
#-------------------------------------------------------------------------------------------------------------------------
"""
Schmidt decomposition specialized to the full [`TensorBasis`](@ref).

When the subsystem bases are the default tensor-product bases, this reshapes the
state vector into the full product tensor and permutes axes into `(A, B)` order
instead of iterating over product states. Custom `B1` or `B2` arguments keep the
generic onsite assembly path.
"""
function schmidt(v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::TensorBasis; B1=nothing, B2=nothing)
    if !isnothing(B1) || !isnothing(B2)
        return _generic_onsite_schmidt(v, Ainds, b, B1, B2)
    end

    L = length(b)
    base = b.B
    if _is_contiguous_site_block(Ainds)
        nA = length(Ainds)
        firstA = isempty(Ainds) ? 1 : first(Ainds)
        lastA = isempty(Ainds) ? 0 : last(Ainds)
        nright = L - lastA
        nleft = firstA - 1
        iszero(nleft) && return permutedims(reshape(v, base^nright, base^nA))
        iszero(nright) && return reshape(copy(v), base^nA, base^nleft)
        tensor = reshape(v, base^nright, base^nA, base^nleft)
        return reshape(permutedims(tensor, (2, 1, 3)), base^nA, base^(L - nA))
    end

    Binds = [i for i in 1:L if !in(i, Ainds)]
    site_to_dim(site) = L - site + 1
    perm = vcat(site_to_dim.(reverse(Ainds)), site_to_dim.(reverse(Binds)))
    tensor = reshape(v, ntuple(_ -> base, L))
    reshape(permutedims(tensor, perm), base^length(Ainds), base^length(Binds))
end

function _is_contiguous_site_block(Ainds::AbstractVector{<:Integer})
    isempty(Ainds) && return true
    start = first(Ainds)
    for (offset, site) in enumerate(Ainds)
        site == start + offset - 1 || return false
    end
    true
end
#-------------------------------------------------------------------------------------------------------------------------
"""
Schmidt decomposition specialized to [`TranslationalBasis`](@ref).

Each reduced-basis coefficient is expanded across the full translation orbit
with the appropriate momentum phase before contributing to the bipartite matrix.
"""
function schmidt(v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::TranslationalBasis; B1=nothing, B2=nothing)
    dgt = similar(b.dgt)
    R, phase = b.R, b.C[2]
    S = schmidtmatrix(promote_type(eltype(v), eltype(b)), b, Ainds, B1, B2; dgt)
    for i = 1:length(v)
        change!(b, i, dgt)
        val = v[i] / R[i]
        for j in 1:length(dgt)÷b.A
            addto!(S, val)
            circshift!(dgt, b.A)
            val *= phase
        end
    end
    S.M
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    spinflip(v::AbstractVector{<:Integer}, base::Integer)

Flip spins Sz on each site.
"""
function spinflip(v::AbstractVector{<:Integer}, base::Integer)
    vf = Vector{eltype(v)}(undef, length(v))
    base -= 1
    for i = 1:length(vf)
        vf[i] = base - v[i]
    end
    vf
end

"""
    spinflip!(v::AbstractVector{<:Integer}, base::Integer)

In-place version of [`spinflip`](@ref).
"""
function spinflip!(v::AbstractVector{<:Integer}, base::Integer)
    base -= 1
    for i in eachindex(v)
        v[i] = base - v[i]
    end
    v
end
#-------------------------------------------------------------------------------------------------------------------------
"""
Internal helper for Schmidt decomposition in bases that combine translation with
an involutive discrete symmetry such as parity or spin flip.
"""
function parity_schmidt(parity!, v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::AbstractTranslationalParityBasis;B1=nothing, B2=nothing)
    dgt = similar(b.dgt)
    R, phase = b.R, b.C[2]
    S = schmidtmatrix(promote_type(eltype(v), eltype(b)), b, Ainds, B1, B2; dgt)
    for i = 1:length(v)
        change!(b, i, dgt)
        val = v[i] / R[i]
        for j in 1:length(dgt)÷b.A
            addto!(S, val)
            circshift!(dgt, b.A)
            val *= phase
        end
        parity!(dgt)
        val *= b.P
        for j in 1:length(dgt)÷b.A
            addto!(S, val)
            circshift!(dgt, b.A)
            val *= phase
        end
    end
    S.M
end
#-------------------------------------------------------------------------------------------------------------------------
schmidt(v, Ainds, b::TranslationParityBasis;B1=nothing, B2=nothing) = parity_schmidt(reverse!, v, Ainds, b; B1, B2)
schmidt(v, Ainds, b::TranslationFlipBasis;B1=nothing, B2=nothing) = parity_schmidt(x -> spinflip!(x, b.B), v, Ainds, b; B1, B2)

#-------------------------------------------------------------------------------------------------------------------------
"""
Schmidt decomposition specialized to [`FlipBasis`](@ref).
"""
function schmidt(v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::FlipBasis; B1=nothing, B2=nothing)
    dgt = similar(b.dgt)
    R, phase = b.R, b.P
    S = schmidtmatrix(promote_type(eltype(v), eltype(b)), b, Ainds, B1, B2; dgt)
    for i = 1:length(v)
        change!(b, i, dgt)
        val = v[i] / R[i]
        addto!(S, val)
        spinflip!(dgt, b.B)
        addto!(S, phase * val)
    end
    S.M
end
#-------------------------------------------------------------------------------------------------------------------------
"""
Schmidt decomposition specialized to [`ParityBasis`](@ref).
"""
function schmidt(v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::ParityBasis; B1=nothing, B2=nothing)
    dgt = similar(b.dgt)
    R, phase = b.R, b.P
    S = schmidtmatrix(promote_type(eltype(v), eltype(b)), b, Ainds, B1, B2; dgt)
    for i = 1:length(v)
        change!(b, i, dgt)
        val = v[i] / R[i]
        addto!(S, val)
        reverse!(dgt)
        addto!(S, phase * val)
    end
    S.M
end
#-------------------------------------------------------------------------------------------------------------------------
"""
Schmidt decomposition specialized to [`ParityFlipBasis`](@ref).
"""
function schmidt(v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::ParityFlipBasis; B1=nothing, B2=nothing)
    dgt = similar(b.dgt)
    R, p1, p2 = b.R, b.P, b.Z
    S = schmidtmatrix(promote_type(eltype(v), eltype(b)), b, Ainds, B1, B2; dgt)
    for i = 1:length(v)
        # (P,Z) = (0,0)
        change!(b, i, dgt)
        val = v[i] / R[i]
        addto!(S, val)
        # (P,Z) = (1,0)
        reverse!(dgt)
        val *= p1
        addto!(S, val)
        # (P,Z) = (1,1)
        spinflip!(dgt, b.B)
        val *= p2
        addto!(S, val)
        # (P,Z) = (0,1)
        reverse!(dgt)
        val *= p1
        addto!(S, val)
    end
    S.M
end
