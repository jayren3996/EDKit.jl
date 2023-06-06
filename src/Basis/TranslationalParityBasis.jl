#-------------------------------------------------------------------------------------------------------------------------
# Translation + Parity
#-------------------------------------------------------------------------------------------------------------------------
export TranslationParityBasis
"""
    TranslationParityBasis

Basis with translational and reflection symmetries.

Properties:
-----------
- dgt: Digits.
- I  : Representing states.
- R  : Normalization.
- C  : {±1}, momentum phase.
- P  : {±1}, parity.
- A  : Length of unit cell.
- B  : Base.
"""
struct TranslationParityBasis{Ti <: Integer} <: AbstractTranslationalParityBasis
    dgt::Vector{Int64}    # Digits
    I::Vector{Ti}         # Representing states
    R::Vector{Float64}    # Normalization
    C::Vector{Int}        # Momentum phase exp(1im * 0/π) = {±1}
    P::Int                # {±1}, parity
    A::Int                # Length of unit cell
    B::Int64              # Base
end
#-------------------------------------------------------------------------------------------------------------------------
eltype(::TranslationParityBasis) = Float64
copy(b::TranslationParityBasis) = TranslationParityBasis(deepcopy(b.dgt), b.I, b.R, b.C, b.P, b.B)


#-------------------------------------------------------------------------------------------------------------------------
# Functions selecting the basis
#-------------------------------------------------------------------------------------------------------------------------
struct TranslationParityJudge
    F                   # Projective selection
    parity              # Parity operation acting on digits
    K::Int64            # Momentum
    A::Int64            # Length of unit cell
    P::Int64            # Parity eigenvalue {±1}
    L::Int64            # Length of translation
    B::Int64            # Base
    C::Vector{Float64}  # Normalization coefficients
end
#-------------------------------------------------------------------------------------------------------------------------
"""
Double check the symmetry of the state under parity.

Outputs:
--------
- Qs: Whether there is smaller index.
- Qp: Thether the state is symmetric under parity.
- M : (If Qp,) minimal value such that Tᵐ⋅P|dgt⟩ = |dgt⟩.
"""
function translation_double_check!(
    dgt::AbstractVector{<:Integer}, R::Integer, I0::Integer, base::Integer; a::Integer=1
)
    # Test first index
    In = index(dgt, base=base, dtype=typeof(I0))
    In < I0 && return (false, false, 0)
    isequal(In, I0) && return (true, true, 0)

    # Test shifted indices
    for i=1:R-1
        cyclebits!(dgt, a)
        In = index(dgt, base=base, dtype=typeof(I0))
        In < I0 && return (false, false, 0)
        isequal(In, I0) && return (true, true, i)
    end

    # Find no parity symmetry
    true, false, 0
end
#-------------------------------------------------------------------------------------------------------------------------
"""
Check whether a state is a valid representing state, and whether it is reflect-translation-invariant.

Inputs:
-------
- `parity`: Parity function, can be spatio-reflection or spin-reflection.
- `dgt`   : Digits.
- `I0`    : Index of the `dgt`.
- `k`     : Momentum.
- `p`     : Parity.
- `base`  : Base.
- `a`     : Length of unit cell

Outputs:
--------
- `Qm`: Whether `dgt` is a representing vector.
- `Qp`: Whether `dgt` is reflect-translation-invariant.
- `R` : Minimum R that Tᴿ⋅|dgt⟩ = |dgt⟩.
"""
function translation_parity_check!(
    parity, dgt::AbstractVector{<:Integer}, I0::Integer, 
    k::Integer, p::Integer, base::Integer; a::Integer=1
)
    pdgt = parity(dgt)
    # First going through `translation_check` routine.
    Q, R = translation_check!(dgt, I0, k, base, a=a)
    Q || return (false, false, R)
    
    # Then going through `translation_double_check` routine.
    Qs, Qp, M = translation_double_check!(pdgt, R, I0, base, a=a)
    Qs || return (false, false, R)

    # Check parity 
    Qp && isequal(isone(p), isodd(div(k*M, length(dgt)÷a÷2))) && return (false, false, R)
    true, Qp, R
end
#-------------------------------------------------------------------------------------------------------------------------
function (judge::TranslationParityJudge)(dgt::AbstractVector{<:Integer}, i::Integer)
    # First check projective selection rule
    isnothing(judge.F) || judge.F(dgt) || return false, 0.0
    
    # Going through `translation_parity_check` routine:
    Qm, Qp, R = translation_parity_check!(judge.parity, dgt, i, judge.K, judge.P, judge.B, a=judge.A)
    Qm || return false, 0.0

    # Calculating norm
    N = Qp ? judge.C[judge.L+R] : judge.C[R]
    true, N
end

#-------------------------------------------------------------------------------------------------------------------------
# Construction
#-------------------------------------------------------------------------------------------------------------------------
"""
    TranslationParityBasis(f, k, p, L; base=2, alloc=1000, threaded=true)

Construction for `TranslationParityBasis`.

Inputs:
-------
- `f`       : Selection function for the basis contents.
- `k`       : Momentum number from 0 to L-1.
- `p`       : Parity number ±1.
- `L`       : Length of the system.
- `N`       : Particle numper.
- `a`       : Length of unit cell.
- `base`    : Base, default = 2.
- `alloc`   : Size of the prealloc memory for the basis content, used only in multithreading, default = 1000.
- `threaded`: Whether use the multithreading, default = true.

Outputs:
--------
- `b`: TranslationParityBasis.
"""
function TranslationParityBasis(
    dtype::DataType=Int64
    ;L::Integer, f=nothing, k::Integer=0, p::Integer=1, N::Union{Nothing, Integer}=nothing, a::Integer=1,
    base::Integer=2, alloc::Integer=1000, threaded::Bool=true, small_N::Bool=true
)
    len, check_a = divrem(L, a)
    k = mod(k, len)
    @assert iszero(check_a) "Length of unit-cell $a incompatible with L=$L"
    @assert iszero(k) || isequal(2k, len) "Momentum incompatible with parity."
    @assert isone(p) || isone(-p) "Invalid parity"
    
    I, R = begin
        N2 = [2len/ sqrt(i) for i = 1:len]
        N1 = N2 ./ sqrt(2)
        judge = if small_N || isnothing(N)
            TranslationParityJudge(f, reverse, k, a, p, len, base, vcat(N1, N2))
        else
            num = L*(base-1)-N
            g = isnothing(f) ? x -> (sum(x) == num) : x -> (sum(x) == num && f(x))
            TranslationParityJudge(g, reverse, k, a, p, len, base, vcat(N1, N2))
        end
        if small_N && !isnothing(N)
            selectindexnorm_N(judge, L, N, base=base, dtype=dtype)
        else
            threaded ? selectindexnorm_threaded(judge, L, base=base, alloc=alloc) : selectindexnorm(judge, L, 1:base^L, base=base, alloc=alloc)
        end
    end
    C = iszero(k) ? fill(1, len) : [iseven(i) ? 1 : -1 for i=0:len-1]
    TranslationParityBasis(zeros(Int, L), I, R, C, p, a, base)
end

#-------------------------------------------------------------------------------------------------------------------------
# Indexing
#-------------------------------------------------------------------------------------------------------------------------
"""
translation_parity_index(parity, b::AbstractTranslationalParityBasis)

Return normalization and index.
"""
function translation_parity_index(parity, b::AbstractTranslationalParityBasis)
    # Check the following:
    # - reflectQ : Whether the minimum index is in the other parity sector.
    # - min_ind  : Minimum index.
    # - min_trans: Minimum M that Tᴹ⋅|dgt⟩ = |Im⟩ or Tᴹ⋅P⋅|dgt⟩ = |Im⟩.
    reflectQ, min_ind, min_trans = begin
        Im1, T1 = translation_index(b.dgt, b.B, dtype=int_type(b), a=b.A)
        Im2, T2 = translation_index(parity(b.dgt), b.B, dtype=int_type(b), a=b.A)
        Im2 < Im1 ? (true, Im2, T2) : (false, Im1, T1)
    end

    # Search for minimum index
    ind = binary_search(b.I, min_ind)
    iszero(ind) && return 0.0, 1

    # Calculate norm, takin into account the phase from parity.
    n = reflectQ ? b.P * b.C[min_trans+1] : b.C[min_trans+1]
    return n * b.R[ind], ind
end
#-------------------------------------------------------------------------------------------------------------------------
index(b::TranslationParityBasis) = translation_parity_index(reverse, b)



