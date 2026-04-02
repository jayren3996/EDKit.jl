#-------------------------------------------------------------------------------------------------------------------------
# Translation + Parity
#-------------------------------------------------------------------------------------------------------------------------
export TranslationParityBasis
"""
    TranslationParityBasis

Basis reduced simultaneously by translation and reflection parity.

This basis is only valid in momentum sectors compatible with reflection, namely
the sectors where the translation phase is real (`k = 0` or zone-boundary
equivalents). The resulting basis carries both a translation phase table `C`
and a reflection eigenvalue `P`.
"""
struct TranslationParityBasis{Ti <: Integer} <: AbstractTranslationalParityBasis
    dgt::Vector{Ti}       # Digits
    I::Vector{Ti}         # Representing states
    R::Vector{Float64}    # Normalization
    C::Vector{Int}        # Momentum phase exp(1im * 0/π) = {±1}
    P::Int                # {±1}, parity
    A::Int                # Length of unit cell
    B::Ti                 # Base
end
#-------------------------------------------------------------------------------------------------------------------------
eltype(::TranslationParityBasis) = Float64
copy(b::TranslationParityBasis) = TranslationParityBasis(deepcopy(b.dgt), b.I, b.R, b.C, b.P, b.A, b.B)


#-------------------------------------------------------------------------------------------------------------------------
# Functions selecting the basis
#-------------------------------------------------------------------------------------------------------------------------
"""
    rindex(dgt; base)

Return the 1-based index of the reversed digit string.

This is the reflection counterpart of [`index(dgt; base=...)`](@ref) and is used
throughout reflection-based reduced bases.
"""
function rindex(dgt::AbstractVector{<:Integer}; base::T) where T <: Integer
    evalpoly(base, dgt) + one(T)
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    TranslationParityJudge

Internal selector for [`TranslationParityBasis`](@ref).

It simultaneously checks canonicality under translation and reflection,
momentum compatibility, and parity-sector compatibility.
"""
struct TranslationParityJudge{T <: Integer}
    F                   # Projective selection
    K::Int64            # Momentum phase exp(1im * 0/π) = {±1}
    A::Int64            # Length of unit cell
    P::Int64            # Parity eigenvalue {±1}
    L::Int64            # Length of translation
    B::T                # Base
    C::Vector{Float64}  # Normalization coefficients
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    (judge::TranslationParityJudge)(dgt, i)

Return whether state `i` should be kept in the translation-parity basis and its
normalization factor if accepted.
"""
function (judge::TranslationParityJudge)(dgt::AbstractVector{<:Integer}, i::Integer)
    # First check projective selection rule
    isnothing(judge.F) || judge.F(dgt) || return (false, 0.0)

    R = judge.L
    Qp = false

    # Check reflection
    In = rindex(dgt, base=judge.B)
    if In < i 
        return (false, 0.0)
    elseif isequal(In, i)
        Qp = true
        # Check parity 
        isone(judge.P) || return (false, 0.0)
    end

    # Cycle digits
    for n in 1:judge.L-1
        circshift!(dgt, judge.A)
        In1 = index(dgt, base=judge.B)
        In2 = rindex(dgt, base=judge.B)

        # Check cycled index
        In1 < i && return (false, 0.0)
        In2 < i && return (false, 0.0)

        if isequal(In2, i)
            # Check Parity
            isone(judge.P * judge.K ^ n) || return (false, 0.0)
            Qp = true
        end

        if isequal(In1, i)  
            # Find periodicity; check Momentum 
            isone(judge.K) || iseven(n) || return (false, 0.0)
            R = n 
            break
        end
    end

    # Calculating norm
    N = Qp ? judge.C[judge.L+R] : judge.C[R]
    true, N
end

#-------------------------------------------------------------------------------------------------------------------------
# Construction
#-------------------------------------------------------------------------------------------------------------------------
"""
    TranslationParityBasis(f, k, p, L; base=2, alloc=1000, threaded=true)

Construct a translation-plus-parity basis.

Arguments:
- `f`       : Selection function for the basis contents.
- `k`       : Momentum number, restricted to sectors compatible with parity.
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

Notes:
- This constructor only accepts parity-compatible momentum sectors.
- Internally the translation phase is stored as `±1`.
"""
function TranslationParityBasis(
    dtype::DataType=Int64; L::Integer, f=nothing, k::Integer=0, p::Integer=1, N::Union{Nothing, Integer}=nothing, 
    a::Integer=1, base::Integer=2, alloc::Integer=1000, threaded::Bool=true, small_N::Bool=false
)
    len, check_a = divrem(L, a)
    @assert iszero(check_a) "Length of unit-cell $a incompatible with L=$L"
    k = if iszero(mod(k, len))
        1
    elseif iszero(mod(2k, len))
        -1
    else
        error("Momentum $k incompatible with parity.")
    end
    @assert isone(p) || isone(-p) "Invalid parity"

    base = convert(dtype, base)
    I, R = begin
        N2 = [2len/ sqrt(i) for i = 1:len]
        N1 = N2 ./ sqrt(2)
        judge = if small_N || isnothing(N)
            TranslationParityJudge(f, k, a, p, len, base, vcat(N1, N2))
        else
            num = L*(base-1)-N
            g = isnothing(f) ? x -> (sum(x) == num) : x -> (sum(x) == num && f(x))
            TranslationParityJudge(g, k, a, p, len, base, vcat(N1, N2))
        end
        if small_N && !isnothing(N)
            selectindexnorm_N(judge, L, N, base=base)
        else
            threaded ? selectindexnorm_threaded(judge, L, base=base, alloc=alloc) : selectindexnorm(judge, L, 1:base^L, base=base, alloc=alloc)
        end
    end
    C = isone(k) ? fill(1, len) : [iseven(i) ? 1 : -1 for i=0:len-1]
    TranslationParityBasis(zeros(dtype, L), I, R, C, p, a, base)
end

#-------------------------------------------------------------------------------------------------------------------------
# Indexing
#-------------------------------------------------------------------------------------------------------------------------
"""
    index(b::TranslationParityBasis)

Return the coefficient/index pair for the current digit buffer in the
translation-parity basis.

The returned coefficient includes both translation phase and reflection parity
contributions relative to the stored representative.
"""
function index(b::TranslationParityBasis)
    Ia0 = index(b.dgt, base=b.B)
    Ib0 = rindex(b.dgt, base=b.B)
    Iam, Ibm = Ia0, Ib0
    R, M = 0, 0

    # Cycling and find monimum
    resetQ = true
    for i in 1:ncycle(b)-1
        circshift!(b.dgt, b.A)
        Ia1 = index(b.dgt, base=b.B)
        if isequal(Ia1, Ia0)
            resetQ = false
            break
        elseif Ia1 < Iam
            Iam, R = Ia1, i
        end
        Ib1 = rindex(b.dgt, base=b.B)
        if Ib1 < Ibm
            Ibm, M = Ib1, i
        end
    end
    resetQ && circshift!(b.dgt, b.A)

    # Find index and normalization
    i, n = if Ibm < Iam
        binary_search(b.I, Ibm), b.P * b.C[M+1]
    else
        binary_search(b.I, Iam), b.C[R+1]
    end
    iszero(i) && return (zero(eltype(b)), one(eltype(b.I)))
    n * b.R[i], i
end

@inline function _bitreverse(x::T, L::Int) where T <: Integer
    r = zero(T)
    v = x
    for _ in 1:L
        r = (r << 1) | (v & one(T))
        v >>= 1
    end
    r
end

@inline function _index_bits(b::TranslationParityBasis)
    Ia0 = index(b.dgt, base=b.B)
    state = Ia0 - one(eltype(b.I))
    L = length(b.dgt)
    A = b.A
    TI = eltype(b.I)
    mask = (one(TI) << L) - one(TI)

    Ib0 = _bitreverse(state, L) + one(TI)
    Iam, Ibm = Ia0, Ib0
    R, M = 0, 0

    for i in 1:ncycle(b)-1
        state = ((state >> A) | (state << (L - A))) & mask
        Ia1 = state + one(TI)
        if isequal(Ia1, Ia0)
            break
        elseif Ia1 < Iam
            Iam, R = Ia1, i
        end
        Ib1 = _bitreverse(state, L) + one(TI)
        if Ib1 < Ibm
            Ibm, M = Ib1, i
        end
    end

    i, n = if Ibm < Iam
        binary_search(b.I, Ibm), b.P * b.C[M+1]
    else
        binary_search(b.I, Iam), b.C[R+1]
    end
    iszero(i) && return (zero(eltype(b)), one(eltype(b.I)))
    @inbounds n * b.R[i], i
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    order(b::TranslationParityBasis)

Return the maximum orbit size generated by translation and reflection.
"""
order(b::TranslationParityBasis) = 2*length(b.dgt)
