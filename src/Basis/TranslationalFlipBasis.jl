export TranslationFlipBasis
"""
    TranslationFlipBasis

Basis with translational and spin-flip symmetries.
"""
struct TranslationFlipBasis{Ti, T} <: AbstractTranslationalParityBasis
    dgt::Vector{Ti}         # Digits
    I::Vector{Ti}           # Representing states
    R::Vector{Float64}      # Normalization
    C::Vector{T}            # Momentum phase exp(-1im * 2π/L * K)
    P::Int                  # {±1}, parity
    A::Int                  # Length of unit cell
    M::Int                  # MAX = base^length + 1
    B::Ti                   # Base
end
#-------------------------------------------------------------------------------------------------------------------------
eltype(::TranslationFlipBasis{Ti, T}) where {Ti, T} = T
copy(b::TranslationFlipBasis) = TranslationFlipBasis(deepcopy(b.dgt), b.I, b.R, b.C, b.P, b.A, b.B)

#-------------------------------------------------------------------------------------------------------------------------
# Functions selecting the basis
#-------------------------------------------------------------------------------------------------------------------------
struct TranslationFlipJudge{T <: Integer}
    F                   # Projective selection
    K::Int64            # Momentum
    A::Int64            # Length of unit cell
    P::Int64            # Parity eigenvalue {±1}
    L::Int64            # Length of translation
    B::T                # Base
    C::Vector{Float64}  # Normalization coefficients
    MAX::Int            # Base^length + 1
end
#-------------------------------------------------------------------------------------------------------------------------
function (judge::TranslationFlipJudge)(dgt::AbstractVector{<:Integer}, i::Integer)
    # First check projective selection rule
    isnothing(judge.F) || judge.F(dgt) || return (false, 0.0)

    R = judge.L
    Qp = false

    # Check Flip
    In = judge.MAX - i
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
        In2 = judge.MAX - In1

        # Check cycled index
        In1 < i && return (false, 0.0)
        In2 < i && return (false, 0.0)

        if isequal(In2, i)
            # Check flip parity
            isequal(isone(judge.P), iseven(div(n * judge.K, judge.L÷2))) || return (false, 0.0)
            Qp = true
        end

        if isequal(In1, i)  
            # Find periodicity; check Momentum 
            iszero(mod(n * judge.K , judge.L)) || return (false, 0.0)
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
TranslationFlipBasis(f, k, p, L; base=2, alloc=1000, threaded=true)

Construction for `TranslationFlipBasis`.

Inputs:
-------
- `f`       : Selection function for the basis contents.
- `k`       : Momentum number from 0 to L-1.
- `p`       : Parity number (under spin flip) ±1.
- `L`       : Length of the system.
- `base`    : Base, default = 2.
- `alloc`   : Size of the prealloc memory for the basis content, used only in multithreading, default = 1000.
- `threaded`: Whether use the multithreading, default = true.

Outputs:
--------
- `b`: TranslationFlipBasis.
"""
function TranslationFlipBasis(
    dtype::DataType=Int64; f=x->true, k::Integer=0, p::Integer=1, L::Integer, N::Union{Nothing, Integer}=nothing,
    a::Integer=1, base::Integer=2, alloc::Integer=1000, threaded::Bool=false, small_N::Bool=false
)
    len, check_a = divrem(L, a)
    @assert iszero(check_a) "Length of unit-cell $a incompatible with L=$L"
    @assert isone(p) || isone(-p) "Invalid parity"
    @assert isnothing(N) || (isodd(base) && isequal(2N, L*(base-1))) "N = $N not compatible."
    k = mod(k, len)
    base = convert(dtype, base)
    MAX = base ^ L + 1

    I, R = begin
        N2 = [2len/ sqrt(i) for i = 1:len]
        N1 = N2 ./ sqrt(2)
        judge = if small_N || isnothing(N)
            TranslationFlipJudge(f, k, a, p, len, base, vcat(N1, N2), MAX)
        else
            num = L*(base-1)÷2
            g = isnothing(f) ? x -> (sum(x) == num) : x -> (sum(x) == num && f(x))
            TranslationFlipJudge(g, k, a, p, len, base, vcat(N1, N2), MAX)
        end
        if small_N && !isnothing(N)
            selectindexnorm_N(judge, L, N, base=base)
        else
            threaded ? selectindexnorm_threaded(judge, L, base=base, alloc=alloc) : selectindexnorm(judge, L, 1:base^L, base=base, alloc=alloc)
        end
    end
    C = phase_factor(k, len)
    TranslationFlipBasis(zeros(dtype, L), I, R, C, p, a, MAX, base)
end

#-------------------------------------------------------------------------------------------------------------------------
# Indexing
#-------------------------------------------------------------------------------------------------------------------------
"""
Return normalization and index.
"""
function index(b::TranslationFlipBasis)
    Ia0 = index(b.dgt, base=b.B)
    Ib0 = b.M - Ia0
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
        Ib1 = b.M - Ia1
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

