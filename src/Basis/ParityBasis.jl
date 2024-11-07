#-------------------------------------------------------------------------------------------------------------------------
# Parity
#-------------------------------------------------------------------------------------------------------------------------
export ParityBasis
"""
    ParityBasis

Basis with translational and reflection symmetries.

Properties:
-----------
- dgt: Digits.
- I  : Representing states.
- R  : Normalization.
- P  : {±1}, parity.
- B  : Base.
"""
struct ParityBasis{Ti <: Integer} <: AbstractPermuteBasis
    dgt::Vector{Ti}       # Digits
    I::Vector{Ti}         # Representing states
    R::Vector{Float64}    # Normalization
    P::Int                # {±1}, parity
    B::Ti                 # Base
end
#-------------------------------------------------------------------------------------------------------------------------
struct ParityJudge{T}
    F                   # Projective selection
    P::Int64            # Parity
    B::T                # Base
    C::Float64          # Value sqrt(2).
end
#-------------------------------------------------------------------------------------------------------------------------
function (judge::ParityJudge)(dgt::AbstractVector{<:Integer}, i::Integer)
    # If there is a projective selection function F, check `F(dgt)` first
    isnothing(judge.F) || judge.F(dgt) || return (false, 0.0)
    
    # Check parity
    In = rindex(dgt, base=judge.B)
    In < i && return (false, 0.0)
    if isequal(In, i)
        isone(judge.P) || return (false, 0.0)
        true, 2.0
    else
        true, judge.C
    end
end
#-------------------------------------------------------------------------------------------------------------------------
function ParityBasis(
    dtype::DataType=Int64; L::Integer, f=nothing, p::Integer, N::Union{Nothing, Integer}=nothing,
    base::Integer=2, alloc::Integer=1000, threaded::Bool=true, small_N::Bool=true
)
    @assert isone(p) || isone(-p) "Invalid parity"
    base = convert(dtype, base)
    I, R = begin
        judge = if small_N || isnothing(N)
            ParityJudge(f, p, base, sqrt(2))
        else
            num = L*(base-1)-N
            g = isnothing(f) ? x -> (sum(x) == num) : x -> (sum(x) == num && f(x))
            ParityJudge(g, p, base, sqrt(2))
        end
        if small_N && !isnothing(N)
            selectindexnorm_N(judge, L, N, base=base)
        else
            threaded ? selectindexnorm_threaded(judge, L, base=base, alloc=alloc) : selectindexnorm(judge, L, 1:base^L, base=base, alloc=alloc)
        end
    end
    ParityBasis(zeros(dtype, L), I, R, p, base)
end
#-------------------------------------------------------------------------------------------------------------------------
function index(b::ParityBasis)
    Ia = index(b.dgt, base=b.B)
    Ib = rindex(b.dgt, base=b.B)
    i, n = if Ib < Ia
        binary_search(b.I, Ib), b.P
    else
        binary_search(b.I, Ia), 1
    end
    iszero(i) && return (zero(eltype(b)), one(i))
    n * b.R[i], i
end

