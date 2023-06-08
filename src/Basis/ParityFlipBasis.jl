#-------------------------------------------------------------------------------------------------------------------------
# Parity + Flip
#-------------------------------------------------------------------------------------------------------------------------
export ParityFlipBasis
"""
    FlipBasis

Basis with translational and reflection symmetries.
"""
struct ParityFlipBasis{Ti <: Integer} <: AbstractPermuteBasis
    dgt::Vector{Ti}       # Digits
    I::Vector{Ti}         # Representing states
    R::Vector{Float64}    # Normalization
    P::Int                # {±1}, parity
    Z::Int                # {±1}, spin-flip parity
    M::Int                # base^length + 1
    B::Ti                 # Base
end
#-------------------------------------------------------------------------------------------------------------------------
struct ParityFlipJudge{T}
    F                     # Projective selection
    P::Int64              # Parity
    Z::Int64              # Spin-flip parity
    B::T                  # Base
    MAX::Int              # base^length + 1
    C::Vector{Float64}    # [2, 2sqrt(2), 0, 4]
end
#-------------------------------------------------------------------------------------------------------------------------
function double_parity_check(dgt, i, base, max, p, z)
    (Ip = rindex(dgt, base=base)) < i && return (false, 1)
    (Iz = max - i) < i && return (false, 3)
    (Ipz = max - Ip) < i && return (false, 3)
    N = 1
    if isequal(Ip, i)
        isone(p) || return (false, 3)
        N += 1
    end
    if isequal(Iz, i)
        isone(z) || return (false, 3)
        N += 1
    end
    if isequal(Ipz, i)
        isone(p*z) || return (false, 3)
        N += 1
    end            
    true, N
end
#-------------------------------------------------------------------------------------------------------------------------
function (judge::ParityFlipJudge)(dgt::AbstractVector{<:Integer}, i::Integer)
    # If there is a projective selection function F, check `F(dgt)` first
    isnothing(judge.F) || judge.F(dgt) || return (false, 0.0)
    
    Q, n = double_parity_check(dgt, i, judge.B, judge.MAX, judge.P, judge.Z)
    Q, judge.C[n]
end
#-------------------------------------------------------------------------------------------------------------------------
function ParityFlipBasis(
    dtype::DataType=Int64; L::Integer, f=nothing, p::Integer, z::Integer, N::Union{Nothing, Integer}=nothing,
    base::Integer=2, alloc::Integer=1000, threaded::Bool=true, small_N::Bool=true
)
    @assert isone(p) || isone(-p) "Invalid parity"
    @assert isnothing(N) || (isodd(base) && isequal(2N, L*(base-1))) "N = $N not compatible."
    base = convert(dtype, base)
    MAX = base ^ L + 1
    I, R = begin
        C = [2.0, 2*sqrt(2), 0.0, 4.0]
        judge = if small_N || isnothing(N)
            ParityFlipJudge(f, p, z, base, MAX, C)
        else
            num = L*(base-1)-N
            g = isnothing(f) ? x -> (sum(x) == num) : x -> (sum(x) == num && f(x))
            ParityFlipJudge(g, p, z, base, MAX, C)
        end
        if small_N && !isnothing(N)
            selectindexnorm_N(judge, L, N, base=base)
        else
            threaded ? selectindexnorm_threaded(judge, L, base=base, alloc=alloc) : selectindexnorm(judge, L, 1:base^L, base=base, alloc=alloc)
        end
    end
    ParityFlipBasis(zeros(dtype, L), I, R, p, z, MAX, base)
end
#-------------------------------------------------------------------------------------------------------------------------
function index(b::ParityFlipBasis)
    I0 = index(b.dgt, base=b.B)
    Ip = rindex(b.dgt, base=b.B)
    Iz = b.M - I0
    Ipz = b.M - Ip
    Ir = min(I0, Iz, Ip, Ipz)
    i = binary_search(b.I, Ir)
    iszero(i) && return (0.0, one(b.B))
    N = if isequal(Ir, I0)
        b.R[i]
    elseif isequal(Ir, Ip) 
        b.P * b.R[i]
    elseif isequal(Ir, Iz)
        b.Z * b.R[i]
    else
        b.P * b.Z * b.R[i]
    end
    N, i
end

