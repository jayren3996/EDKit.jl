export TranslationFlipBasis
"""
    TranslationFlipBasis

Basis with translational and spin-flip symmetries.

Properties:
-----------
- `dgt`: Digits.
- `I`  : Representing states.
- `R`  : Normalization.
- `C`  : Momentum phase.
- `P`  : {±1}, parity.
- `A`  : Length of unit cell.
- `B`  : Base.
"""
struct TranslationFlipBasis{T} <: AbstractTranslationalParityBasis
    dgt::Vector{Int64}      # Digits
    I::Vector{Int}          # Representing states
    R::Vector{Float64}      # Normalization
    C::Vector{T}            # Momentum phase exp(-1im * 2π/L * K)
    P::Int                  # {±1}, parity
    A::Int                  # Length of unit cell
    B::Int                  # Base
end
#-------------------------------------------------------------------------------------------------------------------------
eltype(::TranslationFlipBasis{T}) where T = T
copy(b::TranslationFlipBasis) = TranslationFlipBasis(deepcopy(b.dgt), b.I, b.R, b.C, b.P, b.B)

#-------------------------------------------------------------------------------------------------------------------------
# Construction
#-------------------------------------------------------------------------------------------------------------------------
"""
    spinflip(v::AbstractVector{<:Integer}, base::Integer)

Flip spins Sz on each site.
"""
function spinflip(v::AbstractVector{<:Integer}, base::Integer)
    vf = Vector{eltype(v)}(undef, length(v))
    for i = 1:length(vf)
        vf[i] = base - v[i] - 1
    end
    vf
end
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
    ;f=x->true, k::Integer=0, p::Integer=1, L::Integer, a::Integer=1,
    base::Integer=2, alloc::Integer=1000, threaded::Bool=false
)
    len, check_a = divrem(L, a)
    k = mod(k, len)
    @assert iszero(check_a) "Length of unit-cell $a incompatible with L=$L"
    @assert isone(p) || isone(-p) "Invalid parity"

    I, R = begin
        N2 = [2len/ sqrt(i) for i = 1:len]
        N1 = N2 ./ sqrt(2)
        judge = TranslationParityJudge(f, x -> spinflip(x, base), k, a, p, len, base, vcat(N1, N2))
        threaded ? selectindexnorm_threaded(judge, L, base=base, alloc=alloc) : selectindexnorm(judge, L, 1:base^L, base=base, alloc=alloc)
    end
    C = if iszero(k)
        fill(1.0, len)
    elseif isequal(2k, L)
        [iseven(i) ? 1.0 : -1.0 for i=0:len-1]
    else
        [exp(-1im * 2π/len * k * i) for i=0:len-1]
    end
    TranslationFlipBasis(zeros(Int64, L), I, R, C, p, a, base)
end

#-------------------------------------------------------------------------------------------------------------------------
# Indexing
#-------------------------------------------------------------------------------------------------------------------------
index(b::TranslationFlipBasis) = translation_parity_index(x -> spinflip(x, b.B), b)

