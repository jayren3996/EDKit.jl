export TranslationalFermionBasis
"""
    TranslationalFermionBasis{Ti<:Integer, T<:Number}

Momentum-resolved basis for **spinless fermions** at fixed particle number `N`.

This is the fermionic analogue of [`TranslationalBasis`](@ref). The many-body
translation operator on the Jordan-Wigner occupation representation differs from
the bosonic site shift by a fermion sign: translating by one site moves the
occupied boundary site past the other `N−1` particles, contributing

    (−1)^((N−1)·n_wrap)

per unit translation, where `n_wrap` is the occupation of the wrapping site. For
**odd** `N` this sign is always `+1` (periodic); for **even** `N` it is `−1`
whenever the boundary site is occupied (antiperiodic), which shifts the
compatible momenta. The sign is baked into the orbit phase returned by
[`index`](@ref), so Jordan-Wigner operators built with
[`trans_inv_fermion_operator`](@ref) act correctly within a momentum sector.

Restrictions (MVP): single-site translation (`a = 1`), `base = 2`, fixed `N`.
"""
struct TranslationalFermionBasis{Ti <: Integer, T <: Number} <: AbstractPermuteBasis
    dgt::Vector{Ti}
    I::Vector{Ti}
    R::Vector{Float64}
    C::Vector{T}
    A::Int64
    B::Ti
    evenN::Bool          # whether the fermion translation sign is active (N even)
end
#-------------------------------------------------------------------------------------------------------------------------
eltype(::TranslationalFermionBasis{Ti, T}) where {Ti, T} = T
copy(b::TranslationalFermionBasis) =
    TranslationalFermionBasis(deepcopy(b.dgt), b.I, b.R, b.C, b.A, b.B, b.evenN)
order(b::TranslationalFermionBasis) = ncycle(b)

#-------------------------------------------------------------------------------------------------------------------------
# Selector
#-------------------------------------------------------------------------------------------------------------------------
"""
    TranslationFermionJudge{TF, T}

Internal selector for [`TranslationalFermionBasis`](@ref): like
`TranslationJudge` but accumulates the fermion translation sign while cycling the
orbit and applies the sign-corrected momentum-compatibility condition.
"""
struct TranslationFermionJudge{TF, T}
    F::TF
    K::Int
    L::Int
    B::T
    C::Vector{Float64}
    evenN::Bool
end

# A representative with period `p` survives in sector `k` iff the combined
# translation phase φ_p = exp(-i 2π k p / L) times the accumulated fermion sign σ
# equals 1. Since σ = ±1, φ_p must itself be ±1 and match σ.
@inline function _fermion_check_momentum(p::Integer, k::Integer, L::Integer, σ::Int)
    r = mod(2 * k * p, 2 * L)
    r == 0 && return σ == 1
    r == L && return σ == -1
    false
end

function (judge::TranslationFermionJudge)(dgt::AbstractVector{<:Integer}, i::Integer)
    isnothing(judge.F) || judge.F(dgt) || return (false, 0.0)
    L = length(dgt)
    σ = 1
    for n in 1:judge.L-1
        wrap = dgt[L]                        # boundary site that wraps this step
        circshift!(dgt, 1)
        judge.evenN && isone(wrap) && (σ = -σ)
        In = index(dgt, base=judge.B)
        if In < i
            return (false, 0.0)
        elseif isequal(In, i)
            _fermion_check_momentum(n, judge.K, judge.L, σ) || return (false, 0.0)
            return (true, judge.C[n])
        end
    end
    # Full period L: σ = +1 (N(N-1) even) and φ_L = 1, so always compatible.
    true, judge.C[judge.L]
end

#-------------------------------------------------------------------------------------------------------------------------
# Construction
#-------------------------------------------------------------------------------------------------------------------------
"""
    TranslationalFermionBasis(dtype::DataType=Int64; L, N, k=0, f=nothing,
                              base=2, alloc=1000, threaded=true, small_N=true)

Construct a momentum-resolved spinless-fermion basis on `L` sites at fixed
particle number `N` (occupied-site count) and momentum label `k` (physical
momentum `2πk/L`).

Restrictions: `base = 2`, single-site translation. `N` is required.

# Example
```julia
L, N = 8, 4
Hk = sum(eigvals(Hermitian(Array(
        let B = TranslationalFermionBasis(L=L, N=N, k=k)
            Hhop = trans_inv_fermion_operator("+-", [1, 2], B)
            -(Hhop + adjoint(Hhop))
        end))) for k in 0:L-1)   # union of k-sector spectra = full N-sector spectrum
```
"""
function TranslationalFermionBasis(dtype::DataType=Int64;
    L::Integer, N::Integer, k::Integer=0, f=nothing,
    base::Integer=2, alloc::Integer=1000, threaded::Bool=true, small_N::Bool=true,
)
    base == 2 || error("TranslationalFermionBasis currently supports base=2 (spinless fermions); got base=$base.")
    0 ≤ N ≤ L || error("N=$N is out of range 0:$L.")
    len = L                                    # a = 1
    k = mod(-k, len)                           # EDKit convention T|k⟩ = exp(-ik)|k⟩
    _check_index_capacity(dtype, base, L)
    base = convert(dtype, base)
    norm = [len / sqrt(i) for i in 1:len]
    evenN = iseven(N)
    judge = TranslationFermionJudge(f, k, len, base, norm, evenN)
    # Enumerate fermion-N (occupied) states. The fixed-weight enumerator counts
    # in the projected (dgt=0) convention, so pass `L - N` (mirrors
    # SpinlessFermionBasis).
    I, R = _run_selectindexnorm(judge, L, L - N, base, alloc, threaded, small_N)
    C = phase_factor(k, len)
    TranslationalFermionBasis(zeros(dtype, L), I, R, C, 1, base, evenN)
end

#-------------------------------------------------------------------------------------------------------------------------
# Indexing
#-------------------------------------------------------------------------------------------------------------------------
"""
    index(b::TranslationalFermionBasis, dgt::AbstractVector)

Return `(N, i)`: the orbit phase/normalization `N` (momentum phase × fermion
sign × `R[i]`) and the representative index `i` of the digit string `dgt`. A
state with zero normalization (incompatible momentum) returns `(0, 1)`.
"""
function index(b::TranslationalFermionBasis, dgt::AbstractVector)
    I0 = index(dgt, base=b.B)
    L = length(dgt)
    TI = eltype(b.I)
    state = I0 - one(TI)
    mask = (one(TI) << L) - one(TI)
    Im, M, σM, σ = I0, 0, 1, 1
    for n in 1:ncycle(b)-1
        wrap = state & one(TI)                 # dgt[L] before the shift
        state = ((state >> 1) | (state << (L - 1))) & mask
        b.evenN && isone(wrap) && (σ = -σ)
        In = state + one(TI)
        isequal(In, I0) && break
        if In < Im
            Im, M, σM = In, n, σ
        end
    end
    i = binary_search(b.I, Im)
    iszero(i) && return (zero(eltype(b)), one(b.B))
    @inbounds N = b.C[M+1] * σM * b.R[i]
    N, i
end
