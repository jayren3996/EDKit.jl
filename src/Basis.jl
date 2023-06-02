export AbstractBasis, content, base, index, change!
"""
    AbstractBasis

Abstract type for all concrete `Basis` types. 
A typical concrete `Basis` type usually has the field:

- `dgt`: Digits representation of a state. For a symmetric basis like translational-invariant basis, `dgt` is just a representation of the state.
- `I`  : List of indices in this basis.
- `R`  : List of normalization of each state in the basis.
- `B`  : Base of the basis.

Correspondingly, there are property functions for a `Basis` type:

- `digits` : Return digits representation.
- `content`: Return list of indices in this basis.
- `norm`   : Return list of norm.
- `base`   : Return the base of the basis.

In addition, we have the following property functions:

- `eltype` : Used for type promotion when constructing matrix representation of an `Operator`.
- `length` : Return the length of the digits of the basis.
- `size`   : Return the size of the matrix that can be constructed for an `Operator`.
- `copy`   : Copy basis with a deep copied `dgt`.

The most inportant function for a `Basis` type is the index of a given digits and the digits representation of its ith content:

- `index`  : Return the norm and the index of the digits of the given basis.
- `change!`: Change the digits to the target index, and return the new normalization.

If we want to calculate the schmidt decomposition (in real space), we can use the `schmidt` function.
"""
abstract type AbstractBasis end
abstract type AbstractOnsiteBasis <: AbstractBasis end
abstract type AbstractPermuteBasis <: AbstractBasis end
abstract type AbstractTranslationalParityBasis <: AbstractPermuteBasis end
#-------------------------------------------------------------------------------------------------------------------------
# Default function definitions:
@inline digits(b::AbstractBasis) = b.dgt
@inline content(b::AbstractBasis) = b.I
@inline content(b::AbstractBasis, i::Integer) = b.I[i]
@inline norm(b::AbstractBasis) = b.R
@inline norm(b::AbstractBasis, i::Integer) = b.R[i]
@inline base(b::AbstractBasis) = b.B
@inline eltype(::AbstractBasis) = ComplexF64
@inline length(b::AbstractBasis) = length(b.dgt)
@inline size(b::AbstractBasis, i::Integer) = isone(i) || isequal(i, 2) ? length(b.I) : 1
@inline size(b::AbstractBasis) = (size(b, 1), size(b, 2))
@inline change!(b::AbstractBasis, i::Integer) = (change!(b.dgt, content(b, i), base=b.B); norm(b, i))

#-----------------------------------------------------------------------------------------------------
# Onsite Basis
#-----------------------------------------------------------------------------------------------------
export TensorBasis
"""
    TensorBasis

Basis without any symmetries. The index is computed by:

`I(i₁, i₂, ⋯, iₙ) = i₁ B^(n-1) + i₂ B^(n-2) + ⋯ + iₙ`
"""
struct TensorBasis <: AbstractOnsiteBasis
    dgt::Vector{Int64}
    B::Int64
    TensorBasis(dgt::Vector{Int64}, B::Integer) = new(dgt, Int64(B))
    TensorBasis(dgt::AbstractVector{<:Integer}, B::Integer) = new(collect(Int64, dgt), Int64(B))
    TensorBasis(L::Integer; base::Integer=2) = new(zeros(Int64, L), Int64(base))
end
TensorBasis(;L::Integer, base::Integer=2) = TensorBasis(L, base=base)
#-------------------------------------------------------------------------------------------------------------------------
@inline content(b::TensorBasis) = 1:b.B^length(b.dgt)
@inline content(::TensorBasis, i::Integer) = i
@inline norm(b::TensorBasis) = ones(Int64, B^length(b.dgt))
@inline norm(b::TensorBasis, ::Integer) = 1
@inline eltype(::TensorBasis) = Int64
@inline size(b::TensorBasis, i::Integer) = isone(i) || isequal(i, 2) ? b.B^length(b.dgt) : 1
@inline size(b::TensorBasis) = (l=b.B^length(b.dgt); (l, l))
@inline index(b::TensorBasis) = 1, index(b.dgt, base=b.B)
@inline copy(b::TensorBasis) = TensorBasis(deepcopy(b.dgt), b.B)

#-------------------------------------------------------------------------------------------------------------------------
# ProjectBasis
#-------------------------------------------------------------------------------------------------------------------------
export ProjectedBasis
"""
    ProjectedBasis

Basis for subspace that is spanned only by product states.
It is basically like the `TensorBasis`, but contains only a subset of all basis vectors, selected by a given function.

Properties:
-----------
- `dgt`: Digits.
- `I`  : List of indicies.
- `B`  : Base.
"""
struct ProjectedBasis <: AbstractOnsiteBasis
    dgt::Vector{Int64}
    I::Vector{Int64}
    B::Int64
    ProjectedBasis(dgt::Vector{Int64}, I::Vector{Int64}, B::Integer) = new(dgt, I, Int64(B))
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    ProjectedBasis(;L, N, base=2, alloc=1000, threaded=true)

Construction method for `ProjectedBasis` with fixed particle number (or total U(1) charge).

Inputs:
-------
- `L`       : Length of the system.
- `N`       : Quantum number of up spins / particle number / U(1) charge.
- `base`    : Base, default = 2.
- `alloc`   : Size of the prealloc memory for the basis content, used only in multithreading, default = 1000.
- `threaded`: Whether use the multithreading, default = true.

Outputs:
--------
- `b` : ProjectedBasis.
"""
function ProjectedBasis(
    ;L::Integer, f=nothing, N::Union{Nothing, Integer}=nothing,
    base::Integer=2, alloc::Integer=1000, 
    threaded::Bool=true, small_N::Bool=true
)
    I = if isnothing(N)
        threaded ? selectindex_threaded(f, L, base=base, alloc=alloc) : selectindex(f, L, 1:base^L, base=base, alloc=alloc)
    elseif small_N
        selectindex_N(f, L, N, base=base)
    else
        num = L*(base-1)-N
        g = isnothing(f) ? x -> sum(x) == num : x -> (sum(x) == num && f(x))
        threaded ? selectindex_threaded(g, L, base=base, alloc=alloc) : selectindex(g, L, 1:base^L, base=base, alloc=alloc)
    end
    ProjectedBasis(zeros(Int64, L), I, base)
end
#-------------------------------------------------------------------------------------------------------------------------
eltype(::ProjectedBasis) = Int
copy(b::ProjectedBasis) = ProjectedBasis(deepcopy(b.dgt), b.I, b.B)
change!(b::ProjectedBasis, i::Integer) = (change!(b.dgt, b.I[i], base=b.B); 1)
#-------------------------------------------------------------------------------------------------------------------------
function index(b::ProjectedBasis; check::Bool=true)
    i = index(b.dgt, base=b.B)
    ind = binary_search(b.I, i)
    ind > 0 && return 1, ind
    check ? error("No such symmetry.") : return zero(eltype(b)), 1
end

#-------------------------------------------------------------------------------------------------------------------------
# TranslationalBasis
#-------------------------------------------------------------------------------------------------------------------------
export TranslationalBasis
"""
    TranslationalBasis

Basis for subspace that is spanned by momentum states, which can also incorporate projected restriction.

Properties:
-----------
- `dgt`: Digits.
- `I`  : List of indicies.
- `R`  : List of normalization.
- `C`  : Unit phase factor.
- `B`  : Base.
"""
struct TranslationalBasis{T <: Number} <: AbstractPermuteBasis
    dgt::Vector{Int64}
    I::Vector{Int}
    R::Vector{Float64}
    C::Vector{T}
    B::Int64
    TranslationalBasis(dgt::Vector{Int64}, I, R, C::Vector{T}, B::Integer) where T = new{T}(dgt, I, R, C, Int64(B))
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    TranslationJudge

Structure used for selecting the basis contents.
"""
struct TranslationJudge
    F                   # Projective selection
    K::Int              # Momentum
    B::Int              # Base
    C::Vector{Float64}  # Normalization coefficient
end
#-------------------------------------------------------------------------------------------------------------------------
function (judge::TranslationJudge)(dgt::AbstractVector{<:Integer}, i::Integer)
    isnothing(judge.F) || judge.F(dgt) || return false, 0.0
    c, r = translation_check(dgt, i, judge.B)
    c && iszero(mod(r * judge.K, length(dgt))) && return true, judge.C[r]
    return false, 0.0
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    TranslationalBasis(f, k, L; base=2, alloc=1000, threaded=true)

Construction for `TranslationalBasis`.

Inputs:
-------
- `L`       : Length of the system.
- `f`       : Selection function for the basis contents.
- `k`       : Momentum number from 0 to L-1.
- `N`       : Particle number.
- `base`    : Base, default = 2.
- `alloc`   : Size of the prealloc memory for the basis content, used only in multithreading, default = 1000.
- `threaded`: Whether use the multithreading, default = true.

Outputs:
--------
- `b`: TranslationalBasis.
"""
function TranslationalBasis(
    ;L::Integer, f=nothing, k::Integer=0, N::Union{Nothing, Integer}=nothing,
    base::Integer=2, alloc::Integer=1000, threaded::Bool=true, small_N::Bool=true
)
    judge = if small_N || isnothing(N)
        TranslationJudge(f, k, base, [L/sqrt(i) for i = 1:L])
    else
        num = L*(base-1)-N
        g = isnothing(f) ? x -> (sum(x) == num) : x -> (sum(x) == num && f(x))
        TranslationJudge(g, k, base, [L/sqrt(i) for i = 1:L])
    end
    I, R = if small_N
        selectindexnorm_N(judge, L, N, base=base)
    else
        threaded ? selectindexnorm_threaded(judge, L, base=base, alloc=alloc) : selectindexnorm(judge, L, 1:base^L, base=base, alloc=alloc)        
    end
    C = if iszero(k)
        fill(1.0, L)
    elseif isequal(2k, L)
        [iseven(i) ? 1.0 : -1.0 for i=0:L-1] 
    else
        [exp(-1im*2π/L*k*i) for i=0:L-1]
    end
    TranslationalBasis(zeros(Int64, L), I, R, C, base)
end
#-------------------------------------------------------------------------------------------------------------------------
eltype(::TranslationalBasis{T}) where T = T
copy(b::TranslationalBasis) = TranslationalBasis(deepcopy(b.dgt), b.I, b.R, b.C, b.B)
#-------------------------------------------------------------------------------------------------------------------------
"""
    index(b::TranslationalBasis)

Index of the state in `TranslationalBasis`, together with the normalization.

Inputs:
-------
- `b`: TranslationalBasis.

Outputs:
--------
- `N`  : Normalization for the state `b.dgt`.
- `ind`: Index for the state `b.dgt`.

Notes:
------
To calculate the index of a given digits, we first shift the digits to that of the minimum index, then search for the place in the basis.
Because of the restriction from the momentum, some state with zero normalization (which is not included in the basis) will appear.
To avoid exception in the matrix constructon of `Operation`, we allow the index to not in the basis content.
When this happend, we return index 1, and normalization 0, so it has no effect on the matrix being filled.
"""
function index(b::TranslationalBasis)
    Im, T = translation_index(b.dgt, b.B)
    i = binary_search(b.I, Im)
    # We allow the case where i is not in the basis.
    # In such case we add 0 to index 1.
    iszero(i) && return zero(eltype(b)), 1
    return b.C[T+1] * b.R[i], i
end

#-------------------------------------------------------------------------------------------------------------------------
# Translational + Parity Basis
#-------------------------------------------------------------------------------------------------------------------------
export TranslationParityBasis, TranslationFlipBasis
"""
    TranslationParityBasis

Basis with translational and reflection symmetries.

Properties:
-----------
- `dgt`: Digits.
- `I`  : Representing states.
- `R`  : Normalization.
- `K`  : {±1}, momentum phase.
- `P`  : {±1}, parity.
- `B`  : Base.
"""
struct TranslationParityBasis <: AbstractTranslationalParityBasis
    dgt::Vector{Int64}    # Digits
    I::Vector{Int}          # Representing states
    R::Vector{Float64}      # Normalization
    C::Vector{Int}          # Momentum phase exp(1im * 0/π) = {±1}
    P::Int                  # {±1}, parity
    B::Int64                # Base
    TranslationParityBasis(dgt, I, R, C, P, B::Integer) = new(dgt, I, R, C, P, Int64(B))
end
#-------------------------------------------------------------------------------------------------------------------------
struct TranslationParityJudge
    F                   # Projective selection
    parity              # Parity operation acting on digits
    K::Int64            # Momentum
    P::Int64            # Parity eigenvalue {±1}
    L::Int64            # Length
    B::Int64            # Base
    C::Vector{Float64}  # Normalization coefficients
end
#-------------------------------------------------------------------------------------------------------------------------
function (judge::TranslationParityJudge)(dgt::AbstractVector{<:Integer}, i::Integer)
    isnothing(judge.F) || judge.F(dgt) || return false, 0.0
    
    isrep, reflect, r, m = translation_parity_check(judge.parity, dgt, i, judge.B)
    isrep || return false, 0.0
    iszero(mod(judge.K * r, judge.L)) || return false, 0.0
    reflect || return true, judge.C[r]

    trans = mod(judge.K * m, judge.L)
    isequal(2 * trans, judge.L) && isequal(judge.P, -1) && return true, judge.C[judge.L+r]
    iszero(trans) && isone(judge.P) && return true, judge.C[judge.L+r]
    return false, 0.0
end
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
- `base`    : Base, default = 2.
- `alloc`   : Size of the prealloc memory for the basis content, used only in multithreading, default = 1000.
- `threaded`: Whether use the multithreading, default = true.

Outputs:
--------
- `b`: TranslationParityBasis.
"""
function TranslationParityBasis(
    ;L::Integer, f=nothing, k::Integer=0, p::Integer=1, N::Union{Nothing, Integer}=nothing,
    base::Integer=2, alloc::Integer=1000, threaded::Bool=true, small_N::Bool=true
) 
    @assert iszero(k) || isequal(2k, L) "Momentum incompatible with parity."
    @assert isone(p) || isequal(p, -1) "Invalid parity"
    N2 = [2L / sqrt(i) for i = 1:L]
    N1 = N2 ./ sqrt(2)
    judge = if small_N || isnothing(N)
        TranslationParityJudge(f, reverse, k, p, L, base, vcat(N1, N2))
    else
        num = L*(base-1)-N
        g = isnothing(f) ? x -> (sum(x) == num) : x -> (sum(x) == num && f(x))
        TranslationParityJudge(g, reverse, k, p, L, base, vcat(N1, N2))
    end
    I, R = if small_N
        selectindexnorm_N(judge, L, N, base=base)
    else
        threaded ? selectindexnorm_threaded(judge, L, base=base, alloc=alloc) : selectindexnorm(judge, L, 1:base^L, base=base, alloc=alloc)
    end
    C = iszero(k) ? fill(1.0, L) : [iseven(i) ? 1.0 : -1.0 for i=0:L-1]
    TranslationParityBasis(zeros(Int, L), I, R, C, p, base)
end
#-------------------------------------------------------------------------------------------------------------------------
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
- `B`  : Base.
"""
struct TranslationFlipBasis{T} <: AbstractTranslationalParityBasis
    dgt::Vector{Int64}    # Digits
    I::Vector{Int}          # Representing states
    R::Vector{Float64}      # Normalization
    C::Vector{T}            # Momentum phase exp(-1im * 2π/L * K)
    P::Int                  # {±1}, parity
    B::Int                  # Base
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
    ;f=x->true, k::Integer=0, p::Integer=1, L::Integer,
    base::Integer=2, alloc::Integer=1000, threaded::Bool=false
)
    @assert isone(p) || isequal(p, -1) "Invalid parity"
    judge = begin
        N2 = [2L / sqrt(i) for i = 1:L]
        N1 = N2 ./ sqrt(2)
        TranslationParityJudge(f, x -> spinflip(x, base), k, p, L, base, vcat(N1, N2))
    end
    I, R = threaded ? selectindexnorm_threaded(judge, L, base=base, alloc=alloc) : selectindexnorm(judge, L, 1:base^L, base=base, alloc=alloc)
    C = if iszero(k)
        fill(1.0, L)
    elseif isequal(2k, L)
        [iseven(i) ? 1.0 : -1.0 for i=0:L-1]
    else
        [exp(-1im * 2π/L * k * i) for i=0:L-1]
    end
    TranslationFlipBasis(zeros(Int64, L), I, R, C, p, base)
end
#-------------------------------------------------------------------------------------------------------------------------
eltype(::TranslationParityBasis) = Float64
eltype(::TranslationFlipBasis{T}) where T = T
copy(b::TranslationParityBasis) = TranslationParityBasis(deepcopy(b.dgt), b.I, b.R, b.C, b.P, b.B)
copy(b::TranslationFlipBasis) = TranslationFlipBasis(deepcopy(b.dgt), b.I, b.R, b.C, b.P, b.B)
#-------------------------------------------------------------------------------------------------------------------------
function translation_parity_index(parity, b::AbstractTranslationalParityBasis)
    reflect, i, t = translation_parity_index(parity, b.dgt, b.B)
    ind = binary_search(b.I, i)
    iszero(ind) && return 0.0, 1
    n = reflect ? b.P * b.C[t+1] : b.C[t+1]
    return n * b.R[ind], ind
end
#-------------------------------------------------------------------------------------------------------------------------
index(b::TranslationParityBasis) = translation_parity_index(reverse, b)
index(b::TranslationFlipBasis) = translation_parity_index(x -> spinflip(x, b.B), b)

