#-----------------------------------------------------------------------------------------------------
# Tensor Basis
#-----------------------------------------------------------------------------------------------------
export TensorBasis
"""
    TensorBasis

Basis without any symmetries. The index is computed by:

`I(i₁, i₂, ⋯, iₙ) = i₁ B^(n-1) + i₂ B^(n-2) + ⋯ + iₙ`
"""
struct TensorBasis <: AbstractBasis
    dgt::Vector{Int64}
    B::Int64
    TensorBasis(dgt::Vector{Int64}, B::Integer) = new(dgt, Int64(B))
    TensorBasis(dgt::AbstractVector{<:Integer}, B::Integer) = new(collect(Int64, dgt), Int64(B))
    TensorBasis(L::Integer; base::Integer=2) = new(zeros(Int64, L), Int64(base))
end

@inline content(b::TensorBasis) = 1:b.B^length(b.dgt)
@inline content(::TensorBasis, i::Integer) = i
@inline norm(b::TensorBasis) = ones(Int64, B^length(b.dgt))
@inline norm(b::TensorBasis, ::Integer) = 1
@inline eltype(::TensorBasis) = Int64
@inline size(b::TensorBasis, i::Integer) = isone(i) || isequal(i, 2) ? b.B^length(b.dgt) : 1
@inline size(b::TensorBasis) = (l=b.B^length(b.dgt); (l, l))
@inline index(b::TensorBasis) = 1, index(b.dgt, base=b.B)

function schmidt(v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::TensorBasis)
    S = SchmidtMatrix(eltype(v), b, Ainds)
    for i = 1:length(v)
        change!(b, i)
        addto!(S, v[i])
    end
    Matrix(S)
end

#-------------------------------------------------------------------------------------------------------------------------
# Projected Basis
#-------------------------------------------------------------------------------------------------------------------------
export ProjectedBasis
"""
    ProjectedBasis

Basis for subspace that is spanned only by product states.
"""
struct ProjectedBasis <: AbstractBasis
    dgt::Vector{Int64}
    I::Vector{Int64}
    B::Int64
    ProjectedBasis(dgt::Vector{Int64}, I::Vector{Int64}, B::Integer) = new(dgt, I, Int64(B))
end

@inline eltype(::ProjectedBasis) = Int
@inline change!(b::ProjectedBasis, i::Integer) = (change!(b.dgt, b.I[i], base=b.B); 1)
@inline function index(b::ProjectedBasis)
    i = index(b.dgt, base=b.B)
    ind = binary_search(b.I, i)
    ind > 0 ? (1, ind) : error("No such symmetry.")
end

function ProjectedBasis(f, L::Integer; base::Integer=2, alloc::Integer=1000, threaded::Bool=false)
    dgt = zeros(Int64, L)
    I = threaded ? selectindex_threaded(f, L, base=base, alloc=alloc) : selectindex(f, L, 1:base^L, base=base, alloc=alloc)
    ProjectedBasis(dgt, I, base)
end

function schmidt(v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::ProjectedBasis)
    S = SchmidtMatrix(eltype(v), b, Ainds)
    for i = 1:length(v)
        change!(b, i)
        addto!(S, v[i])
    end
    Matrix(S)
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
struct TranslationalBasis{T <: Number} <: AbstractBasis
    dgt::Vector{Int64}
    I::Vector{Int}
    R::Vector{Float64}
    C::Vector{T}
    B::Int64
    TranslationalBasis(dgt::Vector{Int64}, I, R, C::Vector{T}, B::Integer) where T = new{T}(dgt, I, R, C, Int64(B))
end

@inline eltype(::TranslationalBasis{T}) where T = T

"""
    index(b::TranslationalBasis)

To calculate the index of a given digits, we first shift the digits to that of the minimum index, then search for the place in the basis.
Because of the restriction from the momentum, some state with zero normalization (which is not included in the basis) will appear.
To avoid exception in the matrix constructon of `Operation`, we allow the index to not in the basis content.
When this happend, we return index 1, and normalization 0, so it has no effect on the matrix being filled.
"""
@inline function index(b::TranslationalBasis)
    Im, T = translation_index(b.dgt, b.B)
    i = binary_search(b.I, Im)
    N, ind = if iszero(i)
        # Here we allow the case where i is not in the basis.
        # In such case we add 0 to index 1.
        zero(eltype(b)), 1
    else
        b.C[T+1] * b.R[i], i
    end
    N, ind
end

struct TranslationJudge
    F                   # Projective selection
    K::Int              # Momentum
    B::Int              # Base
    C::Vector{Float64}  # Normalization coefficient
end

function (judge::TranslationJudge)(dgt::AbstractVector{<:Integer}, i::Integer)
    Q, N = if judge.F(dgt)
        c, r = translation_check(dgt, i, judge.B)
        if c && iszero( mod(r * judge.K, length(dgt) ) )
            true, judge.C[r]
        else
            false, 0.0
        end
    else
        false, 0.0
    end
    Q, N
end

function TranslationalBasis(f, k::Integer, L::Integer; base::Integer=2, alloc::Integer=1000, threaded::Bool=false)
    dgt = zeros(Int64, L)
    judge = TranslationJudge(f, k, base, [L/sqrt(i) for i = 1:L])
    I, R = if threaded
        selectindexnorm_threaded(judge, L, base=base, alloc=alloc)
    else
        selectindexnorm(judge, L, 1:base^L, base=base, alloc=alloc)
    end
    C = if iszero(k)
        fill(1.0, L)
    elseif isequal(2k, L)
        [iseven(i) ? 1.0 : -1.0 for i=0:L-1] 
    else
        [exp(-1im*2π/L*k*i) for i=0:L-1]
    end
    TranslationalBasis(dgt, I, R, C, base)
end

function TranslationalBasis(k::Integer, L::Integer; base::Integer=2, alloc::Integer=1000, threaded::Bool=false)
    TranslationalBasis(x -> true, k, L, base=base, alloc=alloc, threaded=threaded)
end

function schmidt(v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::TranslationalBasis)
    dgt, R, phase = b.dgt, b.R, b.C[2]
    S = SchmidtMatrix(promote_type(eltype(v), eltype(b)), b, Ainds)
    for i = 1:length(v)
        change!(b, i)
        val = v[i] / R[i]
        for j in 1:length(dgt)
            addto!(S, val)
            cyclebits!(dgt)
            val *= phase
        end
    end
    Matrix(S)
end

#-------------------------------------------------------------------------------------------------------------------------
# Translational + Parity Basis
#-------------------------------------------------------------------------------------------------------------------------
export TranslationParityBasis, TranslationFlipBasis
abstract type AbstractTranslationalParityBasis <: AbstractBasis end
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

@inline eltype(::TranslationParityBasis) = Float64
@inline eltype(::TranslationFlipBasis{T}) where T = T

function translation_parity_index(parity, b::AbstractTranslationalParityBasis)
    reflect, i, t = translation_parity_index(parity, b.dgt, b.B)
    ind = binary_search(b.I, i)
    N, I = if iszero(ind)
        0.0, 1
    else
        n = reflect ? b.P * b.C[t+1] : b.C[t+1]
        n * b.R[ind], ind
    end
    N, I
end

@inline index(b::TranslationParityBasis) = translation_parity_index(reverse, b)
@inline index(b::TranslationFlipBasis) = translation_parity_index(x -> spinflip(x, b.B), b)

struct TranslationParityJudge
    F                   # Projective selection
    parity              # Parity operation acting on digits
    K::Int64            # Momentum
    P::Int64            # Parity eigenvalue {±1}
    L::Int64            # Length
    B::Int64            # Base
    C::Vector{Float64}  # Normalization coefficients
end

function (judge::TranslationParityJudge)(dgt::AbstractVector{<:Integer}, i::Integer)
    Q, N = if judge.F(dgt)
        isrep, reflect, r, m = translation_parity_check(judge.parity, dgt, i, judge.B)
        if isrep && iszero(mod(judge.K * r, judge.L))
            if reflect
                trans = mod(judge.K * m, judge.L)
                if ( isequal(2 * trans, judge.L) && isequal(judge.P, -1) ) || ( iszero(trans) && isone(judge.P) )
                    true, judge.C[judge.L+r]
                else
                    false, 0.0
                end
            else
                true, judge.C[r]
            end
        else
            false, 0.0
        end
    else
        false, 0.0
    end
    Q, N
end

function TranslationParityBasis(
    f, K::Integer, P::Integer, L::Integer; 
    base::Integer=2, alloc::Integer=1000, threaded::Bool=false
)
    @assert iszero(K) || isequal(2K, L) "Momentum incompatible with parity."
    @assert isone(P) || isequal(P, -1) "Invalid parity"
    judge = begin
        N2 = [2L / sqrt(i) for i = 1:L]
        N1 = N2 ./ sqrt(2)
        TranslationParityJudge(f, reverse, K, P, L, base, vcat(N1, N2))
    end

    I, R = if threaded
        selectindexnorm_threaded(judge, L, base=base, alloc=alloc)
    else
        selectindexnorm(judge, L, 1:base^L, base=base, alloc=alloc)
    end

    C = if iszero(K)
        fill(1.0, L)
    else
        [iseven(i) ? 1.0 : -1.0 for i=0:L-1]
    end

    TranslationParityBasis(zeros(Int, L), I, R, C, P, base)
end

TranslationParityBasis(K, P, L; base=2, alloc=1000, threaded=false) = TranslationParityBasis(x -> true, K, P, L, base=base, alloc=alloc, threaded=threaded)

function TranslationFlipBasis(
    f, K::Integer, P::Integer, L::Integer; 
    base::Integer=2, alloc::Integer=1000, threaded::Bool=false
)
    @assert isone(P) || isequal(P, -1) "Invalid parity"
    judge = begin
        N2 = [2L / sqrt(i) for i = 1:L]
        N1 = N2 ./ sqrt(2)
        TranslationParityJudge(f, x -> spinflip(x, base), K, P, L, base, vcat(N1, N2))
    end

    I, R = if threaded
        selectindexnorm_threaded(judge, L, base=base, alloc=alloc)
    else
        selectindexnorm(judge, L, 1:base^L, base=base, alloc=alloc)
    end

    C = if iszero(K)
        fill(1.0, L)
    elseif isequal(2K, L)
        [iseven(i) ? 1.0 : -1.0 for i=0:L-1]
    else
        [exp(-1im * 2π/L * K * i) for i=0:L-1]
    end

    TranslationFlipBasis(zeros(Int64, L), I, R, C, P, base)
end

TranslationFlipBasis(K, P, L; base=2, alloc=1000, threaded=false) = TranslationFlipBasis(x -> true, K, P, L, base=base, alloc=alloc, threaded=threaded)

function parity_schmidt(parity, v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::AbstractTranslationalParityBasis)
    dgt, R, phase = b.dgt, b.R, b.C[2]
    S = SchmidtMatrix(promote_type(eltype(v), eltype(b)), b, Ainds)
    for i = 1:length(v)
        change!(b, i)
        val = v[i] / R[i]
        for j in 1:length(dgt)
            addto!(S, val)
            cyclebits!(dgt)
            val *= phase
            
        end
        dgt .= parity(dgt)
        val *= b.P
        for j in 1:length(dgt)
            addto!(S, val)
            cyclebits!(dgt)
            val *= phase
        end
    end
    Matrix(S)
end

schmidt(v, Ainds, b::TranslationParityBasis) = parity_schmidt(reverse, v, Ainds, b)
schmidt(v, Ainds, b::TranslationFlipBasis) = parity_schmidt(x -> spinflip(x, b.B), v, Ainds, b)

#-----------------------------------------------------------------------------------------------------
# Double Basis
#-----------------------------------------------------------------------------------------------------
export DoubleBasis
"""
    DoubleBasis{Tb1<:AbstractBasis, Tb2<:AbstractBasis}

Basis for constructing transition matrix from one symmetry sector to another.

Properties:
-----------
- `dgt`: Digits.
- `B1` : Basis of the target symmetry sector.
- `B2` : Basis of the starting symmetry sector.
- `B`  : Base.
"""
struct DoubleBasis{Tb1<:AbstractBasis, Tb2<:AbstractBasis} <: AbstractBasis
    dgt::Vector{Int64}
    B1::Tb1
    B2::Tb2
    B::Int64
    DoubleBasis(B1::AbstractBasis, B2::AbstractBasis) = new{typeof(B1), typeof(B2)}(B1.dgt, B1, B2, B2.B)
end

@inline eltype(b::DoubleBasis) = promote_type(eltype(b.B1), eltype(b.B2))
@inline function change!(b::DoubleBasis, i::Integer)
    N = change!(b.B2, i)
    b.dgt .= b.B2.dgt
    N
end
@inline index(b::DoubleBasis) = index(b.B1)
@inline size(b::DoubleBasis) = size(b.B1, 1), size(b.B2, 2)
@inline size(b::DoubleBasis, i::Integer) = isone(i) ? size(b.B1, 1) : isequal(i, 2) ? size(b.B2, 2) : 1
