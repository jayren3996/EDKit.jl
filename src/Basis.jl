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

The most inportant function for a `Basis` type is the index of a given digits and the digits representation of its ith content:

- `index`  : Return the norm and the index of the digits of the given basis.
- `change!`: Change the digits to the target index, and return the new normalization.

If we want to calculate the schmidt decomposition (in real space), we can use the `schmidt` function.
"""
abstract type AbstractBasis end

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

TensorBasis(;L::Integer, base::Integer=2) = TensorBasis(L, base=base)

@inline content(b::TensorBasis) = 1:b.B^length(b.dgt)
@inline content(::TensorBasis, i::Integer) = i
@inline norm(b::TensorBasis) = ones(Int64, B^length(b.dgt))
@inline norm(b::TensorBasis, ::Integer) = 1
@inline eltype(::TensorBasis) = Int64
@inline size(b::TensorBasis, i::Integer) = isone(i) || isequal(i, 2) ? b.B^length(b.dgt) : 1
@inline size(b::TensorBasis) = (l=b.B^length(b.dgt); (l, l))
@inline index(b::TensorBasis) = 1, index(b.dgt, base=b.B)

"""
    schmidt(v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::TensorBasis)

Schmidt decomposition of state `v`, with respect to given lattice bipartition.

Inputs:
-------
- `v`    : State represented by a (abstract) vector. 
- `Ainds`: List of indices in subsystem `A`, the remaining indices are regarded as subsystem `B`.
- `b`    : Basis.

Outputs:
- `S`: Matrix S in the decomposition: |v⟩ = Sᵢⱼ |Aᵢ⟩|Bⱼ⟩.
"""
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
It is basically like the `TensorBasis`, but contains only a subset of all basis vectors, selected by a given function.

Properties:
-----------
- `dgt`: Digits.
- `I`  : List of indicies.
- `B`  : Base.
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

"""
    ProjectedBasis(f, L; base=2, alloc=1000, threaded=false)

Construction method for `ProjectedBasis`.

Inputs:
-------
- `f`       : Selection function for the basis contents.
- `L`       : Length of the system.
- `base`    : Base, default = 2.
- `alloc`   : Size of the prealloc memory for the basis content, used only in multithreading, default = 1000.
- `threaded`: Whether use the multithreading, default = true.

Outputs:
--------
- `b`: ProjectedBasis.
"""
function ProjectedBasis(
    f, L::Integer; 
    base::Integer=2, 
    alloc::Integer=1000, 
    threaded::Bool=true
)
    dgt = zeros(Int64, L)
    I = threaded ? selectindex_threaded(f, L, base=base, alloc=alloc) : selectindex(f, L, 1:base^L, base=base, alloc=alloc)
    ProjectedBasis(dgt, I, base)
end

function ProjectedBasis(
    f; L::Integer, base::Integer=2, 
    alloc::Integer=1000, threaded::Bool=true
)
    ProjectedBasis(f, L, base=base, alloc=alloc, threaded=threaded)
end

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
    ;L::Integer, N::Integer, base::Integer=2, 
    alloc::Integer=1000, threaded::Bool=true
)
    num = L*(base-1)-N
    ProjectedBasis(x->sum(x)==num, L, base=base, alloc=alloc, threaded=threaded)
end

"""
    schmidt(v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::ProjectedBasis)

Schmidt decomposition for `ProjectedBasis`.
"""
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

"""
    TranslationalBasis(f, k, L; base=2, alloc=1000, threaded=true)

Construction for `TranslationalBasis`.

Inputs:
-------
- `f`       : Selection function for the basis contents.
- `k`       : Momentum number from 0 to L-1.
- `L`       : Length of the system.
- `base`    : Base, default = 2.
- `alloc`   : Size of the prealloc memory for the basis content, used only in multithreading, default = 1000.
- `threaded`: Whether use the multithreading, default = true.

Outputs:
--------
- `b`: TranslationalBasis.
"""
function TranslationalBasis(
    f, k::Integer, L::Integer; 
    base::Integer=2, alloc::Integer=1000, threaded::Bool=true
)
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

function TranslationalBasis(
    k::Integer, L::Integer; 
    base::Integer=2, alloc::Integer=1000, threaded::Bool=true
)
    TranslationalBasis(x -> true, k, L, base=base, alloc=alloc, threaded=threaded)
end

function TranslationalBasis(
    ;L::Integer, f=x->true, k::Integer=0, N::Union{Nothing, Integer}=nothing,
    base::Integer=2, alloc::Integer=1000, threaded::Bool=true
)
    g = if !isnothing(N)
        num = L*(base-1)-N
        x -> (sum(x) == num && f(x))
    else 
        f
    end
    TranslationalBasis(g, k, L, base=base, alloc=alloc, threaded=threaded)
end

"""
Schmidt decomposition for `TranslationalBasis`.
"""
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
    f, K::Integer, P::Integer, L::Integer; 
    base::Integer=2, alloc::Integer=1000, threaded::Bool=true
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

function TranslationParityBasis(
    K::Integer, P::Integer, L::Integer; 
    base::Integer=2, alloc::Integer=1000, threaded::Bool=true
) 
    TranslationParityBasis(x -> true, K, P, L, base=base, alloc=alloc, threaded=threaded)
end

function TranslationParityBasis(
    ;L::Integer, f=x->true, k::Integer=0, p::Integer=1, N::Union{Nothing, Integer}=nothing,
    base::Integer=2, alloc::Integer=1000, threaded::Bool=true
) 
    g = if !isnothing(N)
        num = L*(base-1)-N
        x -> (sum(x) == num && f(x))
    else 
        f
    end
    TranslationParityBasis(g, k, p, L, base=base, alloc=alloc, threaded=threaded)
end



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

function TranslationFlipBasis(
    K::Integer, P::Integer, L::Integer; 
    base::Integer=2, alloc::Integer=1000, threaded::Bool=false
) 
    TranslationFlipBasis(x -> true, K, P, L, base=base, alloc=alloc, threaded=threaded)
end

function TranslationFlipBasis(
    ;f=x->true, k::Integer=0, p::Integer=1, L::Integer,
    base::Integer=2, alloc::Integer=1000, threaded::Bool=false
) 
    TranslationFlipBasis(f, k, p, L, base=base, alloc=alloc, threaded=threaded)
end

"""
Helper function for `schmidt` on parity basis.
"""
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

#-----------------------------------------------------------------------------------------------------
# State
#-----------------------------------------------------------------------------------------------------
export productstate
"""
    productstate(v::AbstractVector{<:Integer}, B::AbstractBasis)

Construction for product state.

Inputs:
-------
- `v`: Vector representing the product state.
- `B`: Basis.

Outputs:
- `s`: Vector representing the many-body state.
"""
function productstate(v::AbstractVector{<:Integer}, B::AbstractBasis)
    s = zeros(size(B, 1))
    B.dgt .= v 
    I = index(B)[2]
    s[I] = 1 
    s
end
