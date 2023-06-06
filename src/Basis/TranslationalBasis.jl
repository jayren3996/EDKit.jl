export TranslationalBasis
"""
    TranslationalBasis

Basis for subspace that is spanned by momentum states, which can also incorporate projected restriction.
For each basis vector represented by `dgt`, the vector is
    |aₙ⟩ = ∑ᵢ Tⁱ |dgt⟩.
The normalized basis is
    |n⟩ = Rₙ⁻¹ |aₙ⟩,
where Rₙ is the normalization.

Properties:
-----------
- `dgt`: Digits.
- `I`  : List of indicies.
- `R`  : List of normalization.
- `C`  : Unit phase factor.
- `B`  : Base.
"""
struct TranslationalBasis{Ti <: Integer, T <: Number} <: AbstractPermuteBasis
    dgt::Vector{Int64}
    I::Vector{Ti}
    R::Vector{Float64}
    C::Vector{T}
    A::Int64
    B::Int64
end
#-------------------------------------------------------------------------------------------------------------------------
eltype(::TranslationalBasis{Ti, T}) where {Ti, T} = T
copy(b::TranslationalBasis) = TranslationalBasis(deepcopy(b.dgt), b.I, b.R, b.C, b.B)


#-------------------------------------------------------------------------------------------------------------------------
# Functions selecting the basis
#-------------------------------------------------------------------------------------------------------------------------
"""
TranslationJudge

Structure used for selecting the basis contents.

Properties:
-----------
- F: Projective selection
- K: Momentum
- A: Length of unit cell
- B: Base
- C: Normalization coefficient
"""
struct TranslationJudge
    F                   # Projective selection
    K::Int              # Momentum
    A::Int              # Length of unit cell
    B::Int              # Base
    C::Vector{Float64}  # Normalization coefficient
end
#-------------------------------------------------------------------------------------------------------------------------
"""
translation_check(dgt, I0, k, base; a=1)

Check if the digits is of minimum index under cycling.

Outputs:
--------
- `Q`: Whether `dgt` is minimum.
- `R`: Periodicity of the `dgt`, return 0 if !Q.
"""
function translation_check!(dgt::AbstractVector{<:Integer}, I0::Integer, k::Integer, base::Integer; a::Integer=1)
    len = length(dgt)÷a  # Maximum periodicity
    R = len
    for i=1:len-1
        circshift!(dgt, a)
        In = index(dgt, base=base, dtype=typeof(I0))
        In < I0 && return (false, 0)  # Find smaller indices, return false.
        if isequal(In, I0)            # Find periodicity.
            R = i
            break
        end
    end
    # Check if momentum is compatible
    iszero(mod(R * k , len)) ? (true, R) : (false, 0)
end
#-------------------------------------------------------------------------------------------------------------------------
"""
(::TranslationJudge)(dgt, i)

Select basis and return normalization.
"""
function (judge::TranslationJudge)(dgt::AbstractVector{<:Integer}, i::Integer)
    # If there is a projective selection function F, check `F(dgt)` first.
    isnothing(judge.F) || judge.F(dgt) || return false, 0.0
    
    # Going through `translation_check` routine. 
    c, r = translation_check!(dgt, i, judge.K, judge.B, a=judge.A)
    c ? (true, judge.C[r]) : (false, 0.0)
end

#-------------------------------------------------------------------------------------------------------------------------
# Construction
#-------------------------------------------------------------------------------------------------------------------------
"""
    selectindexnorm(f, L, rg; base=2, alloc=1000)

Select the legal indices and corresponding norm for a basis.

Inputs:
-------
- `f`    : Function for digits that tells whether a digits is valid as well as its normalization.
- `L`    : Length of the system.
- `rg`   : Range of interation.
- `base` : (Optional) Base, `base=2` by default.
- `alloc`: (Optional) Pre-allocation memory for the list of indices, `alloc=1000` by default.

Outputs:
--------
- `I`: List of indices in a basis.
- `R`: List of normalization for each states.
"""
function selectindexnorm(
    f, L::Integer, rg::UnitRange; 
    base::Integer=2, alloc::Integer=1000
)
    dgt = zeros(Int64, L)
    I, R = Int[], Float64[]
    sizehint!(I, alloc)
    sizehint!(R, alloc)
    for i in rg
        change!(dgt, i, base=base)
        Q, N = f(dgt, i)
        Q || continue
        append!(I, i)
        append!(R, N)
    end
    I, R
end
#-------------------------------------------------------------------------------------------------------------------------
function selectindexnorm_N(
    f, L::Integer, N::Integer;
    base::Integer=2, dtype::DataType=Int64,
    alloc::Integer=1000, sorted::Bool=true
)
    I, R = dtype[], Float64[]
    sizehint!(I, alloc)
    sizehint!(R, alloc)
    for dgt in multiexponents(L, N)
        all(b < base for b in dgt) || continue
        i = index(dgt, base=base, dtype=dtype)
        Q, N = f(dgt, i)
        Q || continue
        append!(I, i)
        append!(R, N)
    end
    sorted || return I, R
    sperm = sortperm(I)
    I[sperm], R[sperm]
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    selectindexnorm_threaded(f, L, rg; base=2, alloc=1000)

Select the legal indices and corresponding norm for a basis with multi-threads.

Inputs:
-------
- `f`    : Function for digits that tells whether a digits is valid as well as its normalization.
- `L`    : Length of the system.
- `base` : (Optional) Base, `base=2` by default.
- `alloc`: (Optional) Pre-allocation memory for the list of indices, `alloc=1000` by default.

Outputs:
--------
- `I`: List of indices in a basis.
- `R`: List of normalization for each states.
"""
function selectindexnorm_threaded(
    f, L::Integer; 
    base::Integer=2, alloc::Integer=1000
)
    nt = Threads.nthreads()
    ni = dividerange(base^L, nt)
    nI = Vector{Vector{Int}}(undef, nt)
    nR = Vector{Vector{Float64}}(undef, nt)
    Threads.@threads for ti in 1:nt
        nI[ti], nR[ti] = selectindexnorm(f, L, ni[ti], base=base, alloc=alloc)
    end
    I, R = vcat(nI...), vcat(nR...)
    I, R
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    TranslationalBasis(f, k, L; base=2, alloc=1000, threaded=true)

Construction for `TranslationalBasis`.

Inputs:
-------
- `dtype`   : Data type of index.
- `L`       : Length of the system.
- `f`       : Selection function for the basis contents.
- `k`       : Momentum number from 0 to L-1.
- `N`       : Particle number.
- `a`       : Length of unit cell.
- `base`    : Base, default = 2.
- `alloc`   : Size of the prealloc memory for the basis content, used only in multithreading, default = 1000.
- `threaded`: Whether use the multithreading, default = true.

Outputs:
--------
- `b`: TranslationalBasis.
"""
function TranslationalBasis(
    dtype::DataType=Int64
    ;L::Integer, f=nothing, k::Integer=0, N::Union{Nothing, Integer}=nothing, a::Integer=1,
    base::Integer=2, alloc::Integer=1000, threaded::Bool=true, small_N::Bool=true
)
    len, check_a = divrem(L, a)
    @assert iszero(check_a) "Length of unit-cell $a incompatible with L=$L"
    k = mod(k, len)
    I, R = begin
        norm = [len/sqrt(i) for i = 1:len]
        judge = if small_N || isnothing(N)
            TranslationJudge(f, k, a, base, norm)
        else
            num = L*(base-1)-N
            g = isnothing(f) ? x -> (sum(x) == num) : x -> (sum(x) == num && f(x))
            TranslationJudge(g, k, a, base, norm)
        end
        if small_N && !isnothing(N)
            selectindexnorm_N(judge, L, N, base=base, dtype=dtype)
        else
            threaded ? selectindexnorm_threaded(judge, L, base=base, alloc=alloc) : selectindexnorm(judge, L, 1:base^L, base=base, alloc=alloc)        
        end
    end
    C = if iszero(k)
        fill(1.0, len)
    elseif isequal(2k, len)
        [iseven(i) ? 1.0 : -1.0 for i=0:len-1] 
    else
        [exp(-1im * 2π/len * k*i) for i=0:len-1]
    end
    TranslationalBasis(zeros(Int64, L), I, R, C, a, base)
end

#-------------------------------------------------------------------------------------------------------------------------
# Indexing
#-------------------------------------------------------------------------------------------------------------------------
"""
translation_index(dgt::AbstractVector{<:Integer}, base::Integer)

Given a digits, return the minimum index under shifting.

Inputs:
-------
- `dgt   : Digits.
- `base` : Base.
- `dtype`: Specify data type of index.
- `a`    : Length of unit cell

Outputs:
--------
- `Im`: Minimum index.
- `M` : Cycling number by which `dgt` is shifted to have index `Im`.
"""
function translation_index(
    dgt::AbstractVector{<:Integer}, base::Integer; 
    dtype::DataType=Int64, a::Integer=1
)
    I0 = index(dgt, base=base, dtype=dtype)
    Im, M = I0, 0
    len = length(dgt)÷a
    for i=1:len-1
        circshift!(dgt, a)
        In = index(dgt, base=base, dtype=dtype)
        isequal(In, I0) && return Im, M  # Reach one period, return result immediately.
        if In < Im                       # Find smaller index.
            Im, M = In, i
        end
    end
    circshift!(dgt, a)  # Reset digit.
    Im, M
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    index(b::TranslationalBasis)

Index of the state in `TranslationalBasis`, together with the normalization.

Outputs:
--------
- `N`: Normalization for the state `b.dgt`.
- `i`: Index for the state `b.dgt`.

Notes:
------
To calculate the index of a given digits, we first shift the digits to that of the minimum index, then search for the place in the basis.
Because of the restriction from the momentum, some state with zero normalization (which is not included in the basis) will appear.
To avoid exception in the matrix constructon of `Operation`, we allow the index to not in the basis content.
When this happend, we return index 1, and normalization 0, so it has no effect on the matrix being filled.
"""
function index(b::TranslationalBasis)
    Im, T = translation_index(b.dgt, b.B, dtype=int_type(b), a=b.A)

    # Search for minimal index 
    i = binary_search(b.I, Im)
    iszero(i) && return zero(eltype(b)), 1

    # Find normalization and return result
    N = b.C[T+1] * b.R[i]
    N, i
end
