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
    dgt::Vector{Ti}
    I::Vector{Ti}
    R::Vector{Float64}
    C::Vector{T}
    A::Int64
    B::Ti
end
#-------------------------------------------------------------------------------------------------------------------------
eltype(::TranslationalBasis{Ti, T}) where {Ti, T} = T
copy(b::TranslationalBasis) = TranslationalBasis(deepcopy(b.dgt), b.I, b.R, b.C, b.A, b.B)
ncycle(b) = length(b.dgt) ÷ b.A

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
struct TranslationJudge{T}
    F                   # Projective selection
    K::Int              # Momentum
    A::Int              # Length of unit cell
    L::Int              # Length of translation
    B::T                # Base
    C::Vector{Float64}  # Normalization coefficient
end
#-------------------------------------------------------------------------------------------------------------------------
@inline check_momentum(n::Integer, k::Integer, L::Integer) = iszero(mod(n * k , L))
"""
(::TranslationJudge)(dgt, i)

Select basis and return normalization.
"""
function (judge::TranslationJudge)(dgt::AbstractVector{<:Integer}, i::Integer)
    # If there is a projective selection function F, check `F(dgt)` first
    isnothing(judge.F) || judge.F(dgt) || return (false, 0.0)
    
    # Check translation 
    for n in 1:judge.L-1
        circshift!(dgt, judge.A)
        In = index(dgt, base=judge.B)
        if In < i 
            # Find smaller indices, return false.
            return (false, 0.0)  
        elseif isequal(In, i)
            # Find periodicity; check if momentum is compatible
            check_momentum(n, judge.K, judge.L) || return (false, 0.0)
            return (true, judge.C[n])
        end
    end

    # Finish a cycle
    true, judge.C[judge.L]
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
function selectindexnorm(f, L::Integer, rg::UnitRange{T}; base::Integer=2, alloc::Integer=1000) where T <: Integer
    dgt = zeros(T, L)
    I, R = T[], Float64[]
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
function selectindexnorm_N(f, L::Integer, N::Integer; base::T=2, alloc::Integer=1000, sorted::Bool=true) where T <: Integer
    I, R = T[], Float64[]
    sizehint!(I, alloc)
    sizehint!(R, alloc)
    for fdgt in multiexponents(L, N)
        all(b < base for b in fdgt) || continue
        dgt = (base-1) .- fdgt
        i = index(dgt, base=base)
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
function selectindexnorm_threaded(f, L::Integer; base::T=2, alloc::Integer=1000) where T <: Integer
    nt = Threads.nthreads()
    ni = dividerange(base^L, nt)
    nI = Vector{Vector{T}}(undef, nt)
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
function TranslationalBasis(dtype::DataType=Int64;
    L::Integer, f=nothing, k::Integer=0, N::Union{Nothing, Integer}=nothing, a::Integer=1,
    base::Integer=2, alloc::Integer=1000, threaded::Bool=true, small_N::Bool=false
)
    len, check_a = divrem(L, a)
    @assert iszero(check_a) "Length of unit-cell $a incompatible with L=$L"
    #=
    Change of definition: 
        old: T|k⟩ = exp(+ik)|k⟩,
        new: T|k⟩ = exp(-ik)|k⟩.
    In the new definition, cₖ⁺|VAC⟩ = |k⟩.
    =#
    k = mod(-k, len) 
    base = convert(dtype, base)
    I, R = begin
        norm = [len/sqrt(i) for i = 1:len]
        judge = if small_N || isnothing(N)
            TranslationJudge(f, k, a, len, base, norm)
        else
            num = L*(base-1)-N
            g = isnothing(f) ? x -> (sum(x) == num) : x -> (sum(x) == num && f(x))
            TranslationJudge(g, k, a, len, base, norm)
        end
        if small_N && !isnothing(N)
            selectindexnorm_N(judge, L, N, base=base)
        else
            threaded ? selectindexnorm_threaded(judge, L, base=base, alloc=alloc) : selectindexnorm(judge, L, 1:base^L, base=base, alloc=alloc)        
        end
    end
    C = phase_factor(k, len)
    TranslationalBasis(zeros(dtype, L), I, R, C, a, base)
end
#-------------------------------------------------------------------------------------------------------------------------
function phase_factor(k::Integer, len::Integer)
    if iszero(k)
        fill(1.0, len)
    elseif isequal(2k, len)
        [iseven(i) ? 1.0 : -1.0 for i=0:len-1] 
    else
        [exp(-1im * 2π/len * k*i) for i=0:len-1]
    end
end

#-------------------------------------------------------------------------------------------------------------------------
# Indexing
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
    I0 = index(b.dgt, base=b.B)
    Im, M = I0, 0

    # Cycle digits
    resetQ = true
    for i in 1:ncycle(b)-1
        circshift!(b.dgt, b.A)
        In = index(b.dgt, base=b.B)
        if isequal(In, I0) 
            resetQ = false
            break
        elseif In < Im
            Im, M = In, i
        end
    end
    resetQ && circshift!(b.dgt, b.A)

    # Search for minimal index 
    i = binary_search(b.I, Im)
    iszero(i) && return (zero(eltype(b)), one(b.B))

    # Find normalization and return result
    N = b.C[M+1] * b.R[i]
    N, i
end
