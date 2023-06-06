export ProjectedBasis
"""
    ProjectedBasis{Ti}

Basis for subspace that is spanned only by product states.
It is basically like the `TensorBasis`, but contains only a subset of all basis vectors, selected by a given function.

Properties:
-----------
- `dgt`: Digits.
- `I`  : List of indicies.
- `B`  : Base.
"""
struct ProjectedBasis{Ti <: Integer} <: AbstractOnsiteBasis
    dgt::Vector{Int64}
    I::Vector{Ti}
    B::Int64
end
#-------------------------------------------------------------------------------------------------------------------------
"""
Deep copy digits.
"""
copy(b::ProjectedBasis) = ProjectedBasis(deepcopy(b.dgt), b.I, b.B)
#-------------------------------------------------------------------------------------------------------------------------
function index(b::ProjectedBasis; check::Bool=true)
    i = index(b.dgt, base=b.B, dtype=int_type(b))
    ind = binary_search(b.I, i)
    ind > 0 && return 1, ind
    check ? error("No such symmetry.") : return zero(eltype(b)), 1
end


#-------------------------------------------------------------------------------------------------------------------------
# Select basis vectors
#-------------------------------------------------------------------------------------------------------------------------
"""
binary_search(list::AbstractVector{<:Integer}, i<:Integer)

Return the position of i in a sorted list using binary search algorithm.
"""
function binary_search(list::AbstractVector{<:Integer}, i::Integer)
    l::Int = 1
    r::Int = length(list)
    c::Int = (l + r) ÷ 2
    while true
        t = list[c]
        (i < t) ? (r = c - 1) : (i > t) ? (l = c + 1) : break
        (l > r) ? (c = 0; break) : (c = (l + r) ÷ 2)
    end
    c
end



#-------------------------------------------------------------------------------------------------------------------------
# Construction
#-------------------------------------------------------------------------------------------------------------------------
"""
dividerange(maxnum::Integer, nthreads::Integer)

Divide the interation range equally according to the number of threads.

Inputs:
-------
- `maxnum`  : Maximum number of range.
- `nthreads`: Number of threads.

Outputs:
--------
- `list`: List of ranges.
"""
function dividerange(maxnum::Integer, nthreads::Integer)
    list = Vector{UnitRange{Int64}}(undef, nthreads)
    eachthreads, left = divrem(maxnum, nthreads)
    start = 1
    for i = 1:left
        stop  = start + eachthreads
        list[i] = start:stop
        start = stop+1
    end
    for i = left+1:nthreads-1
        stop  = start + eachthreads - 1
        list[i] = start:stop
        start = stop+1
    end
    list[nthreads] = start:maxnum
    list
end
#-------------------------------------------------------------------------------------------------------------------------
"""
selectindex(f, L, rg; base=2, alloc=1000)

Select the legal indices for a basis.

Inputs:
-------
- `f`    : Function for digits that tells whether a digits is valid.
- `L`    : Length of the system.
- `rg`   : Range of interation.
- `base` : (Optional) Base, `base=2` by default.
- `alloc`: (Optional) Pre-allocation memory for the list of indices, `alloc=1000` by default.

Outputs:
--------
- `I`: List of indices in a basis.
"""
function selectindex(
    f, L::Integer, rg::UnitRange;
    base::Integer=2, alloc::Integer=1000
)
    dgt = zeros(Int64, L)
    I = Int[]
    sizehint!(I, alloc)
    for i in rg
        change!(dgt, i, base=base)
        f(dgt) && append!(I, i)
    end
    I
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    selectindex_threaded(f, L, rg; base=2, alloc=1000)

Select the legal indices with multi-threads.

Inputs:
-------
- `f`    : Function for digits that tells whether a digits is valid as well as its normalization.
- `L`    : Length of the system.
- `base` : (Optional) Base, `base=2` by default.
- `alloc`: (Optional) Pre-allocation memory for the list of indices, `alloc=1000` by default.

Outputs:
--------
- `I`: List of indices in a basis.
"""
function selectindex_threaded(
    f, L::Integer; 
    base::Integer=2, alloc::Integer=1000
)
    nt = Threads.nthreads()
    ni = dividerange(base^L, nt)
    nI = Vector{Vector{Int}}(undef, nt)
    Threads.@threads for ti in 1:nt
        nI[ti] = selectindex(f, L, ni[ti], base=base, alloc=alloc)
    end
    I = vcat(nI...)
    I
end
#-------------------------------------------------------------------------------------------------------------------------
function selectindex_N(
    f, L::Integer, N::Integer;
    base::Integer=2, dtype::DataType=Int64, 
    alloc::Integer=1000, sorted::Bool=true
)
    I = dtype[]
    sizehint!(I, alloc)
    for dgt in multiexponents(L, N)
        isnothing(f) || f(dgt) || continue
        all(b < base for b in dgt) || continue
        ind = index(dgt, base=base, dtype=dtype)
        append!(I, ind)
    end
    sorted ? sort(I) : I
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    ProjectedBasis(;L, N, base=2, alloc=1000, threaded=true)

Construction method for `ProjectedBasis` with fixed particle number (or total U(1) charge).

Inputs:
-------
- `dtype`   : Data type for indices.
- `L`       : Length of the system.
- `N`       : Quantum number of up spins / particle number / U(1) charge.
- `base`    : Base, default = 2.
- `alloc`   : Size of the prealloc memory for the basis content, used only in multithreading, default = 1000.
- `threaded`: Whether use the multithreading, default = true.
- `small_N` : Whether use the small-N algorithm.

Outputs:
--------
- `b` : ProjectedBasis.
"""
function ProjectedBasis(
    dtype::DataType=Int64
    ;L::Integer, f=nothing, N::Union{Nothing, Integer}=nothing,
    base::Integer=2, alloc::Integer=1000, 
    threaded::Bool=true, small_N::Bool=true
)
    I = if isnothing(N)
        threaded ? selectindex_threaded(f, L, base=base, alloc=alloc) : selectindex(f, L, 1:base^L, base=base, alloc=alloc)
    elseif small_N
        selectindex_N(f, L, N, base=base, dtype=dtype)
    else
        num = L*(base-1)-N
        g = isnothing(f) ? x -> sum(x) == num : x -> (sum(x) == num && f(x))
        threaded ? selectindex_threaded(g, L, base=base, alloc=alloc) : selectindex(g, L, 1:base^L, base=base, alloc=alloc)
    end
    ProjectedBasis(zeros(Int64, L), I, base)
end

