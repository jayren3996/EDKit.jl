export AbstractBasis, content, base, index, change!
"""
    AbstractBasis

Abstract type for all concrete `Basis` types. 
The most important function for a `Basis` type is the index of a given digits and the digits representation of its ith content:

- `index`  : Return the norm and the index of the digits of the given basis.
- `change!`: Change the digits to the target index, and return the new normalization.

If we want to calculate the schmidt decomposition (in real space), we can use the `schmidt` function.
"""
abstract type AbstractBasis end
abstract type AbstractOnsiteBasis <: AbstractBasis end
abstract type AbstractPermuteBasis <: AbstractBasis end
abstract type AbstractTranslationalParityBasis <: AbstractPermuteBasis end

#-------------------------------------------------------------------------------------------------------------------------
# Default function definitions
#-------------------------------------------------------------------------------------------------------------------------
"""
content(b::AbstractBasis, i::Integer)

Return index of ith vector in this basis.
"""
content(b::AbstractBasis, i::Integer) = b.I[i]
#-------------------------------------------------------------------------------------------------------------------------
"""
norm(b::AbstractBasis, i::Integer)

Return norm of ith vector in this basis.
"""
norm(b::AbstractBasis, i::Integer) = b.R[i]
norm(::AbstractOnsiteBasis, ::Integer) = 1
#-------------------------------------------------------------------------------------------------------------------------
"""
eltype(b::AbstractBasis)

Element data type of the basis. 
- Return `Int64` for onsite basis.
- Return `ComplexF64` by default otherwise.

Used for type promotion when constructing matrix representation of an `Operator`.
"""
eltype(::AbstractBasis) = ComplexF64
eltype(::AbstractOnsiteBasis) = Int64
#-------------------------------------------------------------------------------------------------------------------------
"""
length(b::AbstractBasis) 

Length of the digits. 
"""
length(b::AbstractBasis) = length(b.dgt)
#-------------------------------------------------------------------------------------------------------------------------
"""
size(b::AbstractBasis, [dim])

Dimension of the `Operator` object construct by the basis.
"""
function size(b::AbstractBasis, dim::Integer)
    isone(dim) && return length(b.I)
    isequal(dim, 2) ? length(b.I) : 1
end
size(b::AbstractBasis) = (size(b, 1), size(b, 2))
#-------------------------------------------------------------------------------------------------------------------------
"""
change!(b::AbstractBasis, i::Integer) 

Change the digits to the target index, and return the new normalization.
"""
function change!(b::AbstractBasis, i::Integer) 
    change!(b.dgt, content(b, i), base=b.B)
    norm(b, i)
end
#-------------------------------------------------------------------------------------------------------------------------
"""
int_type(b::AbstractBasis)

Return index data type.
This is for basis with large indices exceeding the upper bound of `Int64`.
"""
int_type(b::AbstractBasis) = eltype(b.I)

#-----------------------------------------------------------------------------------------------------
# TensorBasis
#-----------------------------------------------------------------------------------------------------
export TensorBasis
"""
    TensorBasis

Basis without any symmetries. The index is computed by:

    I(i₁, i₂, ⋯, iₙ) = i₁ B^(n-1) + i₂ B^(n-2) + ⋯ + iₙ

In many bases this basis need not be constructed explicitly.
"""
struct TensorBasis <: AbstractOnsiteBasis
    dgt::Vector{Int64}
    B::Int64
end

function TensorBasis(;L::Integer, base::Integer=2)
    dgt = zeros(Int64, L)
    B = Int64(base)
    TensorBasis(dgt, B)
end
#-------------------------------------------------------------------------------------------------------------------------
content(::TensorBasis, i::Integer) = i
norm(b::TensorBasis, ::Integer) = 1
eltype(::TensorBasis) = Int64
size(b::TensorBasis, i::Integer) = isone(i) || isequal(i, 2) ? b.B^length(b.dgt) : 1
size(b::TensorBasis) = (l=b.B^length(b.dgt); (l, l))
index(b::TensorBasis) = 1, index(b.dgt, base=b.B)
int_type(::TensorBasis) = Int64
#-------------------------------------------------------------------------------------------------------------------------
"""
copy(b::TensorBasis)

Copy basis with a deep copied `dgt`.
"""
function copy(b::TensorBasis) 
    TensorBasis(deepcopy(b.dgt), b.B)
end

#-------------------------------------------------------------------------------------------------------------------------
# Index functions
#-------------------------------------------------------------------------------------------------------------------------
"""
    index(dgt::AbstractVector{T}; base::Integer=2, dtype::DataType=T) where T <: Integer

Convert a digits to the integer index using the relation: 

    number = ∑ᵢ bits[i] * base^(L-i) + 1

The summaton is evaluate using the efficieint polynomial evaluation method.

Inputs:
-------
- `dgt`  : Digits.
- `base` : Base.
- `dtype`: Data type of the index.
"""
@inline function index(dgt::AbstractVector{T}; base::Integer=2, dtype::DataType=T) where T <: Integer
    N = zero(dtype)
    base = convert(dtype, base)
    for i = 1:length(dgt)
        N *= base
        N += dgt[i]
    end
    N + one(dtype)
end

#-------------------------------------------------------------------------------------------------------------------------
"""
    index(dgt::AbstractVector{T}, sites::AbstractVector{<:Integer}; base::Integer=2) where T <: Integer

Convert a sub-digits (subarray of `dgt`) to the integer index.
"""
@inline function index(dgt::AbstractVector{T}, sites::AbstractVector{<:Integer}; base::Integer=2, dtype::DataType=T) where T <: Integer
    N = zero(dtype)
    base = convert(dtype, base)
    for i in sites
        N *= base
        N += dgt[i]
    end
    N + one(dtype)
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    change!(dgt::AbstractVector{<:Integer}, ind::Integer; base::Integer=2) 

Change the digits to that with the target index.  
This method is the inverse of `index`.
"""
@inline function change!(dgt::AbstractVector{<:Integer}, ind::Integer; base::Integer=2)
    T = promote_type(eltype(dgt), typeof(ind))
    N = convert(T, ind) - one(T)
    for i = length(dgt):-1:1
        N, dgt[i] = divrem(N, base)
    end
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    change!(dgt::AbstractVector{<:Integer}, sites::AbstractVector{<:Integer}, ind::Integer; base::Integer=2) 

Change the sub-digits to that with the target index.  
This method is the inverse of `index`.
"""
@inline function change!(dgt::AbstractVector{<:Integer}, sites::AbstractVector{<:Integer}, ind::Integer; base::Integer=2)
    T = promote_type(eltype(dgt), typeof(ind))
    N = convert(T, ind) - one(T)
    for i = length(sites):-1:1
        N, dgt[sites[i]] = divrem(N, base)
    end
end
