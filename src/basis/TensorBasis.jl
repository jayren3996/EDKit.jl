#-----------------------------------------------------------------------------------------------------
# Abstract Basis Type
#-----------------------------------------------------------------------------------------------------
abstract type AbstractBasis end
eltype(::AbstractBasis) = ComplexF64
#-----------------------------------------------------------------------------------------------------
change!(b::AbstractBasis, i::Integer) = (change!(b.dgt, b.I[i], base=b.B); b.R[i])
size(b::AbstractBasis, i::Integer) = isone(i) || isequal(i, 2) ? length(b.I) : 1
size(b::AbstractBasis) = (l=length(b.I); (l,l))
length(b::AbstractBasis) = length(b.dgt)
base(b::AbstractBasis) = b.B
#-----------------------------------------------------------------------------------------------------
# Tensor Basis
#-----------------------------------------------------------------------------------------------------
"""
    TensorBasis

Basis without any symmetries.

Properties:
-----------
- dgt : Vector{BITTYPE}, Digits.
- B   : Int, levels on each site.
"""
struct TensorBasis <: AbstractBasis
    dgt::Vector{BITTYPE}
    B::Int
end
eltype(::TensorBasis) = Int
#-----------------------------------------------------------------------------------------------------
index(b::TensorBasis) = 1, index(b.dgt, base=b.B)
change!(b::TensorBasis, i::Integer) = (change!(b.dgt, i, base=b.B); 1)
size(b::TensorBasis, i::Integer) = isone(i) || isequal(i, 2) ? b.B ^ length(b.dgt) : 1
size(b::TensorBasis) = (l=b.B^length(b.dgt); (l, l))
#-----------------------------------------------------------------------------------------------------
# Construction
#-----------------------------------------------------------------------------------------------------
export tensorbasis
tensorbasis(L::Integer; base::Integer=2) = TensorBasis(zeros(Int, L), base)
function schmidt!(target::AbstractMatrix, v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::TensorBasis)
    dgt = b.dgt
    Binds = complement(Ainds, length(dgt))
    for i = 1:length(v)
        change!(dgt, i, base=b.B)
        perm_element!(target, dgt, v[i], Ainds, Binds, b.B)
    end
    target
end
