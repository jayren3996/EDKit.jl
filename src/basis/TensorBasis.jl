#-----------------------------------------------------------------------------------------------------
# Tensor Basis
#-----------------------------------------------------------------------------------------------------
"""
    TensorBasis

Basis without any symmetries.

Properties:
-----------
- dgt : Vector{Int}, Digits.
- B   : Int, levels on each site.
"""
struct TensorBasis <: AbstractBasis
    dgt::Vector{Int}
    B::Int
end
eltype(::TensorBasis) = Int
#-----------------------------------------------------------------------------------------------------
index(b::TensorBasis) = 1, index(b.dgt, base=b.B)
change!(b::TensorBasis, i::Integer) = (change!(b.dgt, i, base=b.B); 1)
size(b::TensorBasis, i::Integer) = (i == 2 || i == 1) ? b.B ^ length(b.dgt) : 1
size(b::TensorBasis) = (l=b.B^length(b.dgt); (l, l))
#-----------------------------------------------------------------------------------------------------
# Construction
#-----------------------------------------------------------------------------------------------------
export tensorbasis
tensorbasis(L::Integer; base::Integer=2) = TensorBasis(zeros(Int, L), base)
function schmidt!(target::AbstractMatrix, v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::TensorBasis)
    Binds = Int[i for i = 1:length(b.dgt) if !in(i, Ainds)]
    dgt = b.dgt
    for i = 1:length(v)
        change!(dgt, i, base=b.B)
        ia = index(dgt, Ainds, base=b.B)
        ib = index(dgt, Binds, base=b.B)
        target[ia, ib] = v[i]
    end
    target
end
