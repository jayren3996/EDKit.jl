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
#-----------------------------------------------------------------------------------------------------
eltype(::TensorBasis) = Int
change!(b::TensorBasis, i::Integer) = (change!(b.dgt, i, base=b.B); 1)
index(b::TensorBasis) = 1, index(b.dgt, base=b.B)
size(b::TensorBasis, i::Integer) = (i == 2 || i == 1) ? b.B ^ length(b.dgt) : 1
size(b::TensorBasis) = (l=b.B^length(b.dgt); (l, l))
#-----------------------------------------------------------------------------------------------------
# Construction
#-----------------------------------------------------------------------------------------------------
export tensorbasis
tensorbasis(L::Integer; base::Integer=2) = TensorBasis(zeros(Int, L), base)

#-----------------------------------------------------------------------------------------------------
# Helper Functions
#-----------------------------------------------------------------------------------------------------
"""
    index(dgt::AbstractVector{<:Integer})

Digit ⟹ index. 

The method compute the polynomial

number = bits[i] * base^(L-i) + 1

in the most efficient way.
"""
function index(dgt::AbstractVector{Int}; base::Int=2)::Int
    N = 0
    for i = 1:length(dgt)
        N = muladd(base, N, dgt[i])
    end
    N + 1
end
#-----------------------------------------------------------------------------------------------------
function index(dgt::AbstractVector{Int}, sites::Vector{Int}; base::Int=2)::Int
    N = 0
    for i in sites
        N = muladd(base, N, dgt[i])
    end
    N + 1
end
#-----------------------------------------------------------------------------------------------------
"""
    change!(dgt::Digit, index::Integer)

Index ⟹ Digit. 

The method compute the bits vector and write to bits.
"""
function change!(dgt::Vector{Int}, ind::Int; base::Int=2)
    N = ind - 1
    for i = length(dgt):-1:1
        N, dgt[i] = divrem(N, base)
    end
end
#-----------------------------------------------------------------------------------------------------
function change!(dgt::Vector{Int}, sites::Vector{Int}, ind::Int; base::Int=2)
    N = ind - 1
    for i = length(sites):-1:1
        N, dgt[sites[i]] = divrem(N, base)
    end
end
