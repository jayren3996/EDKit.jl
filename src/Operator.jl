#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Operator
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
    Operator{Tv}

Object that store the local operator.

Properties:
-----------
- M : Vector{SparseMatrixCSC{Tv, Int}}, local operators represented by CSC sparse matrices.
- I : Vector{Vector{Int}}, indices of sites that the operators act on.
- B : Basis.
"""
struct Operator{Tv<:Number, Tb<:AbstractBasis}
    M::Vector{SparseMatrixCSC{Tv, Int}}
    I::Vector{Vector{Int}}
    B::Tb
end
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
eltype(opt::Operator{Tv, Tb}) where Tv where Tb = promote_type(Tv, eltype(opt.B))
length(opt::Operator) = length(opt.M)
size(opt::Operator) = size(opt.B)
size(opt::Operator, i::Integer) = size(opt.B, i)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Construction
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function operator(mats::AbstractVector{<:AbstractMatrix}, inds::AbstractVector{<:AbstractVector}, B::AbstractBasis)
    M, I = standard_format(mats, inds)
    Operator(M, I, B)
end
function operator(mats::AbstractVector{<:AbstractMatrix}, inds::AbstractVector{<:AbstractVector}, L::Integer)
    M, I = standard_format(mats, inds)
    B = tensorbasis(L, base=find_base(size(mats[1], 1), length(inds[1])))
    Operator(M, I, B)
end
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
operator(mats::AbstractVector{<:AbstractMatrix}, inds::AbstractVector{<:Integer}, C) = operator(mats, [[i] for i in inds], C)
operator(mat::AbstractMatrix, ind::AbstractVector{<:Integer}, C) = operator([mat], [ind], C)
operator(mats::AbstractVector{<:AbstractMatrix}, L::Integer) = operator(mats, [[i] for i in 1:L], L)
operator(mats::AbstractVector{<:AbstractMatrix}, B::AbstractBasis) = operator(mats, [[i] for i in 1:length(B.dgt)], B)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
export trans_inv_operator
function trans_inv_operator(mat::AbstractMatrix, ind::AbstractVector{<:Integer}, B::AbstractBasis)
    L = length(B.dgt)
    smat = sparse(mat)
    mats = fill(smat, L)
    inds = [mod.(ind .+ i, L) .+ 1 for i = -1:L-2]
    Operator(mats, inds, B)
end
function trans_inv_operator(mat::AbstractMatrix, ind::AbstractVector{<:Integer}, L::Integer)
    B = tensorbasis(L, base=find_base(size(mat, 1), length(ind)))
    trans_inv_operator(mat, ind, B)
end
trans_inv_operator(mat::AbstractMatrix, M::Integer, C) = trans_inv_operator(mat, 1:M, C)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Basic Math
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
*(c::Number, opt::Operator) = Operator(c .* opt.M, opt.I, opt.B)
*(opt::Operator, c::Number) = c * opt
/(opt::Operator, c::Number) = Operator(opt.M ./ c, opt.I, opt.B)
function +(opt1::Operator, opt2::Operator)
    n1, n2 = length(opt1), length(opt2)
    Tv = promote_type(eltype(opt1), eltype(opt2))
    M = Vector{SparseMatrixCSC{Tv, Int}}(undef, n1 + n2)
    I = Vector{Vector{Int}}(undef, n1+n2)
    for i = 1:n1
        M[i], I[i] = opt1.M[i], opt1.I[i]
    end
    P = n1
    for i = 1:n2
        if (p = findposition(opt1.I, opt2.I[i])) > 0
            M[p] += opt2.M[i]
        else
            P += 1
            M[n1+i], I[n1+i] = opt2.M[i], opt2.I[i]
        end
    end
    deleteat!(M, P+1:n1+n2)
    deleteat!(I, P+1:n1+n2)
    Operator(M, I, opt1.B)
end
-(opt::Operator) = Operator(-opt.M, opt.I, opt.B)
-(opt1::Operator, opt2::Operator) = opt1 + (-opt2)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Helper Functions
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
standard_mat(Tv::DataType, mat::AbstractMatrix) = Tv.(mat) |> sparse
standard_mat(mats::Vector{SparseMatrixCSC{Tv, Int}}) where Tv = mats
standard_mat(mats::AbstractVector{<:AbstractMatrix}) = standard_mat.(promote_type(eltype.(mats)...), mats)
standard_ind(ind::Vector{Int}) = ind
standard_ind(ind::AbstractVector{<:Integer}) = Int.(ind)
standard_ind(inds::Vector{Vector{Int}}) = inds
standard_ind(inds::AbstractVector{<:AbstractVector{<:Integer}}) = standard_ind.(inds)
standard_format(mats, inds) = standard_mat(mats), standard_ind(inds)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function find_base(a::Integer, b::Integer)::Int
    if b == 1
        return a
    else
        for i=2:a
            if i^b == a
                return i
            end
        end
        error("Incompatible dimension: ($a, $b)")
    end
end
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function findposition(l::Vector, b)::Int
    for i in length(l)
        if l[i] == b
            return i
        end
    end
    0
end
