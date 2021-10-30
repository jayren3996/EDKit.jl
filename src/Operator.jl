"""
Operator

Construction of `Operator` object.
"""

export operator, trans_inv_operator
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

@inline eltype(opt::Operator{Tv, Tb}) where Tv where Tb = promote_type(Tv, eltype(opt.B))
@inline length(opt::Operator) = length(opt.M)
@inline size(opt::Operator) = size(opt.B)
@inline size(opt::Operator, i::Integer) = size(opt.B, i)

function operator(mats::AbstractVector{<:AbstractMatrix}, inds::AbstractVector{<:AbstractVector}, B::AbstractBasis)
    M, I = standard_format(mats, inds)
    Operator(M, I, B)
end

function operator(mats::AbstractVector{<:AbstractMatrix}, inds::AbstractVector{<:AbstractVector}, L::Integer)
    M, I = standard_format(mats, inds)
    B = TensorBasis(L, base=find_base(size(mats[1], 1), length(inds[1])))
    Operator(M, I, B)
end

operator(mats::AbstractVector{<:AbstractMatrix}, inds::AbstractVector{<:Integer}, C) = operator(mats, [[i] for i in inds], C)
operator(mat::AbstractMatrix, ind::AbstractVector{<:Integer}, C) = operator([mat], [ind], C)
operator(mats::AbstractVector{<:AbstractMatrix}, L::Integer) = operator(mats, [[i] for i in 1:L], L)
operator(mats::AbstractVector{<:AbstractMatrix}, B::AbstractBasis) = operator(mats, [[i] for i in 1:length(B.dgt)], B)

function trans_inv_operator(mat::AbstractMatrix, ind::AbstractVector{<:Integer}, B::AbstractBasis)
    L = length(B.dgt)
    smat = sparse(mat)
    mats = fill(smat, L)
    inds = [mod.(ind .+ i, L) .+ 1 for i = -1:L-2]
    Operator(mats, inds, B)
end

function trans_inv_operator(mat::AbstractMatrix, ind::AbstractVector{<:Integer}, L::Integer)
    B = TensorBasis(L, base=find_base(size(mat, 1), length(ind)))
    trans_inv_operator(mat, ind, B)
end

trans_inv_operator(mat::AbstractMatrix, M::Integer, C) = trans_inv_operator(mat, 1:M, C)

*(c::Number, opt::Operator) = Operator(c .* opt.M, opt.I, opt.B)
*(opt::Operator, c::Number) = c * opt
/(opt::Operator, c::Number) = Operator(opt.M ./ c, opt.I, opt.B)

function +(opt1::Operator, opt2::Operator)
    n1, n2 = length(opt1), length(opt2)
    Tv = promote_type(eltype(opt1), eltype(opt2))
    M = Vector{SparseMatrixCSC{Tv, Int}}(undef, n1 + n2)
    I = Vector{Vector{Int}}(undef, n1+n2)
    M[1:n1] .= opt1.M
    I[1:n1] .= opt1.I
    P = n1
    for i = 1:n2
        if (p = findposition(opt1.I, opt2.I[i])) > 0
            M[p] += opt2.M[i]
        else
            P += 1
            M[P], I[P] = opt2.M[i], opt2.I[i]
        end
    end
    deleteat!(M, P+1:n1+n2)
    deleteat!(I, P+1:n1+n2)
    Operator(M, I, opt1.B)
end

-(opt::Operator) = Operator(-opt.M, opt.I, opt.B)
-(opt1::Operator, opt2::Operator) = opt1 + (-opt2)

# Helper Functions
@inline standard_mat(Tv::DataType, mat::AbstractMatrix) = Tv.(mat) |> sparse
@inline standard_mat(mats::Vector{SparseMatrixCSC{Tv, Int}}) where Tv = mats
@inline standard_mat(mats::AbstractVector{<:AbstractMatrix}) = standard_mat.(promote_type(eltype.(mats)...), mats)
@inline standard_ind(ind::Vector{Int}) = ind
@inline standard_ind(ind::AbstractVector{<:Integer}) = Int.(ind)
@inline standard_ind(inds::Vector{Vector{Int}}) = inds
@inline standard_ind(inds::AbstractVector{<:AbstractVector{<:Integer}}) = standard_ind.(inds)
@inline standard_format(mats, inds) = standard_mat(mats), standard_ind(inds)

@inline function find_base(a::Integer, b::Integer)
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

@inline function findposition(l::AbstractVector, b)
    for i in length(l)
        if l[i] == b
            return i
        end
    end
    0
end

# Operator to matrices
export addto!
function addto!(M::AbstractMatrix, opt::Operator)
    for j = 1:size(opt.B, 2)
        colmn!(view(M, :, j), opt, j)
    end
    M
end

function Array(opt::Operator)
    M = zeros(eltype(opt), size(opt))
    if size(M, 1) > 0 && size(M, 2) > 0
        addto!(M, opt)
    end
    M
end

function mul!(target::AbstractVector, opt::Operator, v::AbstractVector)
    for j = 1:length(v)
        colmn!(target, opt, j, v[j])
    end
    target
end

function mul!(target::AbstractMatrix, opt::Operator, m::AbstractMatrix)
    for j = 1:size(m, 1)
        colmn!(target, opt, j, view(m, j, :))
    end
    target
end

function *(opt::Operator, v::AbstractVector)
    ctype = promote_type(eltype(opt), eltype(v))
    M = zeros(ctype, size(opt, 1))
    mul!(M, opt, v)
end

function *(opt::Operator, m::AbstractMatrix)
    ctype = promote_type(eltype(opt), eltype(m))
    M = zeros(ctype, size(opt, 1), size(m, 2))
    mul!(M, opt, m)
end

# Helper functions
function colmn!(target::AbstractVector, M::SparseMatrixCSC, I::Vector{Int}, b::AbstractBasis, coeff::Number=1)
    rows, vals = rowvals(M), nonzeros(M)
    j = index(b.dgt, I, base=b.B)
    change = false
    for i in nzrange(M, j)
        row, val = rows[i], vals[i]
        change!(b.dgt, I, row, base=b.B)
        C, pos = index(b)
        target[pos] += coeff * C * val
        change = true
    end
    change ? change!(b.dgt, I, j, base=b.B) : nothing
end

function colmn!(target::AbstractVector, opt::Operator, j::Integer, coeff::Number=1)
    b, M, I = opt.B, opt.M, opt.I
    C = coeff / change!(b, j)
    for i = 1:length(M)
        colmn!(target, M[i], I[i], b, C)
    end
end

function colmn!(target::AbstractMatrix, M::SparseMatrixCSC, I::Vector{Int}, b::AbstractBasis, coeff::AbstractVector{<:Number})
    rows, vals = rowvals(M), nonzeros(M)
    j = index(b.dgt, I, base=b.B)
    change = false
    for i in nzrange(M, j)
        row, val = rows[i], vals[i]
        change!(b.dgt, I, row, base=b.B)
        C, pos = index(b)
        target[pos, :] .+= (C * val) .* coeff
        change = true
    end
    change ? change!(b.dgt, I, j, base=b.B) : nothing
end

function colmn!(target::AbstractMatrix, opt::Operator, j::Integer, coeff::AbstractVector{<:Number})
    b, M, I = opt.B, opt.M, opt.I
    N = change!(b, j) .* coeff
    for i = 1:length(M)
        colmn!(target, M[i], I[i], b, N)
    end
end


# Spin Matrices
export spin

@inline spin_coeff(D::Integer) = [sqrt(i*(D-i)) for i = 1:D-1]
@inline spin_Sp(D::Integer) = sparse(1:D-1, 2:D, spin_coeff(D), D, D)
@inline spin_Sm(D::Integer) = sparse(2:D, 1:D-1, spin_coeff(D), D, D)
@inline function spin_Sx(D::Integer)
    coeff = spin_coeff(D) / 2
    sp = sparse(1:D-1, 2:D, coeff, D, D)
    sp + sp'
end
@inline function spin_iSy(D::Integer)
    coeff = spin_coeff(D) / 2
    sp = sparse(1:D-1, 2:D, coeff, D, D)
    sp - sp'
end
@inline function spin_Sz(D::Integer)
    J = (D-1) / 2
    sparse(1:D, 1:D, J:-1:-J)
end

function spin_dict(c::Char, D::Integer)
    if     isequal(c, '+') spin_Sp(D)
    elseif isequal(c, '-') spin_Sm(D)
    elseif isequal(c, 'x') spin_Sx(D)
    elseif isequal(c, 'y') spin_iSy(D)
    elseif isequal(c, 'z') spin_Sz(D)
    elseif isequal(c, '1') spdiagm(ones(D))
    elseif isequal(D, 2)
        if     isequal(c, 'X') sparse([0 1; 1  0])
        elseif isequal(c, 'Y') sparse([0 1;-1  0])
        elseif isequal(c, 'Z') sparse([1 0; 0 -1])
        else error("Invalid spin symbol: $c.")
        end
    else error("Invalid spin symbol: $c.")
    end
end

"""
    spin(s::String; D::Integer=2)

Return matrix for spin operators. 
"""
function spin(s::String; D::Integer=2)
    ny = if isequal(D, 2) 
        sum(isequal(si, 'y')||isequal(si, 'Y') for si in s)
    else
        sum(isequal(si, 'y') for si in s)
    end
    mat = isone(length(s)) ? spin_dict(s[1], D) : kron([spin_dict(si, D) for si in s]...)
    sign = iszero(mod(ny, 2)) ? (-1)^(ny÷2) : (-1im)^ny
    sign * mat
end
spin(c::Number, s::String; D::Integer=2) = c * spin(s, D=D)

"""
    spin(spins::Tuple{<:Number, String}...; D::Integer=2)

Return matrix for spin operators. 
The spins should be an iterable onject, each item is of the form (::Number, ::String).
"""
function spin(spins::Tuple{<:Number, String}...; D::Integer=2)
    sum(ci * spin(si, D=D) for (ci, si) in spins)
end

spin(spins::AbstractVector{<:Tuple{<:Number, String}}; D::Integer=2) = spin(spins..., D=D)

function operator(s::String, inds::AbstractVector{<:Integer}, basis::AbstractBasis)
    mat = spin(s, D=base(basis))
    operator(mat, inds, basis)
end

function operator(s::String, inds::AbstractVector{<:Integer}, L::Integer; base::Integer=2)
    basis = tensorbasis(L, base=base)
    operator(s, inds, basis)
end

function trans_inv_operator(s::String, inds::AbstractVector{<:Integer}, basis::AbstractBasis)
    mat = spin(s, D=base(basis))
    trans_inv_operator(mat, inds, basis)
end

function trans_inv_operator(s::String, basis::AbstractBasis)
    mat = spin(s, D=base(basis))
    trans_inv_operator(mat, length(s), basis)
end

function trans_inv_operator(s::String, L::Integer; base::Integer=2)
    basis = TensorBasis(L, base=base)
    trans_inv_operator(s, basis)
end
