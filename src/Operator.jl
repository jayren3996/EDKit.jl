"""
Operator

Construction of `Operator` object.
"""
#---------------------------------------------------------------------------------------------------
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
#---------------------------------------------------------------------------------------------------
eltype(opt::Operator{Tv, Tb}) where Tv where Tb = promote_type(Tv, eltype(opt.B))
length(opt::Operator) = length(opt.M)
size(opt::Operator) = size(opt.B)
size(opt::Operator, i::Integer) = size(opt.B, i)
function Base.display(opt::Operator)
    println("Operator of size $(size(opt)) with $(length(opt)) terms.")
end
#---------------------------------------------------------------------------------------------------
"""
    operator(mats::AbstractVector{<:AbstractMatrix}, inds::AbstractVector{<:AbstractVector}, B::AbstractBasis)

Constructor for `Operator`. 

Inputs:
-------
- `mats`: List of matrices for local operators.
- `inds`: List of sites on which local operators act.
- `B`   : Basis. 

Outputs:
--------
- `O` : Operator.
"""
function operator(mats::AbstractVector{<:AbstractMatrix}, inds::AbstractVector{<:AbstractVector}, B::AbstractBasis)
    num = length(mats)
    @assert num == length(inds) "Numbers mismatch: $num matrices and $(length(inds)) indices."
    dtype = promote_type(eltype.(mats)...)
    M = Vector{SparseMatrixCSC{dtype, Int64}}(undef, num)
    I = Vector{Vector{Int64}}(undef, num)
    N = 0
    for i = 1:num
        iszero(mats[i]) && continue
        ind = inds[i]
        pos = findfirst(x -> isequal(x, ind), view(I, 1:N))
        if isnothing(pos)
            N += 1
            I[N] = ind
            M[N] = sparse(mats[i])
        else
            M[pos] += mats[i]
        end
    end
    deleteat!(M, N+1:num)
    deleteat!(I, N+1:num)
    Operator(M, I, B)
end

function operator(mats::AbstractVector{<:AbstractMatrix}, inds::AbstractVector{<:AbstractVector}, L::Integer)
    base = find_base(size(mats[1], 1), length(inds[1]))
    basis = TensorBasis(L=L, base=base)
    operator(mats, inds, basis)
end
operator(mats::AbstractVector{<:AbstractMatrix}, inds::AbstractVector{<:Integer}, C) = operator(mats, [[i] for i in inds], C)
operator(mats::AbstractVector{<:AbstractMatrix}, L::Integer) = operator(mats, [[i] for i in 1:L], L)
operator(mats::AbstractVector{<:AbstractMatrix}, B::AbstractBasis) = operator(mats, [[i] for i in 1:length(B.dgt)], B)
operator(mat::AbstractMatrix, ind::AbstractVector{<:Integer}, C) = operator([mat], [ind], C)
#---------------------------------------------------------------------------------------------------
function trans_inv_operator(mat::AbstractMatrix, ind::AbstractVector{<:Integer}, B::AbstractBasis)
    L = length(B.dgt)
    smat = sparse(mat)
    mats = fill(smat, L)
    inds = [mod.(ind .+ i, L) .+ 1 for i = -1:L-2]
    operator(mats, inds, B)
end

function trans_inv_operator(mat::AbstractMatrix, ind::AbstractVector{<:Integer}, L::Integer)
    base = find_base(size(mat, 1), length(ind))
    B = TensorBasis(L=L, base=base)
    trans_inv_operator(mat, ind, B)
end

trans_inv_operator(mat::AbstractMatrix, M::Integer, C) = trans_inv_operator(mat, 1:M, C)
#---------------------------------------------------------------------------------------------------
*(c::Number, opt::Operator) = operator(c .* opt.M, opt.I, opt.B)
*(opt::Operator, c::Number) = c * opt
/(opt::Operator, c::Number) = operator(opt.M ./ c, opt.I, opt.B)
function +(opt1::Operator, opt2::Operator)
    mats = vcat(opt1.M, opt2.M)
    inds = vcat(opt1.I, opt2.I)
    operator(mats, inds, opt1.B)
end
+(opt1::Operator, ::Nothing) = opt1
+(::Nothing, opt1::Operator) = opt1
-(opt::Operator) = Operator(-opt.M, opt.I, opt.B)
-(opt1::Operator, opt2::Operator) = opt1 + (-opt2)
#---------------------------------------------------------------------------------------------------
function LinearAlgebra.adjoint(opt::Operator)
    M = adjoint.(opt.M)
    operator(M, opt.I, opt.B)
end
#---------------------------------------------------------------------------------------------------
function find_base(a::Integer, b::Integer)
    isone(b) && return a
    for i in 2:a
        i^b == a && return i
    end
    error("Incompatible dimension: ($a, $b).")
end

#---------------------------------------------------------------------------------------------------
# Operator to matrices
#---------------------------------------------------------------------------------------------------
export addto!
function addto!(M::AbstractMatrix, opt::Operator)
    for j = 1:size(opt.B, 2)
        colmn!(view(M, :, j), opt, j)
    end
    M
end
#---------------------------------------------------------------------------------------------------
function Array(opt::Operator)
    M = zeros(eltype(opt), size(opt))
    if size(M, 1) > 0 && size(M, 2) > 0
        addto!(M, opt)
    end
    M
end
#---------------------------------------------------------------------------------------------------
function SparseArrays.sparse(opt::Operator)
    M = spzeros(eltype(opt), size(opt)...)
    if size(M, 1) > 0 && size(M, 2) > 0
        addto!(M, opt)
    end
    M
end
#---------------------------------------------------------------------------------------------------
LinearAlgebra.Hermitian(opt::Operator) = Array(opt) |> Hermitian
LinearAlgebra.Symmetric(opt::Operator) = Array(opt) |> Symmetric
LinearAlgebra.eigen(opt::Operator) = Array(opt) |> eigen
LinearAlgebra.eigvals(opt::Operator) = Array(opt) |>eigvals
LinearAlgebra.svd(opt::Operator) = Array(opt) |> svd
LinearAlgebra.svdvals(opt::Operator) = Array(opt) |> svdvals
#---------------------------------------------------------------------------------------------------
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

export mul
"""
Multi-threaded operator multiplication.
"""
function mul(opt::Operator, v::AbstractVector)
    ctype = promote_type(eltype(opt), eltype(v))
    nt = Threads.nthreads()
    ni = dividerange(length(v), nt)
    Ms = [zeros(ctype, size(opt, 1)) for i in 1:nt]
    Threads.@threads for i in 1:nt
        opt_c = Operator(opt.M, opt.I, copy(opt.B))
        for j in ni[i]
            colmn!(Ms[i], opt_c, j, v[j])
        end
    end
    sum(m for m in Ms)
end

function mul(opt::Operator, m::AbstractMatrix)
    ctype = promote_type(eltype(opt), eltype(m))
    nt = Threads.nthreads()
    ni = dividerange(size(m,1), nt)
    Ms = [zeros(ctype, size(opt, 1), size(m, 2)) for i in 1:nt]
    Threads.@threads for i in 1:nt
        opt_c = Operator(opt.M, opt.I, copy(opt.B))
        for j in ni[i]
            colmn!(Ms[i], opt_c, j, view(m, j, :))
        end
    end
    sum(m for m in Ms)
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

#---------------------------------------------------------------------------------------------------
# Helper functions
#---------------------------------------------------------------------------------------------------
"""
    colmn!(target::AbstractVecOrMat, M::SparseMatrixCSC, I::Vector{Int}, b::AbstractBasis, coeff=1)

Central helper function for operator multiplication. 

For a local matrix `M` acting on indices `I`, `colmn!` return the j-th colume (given by `b.dgt`) in 
the many-body basis. The result is writen inplace on the vector `target`.
"""
function colmn!(target::AbstractVecOrMat, M::SparseMatrixCSC, I::Vector{Int}, b::AbstractBasis, coeff=1)
    rows, vals = rowvals(M), nonzeros(M)
    j = index(b.dgt, I, base=b.B)
    change = false
    for i in nzrange(M, j)
        row, val = rows[i], vals[i]
        change!(b.dgt, I, row, base=b.B)
        C, pos = index(b)
        isa(target, AbstractVector) ? (target[pos] += coeff * C * val) : (target[pos, :] .+= (C * val) .* coeff)
        change = true
    end
    change && change!(b.dgt, I, j, base=b.B) 
    nothing
end
#---------------------------------------------------------------------------------------------------
function colmn!(target::AbstractVecOrMat, opt::Operator, j::Integer, coeff=1)
    b, M, I = opt.B, opt.M, opt.I
    r = change!(b, j)
    C = isone(r) ? coeff : coeff / r
    for i = 1:length(M)
        colmn!(target, M[i], I[i], b, C)
    end
end

#---------------------------------------------------------------------------------------------------
# Spin Matrices
#---------------------------------------------------------------------------------------------------
const PAULI = sparse.(
    [ I(2), [0 1; 1 0], [0 -1im; 1im 0], [1 0; 0 -1] ]
)
const GELLMANN = sparse.([
    I(3),
    [0 1 0; 1 0 0; 0 0 0], 
    [0 -1im 0; 1im 0 0; 0 0 0], 
    [1 0 0; 0 -1 0; 0 0 0], 
    [0 0 1; 0 0 0; 1 0 0], 
    [0 0 -1im; 0 0 0; 1im 0 0], 
    [0 0 0; 0 0 1; 0 1 0],
    [0 0 0; 0 0 -1im; 0 1im 0], 
    [1 0 0; 0 1 0; 0 0 -2] / sqrt(3),
    [1 0 0; 0 1 0; 0 0 1]
])

σ(i::Integer) = PAULI[i+1]
σ(Is::Integer...) = kron((σ(i) for i in Is)...)
σ(Is::Union{<:Tuple, <:AbstractVector}) = σ(Is...)

λ(i::Integer) = GELLMANN[i+1]
λ(Is::Integer...) = kron((λ(i) for i in Is)...)
λ(Is::Union{<:Tuple, <:AbstractVector}) = λ(Is...)
#---------------------------------------------------------------------------------------------------
export spin
spin_coeff(D::Integer) = [sqrt(i*(D-i)) for i = 1:D-1]
spin_Sp(D::Integer) = sparse(1:D-1, 2:D, spin_coeff(D), D, D)
spin_Sm(D::Integer) = sparse(2:D, 1:D-1, spin_coeff(D), D, D)
function spin_Sx(D::Integer)
    coeff = spin_coeff(D) / 2
    sp = sparse(1:D-1, 2:D, coeff, D, D)
    sp + sp'
end
function spin_iSy(D::Integer)
    coeff = spin_coeff(D) / 2
    sp = sparse(1:D-1, 2:D, coeff, D, D)
    sp - sp'
end
function spin_Sz(D::Integer)
    J = (D-1) / 2
    sparse(1:D, 1:D, J:-1:-J)
end
#---------------------------------------------------------------------------------------------------
function spin_dict(c::Char, D::Integer)
    isequal(c, '+') && return spin_Sp(D)
    isequal(c, '-') && return spin_Sm(D)
    isequal(c, 'x') && return spin_Sx(D)
    isequal(c, 'y') && return spin_iSy(D)
    isequal(c, 'z') && return spin_Sz(D)
    isequal(c, '1') && return spdiagm(ones(D))
    isequal(c, 'I') && return spdiagm(ones(D))

    @assert isequal(D, 2) "Invalid spin symbol: $c."
    isequal(c, 'X') && return sparse([0 1; 1  0])
    isequal(c, 'Y') && return sparse([0 1;-1  0])
    isequal(c, 'Z') && return sparse([1 0; 0 -1])
    error("Invalid spin symbol: $c.")
end
#---------------------------------------------------------------------------------------------------
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
#---------------------------------------------------------------------------------------------------
"""
    spin(spins::Tuple{<:Number, String}...; D::Integer=2)

Return matrix for spin operators. 
The spins should be an iterable onject, each item is of the form (::Number, ::String).
"""
function spin(spins::Tuple{<:Number, String}...; D::Integer=2)
    sum(ci * spin(si, D=D) for (ci, si) in spins)
end

spin(spins::AbstractVector{<:Tuple{<:Number, String}}; D::Integer=2) = spin(spins..., D=D)
#---------------------------------------------------------------------------------------------------
function operator(s::String, inds::AbstractVector{<:Integer}, basis::AbstractBasis)
    mat = spin(s, D=basis.B)
    operator(mat, inds, basis)
end
#---------------------------------------------------------------------------------------------------
function operator(s::String, inds::AbstractVector{<:Integer}, L::Integer; base::Integer=2)
    basis = TensorBasis(L=L, base=base)
    operator(s, inds, basis)
end
#---------------------------------------------------------------------------------------------------
function trans_inv_operator(s::String, inds::AbstractVector{<:Integer}, basis::AbstractBasis)
    mat = spin(s, D=basis.B)
    trans_inv_operator(mat, inds, basis)
end
#---------------------------------------------------------------------------------------------------
function trans_inv_operator(s::String, basis::AbstractBasis)
    mat = spin(s, D=basis.B)
    trans_inv_operator(mat, length(s), basis)
end
#---------------------------------------------------------------------------------------------------
function trans_inv_operator(s::String, L::Integer; base::Integer=2)
    basis = TensorBasis(L=L, base=base)
    trans_inv_operator(s, basis)
end


