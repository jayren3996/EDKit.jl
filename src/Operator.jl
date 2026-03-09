"""
Operator

Construction of `Operator` object.
"""
#---------------------------------------------------------------------------------------------------
export operator, trans_inv_operator
"""
    Operator{Tv}

Many-body operator stored as a sum of local terms on a chosen basis.

Each term is represented by a sparse local matrix `M[n]` acting on lattice sites
`I[n]`, while `B` controls how those local actions are embedded into the full
Hilbert space.
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

Construct an [`Operator`](@ref) from local matrices and the sites they act on.

Arguments:
- `mats`: local operator matrices.
- `inds`: site lists matching `mats`, using 1-based site labels.
- `B`: basis in which the operator should act.

Repeated site patterns are merged automatically by summing their local matrices.

Examples:
```julia
mat = spin((1.0, "xx"), (1.0, "yy"))
H = operator(fill(mat, L), [[i, mod(i, L) + 1] for i in 1:L], TensorBasis(L=L))
```
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
"""
    trans_inv_operator(mat::AbstractMatrix, ind, basis_or_length)

Construct a translationally invariant operator by translating one local term
around the lattice.

`ind` specifies the support of the seed term, for example `1:2` for a nearest-
neighbor bond term. The function then generates all translated copies with
periodic boundary conditions and sums them into a single [`Operator`](@ref).

This accepts either a concrete basis or just a system size `L`.
"""
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
"""
    addto!(M::AbstractMatrix, opt::Operator)

Accumulate the matrix representation of `opt` into the preallocated array `M`.

This is useful when you want control over the output array type or want to reuse
existing storage instead of calling `Array(opt)` or `sparse(opt)`.
"""
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
    mul(opt::Operator, v)
    mul(opt::Operator, m)

Multi-threaded multiplication of an [`Operator`](@ref) with a vector or matrix.

Unlike `opt * v`, this version parallelizes over input columns / basis states
using `Threads.@threads`.
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

Internal helper used to apply one local term to the current basis state in `b`.

`b.dgt` is interpreted as the current many-body basis configuration. The result
of acting with the local matrix `M` on sites `I` is accumulated in `target`.
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
@inline function spin_ny(s::AbstractString, D::Integer)
    n = 0
    if isequal(D, 2)
        for c in s
            n += isequal(c, 'y') || isequal(c, 'Y')
        end
    else
        for c in s
            n += isequal(c, 'y')
        end
    end
    n
end

function spin_product(s::AbstractString, D::Integer)
    isempty(s) && error("Empty spin string.")
    mat = spin_dict(first(s), D)
    for c in Iterators.drop(s, 1)
        mat = kron(mat, spin_dict(c, D))
    end
    mat
end

const SPIN_CACHE = LRU{Tuple{Int, String}, Any}(maxsize=256)
#---------------------------------------------------------------------------------------------------
"""
    spin(s::String; D::Integer=2)

Return the local operator encoded by the string `s`.

For `D = 2`, uppercase Pauli-style labels such as `"X"`, `"Y"`, `"Z"` are
supported, as well as lowercase spin-operator labels such as `"x"`, `"y"`,
`"z"`, `"+"`, `"-"`, and `"1"`. Multi-site strings like `"xx"` or `"x1z"`
build Kronecker products in the given order.
"""
function spin(s::String; D::Integer=2)
    key = (Int(D), s)
    mat = get!(SPIN_CACHE, key) do
        ny = spin_ny(s, D)
        raw = spin_product(s, D)
        sign = iszero(mod(ny, 2)) ? (-1)^(ny÷2) : (-1im)^ny
        sign * raw
    end
    copy(mat)
end
spin(c::Number, s::String; D::Integer=2) = c * spin(s, D=D)
#---------------------------------------------------------------------------------------------------
"""
    spin(spins::Tuple{<:Number, String}...; D::Integer=2)

Build a linear combination of local operators.

Each argument should be a pair `(coefficient, "operator_string")`, for example
`spin((1.0, "xx"), (1.0, "yy"), (0.5, "zz"))`.
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
