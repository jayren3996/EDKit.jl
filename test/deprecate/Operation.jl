#-----------------------------------------------------------------------------------------------------
# Functions for Operations.

# Operation is the collection of "Operator"s, together with a "Basis".
# Functions on "Operation"s are design to mimic the behavior of matrices.
#-----------------------------------------------------------------------------------------------------
"""
    Operator{T1<:Number, T2<:Integer}

Type that contains 2 fields:

1. `mat`: Matrices representation of local operator.
2. `inds`: Indeces of sites the operator acts on.
"""
struct Operator{T1<:Number, T2<:Integer}
    mat::Matrix{T1}
    inds::Vector{T2}
end
#-----------------------------------------------------------------------------------------------------
# Basic functions on "Operator":

# 1. *: Multiplication by a number.
# 2. /: Division by a number.
#-----------------------------------------------------------------------------------------------------
*(c::Number, o::Operator) = Operator(c * o.mat, o.inds)
/(o::Operator, c::Number) = Operator(o.mat / c, o.inds)
#-----------------------------------------------------------------------------------------------------
# Operator to fill vector/matrix.

# To finally get the multiplication function, here we focus on a single row manipulation:

# 1. Start with a given product state(digits), we cut a segment from it;
# 2. Iterate all possible digits in the segment, and read the matrix element of operator;
# 3. Add to the vector/matrix element.
#-----------------------------------------------------------------------------------------------------
"""
    addtovec!(vec::AbstractVector, opt::Operator, basis::Basis, coeff::Number=1)

Elementary multiplication step. Calculate single vector element of an operator * vector.

- `vec`: Vector to fill
- `opt`: Operator
- `basis`: Indicate the specific row
- `coeff`: The vector element of input state
"""
function addtovec!(vec::AbstractVector, opt::Operator, basis::Basis, coeff::Number=1)
    basis_view = view(basis, opt.inds)              # Create a view on the segment of the digits
    index_view = index(basis_view)                  # Get the initial index of the viewed segment
    opt_column = opt.mat[:, index_view] * coeff     # Get the elements in row of index_view
    row_number = length(opt_column)
    for k = 1:row_number
        change!(basis_view, k)
        i = index(basis)
        vec[i] += opt_column[k]
    end                                             # Fill the vector
    change!(basis_view, index_view)                 # Reset the segment
end
#-----------------------------------------------------------------------------------------------------
"""
    addtovecs!(vecs::AbstractMatrix, opt::Operator, basis::Basis, coeff::AbstractVector)

Elementary multiplication step. Calculate single row of matrix element of an operator * matrix.

- vecs : Vectors to fill
- opt  : Operator
- basis: Indicate the specific row
- coeff: The row of matrix elements of input states
"""
function addtovecs!(vecs::AbstractMatrix, opt::Operator, basis::Basis, coeff::AbstractVector{<:Number})
    basis_view = view(basis, opt.inds)              # Create a view on the segment of the digits
    index_view = index(basis_view)                  # Get the initial index of the viewed segment
    opt_column = opt.mat[:, index_view]             # Get the elements in row of index_view
    row_number = length(opt_column)                 # Fill the matrix
    for k = 1:row_number
        change!(basis_view, k)
        i = index(basis)
        vecs[i, :] += opt_column[k] * coeff
    end
    change!(basis_view, index_view)                 # Reset the segment
end
#-----------------------------------------------------------------------------------------------------
# Type: Operation
#-----------------------------------------------------------------------------------------------------
export Operation
"""
    Operation{OptType<:AbstractVector{<:Operator}, BasType<:Basis}

Type that contains 2 fields:

1. `opts`: List of operators
2. `basis`: Basis for the system
"""
struct Operation{T1 <: Number, T2 <: Integer, T3}
    opts::Vector{Operator{T1, T2}}
    basis::Basis{T3}
end
eltype(::Operation{T1, T2, T3}) where T1 where T2 where T3 = T1
function size(opt::Operation)
    basis = opt.basis
    dim = basis.base ^ basis.len
    (dim, dim)
end
#-----------------------------------------------------------------------------------------------------
# Basis Operation initiation
#-----------------------------------------------------------------------------------------------------
export operation
"""
    operation(mats, inds, len=0; base=0)

Canonical construction method for `Operation` object.

- `mats`: List of matrix representation of operators.
- `inds`: List of sites the operastors act on.
- `len`: Specify the size of the system. Default `len=0` will induce auto deduction.
- `base`: Quantum number of each sites. Default `base=0` will induce auto deduction.
"""
function operation(
    mats::AbstractVector{<:AbstractMatrix}, 
    inds::AbstractVector{<:AbstractVector}, 
    len::Integer=0; 
    base::Integer=0
)
    B = begin
        b = base==0 ? round(Int64, size(mats[1], 1)^(1/length(inds[1]))) : Int64(base)
        l = len==0 ? maximum(maximum.(inds)) : Int64(len)
        basis(b, l)
    end
    mat_type = promote_type(eltype.(mats)...)
    ind_type = promote_type(eltype.(inds)...)
    O = [Operator(Array{mat_type}(mats[i]), Array{ind_type}(inds[i])) for i=1:length(mats)]
    Operation(O, B)
end
#-----------------------------------------------------------------------------------------------------
# Basic Functions for type Operation
#-----------------------------------------------------------------------------------------------------
*(c::Number, o::Operation) = Operation(c .* o.opts, o.basis)
/(o::Operation, c::Number) = Operation(o.opts ./ c, o.basis)
+(opt1::Operation, opt2::Operation) = Operation(vcat(opt1.opts, opt2.opts), opt1.basis)
-(o1::Operation, o2::Operation) = o1 + ((-1) * o2)
function sum(ol::AbstractVector{<:Operation})
    basis = ol[1].basis
    opts = vcat([oi.opts for oi in ol]...)
    Operation(opts, basis)
end
#-----------------------------------------------------------------------------------------------------
# Operation to fill vector/matrix:

# The method is basically iterate each operator to the given digits-represented basis
#-----------------------------------------------------------------------------------------------------
"""
    addtovec!(vec::AbstractVector, opt::Operation, coeff::Number=1)

Elementary multiplication step. Calculate single vector element of an operation * vector.
The row information is stored in the basis of operation.

- `vec`: Vector to fill.
- `opt`: Operation.
- `coeff`: The vector element of input state.
"""
function addtovec!(vec::AbstractVector, opt::Operation, coeff::Number=1)
    basis = opt.basis
    opts = opt.opts
    num_of_opts = length(opts)
    for i = 1:num_of_opts
        addtovec!(vec, opts[i], basis, coeff)
    end # Multiply by each operator and add them all to vector
end
#-----------------------------------------------------------------------------------------------------
"""
    addtovecs!(vecs::AbstractMatrix, opt::Operation, coeff::AbstractVector)

Elementary multiplication step. Calculate single row of matrix element of an operation * matrix.
The row information is stored in the basis of operation.

- `vecs`: Vectors to fill
- `opt`: Operation
- `coeff`: The row of matrix elements of input states
"""
function addtovecs!(vecs::AbstractMatrix, opt::Operation, coeff::AbstractVector{<:Number})
    basis = opt.basis
    opts = opt.opts
    num_of_opts = length(opts)
    for i = 1:num_of_opts
        addtovecs!(vecs, opts[i], basis, coeff)
    end # Multiply by each operator and add them all to vectors
end
#-----------------------------------------------------------------------------------------------------
# Multiplication

# The idea is to iterate all product state basis, and get all the vector/matrix elements
#-----------------------------------------------------------------------------------------------------
"""
    mul!(vec::AbstractVector, opt::Operation, state::AbstractVector)

Full multiplication for operation and vector.

- `vec`: Vector to fill.
- `opt`: Operation.
- `state`: Input vector.
"""
function mul!(vec::AbstractVector, opt::Operation, state::AbstractVector)
    basis = opt.basis
    for j = 1:length(state)
        change!(basis, j)
        addtovec!(vec, opt, state[j])
    end
end
#-----------------------------------------------------------------------------------------------------
"""
    mul!(mat::AbstractVector, opt::Operation, state::AbstractVector)

Full multiplication for operation and matrix.

- `mat`: Matrix to fill.
- `opt`: Operation.
- `states`: Input states(matrix).
"""
function mul!(mat::AbstractMatrix, opt::Operation, states::AbstractMatrix)
    basis = opt.basis
    for j = 1:size(states, 1)
        change!(basis, j)
        addtovecs!(mat, opt, states[j, :])
    end
end
#-----------------------------------------------------------------------------------------------------
"""
    *(opt::Operation, vec_or_mat::AbstractVecOrMat)

General multiplication for operation and vector/matrix.

- `opt`: Operation.
- `vec_or_mat`: Input vector/matrix.
"""
function *(opt::Operation, vec_or_mat::AbstractVecOrMat)
    ctype = promote_type(eltype(opt), eltype(vec_or_mat))
    out = zeros(ctype, size(vec_or_mat))
    mul!(out, opt, vec_or_mat)
    out
end
#-----------------------------------------------------------------------------------------------------
# Fill matrix
#-----------------------------------------------------------------------------------------------------
export fillmat!
"""
    fillmat!(mat::AbstractMatrix, opt::Operation)

Fill the zero matrix with matrix element from operation.

- `mat`: Matrix to fill.
- `opt`: Operation.
"""
function fillmat!(mat::AbstractMatrix, opt::Operation)
    basis = opt.basis
    row_number = size(mat, 1)
    for j = 1:row_number
        change!(basis, j)
        addtovec!(view(mat, :, j), opt)
    end
end
#-----------------------------------------------------------------------------------------------------
"""
    Array(opt::Operation)

Return the full matrix form of `opt`.
"""
function Array(opt::Operation)
    mat = zeros(eltype(opt), size(opt))
    fillmat!(mat, opt)
    mat
end