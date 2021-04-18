#-----------------------------------------------------------------------------------------------------
# More construction method for Operation
#-----------------------------------------------------------------------------------------------------
export onsite_operation
"""
    onsite_operation(mats)

Return operation of the product of onsite operators on each site.
"""
function onsite_operation(mats::Vector{<:Matrix})
    base = size(mats[1], 1)
    len = length(mats)
    inds = [[i] for i in 1:len]
    operation(mats, inds, len, base=base)
end
#-----------------------------------------------------------------------------------------------------
"""
    onsite_operation(mats, ind, len)

Return operation of the product of onsite operator on a set of given sites.
The sites that is acted on is given by the input `ind`.
`len` indicate the size of the system.
"""
function onsite_operation(
    mats::Vector{<:Matrix},
    ind::AbstractVector{<:Integer},
    len::Int64
)
    base = size(mats[1], 1)
    inds = [[i] for i in ind]
    operation(mats, inds, len, base=base)
end

#-----------------------------------------------------------------------------------------------------
"""
    onsite_operation(mat, len)

Return operation of the product of onsite operator `mat` on each site.
`len` indicate the size of the system.
"""
function onsite_operation(mat::Matrix, len::Int64)
    mats = fill(mat, len)
    base = size(mat, 1)
    inds = [[i] for i in 1:len]
    operation(mats, inds, len, base=base)
end

#-----------------------------------------------------------------------------------------------------
export trans_inv_operation
"""
    trans_inv_operation(mat, ind, len)

Return a 1d translational invariant operator.
`mat` is the operator acting on sites `ind`, and all `ind`+i sites.
`len` indicate the system length.
"""
function trans_inv_operation(
    mat::Matrix, 
    ind::AbstractVector{<:Integer}, 
    len::Int64
)
    mats = fill(mat, len)
    inds = [mod.(ind .+ (i-1), len) .+ 1 for i=0:len-1]
    operation(mats, inds, len)
end
