#-----------------------------------------------------------------------------------------------------
# Spin Matrices
#-----------------------------------------------------------------------------------------------------
# Dictionary
spin_coeff(D::Integer) = [sqrt(i*(D-i)) for i = 1:D-1]
spin_Sp(D::Integer) = sparse(1:D-1, 2:D, spin_coeff(D), D, D)
spin_Sm(D::Integer) = sparse(2:D, 1:D-1, spin_coeff(D), D, D)

function spin_Sx(D::Integer)
    coeff = spin_coeff(D)
    sp = sparse(1:D-1, 2:D, coeff, D, D)
    sp + sp'
end

function spin_iSy(D::Integer)
    coeff = spin_coeff(D)
    sp = sparse(1:D-1, 2:D, coeff, D, D)
    sp - sp'
end

function spin_Sz(D::Integer)
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
    else error("Invalid spin symbol: $c.")
    end
end

# Spin matrix
export spin
"""
    spin(s::String, D::Integer)

Return matrix for spin operators. 
"""
function spin(s::String, D::Integer)
    ny = sum(isequal(si, 'y') for si in s)
    mat = isone(length(s)) ? spin_dict(s[1], D) : kron([spin_dict(si, D) for si in s]...)
    sign = iszero(mod(ny, 2)) ? (-1)^(ny÷2) : (1im)^ny
    sign * mat
end


"""
    spin(spins::Tuple{<:Number, String}...; D::Integer=2)

Return matrix for spin operators. 
The spins should be an iterable onject, each item is of the form (::Number, ::String).
"""
function spin(spins::Tuple{<:Number, String}...; D::Integer=2)
    sum(ci * spin(si, D) for (ci, si) in spins)
end

spin(spins::AbstractVector{<:Tuple{<:Number, String}}; D::Integer=2) = spin(spins..., D=D)

function operator(s::String, inds::AbstractVector{<:Integer}, basis::AbstractBasis)
    mat = spin(s, base(basis))
    operator(mat, inds, basis)
end

function trans_inv_operator(s::String, inds, basis::AbstractBasis)
    mat = spin(s, base(basis))
    trans_inv_operator(mat, inds, basis)
end

#-----------------------------------------------------------------------------------------------------
# Level statistics
#-----------------------------------------------------------------------------------------------------
export gapratio, meangapratio
function gapratio(E::AbstractVector{<:Real})
    dE = diff(E)
    r = zeros(length(dE)-1)
    for i = 1:length(r)
        if dE[i] < dE[i+1]
            r[i] = dE[i]/dE[i+1]
        elseif dE[i] > dE[i+1]
            r[i] = dE[i+1]/dE[i]
        else
            r[i] = 1.0
        end
    end
    r
end

meangapratio(E::AbstractVector{<:Real}) = sum(gapratio(E)) / (length(E) - 2)
