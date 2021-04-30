#-----------------------------------------------------------------------------------------------------
# Spin Matrices
#-----------------------------------------------------------------------------------------------------
# Dictionary
spin_coeff(D::Integer) = [sqrt(i*(D-i)) for i = 1:D-1]
function spin_dict(D::Integer)
    J = (D-1)/2
    coeff = [sqrt(i*(D-i)) for i = 1:D-1]
    sp = sparse(1:D-1, 2:D, coeff, D, D)
    sm = sp'
    sx = (sp+sm)/2
    sy = (sp-sm)/2
    sz = sparse(1:D, 1:D, J:-1:-J)
    s0 = sparse(1.0I, D, D)
    Dict([('+', sp), ('-', sm), ('x', sx), ('y', sy), ('z', sz), ('1', s0)])
end
#-----------------------------------------------------------------------------------------------------
# Atomic spin matrix
function spin_atom(s::String, dic::Dict)
    n = length(s)
    ny = sum(si == 'y' for si in s)
    temp = n == 1 ? dic[s[1]] : kron([dic[si] for si in s]...)
    P = mod(ny, 2) == 0 ? (-1)^(ny÷2) : (1im)^ny
    P * temp
end
#-----------------------------------------------------------------------------------------------------
export spin
"""
    spin(spins; D=2)

Return matrix for spin operators. 
The spins should be an iterable onject, each item is of the form (::Number, ::String).
"""
function spin(spins...; D::Integer=2)
    dic = spin_dict(D)
    res = sum(ci * spin_atom(si, dic) for (ci, si) in spins)
    Array(res)
end

spin(spins::AbstractVector{<:Tuple{<:Number, String}}; D::Integer=2) = spin(spins..., D=D)

#-----------------------------------------------------------------------------------------------------
# Level statistics
#-----------------------------------------------------------------------------------------------------
function gapratio(dE::AbstractVector{<:Real})
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

function meangapratio(E::AbstractVector{<:Real})
    dE = diff(E)
    r = gapratio(dE)
    sum(r)/length(r)
end
