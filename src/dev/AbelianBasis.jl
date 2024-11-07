include("../EDKit.jl")
using Main.EDKit
using LinearAlgebra
import Base.:+
import EDKit: binary_search
#-------------------------------------------------------------------------------------------------------------------------
# Abelian Operator
#-------------------------------------------------------------------------------------------------------------------------
struct AbelianOperator{Tp <: Number}
    s::Vector{Int}
    g::Vector{Int}
    c::Vector{Vector{Tp}}
    f::Vector
end
#-------------------------------------------------------------------------------------------------------------------------
function AbelianOperator(g::Int, k::Integer, f)
    phase = 1im * 2π * k / g
    c = if iszero(k)
        ones(g)
    elseif 2k == g
        [iseven(j) ? 1 : -1 for j in 0:g-1]
    else
        [exp(phase * j) for j in 0:g-1]
    end
    AbelianOperator([1], [g], [c], [f])
end
#-------------------------------------------------------------------------------------------------------------------------
function +(g1::AbelianOperator, g2::AbelianOperator)
    AbelianOperator([g1.s; g2.s], [g1.g; g2.g], [g1.c; g2.c], [g1.f; g2.f])
end
#-------------------------------------------------------------------------------------------------------------------------
function order(g::AbelianOperator) 
    prod(g.g)
end
#-------------------------------------------------------------------------------------------------------------------------
function phase(g::AbelianOperator)
    prod(g.c[i][g.s[i]] for i in eachindex(g.c))
end
#-------------------------------------------------------------------------------------------------------------------------
function init!(g::AbelianOperator)
    for i in eachindex(g.s)
        g.s[i] = 1 
    end
    g
end
#-------------------------------------------------------------------------------------------------------------------------
function (ag::AbelianOperator)(dgt::Vector)
    ag.f[1](dgt) 
    ag.s[1] += 1 
    for i in 1:length(ag.s)-1 
        ag.s[i] > ag.g[i] || break
        ag.s[i] = 1
        ag.f[i+1](dgt)
        ag.s[i+1] += 1 
    end
    dgt
end
#-------------------------------------------------------------------------------------------------------------------------
function check_min(dgt, g::AbelianOperator; base=2)
    init!(g)
    I0 = index(dgt; base)
    N = 1 
    for _ in 2:order(g)
        g(dgt) 
        In = index(dgt; base)
        if In < I0 
            return false, 0 
        elseif In == I0 
            abs(phase(g)-1) < 1e-14 || return false, 0
            N += 1 
        end
    end
    true, N
end
#-------------------------------------------------------------------------------------------------------------------------
function shift_canonical!(dgt, g::AbelianOperator; base=2)
    init!(g)
    Im = index(dgt; base)
    ms = deepcopy(g.s)
    for i in 1:order(g)
        g(dgt) 
        In = index(dgt; base)
        if In < Im 
            Im = In
            ms = deepcopy(g.s)
        end
    end
    g.s .= ms 
    Im, g
end

#-------------------------------------------------------------------------------------------------------------------------
# Abelian Basis
#-------------------------------------------------------------------------------------------------------------------------
struct AbelianBasis{Ti <: Integer, Tg <: Number} <: EDKit.AbstractPermuteBasis
    dgt::Vector{Ti}         # Digits
    I::Vector{Ti}           # Representing states
    R::Vector{Float64}      # Normalization
    G::AbelianOperator{Tg}  # Generator
    B::Ti                   # Base
end

#-------------------------------------------------------------------------------------------------------------------------
function AbelianBasis(
    dtype::DataType=Int64;
    L::Integer, G::AbelianOperator, base::Integer=2, f=x->true
)
    MAX = base ^ L + 1
    dgt = zeros(dtype, L)
    Ng = order(G)
    C = [sqrt(Ng * i) for i in 1:Ng]
    I = dtype[]
    R = Float64[]
    for i in 1:MAX 
        change!(dgt, i; base)
        f(dgt) || continue
        Q, n = check_min(dgt, G; base)
        if Q 
            push!(I, i)
            push!(R, C[n])
        end
    end
    AbelianBasis(zeros(dtype, L), I, R, G, base)
end

function EDKit.index(B::AbelianBasis)
    Im, g = shift_canonical!(B.dgt, B.G; base=B.B)
    ind = binary_search(B.I, Im)
    if iszero(ind)
        return zero(eltype(B)), one(eltype(B.I))
    else
        return phase(g) * B.R[ind], ind
    end
end

function EDKit.schmidt(v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::AbelianBasis)
    dgt, R, g = b.dgt, b.R, b.G
    S = EDKit.SchmidtMatrix(promote_type(eltype(v), eltype(b)), b, Ainds)
    for i in eachindex(v)
        init!(g)
        change!(b, i)
        val = v[i] / R[i]
        addto!(S, val)
        for j in 2:order(g)
            g(dgt)
            addto!(S, phase(g) * val)
        end
    end
    S.M
end






