#-------------------------------------------------------------------------------------------------------------------------
# Abelian Operator
#-------------------------------------------------------------------------------------------------------------------------
"""
    AbelianOperator{Tp <: Number} 

Operator of Abelian group, which can be product of elementary cycling groups.
The operator keeps tract of the state indices, which gives coefficients.

Properties:
-----------
- s: State indices
- g: Group orders
- c: Coefficients
- f: Permutation functions
"""
struct AbelianOperator{Tp <: Number}
    s::Vector{Int}
    g::Vector{Int}
    c::Vector{Vector{Tp}}
    f::Vector
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    AbelianOperator(g::Int, k::Integer, f)

Create elementary cycling group operator.

Inputs:
-------
- g: Group order 
- k: Momentum 
- f: Permutation function
"""
function AbelianOperator(g::Int, k::Integer, f)
    c = if iszero(k)
        ones(g)
    elseif 2k == g
        [iseven(j) ? 1 : -1 for j in 0:g-1]
    else
        phase = 1im * 2π * k / g
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
"""
Return the phase 

θ(s) = ∏ⱼexp(i*k*sⱼ)
"""
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
"""
Apply Abelian Operator on digits.
"""
function (ag::AbelianOperator)(dgt::Vector)
    for i in eachindex(ag.s)
        ag.f[i](dgt) 
        (ag.s[i] = ag.s[i] + 1) > ag.g[i] ? (ag.s[i] = 1) : break
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
"""
Find the group element that shift |d⟩ to the representing digit |d̃⟩.

|d̃⟩ = ∏ⱼ gⱼ^sⱼ|d⟩ 
"""
function shift_canonical!(dgt, g::AbelianOperator; base=2)
    init!(g)
    Im = index(dgt; base)
    ms = g.s[:]
    for _ in 1:order(g)
        g(dgt) 
        In = index(dgt; base)
        if In < Im 
            Im = In
            ms .= g.s
        end
    end
    g.s .= ms 
    Im, g
end

#-------------------------------------------------------------------------------------------------------------------------
# Abelian Basis
#-------------------------------------------------------------------------------------------------------------------------
struct AbelianBasis{Ti <: Integer, Tg <: Number} <: AbstractPermuteBasis
    dgt::Vector{Ti}         # Digits
    I::Vector{Ti}           # Representing states
    R::Vector{Float64}      # Normalization
    G::AbelianOperator{Tg}  # Generator
    B::Ti                   # Base
end

#-------------------------------------------------------------------------------------------------------------------------
function AbelianBasis(
    dtype::DataType=Int64;
    L::Integer, G::AbelianOperator, base::Integer=2, f=x->true,
    alloc=1000, threaded::Bool=true
)
    Ng = order(G)
    
    C = zeros(Ng)
    for i in eachindex(C)
        iszero(mod(Ng, i)) && (C[i] = sqrt(Ng * i))
    end

    I, R = if threaded
        nt = Threads.nthreads()
        ni = dividerange(base^L, nt)
        nI = Vector{Vector{dtype}}(undef, nt)
        nR = Vector{Vector{Float64}}(undef, nt)
        Threads.@threads for ti in 1:nt
            nI[ti], nR[ti] = _abelian_select(f, deepcopy(G), L, ni[ti], C; base, alloc)
        end
        vcat(nI...), vcat(nR...)
    else
        _abelian_select(f, G, L, 1:base^L, C; base, alloc)
    end
    AbelianBasis(zeros(dtype, L), I, R, G, base)
end

function _abelian_select(
    f, 
    G,
    L::Integer, 
    rg::UnitRange{T},
    C; 
    base::Integer=2, 
    alloc::Integer=1000
) where T <: Integer
    dgt = zeros(T, L)
    Is = T[]
    Rs = Float64[]
    sizehint!(Is, alloc)
    sizehint!(Rs, alloc)
    for i in rg
        change!(dgt, i; base)
        f(dgt) || continue
        Q, n = check_min(dgt, G; base)
        Q || continue 
        push!(Is, i)
        push!(Rs, C[n])
    end
    Is, Rs
end
#-------------------------------------------------------------------------------------------------------------------------
function index(B::AbelianBasis)
    Im, g = shift_canonical!(B.dgt, B.G; base=B.B)
    ind = binary_search(B.I, Im)
    if iszero(ind)
        return zero(eltype(B)), one(eltype(B.I))
    else
        return phase(g) * B.R[ind], ind
    end
end
#-------------------------------------------------------------------------------------------------------------------------
function schmidt(v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::AbelianBasis)
    dgt, R, g = b.dgt, b.R, b.G
    S = SchmidtMatrix(promote_type(eltype(v), eltype(b)), b, Ainds)
    for i in eachindex(v)
        init!(g)
        change!(b, i)
        val = v[i] / R[i]
        addto!(S, val)
        for _ in 2:order(g)
            g(dgt)
            addto!(S, phase(g) * val)
        end
    end
    S.M
end
#-------------------------------------------------------------------------------------------------------------------------
export basis
function basis(
    dtype::DataType=Int64; 
    L::Integer, 
    f=nothing,
    base::Integer=2, 
    N::Union{Nothing, Integer}=nothing, 
    k::Union{Nothing, Integer}=nothing, a::Integer=1, 
    p::Union{Nothing, Integer}=nothing, 
    z::Union{Nothing, Integer}=nothing,
    threaded::Bool=base^L>3000
)
    gs = AbelianOperator[]
    isnothing(k) || push!(gs, AbelianOperator(L÷a, k, x -> circshift!(x, a)))
    isnothing(p) || push!(gs, AbelianOperator(2, isone(-p) ? 1 : 0, reverse!))
    isnothing(z) || push!(gs, AbelianOperator(2, isone(-z) ? 1 : 0, x -> spinflip!(x, base)))
    
    if isempty(gs)
        isnothing(f) && isnothing(N) && return TensorBasis(;L, base)
        return ProjectedBasis(dtype; L, f, N, base, threaded)
    else
        g = if isnothing(N)
            isnothing(f) ? x -> true : f
        else
            num = L*(base-1)-N
            isnothing(f) ? x -> (sum(x) == num) : x -> (sum(x) == num && f(x))
        end
        return AbelianBasis(dtype; L, G=sum(gs), base, f=g, threaded)
    end
end





