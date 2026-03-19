#-------------------------------------------------------------------------------------------------------------------------
# Abelian Operator
#-------------------------------------------------------------------------------------------------------------------------
"""
    AbelianOperator{Tp <: Number} 

Internal representation of a product of commuting cyclic symmetry generators.

An `AbelianOperator` stores enough state to iterate through all elements of a
finite Abelian group action on a digit string while tracking the accumulated
phase associated with the current group element.

This is the symmetry backend used by [`AbelianBasis`](@ref).
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

Construct one elementary cyclic generator.

Arguments:
- `g`: order of the cyclic group.
- `k`: momentum or character label selecting the phase representation.
- `f`: in-place permutation acting on a digit buffer.

Returns:
- An [`AbelianOperator`](@ref) initialized at the identity group element.
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
"""
    +(g1::AbelianOperator, g2::AbelianOperator)

Form the direct-product combination of two commuting Abelian generators.
"""
function +(g1::AbelianOperator, g2::AbelianOperator)
    AbelianOperator([g1.s; g2.s], [g1.g; g2.g], [g1.c; g2.c], [g1.f; g2.f])
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    order(g::AbelianOperator)

Return the total number of group elements represented by `g`.
"""
function order(g::AbelianOperator) 
    prod(g.g)
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    phase(g::AbelianOperator)

Return the current character phase associated with the internal group state of
`g`.

θ(s) = ∏ⱼexp(i*k*sⱼ)
"""
function phase(g::AbelianOperator)
    prod(g.c[i][g.s[i]] for i in eachindex(g.c))
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    init!(g::AbelianOperator)

Reset the internal group-element counters of `g` to the identity element.
"""
function init!(g::AbelianOperator)
    for i in eachindex(g.s)
        g.s[i] = 1 
    end
    g
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    (ag::AbelianOperator)(dgt)

Apply the next group action to `dgt` and advance the internal group state.

This mutates both `dgt` and the internal counters stored in `ag`.
"""
function (ag::AbelianOperator)(dgt::Vector)
    for i in eachindex(ag.s)
        ag.f[i](dgt) 
        (ag.s[i] = ag.s[i] + 1) > ag.g[i] ? (ag.s[i] = 1) : break
    end
    dgt
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    check_min(dgt, g::AbelianOperator; base=2)

Check whether `dgt` is the canonical representative of its full Abelian orbit.

Returns:
- `(true, n)` when `dgt` is canonical and has stabilizer size `n`,
- `(false, 0)` otherwise.

This helper is used during basis construction to decide whether a product state
should be stored as a representative.
"""
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
    shift_canonical!(dgt, g::AbelianOperator; base=2)

Find the canonical representative of the Abelian orbit of `dgt` and record the
group element that maps the current state to that representative.

|d̃⟩ = ∏ⱼ gⱼ^sⱼ|d⟩ 

Returns:
- `Im`: the representative index,
- `g`: the same operator with internal counters set to the canonicalizing group
  element.
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
"""
    AbelianBasis

Basis built from a collection of commuting discrete symmetries represented as
Abelian actions on digit strings.

This is the general symmetry-reduction backend used by the high-level
[`basis`](@ref) constructor when one or more symmetry quantum numbers such as
translation momentum, reflection parity, or spin-flip parity are requested.
"""
struct AbelianBasis{Ti <: Integer, Tg <: Number} <: AbstractPermuteBasis
    dgt::Vector{Ti}         # Digits
    I::Vector{Ti}           # Representing states
    R::Vector{Float64}      # Normalization
    G::AbelianOperator{Tg}  # Generator
    B::Ti                   # Base
end
order(b::AbelianBasis) = order(b.G)
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

"""
    _abelian_select(f, G, L, rg, C; base=2, alloc=1000)

Low-level worker that scans a range of product-state indices and selects
canonical Abelian-orbit representatives.

This helper is used by both the threaded and non-threaded [`AbelianBasis`](@ref)
constructor paths.
"""
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
"""
    index(B::AbelianBasis)

Interpret the current digit buffer as coordinates in the Abelian-reduced basis.

The returned coefficient includes both the stored normalization and the phase
associated with the group element that maps the current digits to the canonical
representative.
"""
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
"""
    schmidt(v, Ainds, b::AbelianBasis; B1=nothing, B2=nothing)

Construct the Schmidt matrix of a state represented in an [`AbelianBasis`](@ref).

Unlike onsite bases, each basis coefficient must be expanded over the entire
Abelian orbit with the correct symmetry phase before it contributes to the
bipartite decomposition.
"""
function schmidt(v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::AbelianBasis; B1=nothing, B2=nothing)
    dgt, R, g = b.dgt, b.R, b.G
    S = schmidtmatrix(promote_type(eltype(v), eltype(b)), b, Ainds, B1, B2)
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
"""
    basis(dtype::DataType=Int64; L, f=nothing, base=2, N=nothing, k=nothing, a=1, p=nothing, z=nothing, threaded=base^L>3000)

High-level basis constructor for the common symmetry combinations supported by
EDKit.

Keywords:
- `L`: system size.
- `f`: optional predicate used to impose additional local constraints.
- `base`: local Hilbert-space dimension.
- `N`: fixed U(1)-like charge sector.
- `k`: translation quantum number.
- `a`: unit-cell length used together with `k`.
- `p`: reflection-parity eigenvalue `±1`.
- `z`: spin-flip eigenvalue `±1`.

Return value:
- `TensorBasis` if no symmetry or constraint is requested.
- `ProjectedBasis` if only `f` and/or `N` is requested.
- `AbelianBasis` when one or more discrete symmetries (`k`, `p`, `z`) are used.

Notes:
- `k` and `p` are only simultaneously valid in the symmetry-compatible momentum
  sectors handled by the underlying basis implementation.
- `N` and `z` are only compatible in the half-filling sector for spin-1/2
  systems, mirroring the restrictions of the dedicated basis types.
- This is the most convenient user-facing entry point when you want to combine
  several commuting symmetries without manually choosing a concrete basis type.
"""
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


