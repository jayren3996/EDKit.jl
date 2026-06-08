export SpinfulFermionBasis, fermionmode
"""
    SpinfulFermionBasis{T<:Integer}

Occupation-number basis for `L` lattice sites carrying `S` fermion species
(e.g. spin-½ uses `S = 2`), optionally restricted to fixed per-species or total
particle-number sectors.

Internally this is a `base = 2` [`SpinlessFermionBasis`](@ref) over `M = S·L`
fermionic *modes* with a fixed **blocked** ordering

    mode(i, σ) = (σ − 1)·L + i        # all species-1 modes, then all species-2, …

so that same-species hopping `c†_{iσ} c_{i+1,σ}` is nearest-neighbour and the
on-site interaction `n_{i↑} n_{i↓}` carries no Jordan-Wigner string. The distinct
type stores the site count `L` and species count `S` so that
[`fermionmode`](@ref), the `(site, spin)` operator overloads, and
[`hubbard`](@ref) can translate physical labels into mode indices.

Because it is an [`AbstractOnsiteBasis`](@ref) (not a permutation/symmetry
basis), Jordan-Wigner-bearing operators (`c†`, `c`, hopping, pairing) are valid
on it via [`fermion_operator`](@ref).
"""
struct SpinfulFermionBasis{T <: Integer} <: AbstractOnsiteBasis
    dgt::Vector{T}
    I::Vector{T}
    B::T
    L::Int
    S::Int
end

"""
    SpinfulFermionBasis(dtype::DataType=Int64; L, S=2, N=nothing, f=nothing,
                        alloc=1000, threaded=true, small_N=true)

Construct a spinful / multi-species fermion occupation basis on `L` sites with
`S` species (default `S = 2`, i.e. spin-½).

Arguments:
- `dtype` : Integer type for stored representatives.
- `L`     : Number of lattice sites.
- `S`     : Number of fermion species.
- `N`     : (Optional) particle-number constraint:
    - `nothing` — all sectors (full `2^{S·L}` Fock space);
    - an `Integer` — fixed **total** particle number across all species;
    - an `NTuple{S}` / length-`S` `Vector` — fixed **per-species** numbers
      `(N₁, …, N_S)` (e.g. `(N↑, N↓)`), each in `0:L`.
- `f`     : (Optional) predicate `f(dgt) -> Bool` on the `M = S·L`-mode
            occupation string, further restricting which states are kept.
- `alloc`, `threaded`, `small_N` : forwarded to [`SpinlessFermionBasis`](@ref).
"""
function SpinfulFermionBasis(dtype::DataType=Int64;
    L::Integer, S::Integer=2, N=nothing, f=nothing, kwargs...,
)
    L ≥ 1 || error("L=$L must be ≥ 1.")
    S ≥ 1 || error("S=$S must be ≥ 1.")
    M = S * L
    totalN, Nvec = _spinful_sector(N, S, L)
    pred = _spinful_predicate(f, Nvec, L, S)
    sub = SpinlessFermionBasis(dtype; L=M, N=totalN, f=pred, kwargs...)
    SpinfulFermionBasis{dtype}(zeros(dtype, M), sub.I, sub.B, Int(L), Int(S))
end

# Interpret the `N` argument into (total particle number or nothing, per-species
# vector or nothing).
function _spinful_sector(N, S::Integer, L::Integer)
    if isnothing(N)
        return nothing, nothing
    elseif N isa Integer
        0 ≤ N ≤ S * L || error("Total N=$N is out of range 0:$(S*L) for $S species on $L sites.")
        return Int(N), nothing
    elseif N isa Union{Tuple, AbstractVector}
        length(N) == S || error("Per-species N must have one entry per species (expected $S, got $(length(N))).")
        all(n -> 0 ≤ n ≤ L, N) || error("Every per-species particle number must lie in 0:$L (got $N).")
        Nvec = Int[Int(n) for n in N]
        return sum(Nvec), Nvec
    else
        error("N must be `nothing`, an Integer (total), or a length-$S tuple/vector (per-species); got $(typeof(N)).")
    end
end

# Build the combined occupation predicate over the M-mode digit string, or
# `nothing` when neither a per-species constraint nor a user predicate applies.
function _spinful_predicate(f, Nvec, L::Integer, S::Integer)
    if isnothing(Nvec)
        return f
    elseif isnothing(f)
        return dgt -> _species_weights_ok(dgt, Nvec, L, S)
    else
        return dgt -> _species_weights_ok(dgt, Nvec, L, S) && f(dgt)
    end
end

# Each blocked species window [(σ-1)L+1 : σL] must have Hamming weight Nvec[σ].
function _species_weights_ok(dgt, Nvec, L::Integer, S::Integer)
    @inbounds for σ in 1:S
        c = 0
        base = (σ - 1) * L
        for i in 1:L
            c += dgt[base + i]
        end
        c == Nvec[σ] || return false
    end
    true
end

copy(b::SpinfulFermionBasis) = SpinfulFermionBasis(deepcopy(b.dgt), b.I, b.B, b.L, b.S)

# Resolve a species/spin label to its integer species index in 1:S.
function _species_index(spin, S::Integer)
    if spin isa Integer
        1 ≤ spin ≤ S || error("species $spin out of range 1:$S.")
        return Int(spin)
    elseif spin === :↑ || spin === :up
        return 1
    elseif spin === :↓ || spin === :down
        S ≥ 2 || error("spin label $spin requires S ≥ 2 (got S=$S).")
        return 2
    else
        error("Unrecognized species/spin label $spin; use an Integer in 1:$S, or :↑/:↓ (:up/:down) for spin-½.")
    end
end

"""
    fermionmode(B::SpinfulFermionBasis, site, spin) -> Int

Mode index of `(site, spin)` under the blocked ordering `mode = (σ−1)·L + site`,
where `σ` is the integer species index. `spin` may be an `Integer` species in
`1:S`, or for `S = 2` one of the symbols `:↑`/`:↓` (`:up`/`:down`).
"""
function fermionmode(B::SpinfulFermionBasis, site::Integer, spin)
    1 ≤ site ≤ B.L || error("site $site out of range 1:$(B.L).")
    σ = _species_index(spin, B.S)
    (σ - 1) * B.L + site
end
