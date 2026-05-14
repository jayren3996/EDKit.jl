#---------------------------------------------------------------------------------------------------
# Adaptive Krylov/Lanczos real-time evolution
#---------------------------------------------------------------------------------------------------
# This file adds a Hermitian real-time propagator based on the following idea:
#
#   Build a Lanczos basis anchored at a state, treat it as a reusable reduced
#   dynamical model, monitor its validity through the Lanczos boundary term, and
#   only extend or restart the basis when the monitor exceeds tolerance.
#
# The main public entry points are `timeevolve`, `timeevolve!`, and the explicit
# `KrylovEvolutionCache` workflow. Everything is matrix-free by default and
# works with EDKit `Operator`, `AbstractMatrix`, `SparseMatrixCSC`, or any other
# object that supports `LinearAlgebra.mul!(y, H, x)`.
#---------------------------------------------------------------------------------------------------

export timeevolve, timeevolve!, KrylovEvolutionCache, KrylovEvolutionDiagnostics

#---------------------------------------------------------------------------------------------------
# Diagnostics
#---------------------------------------------------------------------------------------------------
"""
    KrylovEvolutionDiagnostics

Mutable counter bag describing what the adaptive Krylov propagator did.

Fields:
- `basis_builds`       : how many times a fresh Lanczos basis was constructed
  (initial build + restarts).
- `basis_extensions`   : how many times an existing basis was extended in place.
- `restarts`           : how many times the propagator advanced its anchor and
  built a fresh basis from scratch.
- `matvecs`            : total number of Hamiltonian applications (`H * v`).
- `total_times_served` : number of requested output times returned.
- `max_dim_used`       : largest Lanczos dimension ever reached.
- `accepted_intervals` : list of validated interval lengths `τ_valid`.
"""
mutable struct KrylovEvolutionDiagnostics
    basis_builds::Int
    basis_extensions::Int
    restarts::Int
    matvecs::Int
    total_times_served::Int
    max_dim_used::Int
    accepted_intervals::Vector{Float64}
end

KrylovEvolutionDiagnostics() =
    KrylovEvolutionDiagnostics(0, 0, 0, 0, 0, 0, Float64[])

function Base.show(io::IO, d::KrylovEvolutionDiagnostics)
    println(io, "KrylovEvolutionDiagnostics:")
    println(io, "  basis_builds       = ", d.basis_builds)
    println(io, "  basis_extensions   = ", d.basis_extensions)
    println(io, "  restarts           = ", d.restarts)
    println(io, "  matvecs            = ", d.matvecs)
    println(io, "  total_times_served = ", d.total_times_served)
    println(io, "  max_dim_used       = ", d.max_dim_used)
    print(io,   "  accepted_intervals = ", d.accepted_intervals)
end

#---------------------------------------------------------------------------------------------------
# Cache
#---------------------------------------------------------------------------------------------------
"""
    KrylovEvolutionCache{TH, T, R}

Reusable workspace for adaptive Krylov/Lanczos real-time evolution.

A cache is tied to:
- a Hamiltonian `H`,
- an **anchor time** `t_anchor`,
- and an **anchor state** `ψ_anchor`.

It stores a Lanczos basis built from `ψ_anchor`, the projected tridiagonal
matrix `T_m`, and a cached eigendecomposition of `T_m` so that the reduced
propagator `exp(-i τ T_m) e_1` and the boundary monitor can be evaluated
cheaply for many values of `τ = t - t_anchor`. The full evolved state is
then `ψ(t) = ‖ψ_anchor‖ V_m exp(-i τ T_m) e_1`.

The cache is updated *in place*. Whenever the basis is no longer accurate for a
requested `τ`, it is either extended via additional Lanczos steps or restarted
from a new anchor state. See [`timeevolve`](@ref) for the high-level public API.
"""
mutable struct KrylovEvolutionCache{TH, T<:Number, R<:Real}
    # Hamiltonian and dimensions
    H::TH
    N::Int

    # Anchor state / time
    t_anchor::Float64
    t_latest::Float64
    ψ_anchor::Vector{T}
    norm_anchor::R

    # Lanczos storage
    V::Matrix{T}      # columns hold v_1 … v_{m+1}, last column is the frontier
    α::Vector{R}      # diagonal of T_m
    β::Vector{R}      # sub/super-diagonal; β[m] is the boundary coefficient
    m::Int            # current Lanczos dimension

    # Cached reduced-space eigendecomposition of T_m
    λ::Vector{R}      # eigenvalues of T_m (length m)
    Q::Matrix{R}      # eigenvectors of T_m (m × m)
    Qe1::Vector{R}    # first row of Q  (= Q' * e_1)
    Qem::Vector{R}    # last  row of Q  (= Q' * e_m)

    # Settings
    tol::R
    m_init::Int
    m_max::Int
    extend_step::Int
    nsample::Int
    bisect_iters::Int
    reuse_basis::Bool
    extend_basis::Bool
    normalize_output::Bool

    # Scratch
    w::Vector{T}
    reduced_phase::Vector{Complex{R}}
    reduced_coeffs::Vector{Complex{R}}

    # Diagnostics
    diagnostics::KrylovEvolutionDiagnostics
end

"""
    KrylovEvolutionCache(H, ψ0; kwargs...)

Construct a Krylov evolution cache anchored at `ψ0` with anchor time `0`, build
its initial Lanczos basis, and prepare the reduced-space propagator.

The propagator is **Hermitian-only** — `H` is assumed to satisfy `H' == H`.
This is not checked at runtime; passing a non-Hermitian operator will produce
meaningless output. If you need non-Hermitian generators, use a different
solver.

Supported Hamiltonian types: EDKit [`Operator`](@ref), `AbstractMatrix`,
`SparseMatrixCSC`, `Hermitian`, `Symmetric`, or anything for which
`LinearAlgebra.mul!(y, H, x)` is defined.

# Keyword arguments
- `tol`               : defect-monitor tolerance used to decide whether the
  current basis is still valid. Default `1e-10`.
- `m_init`            : initial Lanczos dimension. Default `30`.
- `m_max`             : largest Lanczos dimension the cache will ever reach
  before restarting from a new anchor. Default `60`.
- `extend_step`       : number of Lanczos vectors added per extension when
  `extend_basis = true`. Default `10`.
- `nsample`           : minimum number of sample points (including both
  endpoints) used to check the boundary monitor on a candidate interval.
  The effective sample count is automatically densified when the spectral
  spread of the projected tridiagonal matrix is large, to avoid missing
  narrow peaks. Default `13`.
- `bisect_iters`      : number of bisection refinements when shrinking a
  candidate interval. Default `12`.
- `reuse_basis`       : allow a single basis to serve multiple output times
  (this is the central optimization). Default `true`.
- `extend_basis`      : allow in-place extension of the current basis before
  falling back to a restart. Default `true`.
- `normalize_output`  : **opt-in post-processing** that renormalizes every
  reconstructed state vector. Default `false`. Setting this to `true` can
  mask accumulated propagation error by hiding norm drift — the accuracy
  guarantee of `timeevolve` is delivered by `tol`, not by renormalization.
  Use this flag only when you explicitly want a normalized output, e.g. for
  plotting, and never rely on it as a substitute for tightening `tol`.
"""
function KrylovEvolutionCache(H, ψ0::AbstractVector;
    tol::Real=1e-10,
    m_init::Integer=30,
    m_max::Integer=60,
    extend_step::Integer=10,
    nsample::Integer=13,
    bisect_iters::Integer=12,
    reuse_basis::Bool=true,
    extend_basis::Bool=true,
    normalize_output::Bool=false,
)
    N = length(ψ0)
    N > 0 || error("Initial state must be non-empty")
    m_init ≥ 1 || error("m_init must be ≥ 1")
    m_max  ≥ m_init || error("m_max must be ≥ m_init")
    m_max  ≤ N || (m_max = N)
    m_init ≤ m_max || (m_init = m_max)
    nsample ≥ 2 || error("nsample must be ≥ 2")
    bisect_iters ≥ 0 || error("bisect_iters must be ≥ 0")
    extend_step ≥ 1 || error("extend_step must be ≥ 1")

    # Promote state eltype so that e^{-iτH}ψ is representable.
    T = promote_type(ComplexF64, eltype(ψ0))
    R = real(T)

    ψ_anchor = convert(Vector{T}, collect(ψ0))
    nrm = norm(ψ_anchor)
    nrm > 0 || error("Initial state has zero norm")
    norm_anchor = R(nrm)

    V = zeros(T, N, m_max + 1)
    V[:, 1] .= ψ_anchor ./ norm_anchor
    α = zeros(R, m_max)
    β = zeros(R, m_max)
    w = zeros(T, N)
    reduced_phase = Vector{Complex{R}}(undef, m_max)
    reduced_coeffs = Vector{Complex{R}}(undef, m_max)

    cache = KrylovEvolutionCache{typeof(H), T, R}(
        H, N, 0.0, 0.0, ψ_anchor, norm_anchor,
        V, α, β, 0,
        R[], Matrix{R}(undef, 0, 0), R[], R[],
        R(tol), Int(m_init), Int(m_max), Int(extend_step),
        Int(nsample), Int(bisect_iters),
        reuse_basis, extend_basis, normalize_output,
        w, reduced_phase, reduced_coeffs,
        KrylovEvolutionDiagnostics(),
    )

    _lanczos_build!(cache, cache.m_init)
    _build_reduced!(cache)
    cache.diagnostics.basis_builds += 1
    cache.diagnostics.max_dim_used = max(cache.diagnostics.max_dim_used, cache.m)
    return cache
end

Base.size(cache::KrylovEvolutionCache) = (cache.N,)
Base.length(cache::KrylovEvolutionCache) = cache.N

function Base.show(io::IO, c::KrylovEvolutionCache{TH,T,R}) where {TH,T,R}
    print(io, "KrylovEvolutionCache{", TH, "} (N=", c.N,
          ", m=", c.m, "/", c.m_max,
          ", tol=", c.tol,
          ", t_anchor=", c.t_anchor, ")")
end

#---------------------------------------------------------------------------------------------------
# Hamiltonian application
#---------------------------------------------------------------------------------------------------
# `mul!(y, H, x)` is LinearAlgebra's overwriting semantics, but EDKit `Operator`
# multiplication accumulates into `y` instead of zeroing it first. This helper
# makes the overwrite contract explicit and works uniformly for `Operator`,
# `AbstractMatrix`, `SparseMatrixCSC`, `Hermitian`, etc.
@inline function _apply!(y::AbstractVector, H, x::AbstractVector)
    fill!(y, zero(eltype(y)))
    mul!(y, H, x)
    return y
end

#---------------------------------------------------------------------------------------------------
# Lanczos basis builder
#---------------------------------------------------------------------------------------------------
# Builds α[m_start+1 .. m_target] and β[m_start+1 .. m_target], filling V[:, j+1]
# with v_{j+1} as it goes. Uses full reorthogonalization (classical Gram–Schmidt
# applied twice) for numerical stability. Returns `true` on lucky breakdown.
#---------------------------------------------------------------------------------------------------
function _lanczos_build!(cache::KrylovEvolutionCache{TH,T,R}, m_target::Integer) where {TH,T,R}
    V, α, β, w = cache.V, cache.α, cache.β, cache.w
    H = cache.H
    m_target = min(Int(m_target), cache.m_max)
    m_start = cache.m

    for j in (m_start + 1):m_target
        vj = @view V[:, j]
        _apply!(w, H, vj)
        cache.diagnostics.matvecs += 1

        if j > 1
            vjm1 = @view V[:, j - 1]
            LinearAlgebra.axpy!(-β[j - 1], vjm1, w)
        end

        αj = real(dot(vj, w))
        α[j] = αj
        LinearAlgebra.axpy!(-αj, vj, w)

        # Full reorthogonalization against V[:, 1:j] (twice for stability).
        for _sweep in 1:2
            for i in 1:j
                vi = @view V[:, i]
                s = dot(vi, w)
                LinearAlgebra.axpy!(-s, vi, w)
            end
        end

        βj = norm(w)
        β[j] = R(βj)

        # Lucky breakdown: the Krylov space is invariant; the anchor state is
        # already exactly representable in the current subspace.
        if βj ≤ sqrt(eps(R)) * max(one(R), abs(αj))
            cache.m = j
            β[j] = zero(R)
            return true
        end

        cache.m = j
        if j + 1 ≤ size(V, 2)
            vjp1 = @view V[:, j + 1]
            invβ = one(T) / T(βj)
            @inbounds @simd for i in eachindex(w)
                vjp1[i] = w[i] * invβ
            end
        end
    end
    return false
end

#---------------------------------------------------------------------------------------------------
# Reduced-space propagator cache (eigendecomposition of the tridiagonal T_m)
#---------------------------------------------------------------------------------------------------
function _build_reduced!(cache::KrylovEvolutionCache{TH,T,R}) where {TH,T,R}
    m = cache.m
    m ≥ 1 || error("Cannot build reduced propagator for empty basis")
    diag   = cache.α[1:m]
    offdiag = m > 1 ? cache.β[1:(m - 1)] : R[]
    Tm = SymTridiagonal(diag, offdiag)
    F  = eigen(Tm)
    cache.λ   = F.values
    cache.Q   = F.vectors
    cache.Qe1 = cache.Q[1, :]
    cache.Qem = cache.Q[m, :]
    resize!(cache.reduced_phase, m)
    resize!(cache.reduced_coeffs, m)
    return cache
end

#---------------------------------------------------------------------------------------------------
# Cheap reduced-space operations
#---------------------------------------------------------------------------------------------------
# c(τ) = e^{-i τ T_m} e_1 = Q diag(e^{-iτλ}) (Q' e_1) = Q (e^{-iτλ} .* Qe1)
function _reduced_coeffs(cache::KrylovEvolutionCache{TH,T,R}, τ::Real) where {TH,T,R}
    m = cache.m
    phase = cache.reduced_phase
    @inbounds @simd for k in 1:m
        phase[k] = cis(-R(τ) * cache.λ[k]) * cache.Qe1[k]
    end
    c = cache.reduced_coeffs
    mul!(c, cache.Q, phase)
    return c
end

# η_m(τ) = |β_m · e_m' c(τ)| where e_m' c(τ) = Qem ⋅ (e^{-iτλ} .* Qe1)
function _defect_monitor(cache::KrylovEvolutionCache{TH,T,R}, τ::Real) where {TH,T,R}
    m = cache.m
    βm = cache.β[m]
    βm == 0 && return 0.0
    acc = zero(Complex{R})
    @inbounds @simd for k in 1:m
        acc += cache.Qem[k] * cis(-R(τ) * cache.λ[k]) * cache.Qe1[k]
    end
    return Float64(abs(βm) * abs(acc))
end

function _reconstruct_state!(ψ_out::AbstractVector, cache::KrylovEvolutionCache{TH,T,R}, τ::Real) where {TH,T,R}
    c = _reduced_coeffs(cache, τ)
    V = @view cache.V[:, 1:cache.m]
    mul!(ψ_out, V, c)
    @. ψ_out *= cache.norm_anchor
    if cache.normalize_output
        nrm = norm(ψ_out)
        nrm > 0 && (@. ψ_out /= nrm)
    end
    return ψ_out
end

#---------------------------------------------------------------------------------------------------
# Adaptive interval selector
#---------------------------------------------------------------------------------------------------
# The defect monitor η(τ) = |β_m · Σ_k Qem[k]·Qe1[k]·exp(-i τ λ_k)| is a
# trigonometric sum with real frequencies {λ_k}, the eigenvalues of the
# projected tridiagonal T_m. Its amplitude is bounded by |β_m| (Cauchy–Schwarz
# on the orthonormal columns of Q), and its oscillation frequency is bounded
# by the spectral spread ω = λ_max − λ_min of T_m. Sampling at spacing
# Δτ ≲ π / ω therefore resolves every oscillation up to the Nyquist rate.
#
# The selector:
#   1. Takes a free early-exit when |β_m| ≤ tol: the monitor cannot exceed tol
#      anywhere, so the full candidate interval is accepted.
#   2. Picks an effective sample count n_eff that is at least the user-supplied
#      `nsample` but also satisfies the Nyquist-style density condition above,
#      with a safety factor of 4 and a hard upper cap to avoid runaway cost.
#   3. Scans the monitor at those samples, finds the first failure, bisects
#      that sub-interval to refine, and returns the largest accepted prefix.
#---------------------------------------------------------------------------------------------------
function _effective_nsample(cache::KrylovEvolutionCache, τ_try::Float64)
    base_n = cache.nsample
    isempty(cache.λ) && return base_n
    m = cache.m
    m ≤ 1 && return base_n
    ω = Float64(cache.λ[end] - cache.λ[1])           # spectral spread of T_m
    ω ≤ 0 && return base_n
    # Nyquist-style: need at least ceil(ω·τ_try / π) + 1 samples; we use a
    # safety factor of 4 for conservative error control, then cap the result.
    n_nyq = ceil(Int, 4 * ω * τ_try / π) + 1
    return clamp(max(base_n, n_nyq), base_n, 513)
end

function _max_valid_interval(cache::KrylovEvolutionCache, τ_try::Real)
    τ_try > 0 || return 0.0
    tol = cache.tol

    # Early exit: the monitor is bounded above by |β_m|. If the Lanczos
    # boundary coefficient itself is already below tolerance, the basis is
    # valid for every τ ≥ 0 and the full candidate interval is safe.
    if abs(cache.β[cache.m]) ≤ tol
        return Float64(τ_try)
    end

    τ_try_f = Float64(τ_try)
    n_eff = _effective_nsample(cache, τ_try_f)
    xs = range(0.0, τ_try_f, length=n_eff)

    for k in 1:n_eff
        η = _defect_monitor(cache, xs[k])
        if !isfinite(η) || η > tol
            k == 1 && return 0.0
            τ_lo = Float64(xs[k - 1])
            τ_hi = Float64(xs[k])
            for _ in 1:cache.bisect_iters
                τ_mid = 0.5 * (τ_lo + τ_hi)
                if _defect_monitor(cache, τ_mid) > tol
                    τ_hi = τ_mid
                else
                    τ_lo = τ_mid
                end
            end
            return τ_lo
        end
    end
    return τ_try_f
end

#---------------------------------------------------------------------------------------------------
# Basis extension and restart
#---------------------------------------------------------------------------------------------------
function _extend_basis!(cache::KrylovEvolutionCache)
    cache.extend_basis || return false
    cache.m < cache.m_max || return false
    old_m = cache.m
    new_m = min(cache.m + cache.extend_step, cache.m_max)
    _lanczos_build!(cache, new_m)
    if cache.m > old_m
        _build_reduced!(cache)
        cache.diagnostics.basis_extensions += 1
        cache.diagnostics.max_dim_used = max(cache.diagnostics.max_dim_used, cache.m)
        return true
    end
    return false
end

function _restart_from_new_anchor!(cache::KrylovEvolutionCache{TH,T,R}, τ_valid::Real) where {TH,T,R}
    # 1. compute current evolved state at τ_valid from existing basis
    ψ_new = similar(cache.ψ_anchor)
    _reconstruct_state!(ψ_new, cache, τ_valid)

    # 2. promote it to the new anchor
    cache.t_anchor += Float64(τ_valid)
    copyto!(cache.ψ_anchor, ψ_new)
    nrm = norm(cache.ψ_anchor)
    nrm > 0 || error("Restart failed: anchor state has zero norm")
    cache.norm_anchor = R(nrm)

    # 3. reset Lanczos storage for a fresh build
    fill!(cache.α, zero(R))
    fill!(cache.β, zero(R))
    cache.m = 0
    V = cache.V
    @inbounds @simd for i in 1:cache.N
        V[i, 1] = cache.ψ_anchor[i] / cache.norm_anchor
    end

    # 4. build a fresh Lanczos basis anchored at the new state
    _lanczos_build!(cache, cache.m_init)
    _build_reduced!(cache)
    cache.diagnostics.basis_builds += 1
    cache.diagnostics.restarts += 1
    cache.diagnostics.max_dim_used = max(cache.diagnostics.max_dim_used, cache.m)
    return cache
end

#---------------------------------------------------------------------------------------------------
# Main propagation loop
#---------------------------------------------------------------------------------------------------
# Advance the cache so that its current basis remains valid for every τ in
# [0, τ_target] *relative to the current anchor*. Any required extensions or
# restarts are performed in place. After this returns, the caller may safely
# reconstruct the state at any τ ∈ [0, τ_target] using `_reconstruct_state!`.
function _ensure_valid_up_to!(cache::KrylovEvolutionCache, τ_target::Real)
    τ_remaining = Float64(τ_target)
    τ_remaining ≤ 0 && return cache

    while true
        τ_valid = _max_valid_interval(cache, τ_remaining)

        if τ_valid ≥ τ_remaining
            return cache
        end

        # Not enough — try to extend the basis first.
        if cache.extend_basis && cache.m < cache.m_max
            if _extend_basis!(cache)
                continue
            end
        end

        # If the current basis is not valid at all (τ_valid ≈ 0), we must not
        # lose progress: try to extend, otherwise error.
        if τ_valid ≤ 0
            if cache.m < cache.m_max && _extend_basis!(cache)
                continue
            end
            error("KrylovEvolutionCache: basis at m=$(cache.m) is not valid at τ=0; " *
                  "consider increasing m_max or loosening tol")
        end

        # Otherwise: accept the valid prefix, restart from the new anchor, and
        # continue consuming the remaining interval.
        push!(cache.diagnostics.accepted_intervals, τ_valid)
        _restart_from_new_anchor!(cache, τ_valid)
        τ_remaining -= τ_valid
    end
end

#---------------------------------------------------------------------------------------------------
# Public API
#---------------------------------------------------------------------------------------------------
"""
    timeevolve(H, ψ0, t::Real; kwargs...) -> Vector
    timeevolve(H, ψ0, ts::AbstractVector; kwargs...) -> Matrix

Adaptive Krylov/Lanczos real-time evolution of `ψ0` under the Hermitian
Hamiltonian `H`.

`H` can be an EDKit [`Operator`](@ref), an `AbstractMatrix`, a `SparseMatrixCSC`,
or any object that supports `LinearAlgebra.mul!(y, H, x)`. The propagator is
**Hermitian-only** and **forward in time only**: the single-time form returns
`exp(-i t H) ψ0` for `t ≥ 0`, and a negative `t` is rejected. The multi-time
form returns a matrix whose `k`-th column is the state at `ts[k]`, in the
order the caller supplied; `ts` is sorted internally before propagation and
the output columns are restored to the input order, so passing an unsorted
`ts` is fine. All times in `ts` must be non-negative. One Lanczos basis is
reused across as many requested times as possible before a restart happens.

Keyword arguments are forwarded to [`KrylovEvolutionCache`](@ref). The
additional keyword `return_diagnostics=true` (multi-time form only) causes the
function to return a tuple `(states, diagnostics)`.

# Example
```julia
L = 10
H = trans_inv_operator(spin((1.0,"xx"),(1.0,"yy"),(1.0,"zz")), 1:2, TensorBasis(L=L, base=2))
ψ0 = productstate([iseven(i) ? 0 : 1 for i in 1:L], H.B) .+ 0im
ψt = timeevolve(H, ψ0, 2.0; tol=1e-10)
ts = range(0, 5; length=51)
ψs = timeevolve(H, ψ0, ts; tol=1e-10)
```
"""
function timeevolve(H, ψ0::AbstractVector, t::Real;
                    return_diagnostics::Bool=false, kwargs...)
    cache = KrylovEvolutionCache(H, ψ0; kwargs...)
    ψt = timeevolve!(cache, t)
    return return_diagnostics ? (ψt, cache.diagnostics) : ψt
end

function timeevolve(H, ψ0::AbstractVector, ts::AbstractVector{<:Real};
                    return_diagnostics::Bool=false, kwargs...)
    any(<(0), ts) && error("timeevolve: all times in ts must be non-negative; backward evolution is not supported")
    cache = KrylovEvolutionCache(H, ψ0; kwargs...)
    T = eltype(cache.ψ_anchor)
    perm = sortperm(ts)
    ts_sorted = collect(Float64.(ts[perm]))
    out_sorted = Matrix{T}(undef, cache.N, length(ts))
    timeevolve!(out_sorted, cache, ts_sorted)
    out = similar(out_sorted)
    @inbounds for i in eachindex(perm)
        out[:, perm[i]] = view(out_sorted, :, i)
    end
    return return_diagnostics ? (out, cache.diagnostics) : out
end

"""
    timeevolve!(cache::KrylovEvolutionCache, t::Real) -> Vector

Advance an existing cache to time `t` (which must be ≥ `cache.t_anchor`) and
return a freshly allocated state vector at that time. The cache's anchor may be
moved forward as a side effect when the current basis cannot cover `t` alone.
"""
function timeevolve!(cache::KrylovEvolutionCache, t::Real)
    tf = Float64(t)
    tf < cache.t_latest && error("timeevolve!: requested time $t precedes the last served time $(cache.t_latest); backward evolution is not supported")
    τ = tf - cache.t_anchor
    τ < 0 && error("timeevolve!: requested time $t precedes anchor time $(cache.t_anchor)")
    _ensure_valid_up_to!(cache, τ)
    ψ_out = similar(cache.ψ_anchor)
    _reconstruct_state!(ψ_out, cache, tf - cache.t_anchor)
    cache.t_latest = tf
    cache.diagnostics.total_times_served += 1
    return ψ_out
end

"""
    timeevolve!(out::AbstractMatrix, cache::KrylovEvolutionCache, ts::AbstractVector{<:Real}) -> out

Fill the columns of `out` with the evolved states at the times `ts`, reusing
the cache's current Lanczos basis for as many requested times as the defect
monitor allows before extension or restart.

`ts` **must be sorted in strictly forward order** (i.e. non-decreasing, and
every entry must be ≥ the cache's last served time). This form deliberately
does not sort internally: the cache is stateful, and reordering the time
stream would either advance the anchor past intermediate requests or require
the caller to know about the internal ordering. Use the stateless
`timeevolve(H, ψ0, ts)` form if you want a convenience sort.
"""
function timeevolve!(out::AbstractMatrix, cache::KrylovEvolutionCache,
                     ts::AbstractVector{<:Real})
    size(out, 1) == cache.N || error("timeevolve!: output has $(size(out,1)) rows, expected $(cache.N)")
    size(out, 2) == length(ts) || error("timeevolve!: output has $(size(out,2)) columns, expected $(length(ts))")
    issorted(ts) || error("timeevolve!: ts must be sorted in ascending order")

    for i in eachindex(ts)
        tf = Float64(ts[i])
        tf < cache.t_latest && error("timeevolve!: requested time $(ts[i]) precedes the last served time $(cache.t_latest); backward evolution is not supported")
        τ = tf - cache.t_anchor
        τ < 0 && error("timeevolve!: requested time $(ts[i]) precedes anchor time $(cache.t_anchor)")
        _ensure_valid_up_to!(cache, τ)
        _reconstruct_state!(view(out, :, i), cache, tf - cache.t_anchor)
        cache.t_latest = tf
        cache.diagnostics.total_times_served += 1
    end
    return out
end

"""
    timeevolve!(out::AbstractVector, H, ψ0, t::Real; kwargs...) -> out

In-place single-time variant. Creates a fresh cache, evolves to `t`, and writes
the resulting state into `out`.
"""
function timeevolve!(out::AbstractVector, H, ψ0::AbstractVector, t::Real; kwargs...)
    length(out) == length(ψ0) || error("timeevolve!: output length mismatch")
    cache = KrylovEvolutionCache(H, ψ0; kwargs...)
    tf = Float64(t)
    τ = tf - cache.t_anchor
    τ < 0 && error("timeevolve!: requested time $t precedes anchor time $(cache.t_anchor)")
    _ensure_valid_up_to!(cache, τ)
    _reconstruct_state!(out, cache, tf - cache.t_anchor)
    cache.t_latest = tf
    cache.diagnostics.total_times_served += 1
    return out
end
