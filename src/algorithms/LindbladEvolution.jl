#---------------------------------------------------------------------------------------------------
# Adaptive Arnoldi real-time evolution for many-body Lindblad dynamics
#---------------------------------------------------------------------------------------------------

export LiouvillianMap, lindblad_timeevolve, lindblad_timeevolve!,
       LindbladArnoldiCache, LindbladArnoldiDiagnostics

"""
    LindbladArnoldiDiagnostics

Diagnostics collected by the adaptive Arnoldi Lindblad propagator.

The vectors `trace_drifts`, `hermiticity_drifts`, and `min_hermitian_eigs`
record raw approximation-quality diagnostics at each returned output time.
They are diagnostics, not guarantees of physicality.
"""
mutable struct LindbladArnoldiDiagnostics
    basis_builds::Int
    basis_extensions::Int
    restarts::Int
    liouvillian_applies::Int
    total_times_served::Int
    max_dim_used::Int
    accepted_intervals::Vector{Float64}
    trace_drifts::Vector{Float64}
    hermiticity_drifts::Vector{Float64}
    min_hermitian_eigs::Vector{Float64}
    max_trace_drift::Float64
    max_hermiticity_drift::Float64
    min_min_hermitian_eig::Float64
end

LindbladArnoldiDiagnostics() = LindbladArnoldiDiagnostics(
    0, 0, 0, 0, 0, 0,
    Float64[], Float64[], Float64[], Float64[],
    0.0, 0.0, Inf,
)

function Base.show(io::IO, d::LindbladArnoldiDiagnostics)
    println(io, "LindbladArnoldiDiagnostics:")
    println(io, "  basis_builds         = ", d.basis_builds)
    println(io, "  basis_extensions     = ", d.basis_extensions)
    println(io, "  restarts             = ", d.restarts)
    println(io, "  liouvillian_applies  = ", d.liouvillian_applies)
    println(io, "  total_times_served   = ", d.total_times_served)
    println(io, "  max_dim_used         = ", d.max_dim_used)
    println(io, "  max_trace_drift      = ", d.max_trace_drift)
    println(io, "  max_hermiticity_drift= ", d.max_hermiticity_drift)
    println(io, "  min_min_hermitian_eig= ", d.min_min_hermitian_eig)
    print(io,   "  accepted_intervals   = ", d.accepted_intervals)
end

"""
    LiouvillianMap(H, jumps)
    LiouvillianMap(lb::Lindblad)

Matrix-free many-body Lindbladian acting on `vec(ρ)` while applying the
generator internally on reshaped density-matrix views.

This is the low-level linear map consumed by [`LindbladArnoldiCache`](@ref).
"""
struct LiouvillianMap{TH, TL, TD, T <: Number}
    H::TH
    L::Vector{TL}
    D::TD
    dim::Int
end

eltype(::LiouvillianMap{TH, TL, TD, T}) where {TH, TL, TD, T} = T
Base.size(A::LiouvillianMap) = (A.dim^2, A.dim^2)
Base.size(A::LiouvillianMap, i::Integer) = i <= 2 ? A.dim^2 : 1
Base.length(A::LiouvillianMap) = A.dim^2

function Base.show(io::IO, A::LiouvillianMap)
    print(io, "LiouvillianMap(dim=", A.dim, ", jumps=", length(A.L), ")")
end

function LiouvillianMap(H::AbstractMatrix, jumps::AbstractVector{<:AbstractMatrix})
    size(H, 1) == size(H, 2) || error("LiouvillianMap: H must be square")
    d = size(H, 1)
    for (j, L) in enumerate(jumps)
        size(L, 1) == d && size(L, 2) == d || error("LiouvillianMap: jump $j has incompatible size")
    end

    T = promote_type(ComplexF64, eltype(H), map(eltype, jumps)...)
    D = if isempty(jumps)
        issparse(H) ? spzeros(T, d, d) : zeros(T, d, d)
    else
        mapreduce(L -> adjoint(L) * L, +, jumps)
    end
    LiouvillianMap{typeof(H), eltype(jumps), typeof(D), T}(H, collect(jumps), D, d)
end

LiouvillianMap(lb::Lindblad) = LiouvillianMap(lb.H, lb.L)

mutable struct _LiouvillianWorkspace{T}
    X::Matrix{T}
    tmp1::Matrix{T}
    tmp2::Matrix{T}
end

function _LiouvillianWorkspace(::Type{T}, d::Integer) where {T}
    _LiouvillianWorkspace(zeros(T, d, d), zeros(T, d, d), zeros(T, d, d))
end

@inline function _row_times_vector(row::AbstractVector, vec::AbstractVector)
    acc = zero(promote_type(eltype(row), eltype(vec)))
    @inbounds @simd for i in eachindex(row, vec)
        acc += row[i] * vec[i]
    end
    acc
end

function _liouvillian_apply!(y::AbstractVector, A::LiouvillianMap{TH, TL, TD, T},
                             x::AbstractVector, ws::_LiouvillianWorkspace{T}) where {TH, TL, TD, T}
    X, tmp1, tmp2 = ws.X, ws.tmp1, ws.tmp2
    copyto!(vec(X), x)
    Y = reshape(y, A.dim, A.dim)
    fill!(Y, zero(T))

    mul!(tmp1, A.H, X)
    mul!(tmp2, X, A.H)
    @. Y = -1im * (tmp1 - tmp2)

    for L in A.L
        mul!(tmp1, L, X)
        mul!(tmp2, tmp1, adjoint(L))
        @. Y += tmp2
    end

    mul!(tmp1, A.D, X)
    @. Y += -0.5 * tmp1
    mul!(tmp2, X, A.D)
    @. Y += -0.5 * tmp2
    return y
end

function Base.:*(A::LiouvillianMap, ρ::AbstractMatrix)
    T = promote_type(eltype(A), eltype(ρ))
    ws = _LiouvillianWorkspace(T, A.dim)
    out = zeros(T, A.dim^2)
    _liouvillian_apply!(out, A, vec(Matrix{T}(ρ)), ws)
    reshape(out, A.dim, A.dim)
end

function Base.:*(A::LiouvillianMap, v::AbstractVector)
    T = promote_type(eltype(A), eltype(v))
    ws = _LiouvillianWorkspace(T, A.dim)
    out = zeros(T, length(v))
    _liouvillian_apply!(out, A, Vector{T}(v), ws)
end

"""
    LindbladArnoldiCache(A, dm::DensityMatrix; kwargs...)
    LindbladArnoldiCache(lb, dm::DensityMatrix; kwargs...)

Adaptive Arnoldi workspace for many-body Lindblad time evolution.

The cache is anchored at a density matrix and time, stores an Arnoldi basis for
the vectorized density matrix, and reuses a cached complex-Schur reduced model
until the Arnoldi defect monitor exceeds tolerance.

# Keyword arguments
- `tol`: defect-monitor tolerance. Default `1e-10`.
- `m_init`: initial Arnoldi dimension. Default `20`.
- `m_max`: largest Arnoldi dimension before restart. Default `50`.
- `extend_step`: Arnoldi vectors added per extension. Default `10`.
- `nsample`: baseline interval-monitor sample count. Default `17`.
- `bisect_iters`: bisection refinements after the first failed sample. Default `12`.
- `extend_basis`: allow in-place Arnoldi basis extension before restart. Default `true`.
- `normalize_trace`: opt-in post-processing that rescales returned density
  matrices to the initial trace. Default `false`.
- `record_min_eig`: also record the minimum eigenvalue of the Hermitian part of
  the returned density matrix. Default `false`.
"""
mutable struct LindbladArnoldiCache{TA, T <: Number, R <: Real}
    A::TA
    dim::Int
    N::Int

    t_anchor::Float64
    t_latest::Float64
    ρ_anchor::Matrix{T}
    norm_anchor::R
    trace_target::T

    V::Matrix{T}
    Hred::Matrix{T}
    m::Int

    schur_T::Matrix{T}
    schur_Z::Matrix{T}
    schur_q1::Vector{T}
    schur_rowm::Vector{T}

    tol::R
    m_init::Int
    m_max::Int
    extend_step::Int
    nsample::Int
    bisect_iters::Int
    extend_basis::Bool
    normalize_trace::Bool
    record_min_eig::Bool

    arnoldi_w::Vector{T}
    reduced_y::Vector{T}
    reduced_c::Vector{T}
    full_vec::Vector{T}
    step_cache::Dict{UInt64, Matrix{T}}
    liou_ws::_LiouvillianWorkspace{T}

    diagnostics::LindbladArnoldiDiagnostics
end

function LindbladArnoldiCache(
    A::Union{Lindblad, LiouvillianMap},
    dm::DensityMatrix;
    tol::Real=1e-10,
    m_init::Integer=20,
    m_max::Integer=50,
    extend_step::Integer=10,
    nsample::Integer=17,
    bisect_iters::Integer=12,
    extend_basis::Bool=true,
    normalize_trace::Bool=false,
    record_min_eig::Bool=false,
)
    map = A isa Lindblad ? LiouvillianMap(A) : A
    d = map.dim
    size(dm.ρ, 1) == d && size(dm.ρ, 2) == d || error("LindbladArnoldiCache: density matrix has incompatible size")
    d > 0 || error("LindbladArnoldiCache: density matrix must be non-empty")
    m_init ≥ 1 || error("m_init must be ≥ 1")
    m_max ≥ m_init || error("m_max must be ≥ m_init")
    extend_step ≥ 1 || error("extend_step must be ≥ 1")
    nsample ≥ 2 || error("nsample must be ≥ 2")
    bisect_iters ≥ 0 || error("bisect_iters must be ≥ 0")

    T = promote_type(ComplexF64, eltype(map), eltype(dm.ρ))
    R = real(T)
    N = d^2
    m_max = min(Int(m_max), N)
    m_init = min(Int(m_init), m_max)

    ρ_anchor = Matrix{T}(dm.ρ)
    nrm = norm(ρ_anchor)
    nrm > 0 || error("LindbladArnoldiCache: initial density matrix has zero Frobenius norm")

    V = zeros(T, N, m_max + 1)
    V[:, 1] .= vec(ρ_anchor) ./ nrm
    Hred = zeros(T, m_max + 1, m_max)

    cache = LindbladArnoldiCache{typeof(map), T, R}(
        map, d, N,
        0.0, 0.0, ρ_anchor, R(nrm), tr(ρ_anchor),
        V, Hred, 0,
        Matrix{T}(undef, 0, 0), Matrix{T}(undef, 0, 0), T[], T[],
        R(tol), Int(m_init), Int(m_max), Int(extend_step), Int(nsample), Int(bisect_iters),
        extend_basis, normalize_trace, record_min_eig,
        zeros(T, N), zeros(T, m_max), zeros(T, m_max), zeros(T, N),
        Dict{UInt64, Matrix{T}}(), _LiouvillianWorkspace(T, d),
        LindbladArnoldiDiagnostics(),
    )

    _arnoldi_build!(cache, cache.m_init)
    _build_reduced!(cache)
    cache.diagnostics.basis_builds += 1
    cache.diagnostics.max_dim_used = max(cache.diagnostics.max_dim_used, cache.m)
    return cache
end

LindbladArnoldiCache(A::Union{Lindblad, LiouvillianMap}, ρ::AbstractMatrix; kwargs...) =
    LindbladArnoldiCache(A, densitymatrix(ρ); kwargs...)

function Base.show(io::IO, cache::LindbladArnoldiCache{TA, T, R}) where {TA, T, R}
    print(io, "LindbladArnoldiCache{", TA, "} (dim=", cache.dim,
          ", m=", cache.m, "/", cache.m_max,
          ", tol=", cache.tol,
          ", t_anchor=", cache.t_anchor, ")")
end

function _arnoldi_build!(cache::LindbladArnoldiCache{TA, T, R}, m_target::Integer) where {TA, T, R}
    V, Hred, w = cache.V, cache.Hred, cache.arnoldi_w
    A, ws = cache.A, cache.liou_ws
    m_target = min(Int(m_target), cache.m_max)
    m_start = cache.m

    for j in (m_start + 1):m_target
        vj = @view V[:, j]
        _liouvillian_apply!(w, A, vj, ws)
        cache.diagnostics.liouvillian_applies += 1

        fill!(@view(Hred[:, j]), zero(T))
        for sweep in 1:2
            for i in 1:j
                vi = @view V[:, i]
                hij = dot(vi, w)
                Hred[i, j] += hij
                @inbounds @simd for k in eachindex(w)
                    w[k] -= hij * vi[k]
                end
            end
        end

        hnext = norm(w)
        Hred[j + 1, j] = T(hnext)
        cache.m = j

        if hnext ≤ sqrt(eps(R)) * max(one(R), norm(@view(Hred[1:j, j])))
            Hred[j + 1, j] = zero(T)
            return true
        end

        vjp1 = @view V[:, j + 1]
        invh = inv(T(hnext))
        @inbounds @simd for k in eachindex(w)
            vjp1[k] = w[k] * invh
        end
    end
    return false
end

function _build_reduced!(cache::LindbladArnoldiCache{TA, T}) where {TA, T}
    m = cache.m
    m ≥ 1 || error("Cannot build reduced Arnoldi model for empty basis")
    Hm = Matrix(@view(cache.Hred[1:m, 1:m]))
    F = schur(Hm)
    cache.schur_T = Matrix(F.T)
    cache.schur_Z = Matrix(F.Z)
    cache.schur_q1 = conj.(cache.schur_Z[1, :])
    cache.schur_rowm = vec(cache.schur_Z[m, :])
    empty!(cache.step_cache)
    return cache
end

@inline _boundary_coeff(cache::LindbladArnoldiCache) = cache.m < 1 ? zero(eltype(cache.Hred)) : cache.Hred[cache.m + 1, cache.m]

function _reduced_step_matrix(cache::LindbladArnoldiCache{TA, T}, τ::Real) where {TA, T}
    τf = Float64(τ)
    key = reinterpret(UInt64, τf)
    get!(cache.step_cache, key) do
        Matrix(exp(UpperTriangular(cache.schur_T .* T(τf))))
    end
end

function _reduced_y_at!(cache::LindbladArnoldiCache{TA, T}, τ::Real) where {TA, T}
    m = cache.m
    y = @view cache.reduced_y[1:m]
    if τ == 0
        copyto!(y, cache.schur_q1)
        return y
    end
    Eτ = _reduced_step_matrix(cache, τ)
    mul!(y, Eτ, cache.schur_q1)
    return y
end

function _reduced_coeffs_at!(cache::LindbladArnoldiCache{TA, T}, τ::Real) where {TA, T}
    m = cache.m
    y = _reduced_y_at!(cache, τ)
    c = @view cache.reduced_c[1:m]
    mul!(c, cache.schur_Z, y)
    return c
end

function _defect_from_y(cache::LindbladArnoldiCache, y::AbstractVector)
    hnext = _boundary_coeff(cache)
    iszero(hnext) && return 0.0
    val = cache.norm_anchor * abs(hnext) * abs(_row_times_vector(cache.schur_rowm, y))
    Float64(val)
end

function _defect_monitor(cache::LindbladArnoldiCache, τ::Real)
    y = _reduced_y_at!(cache, τ)
    _defect_from_y(cache, y)
end

function _effective_nsample(cache::LindbladArnoldiCache, τ_try::Float64)
    base_n = cache.nsample
    cache.m ≤ 1 && return base_n
    γ = opnorm(@view(cache.Hred[1:cache.m, 1:cache.m]), 1)
    γ ≤ 0 && return base_n
    n_eff = ceil(Int, 4 * γ * τ_try) + 1
    clamp(max(base_n, n_eff), base_n, 513)
end

function _max_valid_interval(cache::LindbladArnoldiCache, τ_try::Real)
    τ_try > 0 || return 0.0
    iszero(_boundary_coeff(cache)) && return Float64(τ_try)

    τf = Float64(τ_try)
    n_eff = _effective_nsample(cache, τf)
    xs = range(0.0, τf, length=n_eff)
    y = @view cache.reduced_y[1:cache.m]
    tmp = @view cache.reduced_c[1:cache.m]
    copyto!(y, cache.schur_q1)

    η0 = _defect_from_y(cache, y)
    (!isfinite(η0) || η0 > cache.tol) && return 0.0

    step_mat = n_eff > 1 ? _reduced_step_matrix(cache, step(xs)) : Matrix{eltype(cache.schur_T)}(I, cache.m, cache.m)
    for k in 2:n_eff
        mul!(tmp, step_mat, y)
        copyto!(y, tmp)
        η = _defect_from_y(cache, y)
        if !isfinite(η) || η > cache.tol
            τ_lo = Float64(xs[k - 1])
            τ_hi = Float64(xs[k])
            for _ in 1:cache.bisect_iters
                τ_mid = 0.5 * (τ_lo + τ_hi)
                if _defect_monitor(cache, τ_mid) > cache.tol
                    τ_hi = τ_mid
                else
                    τ_lo = τ_mid
                end
            end
            return τ_lo
        end
    end
    return τf
end

function _reconstruct_density_raw!(ρ_out::AbstractMatrix, cache::LindbladArnoldiCache, τ::Real)
    c = _reduced_coeffs_at!(cache, τ)
    V = @view cache.V[:, 1:cache.m]
    mul!(cache.full_vec, V, c)
    @. cache.full_vec *= cache.norm_anchor
    copyto!(vec(ρ_out), cache.full_vec)
    return ρ_out
end

function _record_output_diagnostics!(cache::LindbladArnoldiCache, ρ_raw::AbstractMatrix)
    diag = cache.diagnostics
    trace_drift = Float64(abs(tr(ρ_raw) - cache.trace_target))
    herm_drift = Float64(norm(ρ_raw - adjoint(ρ_raw)))
    push!(diag.trace_drifts, trace_drift)
    push!(diag.hermiticity_drifts, herm_drift)
    diag.max_trace_drift = max(diag.max_trace_drift, trace_drift)
    diag.max_hermiticity_drift = max(diag.max_hermiticity_drift, herm_drift)

    if cache.record_min_eig
        ρh = (ρ_raw + adjoint(ρ_raw)) / 2
        λmin = eigmin(Hermitian(ρh))
        λminf = Float64(real(λmin))
        push!(diag.min_hermitian_eigs, λminf)
        diag.min_min_hermitian_eig = min(diag.min_min_hermitian_eig, λminf)
    end
    return diag
end

function _copy_output_density(cache::LindbladArnoldiCache, ρ_raw::AbstractMatrix)
    ρ = Matrix(ρ_raw)
    if cache.normalize_trace
        trρ = tr(ρ)
        !iszero(trρ) && (@. ρ *= cache.trace_target / trρ)
    end
    DensityMatrix(ρ)
end

function _extend_basis!(cache::LindbladArnoldiCache)
    cache.extend_basis || return false
    cache.m < cache.m_max || return false
    old_m = cache.m
    new_m = min(cache.m + cache.extend_step, cache.m_max)
    _arnoldi_build!(cache, new_m)
    if cache.m > old_m
        _build_reduced!(cache)
        cache.diagnostics.basis_extensions += 1
        cache.diagnostics.max_dim_used = max(cache.diagnostics.max_dim_used, cache.m)
        return true
    end
    return false
end

function _restart_from_new_anchor!(cache::LindbladArnoldiCache{TA, T, R}, τ_valid::Real) where {TA, T, R}
    ρ_new = similar(cache.ρ_anchor)
    _reconstruct_density_raw!(ρ_new, cache, τ_valid)

    cache.t_anchor += Float64(τ_valid)
    copyto!(cache.ρ_anchor, ρ_new)
    nrm = norm(cache.ρ_anchor)
    nrm > 0 || error("LindbladArnoldiCache: restart failed because anchor density matrix has zero Frobenius norm")
    cache.norm_anchor = R(nrm)

    fill!(cache.Hred, zero(T))
    cache.m = 0
    cache.V[:, 1] .= vec(cache.ρ_anchor) ./ cache.norm_anchor
    _arnoldi_build!(cache, cache.m_init)
    _build_reduced!(cache)
    cache.diagnostics.basis_builds += 1
    cache.diagnostics.restarts += 1
    cache.diagnostics.max_dim_used = max(cache.diagnostics.max_dim_used, cache.m)
    return cache
end

function _ensure_valid_up_to!(cache::LindbladArnoldiCache, τ_target::Real)
    τ_remaining = Float64(τ_target)
    τ_remaining ≤ 0 && return cache

    while true
        τ_valid = _max_valid_interval(cache, τ_remaining)
        if τ_valid ≥ τ_remaining
            return cache
        end

        if cache.extend_basis && cache.m < cache.m_max && _extend_basis!(cache)
            continue
        end

        if τ_valid ≤ 0
            if cache.m < cache.m_max && _extend_basis!(cache)
                continue
            end
            error("LindbladArnoldiCache: basis at m=$(cache.m) is not valid at τ=0; consider increasing m_max or loosening tol")
        end

        push!(cache.diagnostics.accepted_intervals, τ_valid)
        _restart_from_new_anchor!(cache, τ_valid)
        τ_remaining -= τ_valid
    end
end

"""
    lindblad_timeevolve(A, dm, t::Real; kwargs...) -> DensityMatrix
    lindblad_timeevolve(A, dm, ts::AbstractVector; kwargs...) -> Vector{DensityMatrix}

Adaptive Arnoldi propagation for many-body Lindblad dynamics.

The public API accepts either a legacy dense [`Lindblad`](@ref) object or a
[`LiouvillianMap`](@ref). The solver is forward in time only and returns
[`DensityMatrix`](@ref) outputs while storing an internal vectorized Arnoldi
representation.

Pass `return_diagnostics=true` to receive the corresponding
[`LindbladArnoldiDiagnostics`](@ref) alongside the propagated state(s).
"""
function lindblad_timeevolve(
    A::Union{Lindblad, LiouvillianMap},
    dm::DensityMatrix,
    t::Real;
    return_diagnostics::Bool=false,
    kwargs...
)
    tf = Float64(t)
    tf < 0 && error("lindblad_timeevolve: t must be non-negative; backward evolution is not supported")
    cache = LindbladArnoldiCache(A, dm; kwargs...)
    ρt = lindblad_timeevolve!(cache, tf)
    return return_diagnostics ? (ρt, cache.diagnostics) : ρt
end

function lindblad_timeevolve(
    A::Union{Lindblad, LiouvillianMap},
    ρ::AbstractMatrix,
    t::Real;
    kwargs...
)
    lindblad_timeevolve(A, densitymatrix(ρ), t; kwargs...)
end

function lindblad_timeevolve(
    A::Union{Lindblad, LiouvillianMap},
    dm::DensityMatrix,
    ts::AbstractVector{<:Real};
    return_diagnostics::Bool=false,
    kwargs...
)
    any(<(0), ts) && error("lindblad_timeevolve: all times in ts must be non-negative; backward evolution is not supported")
    cache = LindbladArnoldiCache(A, dm; kwargs...)
    perm = sortperm(ts)
    ts_sorted = collect(Float64.(ts[perm]))
    Tρ = DensityMatrix{eltype(cache.ρ_anchor)}
    out_sorted = Vector{Tρ}(undef, length(ts_sorted))
    for (i, t) in enumerate(ts_sorted)
        out_sorted[i] = lindblad_timeevolve!(cache, t)
    end
    out = similar(out_sorted)
    @inbounds for i in eachindex(perm)
        out[perm[i]] = out_sorted[i]
    end
    return return_diagnostics ? (out, cache.diagnostics) : out
end

function lindblad_timeevolve(
    A::Union{Lindblad, LiouvillianMap},
    ρ::AbstractMatrix,
    ts::AbstractVector{<:Real};
    kwargs...
)
    lindblad_timeevolve(A, densitymatrix(ρ), ts; kwargs...)
end

"""
    lindblad_timeevolve!(cache::LindbladArnoldiCache, t::Real) -> DensityMatrix

Advance an existing Lindblad Arnoldi cache to time `t` and return the density
matrix at that time.
"""
function lindblad_timeevolve!(cache::LindbladArnoldiCache, t::Real)
    tf = Float64(t)
    tf < cache.t_latest && error("lindblad_timeevolve!: requested time $t precedes the last served time $(cache.t_latest); backward evolution is not supported")
    τ = tf - cache.t_anchor
    τ < 0 && error("lindblad_timeevolve!: requested time $t precedes anchor time $(cache.t_anchor)")
    _ensure_valid_up_to!(cache, τ)

    ρ_raw = similar(cache.ρ_anchor)
    _reconstruct_density_raw!(ρ_raw, cache, tf - cache.t_anchor)
    _record_output_diagnostics!(cache, ρ_raw)
    cache.t_latest = tf
    cache.diagnostics.total_times_served += 1
    return _copy_output_density(cache, ρ_raw)
end
