using Combinatorics
using EDKit
using ITensors
using ITensorMPS
using LinearAlgebra
using Random
using SparseArrays
using Statistics

const SAMPLES = parse(Int, get(ENV, "EDKIT_BENCH_SAMPLES", "5"))

function measure(f; samples::Int=SAMPLES)
    f()
    GC.gc()
    times = Float64[]
    allocs = Int[]
    for _ in 1:samples
        GC.gc()
        t0 = time_ns()
        alloc = @allocated f()
        push!(times, (time_ns() - t0) / 1e6)
        push!(allocs, alloc)
    end
    median(times), round(Int, median(allocs))
end

function pct_drop(old, new)
    iszero(old) && return iszero(new) ? 0.0 : -Inf
    100 * (old - new) / old
end

function report(label, old, new; samples::Int=SAMPLES)
    old_ms, old_alloc = measure(old; samples)
    new_ms, new_alloc = measure(new; samples)
    speedup = iszero(new_ms) ? Inf : old_ms / new_ms
    println(label)
    println("  old: ", round(old_ms, digits=3), " ms, ", old_alloc, " bytes")
    println("  new: ", round(new_ms, digits=3), " ms, ", new_alloc, " bytes")
    println("  delta: ", round(speedup, digits=2), "x time, ",
            round(pct_drop(old_alloc, new_alloc), digits=1), "% alloc drop")
end

always_true(_) = true

struct OldTranslationJudge{T}
    F
    K::Int
    A::Int
    L::Int
    B::T
    C::Vector{Float64}
end

function (judge::OldTranslationJudge)(dgt::AbstractVector{<:Integer}, i::Integer)
    isnothing(judge.F) || judge.F(dgt) || return (false, 0.0)
    for n in 1:judge.L-1
        EDKit.circshift!(dgt, judge.A)
        In = EDKit.index(dgt, base=judge.B)
        if In < i
            return (false, 0.0)
        elseif isequal(In, i)
            EDKit.check_momentum(n, judge.K, judge.L) || return (false, 0.0)
            return (true, judge.C[n])
        end
    end
    true, judge.C[judge.L]
end

function old_selectindex_append(f, L::Integer, rg::UnitRange{T}; base::Integer=2, alloc::Integer=1000) where T <: Integer
    dgt = zeros(T, L)
    I = T[]
    sizehint!(I, alloc)
    for i in rg
        EDKit.change!(dgt, i, base=base)
        f(dgt) && append!(I, i)
    end
    I
end

function old_selectindexnorm_append(f, L::Integer, rg::UnitRange{T}; base::Integer=2, alloc::Integer=1000) where T <: Integer
    dgt = zeros(T, L)
    I, R = T[], Float64[]
    sizehint!(I, alloc)
    sizehint!(R, alloc)
    for i in rg
        EDKit.change!(dgt, i, base=base)
        Q, N = f(dgt, i)
        Q || continue
        append!(I, i)
        append!(R, N)
    end
    I, R
end

function old_selectindex_threaded(f, L::Integer; base::T=2, alloc::Integer=1000) where T <: Integer
    nt = Threads.nthreads()
    ni = EDKit.dividerange(base^L, nt)
    nI = Vector{Vector{T}}(undef, nt)
    Threads.@threads for ti in 1:nt
        nI[ti] = EDKit.selectindex(f, L, ni[ti], base=base, alloc=alloc)
    end
    vcat(nI...)
end

function old_projected_selectindex_N(f, L::Integer, N::Integer; base::T=2, alloc::Integer=1000, sorted::Bool=true) where T <: Integer
    I = T[]
    sizehint!(I, alloc)
    for fdgt in multiexponents(L, N)
        all(b < base for b in fdgt) || continue
        dgt = (base - 1) .- fdgt
        isnothing(f) || f(dgt) || continue
        push!(I, EDKit.index(dgt, base=base))
    end
    sorted ? sort(I) : I
end

function old_selectindexnorm_N(f, L::Integer, N::Integer; base::T=2, alloc::Integer=1000, sorted::Bool=true) where T <: Integer
    I, R = T[], Float64[]
    sizehint!(I, alloc)
    sizehint!(R, alloc)
    for fdgt in multiexponents(L, N)
        all(b < base for b in fdgt) || continue
        dgt = (base - 1) .- fdgt
        i = EDKit.index(dgt, base=base)
        Q, nrm = f(dgt, i)
        Q || continue
        push!(I, i)
        push!(R, nrm)
    end
    sorted || return I, R
    sperm = sortperm(I)
    I[sperm], R[sperm]
end

function generator_phase(g)
    prod(g.c[i][g.s[i]] for i in eachindex(g.c))
end

function old_check_min(dgt, g::EDKit.AbelianOperator; base=2)
    EDKit.init!(g)
    I0 = EDKit.index(dgt; base)
    N = 1
    tmp = similar(dgt)
    for _ in 2:EDKit.order(g)
        EDKit._apply_group_action!(dgt, g, base, tmp)
        In = EDKit.index(dgt; base)
        In < I0 && return false, N
        In == I0 && return true, N
        N += 1
    end
    true, N
end

function old_abelian_select_gosper(f, G::EDKit.AbelianOperator, L::Int, candidates::Vector{Int}, C; base::Integer=2)
    Is = Int[]
    Rs = Float64[]
    dgt = zeros(Int, L)
    g_local = deepcopy(G)
    for idx in candidates
        EDKit.change!(dgt, idx; base=Int(base))
        f(dgt) || continue
        Q, n = old_check_min(dgt, g_local; base=Int(base))
        Q || continue
        push!(Is, idx)
        push!(Rs, C[n])
    end
    Is, Rs
end

function old_mul_matrix_scaled!(target, opt, m, alpha, beta)
    iszero(beta) ? fill!(target, zero(eltype(target))) : (target .*= beta)
    iszero(alpha) && return target
    dgt = similar(opt.B.dgt)
    for j in 1:size(m, 1)
        old_colmn_operator!(target, opt, j, dgt, alpha .* view(m, j, :))
    end
    target
end

@inline function old_accumulate!(target::AbstractMatrix, pos, C, val, coeff)
    cv = C * val
    @inbounds for k in axes(target, 2)
        target[pos, k] += cv * coeff[k]
    end
end

function old_colmn_local!(target::AbstractMatrix, M, I, b, dgt, coeff=1)
    rows, vals = rowvals(M), nonzeros(M)
    j = EDKit.index(dgt, I, base=b.B)
    changed = false
    @inbounds for i in nzrange(M, j)
        row, val = rows[i], vals[i]
        EDKit.change!(dgt, I, row, base=b.B)
        C, pos = EDKit.index(b, dgt)
        old_accumulate!(target, pos, C, val, coeff)
        changed = true
    end
    changed && EDKit.change!(dgt, I, j, base=b.B)
    nothing
end

function old_colmn_operator!(target::AbstractMatrix, opt, j::Integer, dgt, coeff=1)
    b, M, I = opt.B, opt.M, opt.I
    r = EDKit.change!(b, j, dgt)
    C = isone(r) ? coeff : coeff / r
    for i in eachindex(M)
        old_colmn_local!(target, M[i], I[i], b, dgt, C)
    end
end

function old_threaded_mul_vector(opt, v)
    ctype = promote_type(eltype(opt), eltype(v))
    nt = Threads.nthreads()
    ni = EDKit.dividerange(length(v), nt)
    Ms = [zeros(ctype, size(opt, 1)) for _ in 1:nt]
    Threads.@threads for i in 1:nt
        dgt = similar(opt.B.dgt)
        for j in ni[i]
            EDKit.colmn!(Ms[i], opt, j, dgt, v[j])
        end
    end
    sum(m for m in Ms)
end

function old_threaded_mul_matrix(opt, m)
    ctype = promote_type(eltype(opt), eltype(m))
    nt = Threads.nthreads()
    ni = EDKit.dividerange(size(m, 1), nt)
    Ms = [zeros(ctype, size(opt, 1), size(m, 2)) for _ in 1:nt]
    Threads.@threads for i in 1:nt
        dgt = similar(opt.B.dgt)
        for j in ni[i]
            EDKit.colmn!(Ms[i], opt, j, dgt, view(m, j, :))
        end
    end
    sum(m for m in Ms)
end

function old_reduced_coeffs(cache, tau::Real)
    m = cache.m
    R = eltype(cache.λ)
    Tc = Complex{R}
    phase = Vector{Tc}(undef, m)
    @inbounds @simd for k in 1:m
        phase[k] = cis(-R(tau) * cache.λ[k]) * cache.Qe1[k]
    end
    cache.Q * phase
end

function old_reconstruct_state!(out, cache, tau::Real)
    c = old_reduced_coeffs(cache, tau)
    V = @view cache.V[:, 1:cache.m]
    mul!(out, V, c)
    @. out *= cache.norm_anchor
    if cache.normalize_output
        nrm = norm(out)
        nrm > 0 && (@. out /= nrm)
    end
    out
end

function old_mps2vec(psi::MPS, B::EDKit.AbstractBasis)
    L = length(psi)
    s = siteinds(psi)
    v = Vector{eltype(psi[1])}(undef, size(B, 1))
    for i in eachindex(v)
        R = EDKit.change!(B, i)
        val = ITensor(1.0)
        for j in 1:L
            val *= psi[j] * state(s[j], B.dgt[j] + 1)
        end
        v[i] = scalar(val) * EDKit.order(B) / R
    end
    v
end

function old_parity_schmidt(parity, v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::EDKit.AbstractTranslationalParityBasis; B1=nothing, B2=nothing)
    dgt = similar(b.dgt)
    R, phase = b.R, b.C[2]
    S = EDKit.schmidtmatrix(promote_type(eltype(v), eltype(b)), b, Ainds, B1, B2; dgt)
    for i in eachindex(v)
        EDKit.change!(b, i, dgt)
        val = v[i] / R[i]
        for _ in 1:length(dgt) ÷ b.A
            EDKit.addto!(S, val)
            EDKit.circshift!(dgt, b.A)
            val *= phase
        end
        dgt .= parity(dgt)
        val *= b.P
        for _ in 1:length(dgt) ÷ b.A
            EDKit.addto!(S, val)
            EDKit.circshift!(dgt, b.A)
            val *= phase
        end
    end
    S.M
end

function old_translation_flip_schmidt(v, Ainds, b)
    old_parity_schmidt(x -> EDKit.spinflip(x, b.B), v, Ainds, b)
end

function old_parityflip_schmidt(v::AbstractVector, Ainds::AbstractVector{<:Integer}, b::ParityFlipBasis; B1=nothing, B2=nothing)
    dgt = similar(b.dgt)
    R, p1, p2 = b.R, b.P, b.Z
    S = EDKit.schmidtmatrix(promote_type(eltype(v), eltype(b)), b, Ainds, B1, B2; dgt)
    for i in eachindex(v)
        EDKit.change!(b, i, dgt)
        val = v[i] / R[i]
        EDKit.addto!(S, val)
        reverse!(dgt)
        val *= p1
        EDKit.addto!(S, val)
        dgt .= EDKit.spinflip(dgt, b.B)
        val *= p2
        EDKit.addto!(S, val)
        reverse!(dgt)
        val *= p1
        EDKit.addto!(S, val)
    end
    S.M
end

Random.seed!(7)
println("EDKit performance finding benchmark")
println("samples=", SAMPLES, ", threads=", Threads.nthreads())

report("Basis typed TranslationJudge predicate call x1000",
    let
        L = 10
        norm = [L / sqrt(i) for i in 1:L]
        old_judge = OldTranslationJudge(always_true, 0, 1, L, Int64(2), norm)
        new_judge = EDKit.TranslationJudge(always_true, 0, 1, L, Int64(2), norm)
        old_dgt = zeros(Int64, L)
        new_dgt = zeros(Int64, L)
        idx = EDKit.index(old_dgt, base=2)
        (() -> begin
            acc = 0.0
            for _ in 1:1000
                fill!(old_dgt, 0)
                _, nrm = old_judge(old_dgt, idx)
                acc += nrm
            end
            acc
        end),
        (() -> begin
            acc = 0.0
            for _ in 1:1000
                fill!(new_dgt, 0)
                _, nrm = new_judge(new_dgt, idx)
                acc += nrm
            end
            acc
        end)
    end...)

report("Projected scalar selection append vs push L=12",
    () -> old_selectindex_append(always_true, 12, Int64(1):Int64(2)^12; base=2),
    () -> EDKit.selectindex(always_true, 12, Int64(1):Int64(2)^12; base=2))

report("Translational scalar selection append vs push L=12",
    let
        L = 12
        norm = [L / sqrt(i) for i in 1:L]
        old_judge = EDKit.TranslationJudge(nothing, 0, 1, L, Int64(2), norm)
        new_judge = EDKit.TranslationJudge(nothing, 0, 1, L, Int64(2), norm)
        (() -> old_selectindexnorm_append(old_judge, L, Int64(1):Int64(2)^L; base=2)),
        (() -> EDKit.selectindexnorm(new_judge, L, Int64(1):Int64(2)^L; base=2))
    end...)

report("Projected threaded worker cap L=1 x1000",
    () -> begin
        for _ in 1:1000
            old_selectindex_threaded(always_true, 1; base=Int64(2))
        end
    end,
    () -> begin
        for _ in 1:1000
            EDKit.selectindex_threaded(always_true, 1; base=Int64(2))
        end
    end)

report("Projected fixed-charge selectindex_N L=16 N=8",
    () -> old_projected_selectindex_N(nothing, 16, 8; base=Int64(2)),
    () -> EDKit.selectindex_N(nothing, 16, 8; base=Int64(2)))

report("Translational fixed-charge selectindexnorm_N L=16 N=8",
    let
        L = 16
        norm = [L / sqrt(i) for i in 1:L]
        old_judge = EDKit.TranslationJudge(nothing, 0, 1, L, Int64(2), norm)
        new_judge = EDKit.TranslationJudge(nothing, 0, 1, L, Int64(2), norm)
        (() -> old_selectindexnorm_N(old_judge, L, 8; base=Int64(2))),
        (() -> EDKit.selectindexnorm_N(new_judge, L, 8; base=Int64(2)))
    end...)

report("Projected base-3 fixed-charge reuse/sort! L=8 N=8",
    () -> old_projected_selectindex_N(nothing, 8, 8; base=Int64(3)),
    () -> EDKit.selectindex_N(nothing, 8, 8; base=Int64(3)))

report("TranslationFlipBasis default predicate L=12",
    () -> TranslationFlipBasis(L=12, k=0, p=1, f=_ -> true, threaded=false),
    () -> TranslationFlipBasis(L=12, k=0, p=1, threaded=false))

const PHASE_AG = let
    L = 6
    perm = [mod1(i - 1, L) for i in 1:L]
    EDKit.AbelianOperator(L, 1, perm)
end

report("Abelian phase x10000",
    () -> begin
        acc = zero(ComplexF64)
        for _ in 1:10_000
            acc += generator_phase(PHASE_AG)
        end
        acc
    end,
    () -> begin
        acc = zero(ComplexF64)
        for _ in 1:10_000
            acc += EDKit.phase(PHASE_AG)
        end
        acc
    end)

report("Abelian gosper select L=14 N=7",
    let
        L = 14
        perm = [mod1(i - 1, L) for i in 1:L]
        ag = EDKit.AbelianOperator(L, 0, perm)
        C = zeros(L)
        for i in eachindex(C)
            iszero(mod(L, i)) && (C[i] = sqrt(L * i))
        end
        candidates = EDKit._gosper_enumerate(L, 7)
        (() -> old_abelian_select_gosper(always_true, ag, L, candidates, C; base=2)),
        (() -> EDKit._abelian_select_gosper(always_true, ag, L, candidates, C; base=2, threaded=false))
    end...)

report("Operator mul! matrix alpha/beta L=8 cols=8",
    let
        L = 8
        H = trans_inv_operator(spin((1.0, "xx"), (0.3, "zz")), 2, L)
        M = randn(ComplexF64, size(H, 2), 8)
        target0 = randn(ComplexF64, size(H, 1), size(M, 2))
        old_target = similar(target0)
        new_target = similar(old_target)
        (() -> (copyto!(old_target, target0); old_mul_matrix_scaled!(old_target, H, M, 2, 3))),
        (() -> (copyto!(new_target, target0); mul!(new_target, H, M, 2, 3)))
    end...)

report("Operator reduced-basis mul! matrix alpha/beta L=12 N=6 cols=8",
    let
        L = 12
        B = TranslationalBasis(L=L, k=0, N=6, threaded=false)
        H = trans_inv_operator(spin((1.0, "+-"), (1.0, "-+"), (0.3, "zz")), 2, B)
        M = randn(ComplexF64, size(H, 2), 8)
        target0 = randn(ComplexF64, size(H, 1), size(M, 2))
        old_target = similar(target0)
        new_target = similar(target0)
        (() -> (copyto!(old_target, target0); old_mul_matrix_scaled!(old_target, H, M, 2, 3))),
        (() -> (copyto!(new_target, target0); mul!(new_target, H, M, 2, 3)))
    end...)

report("Operator threaded mul vector L=10",
    let
        L = 10
        H = trans_inv_operator(spin((1.0, "xx"), (0.3, "zz")), 2, L)
        v = randn(ComplexF64, size(H, 2))
        (() -> old_threaded_mul_vector(H, v)),
        (() -> EDKit.mul(H, v))
    end...)

report("Operator threaded mul matrix L=9 cols=8",
    let
        L = 9
        H = trans_inv_operator(spin((1.0, "xx"), (0.3, "zz")), 2, L)
        M = randn(ComplexF64, size(H, 2), 8)
        (() -> old_threaded_mul_matrix(H, M)),
        (() -> EDKit.mul(H, M))
    end...)

report("Krylov reconstruct_state! N=80 m=35",
    let
        A = randn(ComplexF64, 80, 80)
        H = Hermitian((A + A') / 2)
        psi0 = normalize(randn(ComplexF64, 80))
        cache = KrylovEvolutionCache(H, psi0; tol=1e-12, m_init=35, m_max=50)
        old_out = similar(psi0)
        new_out = similar(psi0)
        (() -> old_reconstruct_state!(old_out, cache, 0.7)),
        (() -> EDKit._reconstruct_state!(new_out, cache, 0.7))
    end...)

report("ITensor mps2vec reduced basis L=8 N=4",
    let
        L = 8
        B_old = ProjectedBasis(L=L, N=4)
        B_new = ProjectedBasis(L=L, N=4)
        psi = vec2mps(normalize(randn(ComplexF64, 2^L)), siteinds(2, L))
        (() -> old_mps2vec(psi, B_old)),
        (() -> mps2vec(psi, B_new))
    end...)

report("Schmidt TranslationFlipBasis L=10",
    let
        B = TranslationFlipBasis(L=10, k=0, p=1)
        v = normalize(randn(ComplexF64, size(B, 1)))
        inds = collect(1:5)
        (() -> old_translation_flip_schmidt(v, inds, B)),
        (() -> EDKit.schmidt(v, inds, B))
    end...)

report("Schmidt ParityFlipBasis L=10",
    let
        B = ParityFlipBasis(L=10, p=1, z=1)
        v = normalize(randn(ComplexF64, size(B, 1)))
        inds = collect(1:5)
        (() -> old_parityflip_schmidt(v, inds, B)),
        (() -> EDKit.schmidt(v, inds, B))
    end...)
