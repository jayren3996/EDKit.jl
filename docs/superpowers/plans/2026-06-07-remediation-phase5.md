# EDKit.jl Remediation — Phase 5: long-tail correctness bugs

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Fix the twelve remaining confirmed correctness/quality findings, each with a regression test, grouped into focused commits.

**Branch:** `remediation-0.6.0`. Tests append to `test/remediation_tests.jl` under a "Phase 5" header. Run targeted `julia --project test/remediation_tests.jl` per fix; run the full `Pkg.test()` at the checkpoints noted, and always before declaring done.

**Findings:** `trans-1`, `abelian-3`, `schmidt-1`, `schmidt-2`, `expm-1`, `gapratio-1`, `operator-6`, `bug-4`, `parityflip-4`, `qim-1` (doc), `tebd4` (doc), `parityflip-2`.

---

## Commit C1 — `trans-1`: `order()` returns the orbit size `L/a`

`order` feeds `basis_embedding`/`symmetrizer`/`mps2vec` normalization; returning `L` instead of the orbit size `L/a` mis-normalizes by `a` for `a>1` unit cells (invisible at `a=1`).

**Edits:**
- `src/Basis/TranslationalBasis.jl:353`: `order(b::TranslationalBasis) = length(b.dgt)` → `order(b::TranslationalBasis) = ncycle(b)`
- `src/Basis/TranslationalParityBasis.jl:228`: `order(b::TranslationParityBasis) = 2*length(b.dgt)` → `order(b::TranslationParityBasis) = 2*ncycle(b)`
- `src/Basis/TranslationalFlipBasis.jl:211`: `order(b::TranslationFlipBasis) = 2*length(b.dgt)` → `order(b::TranslationFlipBasis) = 2*ncycle(b)`

**Test (append):**
```julia
# ---------------------------------------------------------------------------
# Phase 5 — long-tail correctness bugs
# ---------------------------------------------------------------------------

@testset "trans-1: order() returns the orbit size L/a" begin
    @test order(TranslationalBasis(L=6, k=0, a=1)) == 6           # unchanged at a=1
    @test order(TranslationalBasis(L=6, k=0, a=2)) == 3           # was 6
    @test order(TranslationalBasis(L=6, k=0, a=3)) == 2           # was 6
    @test order(TranslationParityBasis(L=6, k=0, p=1, a=2)) == 6  # 2*ncycle, was 12
    @test order(TranslationFlipBasis(L=6, k=0, p=1, a=2)) == 6    # 2*ncycle, was 12
end
```
**Checkpoint:** run full `Pkg.test()` after this commit (behavioral). CHANGELOG → Breaking. Commit: `fix(basis): order() returns orbit size L/a, not L (trans-1)`.

---

## Commit C2 — `abelian-3`: thread-safe 2-arg `index`

**Edit** `src/Basis/AbelianBasis.jl:773-775`:
```julia
function index(B::AbelianBasis, dgt::AbstractVector)
    index(B, dgt, B.G)
end
```
→
```julia
function index(B::AbelianBasis, dgt::AbstractVector)
    index(B, dgt, _shallow_workspace(B.G))
end
```

**Test (append):**
```julia
@testset "abelian-3: 2-arg index is thread-safe" begin
    B = basis(L=8, k=1)              # complex sector, exercises the odometer
    cfgs = [rand(0:1, 8) for _ in 1:2000]
    serial = [EDKit.index(B, c) for c in cfgs]
    parallel = Vector{Any}(undef, length(cfgs))
    Threads.@threads for i in eachindex(cfgs)
        parallel[i] = EDKit.index(B, cfgs[i])
    end
    @test parallel == serial
end
```
CHANGELOG → Fixed. Commit: `fix(basis): 2-arg AbelianBasis index uses a private workspace (abelian-3)`.

---

## Commit C3 — `schmidt-1` + `schmidt-2`: Rényi entropy cutoff + normalization

**Edit** `src/Schmidt.jl:142` (`renyi_entropy`):
```julia
renyi_entropy(s::AbstractVector{<:Real}, α::Real) = log(sum(s.^α)) / (1-α)
```
→
```julia
function renyi_entropy(s::AbstractVector{<:Real}, α::Real; cutoff::Real=1e-20)
    Z = 0.0
    for si in s
        si > cutoff && (Z += si)
    end
    iszero(Z) && return 0.0
    acc = 0.0
    for si in s
        si > cutoff && (acc += (si / Z)^α)
    end
    log(acc) / (1 - α)
end
```
**Edit** `src/Schmidt.jl:111` (in `entropy`): `renyi_entropy(s, α)` → `renyi_entropy(s, α, cutoff=cutoff)`

**Test (append):**
```julia
@testset "schmidt-1/2: Renyi entropy respects cutoff and normalizes" begin
    # schmidt-1: a noise eigenvalue below the cutoff is dropped (was ignored entirely).
    s = [0.25, 0.25, 0.25, 0.25, 1e-8]
    @test entropy(s, α=0.5, cutoff=1e-6) ≈ log(4)
    # schmidt-2: unnormalized input is normalized before the Renyi formula.
    @test entropy([0.5, 0.5, 0.5, 0.5], α=2) ≈ log(4)
end
```
**Checkpoint:** run full `Pkg.test()` after this commit. CHANGELOG → Fixed. Commit: `fix(schmidt): Renyi entropy honors cutoff and normalizes input (schmidt-1/2)`.

---

## Commit C4 — `expm-1`: scaling-and-squaring

**Edit** `src/ToolKit.jl:58-68` (`expm`):
```julia
function expm(A; order::Integer=10)
    mat = I + A / order
    order -= 1
    while order > 0
        mat = A * mat
        mat ./= order
        mat += I
        order -= 1
    end
    mat
end
```
→
```julia
function expm(A; order::Integer=10)
    nrm = opnorm(A, 1)
    s = nrm > 0.5 ? ceil(Int, log2(nrm)) + 1 : 0
    B = A / (2.0^s)
    mat = I + B / order
    k = order - 1
    while k > 0
        mat = B * mat
        mat ./= k
        mat += I
        k -= 1
    end
    for _ in 1:s
        mat = mat * mat
    end
    mat
end
```
**Edit** `src/ToolKit.jl:81-91` (`expv`) — route through the now-accurate `expm`:
```julia
function expv(A, v::AbstractVecOrMat; order::Integer=10, λ::Number=1)
    vec = v + λ * A * v / order
    order -= 1
    while order > 0
        vec = λ * A * vec
        vec ./= order
        vec += v
        order -= 1
    end
    vec
end
```
→
```julia
function expv(A, v::AbstractVecOrMat; order::Integer=10, λ::Number=1)
    expm(λ * A; order=order) * v
end
```

**Test (append):**
```julia
@testset "expm-1: scaling-and-squaring is accurate at large norm" begin
    A = [0.0 5.0; -5.0 0.0]            # ‖A‖ ≈ 5, where degree-10 Taylor diverges
    @test expm(A) ≈ exp(A)
    B = 3 .* (randn(4,4) |> M -> M - M')  # antisymmetric, larger norm
    @test expm(B) ≈ exp(B)
    v = randn(4)
    @test expv(B, v) ≈ exp(B) * v
end
```
CHANGELOG → Fixed. Commit: `fix(toolkit): expm/expv use scaling-and-squaring (expm-1)`.

---

## Commit C5 — `gapratio-1`: sorted keyword + degeneracy guard

**Edit** `src/ToolKit.jl:15-29` (`gapratio`) and `:38-41` (`meangapratio`):
```julia
function gapratio(E::AbstractVector{<:Real})
    dE = diff(E)
    length(dE) < 2 && return Float64[]
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
```
→
```julia
function gapratio(E::AbstractVector{<:Real}; sorted::Bool=true)
    Es = sorted ? E : sort(E)
    dE = diff(Es)
    length(dE) < 2 && return Float64[]
    r = zeros(length(dE)-1)
    for i in eachindex(r)
        lo, hi = minmax(dE[i], dE[i+1])
        r[i] = iszero(hi) ? NaN : lo / hi   # 0/0 (full degeneracy) is undefined, not 1
    end
    r
end
```
and
```julia
function meangapratio(E::AbstractVector{<:Real})
    r = gapratio(E)
    isempty(r) ? NaN : sum(r) / length(r)
end
```
→
```julia
function meangapratio(E::AbstractVector{<:Real}; sorted::Bool=true)
    r = filter(!isnan, gapratio(E; sorted))
    isempty(r) ? NaN : sum(r) / length(r)
end
```

**Test (append):**
```julia
@testset "gapratio-1: sorting + degeneracy handling" begin
    @test gapratio([3.0, 1.0, 2.0]; sorted=false) == gapratio([1.0, 2.0, 3.0])
    r = gapratio([1.0, 1.0, 1.0, 2.0])           # 0/0 at the triple-degenerate gap
    @test isnan(r[1]) && r[2] == 0.0             # was 1.0 / 0.0 (wrong)
    @test !isnan(meangapratio([1.0, 1.0, 1.0, 2.0]))  # NaNs filtered out
end
```
CHANGELOG → Fixed. Commit: `fix(toolkit): gapratio sorted keyword and degeneracy guard (gapratio-1)`.

---

## Commit C6 — `operator-6`: O(1) support dedup (behavior-preserving)

**Edit** `src/Operator.jl:62-74` (the dedup loop in `operator(mats, inds, B)`):
```julia
    N = 0
    for i = 1:num
        iszero(mats[i]) && continue
        ind = inds[i]
        pos = findfirst(x -> isequal(x, ind), view(I, 1:N))
        if isnothing(pos)
            N += 1
            I[N] = ind
            M[N] = sparse(mats[i])
        else
            M[pos] += mats[i]
        end
    end
```
→
```julia
    N = 0
    slot = Dict{Vector{Int64}, Int}()
    for i = 1:num
        iszero(mats[i]) && continue
        ind = inds[i]
        key = Vector{Int64}(ind)
        pos = get(slot, key, 0)
        if iszero(pos)
            N += 1
            I[N] = ind
            M[N] = sparse(mats[i])
            slot[key] = N
        else
            M[pos] += mats[i]
        end
    end
```

**Test (append):**
```julia
@testset "operator-6: support dedup sums duplicate terms" begin
    X = [0.0 1.0; 1.0 0.0]
    op = operator([X, 2 .* X], [[1], [1]], TensorBasis(L=4))
    ref = operator([3 .* X], [[1]], TensorBasis(L=4))
    @test Array(op) ≈ Array(ref)
end
```
**Checkpoint:** run full `Pkg.test()` after this commit. CHANGELOG → Performance. Commit: `perf(operator): O(1) support dedup via Dict (operator-6)`.

---

## Commit C7 — `bug-4`: loosen `change!` index/base types

**Edit** `src/Basis/AbstractBasis.jl:360-372` (the `change!(dgt, ind; base)` method):
```julia
@inline function change!(dgt::AbstractVector{T}, ind::T; base::T=2) where T
    N = ind - one(T)
    if base == 2
        @inbounds for i = length(dgt):-1:1
            dgt[i] = N & one(T)
            N >>= 1
        end
    else
        @inbounds for i = length(dgt):-1:1
            N, dgt[i] = divrem(N, base)
        end
    end
end
```
→
```julia
@inline function change!(dgt::AbstractVector{T}, ind::Integer; base::Integer=2) where T
    N = ind - oneunit(ind)
    if base == 2
        @inbounds for i = length(dgt):-1:1
            dgt[i] = N & oneunit(N)
            N >>= 1
        end
    else
        @inbounds for i = length(dgt):-1:1
            N, r = divrem(N, base)
            dgt[i] = r
        end
    end
    dgt
end
```

**Test (append):**
```julia
@testset "bug-4: change! accepts mismatched integer types for ind/base" begin
    dgt = zeros(Int32, 4)
    EDKit.change!(dgt, 6; base=2)     # ind::Int64, dgt::Int32 (was a MethodError)
    @test dgt == [0, 1, 0, 1]
    dgt3 = zeros(Int32, 3)
    EDKit.change!(dgt3, 9; base=3)    # base::Int64
    @test dgt3 == [0, 2, 2]
end
```
CHANGELOG → Fixed. Commit: `fix(basis): loosen change! index/base integer types (bug-4)`.

---

## Commit C8 — `parityflip-4`: add missing `copy` methods

**Edits (append a `copy` method to each file, after the struct/`order`):**
- `src/Basis/ParityBasis.jl` (after `eltype(::ParityBasis) = Float64`):
  ```julia
  copy(b::ParityBasis) = ParityBasis(deepcopy(b.dgt), b.I, b.R, b.P, b.B)
  ```
- `src/Basis/FlipBasis.jl` (after `eltype(::FlipBasis) = Float64`):
  ```julia
  copy(b::FlipBasis) = FlipBasis(deepcopy(b.dgt), b.I, b.R, b.P, b.M, b.B)
  ```
- `src/Basis/ParityFlipBasis.jl` (after `eltype(::ParityFlipBasis) = Float64`):
  ```julia
  copy(b::ParityFlipBasis) = ParityFlipBasis(deepcopy(b.dgt), b.I, b.R, b.P, b.Z, b.M, b.B)
  ```

**Test (append):**
```julia
@testset "parityflip-4: copy methods for parity/flip bases" begin
    for b in (ParityBasis(L=4, p=1), FlipBasis(L=4, p=1), ParityFlipBasis(L=4, p=1, z=1))
        c = copy(b)
        @test c.I === b.I            # representative list shared
        @test c.dgt !== b.dgt        # buffer independent
        @test c.dgt == b.dgt
    end
end
```
**Checkpoint:** run full `Pkg.test()` after this commit. CHANGELOG → Fixed. Commit: `fix(basis): add copy methods for parity/flip bases (parityflip-4)`.

---

## Commit C9 — `qim-1`: document the Hermitian assumption (doc only)

**Edit** `src/algorithms/QIM.jl` (`covmat` docstring) — add after the formula block:
```
Note:
- The `½⟨{hᵢ,hⱼ}⟩` form requires each operator in `ol` to be **Hermitian**. The
  implementation computes `Re⟨hᵢψ|hⱼψ⟩ − ⟨hᵢ⟩⟨hⱼ⟩`, which equals the documented
  anticommutator covariance only for Hermitian `hᵢ`. Pass Hermitian operators.
```
No test (doc only). CHANGELOG → Fixed. Commit: `docs(qim): document covmat's Hermitian-operator assumption (qim-1)`.

---

## Commit C10 — `tebd4`: clarify the time-step convention (doc only)

**Edit** `src/ITensors/TEBD.jl` (`tebd4` docstring, the `τ` argument line):
```
- `τ`: physical time step.
```
→
```
- `τ`: time-step multiplier. The gate built is `exp(τ·h)`, so `τ` multiplies the
  local Hamiltonian directly. For **real-time** evolution by `dt`, pass
  `τ = -im*dt`; passing a real `dt` performs imaginary-time evolution.
```
No test (doc only). CHANGELOG → Fixed. Commit: `docs(tebd): clarify tebd4 exp(τ·h) time-step convention (tebd4)`.

---

## Commit C11 — `parityflip-2`: re-check the predicate on orbit partners

A non-symmetry-invariant `f` must reject orbits whose reflected/flipped partners fail `f` (else the basis over-counts states that cannot form symmetry eigenstates). `F` is `nothing` for the no-predicate and N-sector paths, so this only runs (and only allocates) when a custom `f` is supplied.

**Edit** `src/Basis/FlipBasis.jl` (the `FlipJudge` call operator) — after the initial `judge.F(dgt)` check, before the parity check:
```julia
    isnothing(judge.F) || judge.F(dgt) || return (false, 0.0)
    
    # Check parity
    In = judge.MAX - i
```
→
```julia
    isnothing(judge.F) || judge.F(dgt) || return (false, 0.0)
    # f may not be flip-invariant: the flipped partner must also satisfy f.
    isnothing(judge.F) || judge.F(judge.B .- 1 .- dgt) || return (false, 0.0)

    # Check parity
    In = judge.MAX - i
```

**Edit** `src/Basis/ParityBasis.jl` (the `ParityJudge` call operator) — after the initial `judge.F(dgt)` check:
```julia
    isnothing(judge.F) || judge.F(dgt) || return (false, 0.0)
    
    # Check parity
    In = rindex(dgt, base=judge.B)
```
→
```julia
    isnothing(judge.F) || judge.F(dgt) || return (false, 0.0)
    # f may not be reflection-invariant: the reflected partner must also satisfy f.
    isnothing(judge.F) || judge.F(reverse(dgt)) || return (false, 0.0)

    # Check parity
    In = rindex(dgt, base=judge.B)
```

**Edit** `src/Basis/ParityFlipBasis.jl` (the `ParityFlipJudge` call operator) — after the initial `judge.F(dgt)` check:
```julia
    isnothing(judge.F) || judge.F(dgt) || return (false, 0.0)
```
→
```julia
    isnothing(judge.F) || judge.F(dgt) || return (false, 0.0)
    # f may not be reflection/flip-invariant: every orbit partner must satisfy f.
    if !isnothing(judge.F)
        flp = judge.B .- 1 .- dgt
        (judge.F(reverse(dgt)) && judge.F(flp) && judge.F(reverse(flp))) || return (false, 0.0)
    end
```

**Test (append):**
```julia
@testset "parityflip-2: non-invariant predicate rejects broken orbits" begin
    # f = (dgt[1]==0) is NOT flip/reflection invariant; no orbit is f-closed,
    # so the reduced basis must be empty rather than over-counted.
    @test size(FlipBasis(L=4, p=1, f=dgt->dgt[1]==0), 1) == 0
    @test size(FlipBasis(L=4, p=-1, f=dgt->dgt[1]==0), 1) == 0
    @test size(ParityBasis(L=4, p=1, f=dgt->dgt[1]==0), 1) == 0
    @test size(ParityFlipBasis(L=4, p=1, z=1, f=dgt->dgt[1]==0), 1) == 0
    # An invariant predicate is unaffected (sanity).
    @test size(FlipBasis(L=4, p=1, f=dgt->true), 1) > 0
end
```
**Checkpoint:** run full `Pkg.test()` after this commit. CHANGELOG → Breaking (custom non-invariant `f` now yields a smaller/empty basis). Commit: `fix(basis): re-check predicate on parity/flip orbit partners (parityflip-2)`.

---

## Self-Review

**Spec coverage (Phase 5):** trans-1 (C1), abelian-3 (C2), schmidt-1/2 (C3), expm-1 (C4), gapratio-1 (C5), operator-6 (C6), bug-4 (C7), parityflip-4 (C8), qim-1 (C9), tebd4 (C10), parityflip-2 (C11). All twelve covered.

**Placeholder scan:** none.

**Type/consistency:** `ncycle(b) = length(b.dgt) ÷ b.A` is defined generically and applies to all three translational bases (each has `.A`). `renyi_entropy` keeps its positional `(s, α)` call and adds a `cutoff` keyword; `entropy` forwards it. `gapratio`/`meangapratio` add a `sorted` keyword with the prior default behavior. The `parityflip-2` partner checks are guarded by `isnothing(judge.F)`, so the no-`f` and N-sector hot paths are untouched.

**Risk note:** C1 and C11 are Breaking (change results for `a>1` and for non-invariant custom `f`). If a full-suite failure traces to one of these, confirm the test was asserting the *old buggy* value before adjusting it; do not silence a real regression.
