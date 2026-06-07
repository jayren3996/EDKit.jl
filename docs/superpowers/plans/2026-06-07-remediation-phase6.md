# EDKit.jl Remediation — Phase 6: optimizations

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Two behavior-preserving performance improvements: skip the wasted enumeration in `base>2` fixed-charge construction, and build the basis-embedding/symmetrizer maps sparsely.

**Branch:** `remediation-0.6.0`. These are optimizations, so each test is a behavior-preserving correctness guard (compared against an independent reference), plus an assertion that the new representation is in effect.

**Findings:** `optimization-1`, `linearmap-1`. Parked: `abelian-4`, `o1` (see end).

---

## Commit O1 — `optimization-1`: bounded mixed-radix enumeration

The fixed-`N` constructors for `base>2` use `multiexponents(L,N)` then discard every composition with a digit `≥ base` (measured 70–95% waste). Replace with a generator that emits only in-range digit strings.

**Edit 1** — add a shared helper to `src/Basis/AbstractBasis.jl` (after `_check_index_capacity`):
```julia
#-------------------------------------------------------------------------------------------------------------------------
"""
    _foreach_bounded_digits(body, dgt, target, base)

Call `body(dgt)` for every length-`length(dgt)` digit string with entries in
`0:base-1` summing to `target`. Generates only in-range strings (no filtering),
which is the fixed-charge enumeration for `base>2`.
"""
function _foreach_bounded_digits(body::F, dgt::AbstractVector, target::Integer, base::Integer) where F
    _rec_bounded_digits!(body, dgt, 1, target, base, length(dgt))
end
function _rec_bounded_digits!(body::F, dgt::AbstractVector, pos::Int, remaining::Integer, base::Integer, L::Int) where F
    if pos == L
        (0 <= remaining <= base - 1) || return
        @inbounds dgt[L] = remaining
        body(dgt)
        return
    end
    hi = min(base - 1, remaining)
    lo = max(0, remaining - (L - pos) * (base - 1))
    for d in lo:hi
        @inbounds dgt[pos] = d
        _rec_bounded_digits!(body, dgt, pos + 1, remaining - d, base, L)
    end
    return
end
```

**Edit 2** — `src/Basis/ProjectedBasis.jl`, the `else` (non-binary) branch of `selectindex_N` (currently the `for fdgt in multiexponents(L, N)` loop, lines ~195-205):
```julia
    I = T[]
    sizehint!(I, alloc)
    dgt = Vector{T}(undef, L)
    for fdgt in multiexponents(L, N)
        all(b < base for b in fdgt) || continue
        _complement_digits!(dgt, fdgt, base)
        isnothing(f) || f(dgt) || continue
        ind = index(dgt, base=base)
        push!(I, ind)
    end
    sorted ? sort!(I) : I
```
→
```julia
    I = T[]
    sizehint!(I, alloc)
    dgt = zeros(T, L)
    _foreach_bounded_digits(dgt, L * (base - 1) - N, base) do d
        (isnothing(f) || f(d)) && push!(I, index(d, base=base))
    end
    sorted ? sort!(I) : I
```

**Edit 3** — `src/Basis/TranslationalBasis.jl`, the matching `else` branch of `selectindexnorm_N` (the `for fdgt in multiexponents(L, N)` loop, lines ~176-192):
```julia
    I, R = T[], Float64[]
    sizehint!(I, alloc)
    sizehint!(R, alloc)
    dgt = Vector{T}(undef, L)
    for fdgt in multiexponents(L, N)
        all(b < base for b in fdgt) || continue
        _complement_digits!(dgt, fdgt, base)
        i = index(dgt, base=base)
        Q, nrm = f(dgt, i)
        Q || continue
        push!(I, i)
        push!(R, nrm)
    end
    sorted || return I, R
    sperm = sortperm(I)
    I[sperm], R[sperm]
```
→
```julia
    I, R = T[], Float64[]
    sizehint!(I, alloc)
    sizehint!(R, alloc)
    dgt = Vector{T}(undef, L)
    _foreach_bounded_digits(dgt, L * (base - 1) - N, base) do d
        i = index(d, base=base)
        Q, nrm = f(d, i)
        Q && (push!(I, i); push!(R, nrm))
    end
    sorted || return I, R
    sperm = sortperm(I)
    I[sperm], R[sperm]
```

**Test (append, Phase 6 header):**
```julia
# ---------------------------------------------------------------------------
# Phase 6 — optimizations
# ---------------------------------------------------------------------------

@testset "optimization-1: bounded base>2 enumeration matches full scan" begin
    # The small_N (combinatorial) and full-scan paths must build the identical basis.
    for (L, N, base) in [(4, 4, 3), (5, 6, 3), (4, 5, 4), (3, 4, 5)]
        a = sort(ProjectedBasis(L=L, N=N, base=base, small_N=true).I)
        b = sort(ProjectedBasis(L=L, N=N, base=base, small_N=false).I)
        @test a == b
    end
    # Translational fixed-N path likewise unchanged vs full scan.
    for (L, N, base) in [(4, 4, 3), (6, 6, 3)]
        a = sort(TranslationalBasis(L=L, N=N, k=0, base=base, small_N=true).I)
        b = sort(TranslationalBasis(L=L, N=N, k=0, base=base, small_N=false).I)
        @test a == b
    end
end
```
NOTE: confirm `ProjectedBasis`/`TranslationalBasis` accept `small_N` as a keyword (they do per their signatures). Run the full `Pkg.test()` after this commit. CHANGELOG → Performance. Commit: `perf(basis): bounded mixed-radix enumeration for base>2 fixed-N (optimization-1)`.

---

## Commit O2 — `linearmap-1`: sparse `basis_embedding`

`basis_embedding` allocates a dense `base^L × dim(B)` matrix though it has ≤1 nonzero per row. Build it sparsely with the basis's natural (real or complex) element type; `symmetrizer` then becomes sparse automatically.

**Edit** `src/LinearMap.jl:59-77` (`basis_embedding`):
```julia
function basis_embedding(B::AbstractBasis)
    full = TensorBasis(L = length(B), base = B.B)
    embed = zeros(ComplexF64, size(full, 1), size(B, 1))
    ord = orbit_order(B)
    dgt = similar(B.dgt)
    full_dgt = similar(full.dgt)
    for j in 1:size(full, 1)
        change!(full, j, full_dgt)
        dgt .= full_dgt
        coeff, pos = index_nocheck(B, dgt)
        # ... comment ...
        iszero(coeff) || (embed[j, pos] = conj(coeff) / ord)
    end
    embed
end
```
→
```julia
function basis_embedding(B::AbstractBasis)
    full = TensorBasis(L = length(B), base = B.B)
    T = float(eltype(B))                       # real for Projected/Parity/Flip/real-Abelian sectors
    ord = orbit_order(B)
    dgt = similar(B.dgt)
    full_dgt = similar(full.dgt)
    rows = Int[]; cols = Int[]; vals = T[]
    for j in 1:size(full, 1)
        change!(full, j, full_dgt)
        dgt .= full_dgt
        coeff, pos = index_nocheck(B, dgt)
        # index(B, q) returns the phase needed to shift q INTO the canonical
        # representative s_pos; the matrix element ⟨q|k,j⟩ uses the opposite
        # shift, hence the conjugation. For real-phase bases (Projected, Parity,
        # Flip) the conjugation is a no-op.
        if !iszero(coeff)
            push!(rows, j); push!(cols, pos); push!(vals, conj(coeff) / ord)
        end
    end
    sparse(rows, cols, vals, size(full, 1), size(B, 1))
end
```

**Test (append):**
```julia
@testset "linearmap-1: basis_embedding/symmetrizer are sparse and correct" begin
    B = basis(L=6, k=0)
    E = EDKit.basis_embedding(B)
    @test issparse(E)
    @test size(E) == (2^6, size(B, 1))
    Bfull = TensorBasis(L=6)
    D = DoubleBasis(B, Bfull)
    S = symmetrizer(D)
    @test issparse(S)
    v = randn(ComplexF64, size(Bfull, 1))
    @test S * v ≈ D(v)                         # symmetrizer matrix == matrix-free action
end
```
Run the full `Pkg.test()` after this commit (return type changes dense→sparse; the `DoubleBasis`/symmetrizer suite must stay green). CHANGELOG → Performance. Commit: `perf(maps): build basis_embedding/symmetrizer sparsely (linearmap-1)`.

---

## Self-Review

**Spec coverage (Phase 6):** `optimization-1` (O1), `linearmap-1` (O2). `abelian-4`/`o1` parked.

**Placeholder scan:** none.

**Correctness:** Both are behavior-preserving. `optimization-1`'s generator emits exactly the in-range, fixed-sum digit strings — the same set the old `multiexponents+filter+complement` produced — verified by comparing the `small_N` and full-scan bases. `linearmap-1` keeps the same nonzero pattern/values (`conj(coeff)/ord` at `(j, pos)`); only the container changes to sparse and the element type to `float(eltype(B))`, which `symmetrizer`'s matrix product accepts.

---

## Parked (with rationale)

- `abelian-4` (per-call `ms = g.s[:]` copy in `shift_canonical_int`/`shift_canonical!`): the allocation is a length-`Ng` (1–3) `Int` vector, ~64 B for a single generator (the report's "~570 B/call" headline was overstated). Removing it requires either an `AbelianOperator` struct field or threading a scratch buffer through `colmn!→index→shift_canonical*`, i.e. hot-path/signature surgery disproportionate to the gain. Revisit only if profiling shows GC pressure dominates a real large-build workload.
- `o1` (fermion `2^span` colptr memory): the report's "dense 2^span matrix" framing was refuted — the local term is already built sparse. The residual (`colptr` length `2^span+1` for a single long-range bond) needs an endpoint-only + parity-coefficient reconstruction of `_fermion_string_matrix`; the report rated this directionally valid but with the severity overstated. Deferred as an extension rather than a behavior-preserving optimization.
