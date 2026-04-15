"""
Benchmark micro-optimizations for the matrix-free path.

Tests:
1. base=2 specialization of change!/index using bit operations
2. Replacing custom binary_search with searchsortedlast
3. Integer-based state representation for base=2 (avoid digit array entirely)
4. circshift! + index → direct integer rotation for TranslationalBasis
"""

using EDKit
using LinearAlgebra, SparseArrays
using BenchmarkTools

BenchmarkTools.DEFAULT_PARAMETERS.samples = 20
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 2

# ============================================================================
# 1. Bit-operation specializations for base=2
# ============================================================================

# Current: uses divrem
@inline function change_current!(dgt::AbstractVector{T}, ind::T; base::T=Int64(2)) where T
    N = ind - one(T)
    @inbounds for i = length(dgt):-1:1
        N, dgt[i] = divrem(N, base)
    end
end

# Proposed: bit shifts for base=2
@inline function change_bits!(dgt::AbstractVector{T}, ind::T) where T
    N = ind - one(T)
    @inbounds for i = length(dgt):-1:1
        dgt[i] = N & one(T)
        N >>= 1
    end
end

# Current index
@inline function index_current(dgt::AbstractVector{<:Integer}; base::T=Int64(2)) where T <: Integer
    N = zero(T)
    @inbounds for i = 1:length(dgt)
        N *= base
        N += dgt[i]
    end
    N + one(T)
end

# Proposed: bit shifts for base=2
@inline function index_bits(dgt::AbstractVector{<:Integer})
    N = zero(Int64)
    @inbounds for i = 1:length(dgt)
        N = (N << 1) | dgt[i]
    end
    N + Int64(1)
end

# Proposed: muladd version (works for any base, may help FMA)
@inline function index_muladd(dgt::AbstractVector{<:Integer}; base::T=Int64(2)) where T <: Integer
    N = zero(T)
    @inbounds for i = 1:length(dgt)
        N = muladd(N, base, dgt[i])
    end
    N + one(T)
end

# Partial change! for sites (the most called one in colmn!)
@inline function change_sites_current!(dgt::AbstractVector{T}, sites::AbstractVector{<:Integer}, ind::Integer; base::T=Int64(2)) where T
    N = ind - one(T)
    @inbounds for i = length(sites):-1:1
        N, dgt[sites[i]] = divrem(N, base)
    end
end

@inline function change_sites_bits!(dgt::AbstractVector{T}, sites::AbstractVector{<:Integer}, ind::Integer) where T
    N = ind - one(T)
    @inbounds for i = length(sites):-1:1
        dgt[sites[i]] = N & one(T)
        N >>= 1
    end
end

@inline function index_sites_current(dgt::AbstractVector{T}, sites::AbstractVector{<:Integer}; base::T=Int64(2)) where T <: Integer
    N = zero(T)
    @inbounds for i in sites
        N *= base
        N += dgt[i]
    end
    N + one(T)
end

@inline function index_sites_bits(dgt::AbstractVector{T}, sites::AbstractVector{<:Integer}) where T <: Integer
    N = zero(T)
    @inbounds for i in sites
        N = (N << 1) | dgt[i]
    end
    N + one(T)
end

# ============================================================================
# 2. binary_search vs searchsortedlast
# ============================================================================
function binary_search_current(list::AbstractVector{<:Integer}, i::Integer)
    l::Int = 1
    r::Int = length(list)
    c::Int = (l + r) ÷ 2
    while true
        t = list[c]
        (i < t) ? (r = c - 1) : (i > t) ? (l = c + 1) : break
        (l > r) ? (c = 0; break) : (c = (l + r) ÷ 2)
    end
    c
end

@inline function binary_search_stdlib(list::AbstractVector{<:Integer}, i::Integer)
    idx = searchsortedfirst(list, i)
    @inbounds (idx <= length(list) && list[idx] == i) ? idx : 0
end

# ============================================================================
# 3. Integer-state approach: avoid digit array for index computation
#    For base=2, the entire state is just an integer. circshift! + index
#    can be replaced by bit rotation.
# ============================================================================
@inline function bitrotate_index(state::Int64, L::Int, shift::Int)
    # Circular left shift of an L-bit integer by `shift` positions
    mask = (Int64(1) << L) - Int64(1)
    s = shift % L
    ((state << s) | (state >> (L - s))) & mask
end

# ============================================================================
# Benchmarks
# ============================================================================

println("="^60)
println("  Micro-optimization benchmarks (matrix-free path)")
println("="^60)

# --- index / change! benchmarks ---
for L in [12, 16, 20]
    println("\n--- L = $L ---")
    dgt = zeros(Int64, L)
    sites = [3, 4]  # typical 2-site operator

    # Full index
    change_current!(dgt, Int64(42))
    @assert index_current(dgt) == 42
    @assert index_bits(dgt) == 42
    @assert index_muladd(dgt) == 42

    vals = collect(Int64(1):Int64(2^L))
    print("  index(dgt):         current=")
    b1 = @benchmark for v in $vals; change_current!($dgt, v); index_current($dgt); end
    print("$(BenchmarkTools.prettytime(median(b1).time))  bits=")
    b2 = @benchmark for v in $vals; change_bits!($dgt, v); index_bits($dgt); end
    print("$(BenchmarkTools.prettytime(median(b2).time))  muladd=")
    b3 = @benchmark for v in $vals; change_current!($dgt, v); index_muladd($dgt); end
    println("$(BenchmarkTools.prettytime(median(b3).time))")

    println("  Speedup: bits=$(round(median(b1).time/median(b2).time, digits=2))x  muladd=$(round(median(b1).time/median(b3).time, digits=2))x")

    # Partial index/change (2-site, the hot path in colmn!)
    change_current!(dgt, Int64(42))
    @assert index_sites_current(dgt, sites) == index_sites_bits(dgt, sites)

    print("  index(dgt,sites):   current=")
    b4 = @benchmark for v in $vals; change_sites_current!($dgt, $sites, (v % 4) + 1); index_sites_current($dgt, $sites); end
    print("$(BenchmarkTools.prettytime(median(b4).time))  bits=")
    b5 = @benchmark for v in $vals; change_sites_bits!($dgt, $sites, (v % 4) + 1); index_sites_bits($dgt, $sites); end
    println("$(BenchmarkTools.prettytime(median(b5).time))")
    println("  Speedup: bits=$(round(median(b4).time/median(b5).time, digits=2))x")
end

# --- binary_search benchmarks ---
println("\n--- binary_search ---")
for dim in [246, 810, 3432, 12870]
    sorted = sort(rand(1:100000, dim)) |> unique
    dim_actual = length(sorted)
    targets = rand(sorted, 1000)  # half will be found

    print("  dim≈$dim_actual: current=")
    b1 = @benchmark for t in $targets; binary_search_current($sorted, t); end
    print("$(BenchmarkTools.prettytime(median(b1).time))  stdlib=")
    b2 = @benchmark for t in $targets; binary_search_stdlib($sorted, t); end
    println("$(BenchmarkTools.prettytime(median(b2).time))")
    println("  Speedup: $(round(median(b1).time/median(b2).time, digits=2))x")

    # Verify correctness
    for t in targets
        @assert binary_search_current(sorted, t) == binary_search_stdlib(sorted, t)
    end
end

# --- bit rotation for TranslationalBasis ---
println("\n--- bit rotation vs circshift!+index (TranslationalBasis) ---")
for L in [12, 16, 20]
    dgt = zeros(Int64, L)
    ncycles = L  # single-site unit cell

    # Simulate orbit search: circshift + index L times
    print("  L=$L:  circshift+index=")
    b1 = @benchmark begin
        change_current!($dgt, Int64(42))
        Im = index_current($dgt)
        for _ in 1:$ncycles-1
            circshift!($dgt, 1)
            In = index_current($dgt)
            In < Im && (Im = In)
        end
    end

    # Integer bit rotation approach
    print("$(BenchmarkTools.prettytime(median(b1).time))  bitrotate=")
    b2 = @benchmark begin
        state = Int64(42) - Int64(1)
        Im = state
        mask = (Int64(1) << $L) - Int64(1)
        for _ in 1:$ncycles-1
            state = ((state << 1) | (state >> ($L - 1))) & mask
            state < Im && (Im = state)
        end
    end
    println("$(BenchmarkTools.prettytime(median(b2).time))")
    println("  Speedup: $(round(median(b1).time/median(b2).time, digits=2))x")
end

# --- End-to-end: full Operator * vector ---
println("\n--- End-to-end: Operator * vector (current code with @inbounds) ---")
function heisenberg_operator(L, basis)
    mats = []; inds = []
    Sx = [0 0.5; 0.5 0]; Sy = [0 -0.5im; 0.5im 0]; Sz = [0.5 0; 0 -0.5]
    for i in 1:L
        j = mod1(i+1, L)
        push!(mats, kron(Sx, Sx) + kron(Sy, Sy) + kron(Sz, Sz))
        push!(inds, [i, j])
    end
    operator(sparse.(mats), [collect(ind) for ind in inds], basis)
end

for L in [12, 16]
    for (label, bfn) in [
        ("TensorBasis", () -> TensorBasis(L=L, base=2)),
        ("ProjectedBasis", () -> ProjectedBasis(L=L, base=2, f=x->sum(x)==L÷2)),
        ("TranslationalBasis", () -> TranslationalBasis(L=L, base=2, k=0, f=x->sum(x)==L÷2)),
    ]
        b = bfn()
        opt = heisenberg_operator(L, b)
        dim = size(opt, 1)
        v = randn(ComplexF64, dim)
        bm = @benchmark $opt * $v
        println("  $label L=$L (dim=$dim): $(BenchmarkTools.prettytime(median(bm).time))")
    end
end

println("\nDone!")
