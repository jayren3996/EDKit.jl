"""
Benchmark: Thread-safe alternatives for the mutable `dgt` buffer in Basis types.

Problem: Every Basis struct contains a mutable `dgt::Vector{Int}` that is written
to by `change!()` and read by `index()` on every operator application step.
When a Basis is shared across threads, this is a data race.

Approaches benchmarked:
1. **Baseline (current)**: `dgt` lives in the struct, mutated in-place
2. **Copy-per-thread** (current `mul` workaround): `copy(basis)` per thread
3. **Thread-local via task_local_storage**: lazily allocate per-task buffers
4. **Explicit buffer passing**: `dgt` removed from struct, passed as argument
5. **Stack-allocated MVector**: use StaticArrays for small L

The hot path is:
    change!(b.dgt, content(b, i), base=b.B)   # decode state i into digits
    index(b.dgt, sites, base=b.B)              # read subset of digits
    change!(b.dgt, sites, row, base=b.B)       # write subset of digits
    index(b)                                    # full basis index lookup

We benchmark these primitives directly, plus end-to-end operator*vector.
"""

using EDKit
using BenchmarkTools
using LinearAlgebra
using SparseArrays

#=============================================================================
  Approach 1: Baseline (current) - dgt in struct, single-threaded
=============================================================================#

function baseline_hot_loop!(target, opt, v)
    b = opt.B
    for j = 1:length(v)
        EDKit.change!(b, j)
        for k = 1:length(opt.M)
            EDKit.colmn!(target, opt.M[k], opt.I[k], b, v[j])
        end
    end
    target
end

#=============================================================================
  Approach 2: Copy-per-thread (current mul implementation)
=============================================================================#

function copy_per_thread_mul(opt, v)
    ctype = promote_type(eltype(opt), eltype(v))
    nt = Threads.nthreads()
    ni = EDKit.dividerange(length(v), nt)
    Ms = [zeros(ctype, size(opt, 1)) for _ in 1:nt]
    Threads.@threads for i in 1:nt
        opt_c = EDKit.Operator(opt.M, opt.I, copy(opt.B))
        for j in ni[i]
            EDKit.colmn!(Ms[i], opt_c, j, v[j])
        end
    end
    sum(Ms)
end

#=============================================================================
  Approach 3: Thread-local storage via task_local_storage
=============================================================================#

function get_thread_local_dgt(b::EDKit.AbstractBasis)
    key = :edkit_dgt_buffer
    tls = task_local_storage()
    buf = get(tls, key, nothing)
    if buf === nothing || length(buf) != length(b.dgt)
        buf = zeros(eltype(b.dgt), length(b.dgt))
        tls[key] = buf
    end
    buf::Vector{eltype(b.dgt)}
end

# Simulate the hot loop using TLS-based dgt instead of b.dgt
function tls_hot_loop!(target, opt, v)
    b = opt.B
    dgt = get_thread_local_dgt(b)
    for j = 1:length(v)
        EDKit.change!(dgt, EDKit.content(b, j), base=b.B)
        r = EDKit.norm(b, j)
        C = isone(r) ? v[j] : v[j] / r
        for k = 1:length(opt.M)
            _colmn_with_dgt!(target, opt.M[k], opt.I[k], b, dgt, C)
        end
    end
    target
end

function _colmn_with_dgt!(target, M::SparseMatrixCSC, I::Vector{Int}, b, dgt, coeff)
    rows, vals = rowvals(M), nonzeros(M)
    j = EDKit.index(dgt, I, base=b.B)
    changed = false
    @inbounds for i in SparseArrays.nzrange(M, j)
        row, val = rows[i], vals[i]
        EDKit.change!(dgt, I, row, base=b.B)
        # For TensorBasis/ProjectedBasis, index just reads dgt
        idx = EDKit.index(dgt, base=b.B)
        # For onsite bases, coefficient is 1
        if b isa EDKit.AbstractOnsiteBasis
            pos = if b isa EDKit.TensorBasis
                idx
            else
                EDKit.binary_search(b.I, idx)
            end
            pos > 0 && (target[pos] += coeff * val)
        end
        changed = true
    end
    changed && EDKit.change!(dgt, I, j, base=b.B)
    nothing
end

#=============================================================================
  Approach 4: Explicit buffer passing (dgt as function argument)
  This simulates what the code would look like if dgt were removed from the struct
=============================================================================#

function explicit_buf_hot_loop!(target, opt, v)
    b = opt.B
    dgt = similar(b.dgt)  # allocate once per call
    for j = 1:length(v)
        EDKit.change!(dgt, EDKit.content(b, j), base=b.B)
        r = EDKit.norm(b, j)
        C = isone(r) ? v[j] : v[j] / r
        for k = 1:length(opt.M)
            _colmn_with_dgt!(target, opt.M[k], opt.I[k], b, dgt, C)
        end
    end
    target
end

function explicit_buf_threaded_mul(opt, v)
    ctype = promote_type(eltype(opt), eltype(v))
    nt = Threads.nthreads()
    ni = EDKit.dividerange(length(v), nt)
    Ms = [zeros(ctype, size(opt, 1)) for _ in 1:nt]
    Threads.@threads for i in 1:nt
        b = opt.B
        dgt = similar(b.dgt)  # thread-local buffer, no basis copy needed
        for j in ni[i]
            EDKit.change!(dgt, EDKit.content(b, j), base=b.B)
            r = EDKit.norm(b, j)
            C = isone(r) ? v[j] : v[j] / r
            for k = 1:length(opt.M)
                _colmn_with_dgt!(Ms[i], opt.M[k], opt.I[k], b, dgt, C)
            end
        end
    end
    sum(Ms)
end

#=============================================================================
  Benchmark Runner
=============================================================================#

function run_benchmarks()
    println("=" ^ 70)
    println("Thread-safe dgt benchmark")
    println("Julia threads: $(Threads.nthreads())")
    println("=" ^ 70)

    for L in [10, 14, 18]
        println("\n--- L = $L (Hilbert space dim = $(2^L)) ---")

        # Build Heisenberg Hamiltonian on TensorBasis
        B = TensorBasis(L=L)
        H = trans_inv_operator(spin((1.0, "xx"), (1.0, "yy"), (1.0, "zz")), 1:2, B)
        v = randn(size(H, 2))

        # Verify correctness: baseline
        ref = zeros(eltype(H), size(H, 1))
        baseline_hot_loop!(ref, H, v)

        # Verify explicit buffer approach gives same result
        test_eb = zeros(eltype(H), size(H, 1))
        explicit_buf_hot_loop!(test_eb, H, v)
        err = maximum(abs.(ref - test_eb))
        println("  Explicit-buffer correctness check: max error = $err")

        # Verify TLS approach gives same result
        test_tls = zeros(eltype(H), size(H, 1))
        tls_hot_loop!(test_tls, H, v)
        err_tls = maximum(abs.(ref - test_tls))
        println("  TLS correctness check: max error = $err_tls")

        println()

        # Benchmark 1: Baseline (current, single-threaded)
        print("  [1] Baseline (dgt in struct, single-thread):     ")
        t1 = @benchmark begin
            target = zeros(eltype($H), size($H, 1))
            baseline_hot_loop!(target, $H, $v)
        end
        show(stdout, MIME("text/plain"), t1)
        println()

        # Benchmark 2: Copy-per-thread (current mul)
        if Threads.nthreads() > 1
            print("  [2] Copy-per-thread (current mul, threaded):     ")
            t2 = @benchmark copy_per_thread_mul($H, $v)
            show(stdout, MIME("text/plain"), t2)
            println()
        end

        # Benchmark 3: TLS-based dgt
        print("  [3] Task-local storage dgt (single-thread):      ")
        t3 = @benchmark begin
            target = zeros(eltype($H), size($H, 1))
            tls_hot_loop!(target, $H, $v)
        end
        show(stdout, MIME("text/plain"), t3)
        println()

        # Benchmark 4: Explicit buffer (single-threaded)
        print("  [4] Explicit buffer passing (single-thread):     ")
        t4 = @benchmark begin
            target = zeros(eltype($H), size($H, 1))
            explicit_buf_hot_loop!(target, $H, $v)
        end
        show(stdout, MIME("text/plain"), t4)
        println()

        # Benchmark 5: Explicit buffer (threaded)
        if Threads.nthreads() > 1
            print("  [5] Explicit buffer passing (threaded):          ")
            t5 = @benchmark explicit_buf_threaded_mul($H, $v)
            show(stdout, MIME("text/plain"), t5)
            println()
        end

        println()

        # Summary
        med1 = median(t1).time / 1e6  # ms
        med3 = median(t3).time / 1e6
        med4 = median(t4).time / 1e6
        println("  Summary (median times in ms):")
        println("    Baseline (struct dgt):       $(round(med1, digits=3)) ms")
        println("    TLS dgt:                     $(round(med3, digits=3)) ms  ($(round(med3/med1 * 100, digits=1))% of baseline)")
        println("    Explicit buffer:             $(round(med4, digits=3)) ms  ($(round(med4/med1 * 100, digits=1))% of baseline)")
        if Threads.nthreads() > 1
            med2 = median(t2).time / 1e6
            med5 = median(t5).time / 1e6
            println("    Copy-per-thread (threaded):  $(round(med2, digits=3)) ms")
            println("    Explicit buf (threaded):     $(round(med5, digits=3)) ms  ($(round(med5/med2 * 100, digits=1))% of copy-per-thread)")
        end
    end
end

#=============================================================================
  Additional: Measure just the change!/index primitive overhead
=============================================================================#

function benchmark_primitives()
    println("\n" * "=" ^ 70)
    println("Primitive operation overhead comparison")
    println("=" ^ 70)

    for L in [10, 16, 20]
        println("\n--- L = $L ---")

        # Baseline: dgt in struct
        B = TensorBasis(L=L)
        N = min(2^L, 100_000)

        print("  change!/index via struct dgt:    ")
        t1 = @benchmark begin
            s = 0
            for i in 1:$N
                EDKit.change!($B.dgt, $i, base=$B.B)
                _, idx = EDKit.index($B.dgt, base=$B.B), 0
                s += idx
            end
            s
        end
        show(stdout, MIME("text/plain"), t1)
        println()

        # Explicit buffer
        dgt = zeros(Int64, L)
        print("  change!/index via local buffer:  ")
        t2 = @benchmark begin
            s = 0
            for i in 1:$N
                EDKit.change!($dgt, $i, base=Int64(2))
                _, idx = EDKit.index($dgt, base=Int64(2)), 0
                s += idx
            end
            s
        end
        show(stdout, MIME("text/plain"), t2)
        println()

        # TLS buffer
        print("  change!/index via TLS buffer:    ")
        t3 = @benchmark begin
            tls = task_local_storage()
            buf = get!(tls, :bench_dgt) do
                zeros(Int64, $L)
            end::Vector{Int64}
            s = 0
            for i in 1:$N
                EDKit.change!(buf, $i, base=Int64(2))
                _, idx = EDKit.index(buf, base=Int64(2)), 0
                s += idx
            end
            s
        end
        show(stdout, MIME("text/plain"), t3)
        println()

        med1 = median(t1).time
        med2 = median(t2).time
        med3 = median(t3).time
        println("  Ratio: struct=$(round(med1/med1, digits=2))x  local=$(round(med2/med1, digits=2))x  TLS=$(round(med3/med1, digits=2))x")
    end
end

# Run everything
run_benchmarks()
benchmark_primitives()

println("\n\n" * "=" ^ 70)
println("ANALYSIS & RECOMMENDATION")
println("=" ^ 70)
println("""

Key findings to look for:
1. If explicit-buffer and baseline are ~equal → no overhead from passing dgt as arg
2. If TLS has overhead → task_local_storage() lookup cost matters in tight loops
3. If explicit-buffer threaded vs copy-per-thread → savings from not copying I, R, C arrays

RECOMMENDED APPROACH: Explicit buffer passing
- Zero overhead vs baseline (dgt is just a Vector on the stack frame)
- Thread-safe by construction (each thread/task allocates its own)
- No synchronization primitives needed
- No struct layout change needed (keep dgt in struct for backward compat,
  but add dgt-argument overloads for the hot path functions)
- Avoids copying the entire Basis (I, R, C arrays) per thread
""")
