# Focused benchmark for the basis-consolidation refactor.
# Measures hot-path `index(b, dgt)` and basis construction time.
# Run before and after each step; compare numbers.

using EDKit
using BenchmarkTools
using Random

BenchmarkTools.DEFAULT_PARAMETERS.samples = 200
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 2

const RNG = MersenneTwister(20260521)

function bench_index(b, L; base=2)
    dgt = rand(RNG, 0:base-1, L)
    t = @belapsed EDKit.index($b, $dgt)
    t * 1e9  # ns
end

function bench_index_buf(b, L; base=2)
    # 1-arg variant: uses b.dgt
    dgt = rand(RNG, 0:base-1, L)
    copyto!(b.dgt, dgt)
    t = @belapsed EDKit.index($b)
    t * 1e9
end

function bench_construct(ctor)
    t = @belapsed $ctor()
    t * 1e6  # μs
end

function row(name, t_ns)
    println("  ", rpad(name, 38), lpad(round(t_ns, digits=1), 10), " ns")
end

function rowc(name, t_us)
    println("  ", rpad(name, 38), lpad(round(t_us, digits=2), 10), " μs")
end

println("Julia: ", VERSION)
println("Threads: ", Threads.nthreads())
println()

println("="^60)
println("  index(b, dgt) hot path  [L=14, base=2, N=7 half-filling]")
println("="^60)
L = 14
N = 7

bp = ProjectedBasis(L=L, N=N, threaded=false)
row("ProjectedBasis", bench_index(bp, L))

bf = SpinlessFermionBasis(L=L, N=N, threaded=false)
row("SpinlessFermionBasis", bench_index(bf, L))

bpar = ParityBasis(L=L, p=1, N=N, threaded=false)
row("ParityBasis p=1", bench_index(bpar, L))

bfl = FlipBasis(L=L, p=1, threaded=false)
row("FlipBasis p=1", bench_index(bfl, L))

bpf = ParityFlipBasis(L=L, p=1, z=1, threaded=false)
row("ParityFlipBasis p=1,z=1", bench_index(bpf, L))

bt = TranslationalBasis(L=L, k=0, N=N, threaded=false)
row("TranslationalBasis k=0", bench_index(bt, L))

btp = TranslationParityBasis(L=L, k=0, p=1, N=N, threaded=false)
row("TranslationParityBasis", bench_index(btp, L))

btf = TranslationFlipBasis(L=L, k=0, p=1, threaded=false)
row("TranslationFlipBasis", bench_index(btf, L))

println()
println("="^60)
println("  index(b) 1-arg shim  [L=14, base=2]")
println("="^60)

row("ProjectedBasis", bench_index_buf(bp, L))
row("SpinlessFermionBasis", bench_index_buf(bf, L))
row("ParityBasis", bench_index_buf(bpar, L))
row("FlipBasis", bench_index_buf(bfl, L))
row("ParityFlipBasis", bench_index_buf(bpf, L))
row("TranslationalBasis", bench_index_buf(bt, L))
row("TranslationParityBasis", bench_index_buf(btp, L))
row("TranslationFlipBasis", bench_index_buf(btf, L))

println()
println("="^60)
println("  basis construction  [L=14, base=2, N=7]")
println("="^60)

rowc("ProjectedBasis(L=14,N=7)",
     bench_construct(() -> ProjectedBasis(L=L, N=N, threaded=false)))
rowc("SpinlessFermionBasis(L=14,N=7)",
     bench_construct(() -> SpinlessFermionBasis(L=L, N=N, threaded=false)))
rowc("ParityBasis(L=14,p=1,N=7)",
     bench_construct(() -> ParityBasis(L=L, p=1, N=N, threaded=false)))
rowc("FlipBasis(L=14,p=1)",
     bench_construct(() -> FlipBasis(L=L, p=1, threaded=false)))
rowc("ParityFlipBasis(L=14,p=1,z=1)",
     bench_construct(() -> ParityFlipBasis(L=L, p=1, z=1, threaded=false)))
rowc("TranslationalBasis(L=14,k=0,N=7)",
     bench_construct(() -> TranslationalBasis(L=L, k=0, N=N, threaded=false)))
rowc("TranslationParityBasis(L=14,k=0,p=1,N=7)",
     bench_construct(() -> TranslationParityBasis(L=L, k=0, p=1, N=N, threaded=false)))
rowc("TranslationFlipBasis(L=14,k=0,p=1)",
     bench_construct(() -> TranslationFlipBasis(L=L, k=0, p=1, threaded=false)))

println()
println("Done.")
