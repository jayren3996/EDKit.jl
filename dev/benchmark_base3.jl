# Benchmark: base=2 vs base=3 matrix-free performance
# Tests all basis types to see where the base=2 specializations help
# and how much base=3 (generic path) lags behind.

using EDKit
using LinearAlgebra, SparseArrays
using BenchmarkTools

BenchmarkTools.DEFAULT_PARAMETERS.samples = 10
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 3

function heisenberg_op(L, basis; base=2)
    if base == 2
        Sx = [0 0.5; 0.5 0]; Sy = [0 -0.5im; 0.5im 0]; Sz = [0.5 0; 0 -0.5]
        h = kron(Sx, Sx) + kron(Sy, Sy) + kron(Sz, Sz)
    else
        # Spin-1 Heisenberg: S=1 operators (3x3)
        s1 = 1.0
        Sx = [0 1 0; 1 0 1; 0 1 0] / sqrt(2)
        Sy = [0 -1im 0; 1im 0 -1im; 0 1im 0] / sqrt(2)
        Sz = [1.0 0 0; 0 0 0; 0 0 -1.0]
        h = kron(Sx, Sx) + kron(Sy, Sy) + kron(Sz, Sz)
    end
    mats = [sparse(h) for _ in 1:L]
    inds = [[i, mod1(i+1, L)] for i in 1:L]
    operator(mats, [collect(ind) for ind in inds], basis)
end

function bench(label, opt)
    dim = size(opt, 1)
    dim == 0 && return
    v = randn(ComplexF64, dim)
    m = randn(ComplexF64, dim, 10)

    bv = @benchmark $opt * $v
    bm = @benchmark $opt * $m

    tv = median(bv).time
    tm = median(bm).time
    println("  $label (dim=$dim): vec=$(BenchmarkTools.prettytime(tv))  " *
            "mat10=$(BenchmarkTools.prettytime(tm))")
end

println("="^70)
println("  BASE=2 (spin-1/2) -- with bit optimizations")
println("="^70)

for L in [12, 14, 16]
    println("\nL = $L:")
    bench("TensorBasis", heisenberg_op(L, TensorBasis(L=L, base=2)))

    b = ProjectedBasis(L=L, base=2, f=x->sum(x)==L÷2)
    bench("ProjectedBasis (half-fill)", heisenberg_op(L, b))

    b = TranslationalBasis(L=L, base=2, k=0, f=x->sum(x)==L÷2)
    bench("TranslationalBasis k=0", heisenberg_op(L, b))

    if L <= 14
        b = TranslationParityBasis(L=L, base=2, k=0, p=1, f=x->sum(x)==L÷2)
        bench("TranslationParityBasis k=0,p=1", heisenberg_op(L, b))

        b = TranslationFlipBasis(L=L, base=2, k=0, p=1, f=x->sum(x)==L÷2)
        bench("TranslationFlipBasis k=0,p=1", heisenberg_op(L, b))
    end

    b = ParityBasis(L=L, base=2, p=1, f=x->sum(x)==L÷2)
    bench("ParityBasis p=1", heisenberg_op(L, b))

    b = FlipBasis(L=L, base=2, p=1, f=x->sum(x)==L÷2)
    bench("FlipBasis p=1", heisenberg_op(L, b))
end

println("\n\n")
println("="^70)
println("  BASE=3 (spin-1) -- generic path, no bit optimizations")
println("="^70)

for L in [6, 8, 10]
    println("\nL = $L:")
    bench("TensorBasis", heisenberg_op(L, TensorBasis(L=L, base=3), base=3))

    b = ProjectedBasis(L=L, base=3, f=x->sum(x)==L)
    bench("ProjectedBasis (Sz=0)", heisenberg_op(L, b, base=3))

    b = TranslationalBasis(L=L, base=3, k=0, f=x->sum(x)==L)
    bench("TranslationalBasis k=0", heisenberg_op(L, b, base=3))

    if L <= 8
        b = TranslationParityBasis(L=L, base=3, k=0, p=1, f=x->sum(x)==L)
        bench("TranslationParityBasis k=0,p=1", heisenberg_op(L, b, base=3))
    end

    b = ParityBasis(L=L, base=3, p=1, f=x->sum(x)==L)
    bench("ParityBasis p=1", heisenberg_op(L, b, base=3))
end

println("\n\n")
println("="^70)
println("  INTEGER ORBIT SEARCH TEST: base=3 circshift! vs integer approach")
println("="^70)

# Test: how much would an integer orbit search help for base=3?
# Simulate both approaches manually.
println("\nSimulating orbit search cost (L translations, repeated 10000 times):")

for (L, base) in [(10, 2), (12, 2), (8, 3), (10, 3)]
    dgt = zeros(Int64, L)
    ntrans = L  # single-site unit cell

    # circshift! + index approach
    b_cs = @benchmark begin
        for _ in 1:10000
            # load a random state
            for k in 1:$L; $dgt[k] = rand(0:$base-1); end
            Im = EDKit.index($dgt, base=Int64($base))
            for _ in 1:$ntrans-1
                circshift!($dgt, 1)
                In = EDKit.index($dgt, base=Int64($base))
                In < Im && (Im = In)
            end
        end
    end

    # Integer approach (works for any base)
    # Convert dgt → integer, then "rotate" by recomputing with shifted indices
    b_int = @benchmark begin
        for _ in 1:10000
            for k in 1:$L; $dgt[k] = rand(0:$base-1); end
            state = EDKit.index($dgt, base=Int64($base)) - Int64(1)
            Im = state
            # For general base, "rotation" is: remove last digit, shift others, prepend removed digit
            # state_new = last_digit * base^(L-1) + (state ÷ base)
            pow = Int64($base)^($L-1)
            curr = state
            for _ in 1:$ntrans-1
                last = curr % Int64($base)
                curr = last * pow + curr ÷ Int64($base)
                curr < Im && (Im = curr)
            end
        end
    end

    t_cs = median(b_cs).time
    t_int = median(b_int).time
    println("  L=$L base=$base: circshift=$(BenchmarkTools.prettytime(t_cs))  " *
            "integer=$(BenchmarkTools.prettytime(t_int))  " *
            "speedup=$(round(t_cs/t_int, digits=2))x")
end

println("\nDone!")
