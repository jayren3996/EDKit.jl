# Benchmark AbelianBasis performance and test integer orbit search for general base

using EDKit
using LinearAlgebra, SparseArrays
using BenchmarkTools

BenchmarkTools.DEFAULT_PARAMETERS.samples = 10
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 3

function heisenberg_op(L, basis; base=2)
    if base == 2
        Sx = [0 0.5; 0.5 0]; Sy = [0 -0.5im; 0.5im 0]; Sz = [0.5 0; 0 -0.5]
    else
        Sx = [0 1 0; 1 0 1; 0 1 0] / sqrt(2)
        Sy = [0 -1im 0; 1im 0 -1im; 0 1im 0] / sqrt(2)
        Sz = [1.0 0 0; 0 0 0; 0 0 -1.0]
    end
    h = kron(Sx, Sx) + kron(Sy, Sy) + kron(Sz, Sz)
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
    println("  $label (dim=$dim): vec=$(BenchmarkTools.prettytime(tv))  mat10=$(BenchmarkTools.prettytime(tm))")
end

# ============================================================================
# 1. Current AbelianBasis performance vs specialized bases
# ============================================================================
println("="^70)
println("  AbelianBasis vs specialized bases - base=2")
println("="^70)

for L in [12, 14, 16]
    println("\nL = $L:")

    # Specialized TranslationalBasis
    b_t = TranslationalBasis(L=L, base=2, k=0, f=x->sum(x)==L÷2)
    bench("TranslationalBasis k=0", heisenberg_op(L, b_t))

    # AbelianBasis with same symmetry (translation only)
    b_a = basis(L=L, base=2, k=0, f=x->sum(x)==L÷2)
    bench("AbelianBasis k=0", heisenberg_op(L, b_a))

    if L <= 14
        # Specialized TranslationParityBasis
        b_tp = TranslationParityBasis(L=L, base=2, k=0, p=1, f=x->sum(x)==L÷2)
        bench("TranslationParityBasis k=0,p=1", heisenberg_op(L, b_tp))

        # AbelianBasis with translation + parity
        b_ap = basis(L=L, base=2, k=0, p=1, f=x->sum(x)==L÷2)
        bench("AbelianBasis k=0,p=1", heisenberg_op(L, b_ap))

        # Specialized TranslationFlipBasis
        b_tf = TranslationFlipBasis(L=L, base=2, k=0, p=1, f=x->sum(x)==L÷2)
        bench("TranslationFlipBasis k=0,p=1", heisenberg_op(L, b_tf))

        # AbelianBasis with translation + flip
        b_af = basis(L=L, base=2, k=0, z=1, f=x->sum(x)==L÷2)
        bench("AbelianBasis k=0,z=1", heisenberg_op(L, b_af))
    end
end

# ============================================================================
# 2. AbelianBasis with base=3
# ============================================================================
println("\n\n")
println("="^70)
println("  AbelianBasis - base=3 (spin-1)")
println("="^70)

for L in [6, 8, 10]
    println("\nL = $L:")

    b_t = TranslationalBasis(L=L, base=3, k=0, f=x->sum(x)==L)
    bench("TranslationalBasis k=0", heisenberg_op(L, b_t, base=3))

    b_a = basis(L=L, base=3, k=0, f=x->sum(x)==L)
    bench("AbelianBasis k=0", heisenberg_op(L, b_a, base=3))

    if L <= 8
        b_ap = basis(L=L, base=3, k=0, p=1, f=x->sum(x)==L)
        bench("AbelianBasis k=0,p=1", heisenberg_op(L, b_ap, base=3))
    end
end

# ============================================================================
# 3. Integer orbit search simulation for AbelianBasis
# ============================================================================
println("\n\n")
println("="^70)
println("  Integer orbit search simulation (general base)")
println("="^70)

# For translation: "rotate" an integer by removing last `a` digits and prepending
@inline function int_rotate(state::Int64, L::Int, a::Int, base::Int64)
    pow = base^(L - a)  # precomputed
    # Extract last `a` digits
    tail = state % (base^a)
    # Shift remaining digits right
    head = state ÷ (base^a)
    # Prepend tail digits
    tail * pow + head
end

# For parity (reverse): reverse digit order in the integer
@inline function int_reverse(state::Int64, L::Int, base::Int64)
    r = Int64(0)
    v = state
    for _ in 1:L
        r = r * base + v % base
        v = v ÷ base
    end
    r
end

# For spin-flip (base=2): XOR with all-ones mask
@inline function int_spinflip_b2(state::Int64, L::Int)
    mask = (Int64(1) << L) - Int64(1)
    state ⊻ mask
end

# For spin-flip (general base): complement each digit
@inline function int_spinflip(state::Int64, L::Int, base::Int64)
    bL = base^L
    bL - Int64(1) - state
end

println("\nTranslation orbit search (L translations):")
for (L, base) in [(12, 2), (16, 2), (8, 3), (10, 3)]
    dgt = zeros(Int64, L)
    a = 1  # unit cell
    ntrans = L
    pow_a = Int64(base)^(L - a)
    base_a = Int64(base)^a

    # Current: circshift! + index
    b_cs = @benchmark begin
        s = Int64(0)
        for _ in 1:1000
            for k in 1:$L; $dgt[k] = rand(0:$base-1); end
            Im = EDKit.index($dgt, base=Int64($base))
            for _ in 1:$ntrans-1
                circshift!($dgt, $a)
                In = EDKit.index($dgt, base=Int64($base))
                In < Im && (Im = In)
            end
            s += Im
        end
        s
    end

    # Integer approach
    b_int = @benchmark begin
        s = Int64(0)
        for _ in 1:1000
            for k in 1:$L; $dgt[k] = rand(0:$base-1); end
            state = EDKit.index($dgt, base=Int64($base)) - Int64(1)
            Im = state
            curr = state
            for _ in 1:$ntrans-1
                tail = curr % $base_a
                curr = tail * $pow_a + curr ÷ $base_a
                curr < Im && (Im = curr)
            end
            s += Im
        end
        s
    end

    t_cs = median(b_cs).time
    t_int = median(b_int).time
    println("  L=$L base=$base: circshift=$(BenchmarkTools.prettytime(t_cs))  " *
            "integer=$(BenchmarkTools.prettytime(t_int))  " *
            "speedup=$(round(t_cs/t_int, digits=2))x")
end

println("\nTranslation + parity orbit search:")
for (L, base) in [(12, 2), (14, 2), (8, 3)]
    dgt = zeros(Int64, L)
    a = 1
    ntrans = L
    pow_a = Int64(base)^(L - a)
    base_a = Int64(base)^a
    ngroup = 2 * ntrans

    # Current: circshift! + reverse! + index
    b_cs = @benchmark begin
        s = Int64(0)
        for _ in 1:1000
            for k in 1:$L; $dgt[k] = rand(0:$base-1); end
            Im = EDKit.index($dgt, base=Int64($base))
            for _ in 1:$ntrans-1
                circshift!($dgt, $a)
                In = EDKit.index($dgt, base=Int64($base))
                In < Im && (Im = In)
            end
            reverse!($dgt)
            for _ in 1:$ntrans
                In = EDKit.index($dgt, base=Int64($base))
                In < Im && (Im = In)
                circshift!($dgt, $a)
            end
            s += Im
        end
        s
    end

    # Integer approach
    b_int = @benchmark begin
        s = Int64(0)
        for _ in 1:1000
            for k in 1:$L; $dgt[k] = rand(0:$base-1); end
            state = EDKit.index($dgt, base=Int64($base)) - Int64(1)
            Im = state
            curr = state
            for _ in 1:$ntrans-1
                tail = curr % $base_a
                curr = tail * $pow_a + curr ÷ $base_a
                curr < Im && (Im = curr)
            end
            # Reverse (parity)
            curr_r = int_reverse(state, $L, Int64($base))
            for _ in 1:$ntrans
                curr_r < Im && (Im = curr_r)
                tail = curr_r % $base_a
                curr_r = tail * $pow_a + curr_r ÷ $base_a
            end
            s += Im
        end
        s
    end

    t_cs = median(b_cs).time
    t_int = median(b_int).time
    println("  L=$L base=$base: circshift=$(BenchmarkTools.prettytime(t_cs))  " *
            "integer=$(BenchmarkTools.prettytime(t_int))  " *
            "speedup=$(round(t_cs/t_int, digits=2))x")
end

println("\nDone!")
