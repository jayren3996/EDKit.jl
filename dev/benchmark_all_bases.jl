# Final benchmark: all basis types, base=2 and base=3
# Compare current (post-optimization) vs the pre-optimization baseline

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
    bv = @benchmark $opt * $v
    tv = median(bv).time
    println("  $label (dim=$dim): $(BenchmarkTools.prettytime(tv))")
end

println("="^70)
println("  ALL BASIS TYPES - base=2 (spin-1/2)")
println("="^70)

for L in [12, 14, 16]
    println("\nL = $L:")
    bench("TensorBasis", heisenberg_op(L, TensorBasis(L=L, base=2)))

    b = ProjectedBasis(L=L, base=2, f=x->sum(x)==L÷2)
    bench("ProjectedBasis", heisenberg_op(L, b))

    b = TranslationalBasis(L=L, base=2, k=0, f=x->sum(x)==L÷2)
    bench("TranslationalBasis k=0", heisenberg_op(L, b))

    b = ParityBasis(L=L, base=2, p=1, f=x->sum(x)==L÷2)
    bench("ParityBasis p=1", heisenberg_op(L, b))

    b = FlipBasis(L=L, base=2, p=1, f=x->sum(x)==L÷2)
    bench("FlipBasis p=1", heisenberg_op(L, b))

    if L <= 14
        b = ParityFlipBasis(L=L, base=2, p=1, z=1, f=x->sum(x)==L÷2)
        bench("ParityFlipBasis p=1,z=1", heisenberg_op(L, b))

        b = TranslationParityBasis(L=L, base=2, k=0, p=1, f=x->sum(x)==L÷2)
        bench("TranslationParityBasis k=0,p=1", heisenberg_op(L, b))

        b = TranslationFlipBasis(L=L, base=2, k=0, p=1, f=x->sum(x)==L÷2)
        bench("TranslationFlipBasis k=0,p=1", heisenberg_op(L, b))
    end
end

println("\n\n")
println("="^70)
println("  ALL BASIS TYPES - base=3 (spin-1)")
println("="^70)

for L in [6, 8, 10]
    println("\nL = $L:")
    bench("TensorBasis", heisenberg_op(L, TensorBasis(L=L, base=3), base=3))

    b = ProjectedBasis(L=L, base=3, f=x->sum(x)==L)
    bench("ProjectedBasis Sz=0", heisenberg_op(L, b, base=3))

    b = TranslationalBasis(L=L, base=3, k=0, f=x->sum(x)==L)
    bench("TranslationalBasis k=0", heisenberg_op(L, b, base=3))

    b = ParityBasis(L=L, base=3, p=1, f=x->sum(x)==L)
    bench("ParityBasis p=1", heisenberg_op(L, b, base=3))

    b = FlipBasis(L=L, base=3, p=1, f=x->sum(x)==L)
    bench("FlipBasis p=1", heisenberg_op(L, b, base=3))

    if L <= 8
        b = TranslationParityBasis(L=L, base=3, k=0, p=1, f=x->sum(x)==L)
        bench("TranslationParityBasis k=0,p=1", heisenberg_op(L, b, base=3))

        b = TranslationFlipBasis(L=L, base=3, k=0, p=1, f=x->sum(x)==L)
        bench("TranslationFlipBasis k=0,p=1", heisenberg_op(L, b, base=3))
    end
end

println("\nDone!")
