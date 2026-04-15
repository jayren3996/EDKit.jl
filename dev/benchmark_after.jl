"""
Benchmark after optimization: test sparse! cache acceleration.
"""

using EDKit
using LinearAlgebra, SparseArrays
using BenchmarkTools

BenchmarkTools.DEFAULT_PARAMETERS.samples = 5
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 3

function heisenberg_operator(L, basis)
    mats = []
    inds = []
    Sx = [0 0.5; 0.5 0]
    Sy = [0 -0.5im; 0.5im 0]
    Sz = [0.5 0; 0 -0.5]
    for i in 1:L
        j = mod1(i+1, L)
        push!(mats, kron(Sx, Sx) + kron(Sy, Sy) + kron(Sz, Sz))
        push!(inds, [i, j])
    end
    operator(sparse.(mats), [collect(ind) for ind in inds], basis)
end

function run_benchmark(label, opt, ncols_list)
    dim = size(opt, 1)
    println("\n", "="^70)
    println("  $label  (dim = $dim)")
    println("="^70)

    D = Array(opt)

    for ncols in ncols_list
        m = randn(ComplexF64, dim, ncols)
        ref = D * m

        # Before sparse!: matrix-free
        r1 = opt * m
        @assert norm(r1 - ref) / norm(ref) < 1e-10
        b1 = @benchmark $opt * $m
        t1 = median(b1).time

        println("  ncols=$ncols  before sparse!: $(BenchmarkTools.prettytime(t1))")
    end

    # Now cache with sparse!
    print("  Calling sparse!(opt)... ")
    t_build = @elapsed sparse!(opt)
    println("$(round(t_build, digits=3))s")

    for ncols in ncols_list
        m = randn(ComplexF64, dim, ncols)
        ref = D * m

        # After sparse!: uses cached SpMM
        r2 = opt * m
        @assert norm(r2 - ref) / norm(ref) < 1e-10
        b2 = @benchmark $opt * $m
        t2 = median(b2).time

        println("  ncols=$ncols  after  sparse!: $(BenchmarkTools.prettytime(t2))")
    end

    # Test vector path (should be unchanged, no cache used)
    v = randn(ComplexF64, dim)
    ref_v = D * v
    r_v = opt * v
    @assert norm(r_v - ref_v) / norm(ref_v) < 1e-10
    println("  Vector path: correct (not affected by sparse! cache)")

    clear_sparse_cache!()
    println("  Cache cleared.")
end

println("Julia version: $(VERSION)")
println("Julia threads: $(Threads.nthreads())")

ncols_list = [1, 5, 10, 50]

for L in [10, 12]
    b = TensorBasis(L=L, base=2)
    opt = heisenberg_operator(L, b)
    run_benchmark("TensorBasis L=$L", opt, ncols_list)
end

for L in [14, 16]
    b = ProjectedBasis(L=L, base=2, f = x -> sum(x) == L÷2)
    opt = heisenberg_operator(L, b)
    run_benchmark("ProjectedBasis L=$L (half-filling)", opt, ncols_list)
end

for L in [14, 16]
    b = TranslationalBasis(L=L, base=2, k=0, f = x -> sum(x) == L÷2)
    opt = heisenberg_operator(L, b)
    run_benchmark("TranslationalBasis L=$L k=0", opt, ncols_list)
end

println("\nDone!")
