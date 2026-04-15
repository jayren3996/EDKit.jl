"""
Benchmark: Operator * Matrix performance in EDKit.jl

Tests current performance and validates proposed optimizations:
1. Current `*` (single-threaded matrix-free)
2. Current `mul` (multi-threaded matrix-free)
3. Proposed: sparse(opt) * m  (auto-densification via SparseArrays SpMM)
4. Proposed: column-by-column mul! (cache-friendly)
5. Proposed: Array(opt) * m (dense BLAS, upper bound reference)
"""

using EDKit
using LinearAlgebra, SparseArrays
using BenchmarkTools

BenchmarkTools.DEFAULT_PARAMETERS.samples = 5
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 3

#=============================================================================
  Helper: build a Heisenberg chain operator on a given basis
=============================================================================#
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

#=============================================================================
  Proposed optimization: column-by-column mul!
=============================================================================#
function mul_colwise!(target::AbstractMatrix, opt, m::AbstractMatrix)
    for k = 1:size(m, 2)
        EDKit.mul!(view(target, :, k), opt, view(m, :, k))
    end
    target
end

function star_colwise(opt, m::AbstractMatrix)
    ctype = promote_type(eltype(opt), eltype(m))
    M = zeros(ctype, size(opt, 1), size(m, 2))
    mul_colwise!(M, opt, m)
end

#=============================================================================
  Run benchmarks for one basis type
=============================================================================#
function run_benchmark(label, opt, ncols_list)
    dim = size(opt, 1)
    println("\n", "="^70)
    println("  $label  (dim = $dim)")
    println("="^70)

    # Pre-build sparse and dense representations
    print("  Building sparse(opt)... ")
    t_sparse = @elapsed S = sparse(opt)
    println("$(round(t_sparse, digits=3))s  (nnz = $(nnz(S)))")

    print("  Building Array(opt)... ")
    t_dense = @elapsed D = Array(opt)
    println("$(round(t_dense, digits=3))s")

    for ncols in ncols_list
        println("\n  --- ncols = $ncols ---")
        m = randn(ComplexF64, dim, ncols)

        # Reference: dense BLAS
        ref = D * m

        # 1. Current: opt * m (single-threaded matrix-free)
        r1 = opt * m
        @assert norm(r1 - ref) / norm(ref) < 1e-10 "opt * m incorrect!"
        b1 = @benchmark $opt * $m
        t1 = median(b1).time
        println("  [current]  opt * m          : $(BenchmarkTools.prettytime(t1))  ($(BenchmarkTools.prettymemory(median(b1).memory)))")

        # 2. Current: mul(opt, m) (multi-threaded)
        t2 = nothing
        if Threads.nthreads() > 1
            r2 = mul(opt, m)
            @assert norm(r2 - ref) / norm(ref) < 1e-10 "mul(opt, m) incorrect!"
            b2 = @benchmark mul($opt, $m)
            t2 = median(b2).time
            println("  [current]  mul(opt, m)      : $(BenchmarkTools.prettytime(t2))  ($(BenchmarkTools.prettymemory(median(b2).memory)))")
        end

        # 3. Proposed: sparse(opt) * m
        r3 = S * m
        @assert norm(r3 - ref) / norm(ref) < 1e-10 "sparse * m incorrect!"
        b3 = @benchmark $S * $m
        t3 = median(b3).time
        println("  [proposed] sparse(opt) * m  : $(BenchmarkTools.prettytime(t3))  ($(BenchmarkTools.prettymemory(median(b3).memory)))")

        # 4. Proposed: column-by-column mul!
        r4 = star_colwise(opt, m)
        @assert norm(r4 - ref) / norm(ref) < 1e-10 "colwise mul! incorrect!"
        b4 = @benchmark star_colwise($opt, $m)
        t4 = median(b4).time
        println("  [proposed] colwise mul!     : $(BenchmarkTools.prettytime(t4))  ($(BenchmarkTools.prettymemory(median(b4).memory)))")

        # 5. Reference: dense BLAS
        b5 = @benchmark $D * $m
        t5 = median(b5).time
        println("  [reference] Array(opt) * m  : $(BenchmarkTools.prettytime(t5))  ($(BenchmarkTools.prettymemory(median(b5).memory)))")

        # Summary
        println("  Speedup vs current: sparse=$(round(t1/t3, digits=1))x, " *
                "colwise=$(round(t1/t4, digits=1))x, " *
                "BLAS=$(round(t1/t5, digits=1))x")
    end
end

#=============================================================================
  Main
=============================================================================#
println("Julia version: $(VERSION)")
println("Julia threads: $(Threads.nthreads())")
println("BLAS threads:  $(BLAS.get_num_threads())")

ncols_list = [1, 5, 10, 50]

# --- TensorBasis (full space, no symmetry) ---
for L in [10, 12]
    b = TensorBasis(L=L, base=2)
    opt = heisenberg_operator(L, b)
    run_benchmark("TensorBasis L=$L", opt, ncols_list)
end

# --- ProjectedBasis (half-filling sector) ---
for L in [14, 16]
    b = ProjectedBasis(L=L, base=2, f = x -> sum(x) == L÷2)
    opt = heisenberg_operator(L, b)
    run_benchmark("ProjectedBasis L=$L (half-filling)", opt, ncols_list)
end

# --- TranslationalBasis (momentum sector, most expensive index()) ---
for L in [14, 16]
    b = TranslationalBasis(L=L, base=2, k=0, f = x -> sum(x) == L÷2)
    opt = heisenberg_operator(L, b)
    run_benchmark("TranslationalBasis L=$L k=0", opt, ncols_list)
end

println("\n\nDone!")
