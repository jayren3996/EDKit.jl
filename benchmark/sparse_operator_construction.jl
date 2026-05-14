using EDKit
using SparseArrays

try
    @eval using BenchmarkTools
catch err
    @error "BenchmarkTools is required for this script. Install it in a temporary benchmark environment, not as a package dependency." exception=(err, catch_backtrace())
    exit(1)
end

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 1
BenchmarkTools.DEFAULT_PARAMETERS.samples = 10_000

function old_style_sparse(opt)
    M = spzeros(eltype(opt), size(opt)...)
    if size(M, 1) > 0 && size(M, 2) > 0
        EDKit.addto!(M, opt)
    end
    M
end

function benchmark_case(name, H)
    EDKit.clear_sparse_cache!()
    old_time = @belapsed old_style_sparse($H)
    array_time = @belapsed Array($H)
    sparse_time = @belapsed SparseArrays.sparse($H)
    cached_time = @belapsed begin
        EDKit.clear_sparse_cache!()
        EDKit.sparse!($H)
    end
    EDKit.clear_sparse_cache!()

    println(name)
    println("  size: ", size(H), ", terms: ", length(H))
    println("  old spzeros + addto!: ", round(old_time * 1e3, digits=3), " ms")
    println("  Array(H):             ", round(array_time * 1e3, digits=3), " ms")
    println("  sparse(H):            ", round(sparse_time * 1e3, digits=3), " ms")
    println("  sparse!(H):           ", round(cached_time * 1e3, digits=3), " ms")
    println()
end

function main()
    bond = spin((1.0, "+-"), (1.0, "-+"), (0.2, "zz"))
    shift6 = [2, 3, 4, 5, 6, 1]

    cases = [
        "TensorBasis L=10" => trans_inv_operator(bond, 2, TensorBasis(L = 10, base = 2)),
        "ProjectedBasis L=12 N=6" => trans_inv_operator(bond, 2, ProjectedBasis(L = 12, N = 6, threaded = false)),
        "TranslationalBasis L=12 N=6 k=0" => trans_inv_operator(bond, 2, TranslationalBasis(L = 12, N = 6, k = 0, threaded = false)),
        "Abelian k/p/z L=8 N=4" => trans_inv_operator(bond, 2, basis(L = 8, N = 4, k = 0, p = 1, z = 1, threaded = false)),
        "Custom symmetry L=6 N=3" => trans_inv_operator(bond, 2, basis(L = 6, N = 3, symmetries = [(shift6, 0)], threaded = false)),
    ]

    println("Julia: ", VERSION)
    println("Threads: ", Threads.nthreads())
    println()

    for (name, H) in cases
        benchmark_case(name, H)
    end
end

main()
