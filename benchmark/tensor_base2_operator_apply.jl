using EDKit
using LinearAlgebra

try
    @eval using BenchmarkTools
catch err
    @error "BenchmarkTools is required for this script. Install it in a temporary benchmark environment, not as a package dependency." exception=(err, catch_backtrace())
    exit(1)
end

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 1
BenchmarkTools.DEFAULT_PARAMETERS.samples = 10_000

function old_style_mul_vector(H, v)
    ctype = promote_type(eltype(H), eltype(v))
    target = zeros(ctype, size(H, 1))
    dgt = similar(H.B.dgt)
    for j = 1:length(v)
        EDKit.colmn!(target, H, j, dgt, v[j], 1, nothing)
    end
    target
end

function old_style_threaded_mul_vector(H, v)
    ctype = promote_type(eltype(H), eltype(v))
    nt = min(Threads.nthreads(), max(1, length(v)))
    ni = EDKit.dividerange(length(v), nt)
    Ms = [zeros(ctype, size(H, 1)) for _ in 1:nt]
    Threads.@threads for i in 1:nt
        dgt = similar(H.B.dgt)
        for j in ni[i]
            EDKit.colmn!(Ms[i], H, j, dgt, v[j], 1, nothing)
        end
    end
    target = Ms[1]
    @inbounds for i in 2:nt
        axpy!(one(ctype), Ms[i], target)
    end
    target
end

function old_style_mul_matrix(H, M)
    ctype = promote_type(eltype(H), eltype(M))
    target = zeros(ctype, size(H, 1), size(M, 2))
    dgt = similar(H.B.dgt)
    for j = 1:size(M, 1)
        EDKit._colmn_row!(target, H, j, dgt, M, j, 1, nothing)
    end
    target
end

function old_style_threaded_mul_matrix(H, M)
    ctype = promote_type(eltype(H), eltype(M))
    nt = min(Threads.nthreads(), max(1, size(M, 1)))
    ni = EDKit.dividerange(size(M, 1), nt)
    Ms = [zeros(ctype, size(H, 1), size(M, 2)) for _ in 1:nt]
    Threads.@threads for i in 1:nt
        dgt = similar(H.B.dgt)
        for j in ni[i]
            EDKit._colmn_row!(Ms[i], H, j, dgt, M, j, 1, nothing)
        end
    end
    target = Ms[1]
    @inbounds for i in 2:nt
        axpy!(one(ctype), Ms[i], target)
    end
    target
end

function benchmark_pair(label, old_fn, new_fn)
    old_time = @belapsed $old_fn()
    new_time = @belapsed $new_fn()
    println("  ", rpad(label, 36), round(old_time * 1e3, digits=3), " ms -> ",
        round(new_time * 1e3, digits=3), " ms  (",
        round(old_time / new_time, digits=2), "x)")
end

function benchmark_case(L)
    bond = spin((1.0, "xx"), (1.0, "yy"), (0.2, "zz"))
    H = trans_inv_operator(bond, 2, TensorBasis(L = L, base = 2))
    v = randn(ComplexF64, size(H, 2))
    M = randn(ComplexF64, size(H, 2), 8)

    println("TensorBasis(base=2) L=", L, " size=", size(H), " terms=", length(H))
    EDKit.clear_sparse_cache!()
    benchmark_pair("single-thread H*v", () -> old_style_mul_vector(H, v), () -> H * v)
    benchmark_pair("threaded EDKit.mul(H,v)", () -> old_style_threaded_mul_vector(H, v), () -> EDKit.mul(H, v))

    EDKit.clear_sparse_cache!()
    benchmark_pair("matrix RHS uncached H*M", () -> old_style_mul_matrix(H, M), () -> H * M)
    benchmark_pair("matrix RHS uncached EDKit.mul", () -> old_style_threaded_mul_matrix(H, M), () -> EDKit.mul(H, M))

    EDKit.clear_sparse_cache!()
    EDKit.sparse!(H)
    benchmark_pair("matrix RHS sparse! cache H*M", () -> H * M, () -> H * M)
    benchmark_pair("matrix RHS sparse! cache EDKit.mul", () -> EDKit.mul(H, M), () -> EDKit.mul(H, M))
    EDKit.clear_sparse_cache!()
    println()
end

function main()
    println("Julia: ", VERSION)
    println("Threads: ", Threads.nthreads())
    println("Columns in matrix RHS: 8")
    println("Format: old generic digit-buffer time -> current time (speedup)")
    println()
    for L in (10, 12, 14)
        benchmark_case(L)
    end
end

main()
