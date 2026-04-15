"""
Check memory cost of sparse(opt) for various system sizes.
"""

using EDKit
using LinearAlgebra, SparseArrays

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

function report(label, opt)
    dim = size(opt, 1)
    t = @elapsed S = sparse(opt)
    nnzS = nnz(S)
    # Memory: CSC stores rowval (Int) + nzval (ComplexF64) per nnz, plus colptr (Int) per col+1
    mem_bytes = nnzS * (sizeof(Int) + sizeof(ComplexF64)) + (dim + 1) * sizeof(Int)
    mem_mb = mem_bytes / 1024^2
    sparsity = nnzS / (dim * dim)
    println("  $label: dim=$dim, nnz=$nnzS, sparsity=$(round(sparsity*100, digits=2))%, " *
            "sparse_mem=$(round(mem_mb, digits=2)) MiB, build_time=$(round(t, digits=3))s")
end

println("Memory overhead of sparse(opt) by system size:\n")

# TensorBasis
for L in [10, 12, 14]
    b = TensorBasis(L=L, base=2)
    opt = heisenberg_operator(L, b)
    report("TensorBasis L=$L", opt)
end

println()

# ProjectedBasis half-filling
for L in [14, 16, 18, 20]
    b = ProjectedBasis(L=L, base=2, f = x -> sum(x) == L÷2)
    opt = heisenberg_operator(L, b)
    report("ProjectedBasis L=$L", opt)
end

println()

# TranslationalBasis
for L in [14, 16, 18, 20]
    b = TranslationalBasis(L=L, base=2, k=0, f = x -> sum(x) == L÷2)
    opt = heisenberg_operator(L, b)
    report("TranslationalBasis L=$L k=0", opt)
end

println()

# What about the current allocation cost of opt * m?
println("Comparison: memory allocated by current opt * m:")
for (L, ncols) in [(16, 10), (18, 10), (20, 10)]
    b = ProjectedBasis(L=L, base=2, f = x -> sum(x) == L÷2)
    opt = heisenberg_operator(L, b)
    dim = size(opt, 1)
    m = randn(ComplexF64, dim, ncols)
    # Current * allocates a zeros(ctype, dim, ncols) output matrix
    out_mem = dim * ncols * sizeof(ComplexF64) / 1024^2
    println("  ProjectedBasis L=$L, ncols=$ncols: output=$(round(out_mem, digits=2)) MiB")
end
