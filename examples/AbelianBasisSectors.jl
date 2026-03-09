using LinearAlgebra

const DEV = true

if DEV
    include("../src/EDKit.jl")
    using .EDKit
else
    using EDKit
end

L = 8
target_N = L ÷ 2
target_k = 0
bond = spin((1.0, "xx"), (1.0, "yy"))

full_sector = basis(L = L, N = target_N, k = target_k)
full_vals = eigvals(Hermitian(trans_inv_operator(bond, 2, full_sector))) |> sort

split_vals = Float64[]
split_dims = Dict{Int, Int}()

for parity in (1, -1)
    sector = basis(L = L, N = target_N, k = target_k, p = parity)
    split_dims[parity] = size(sector, 1)
    iszero(size(sector, 1)) && continue
    vals = eigvals(Hermitian(trans_inv_operator(bond, 2, sector)))
    append!(split_vals, vals)
end

split_vals = sort(split_vals)

println("Abelian basis helper for simultaneous symmetries")
println("  L                         = $L")
println("  N                         = $target_N")
println("  k                         = $target_k")
println("  dim(N, k)                 = $(size(full_sector, 1))")
println("  dim(N, k, p=+1)           = $(split_dims[1])")
println("  dim(N, k, p=-1)           = $(split_dims[-1])")
println("  sector recombination err  = $(norm(full_vals - split_vals))")
