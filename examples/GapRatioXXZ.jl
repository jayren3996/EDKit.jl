using LinearAlgebra
using Random

const DEV = true

if DEV
    include("../src/EDKit.jl")
    using .EDKit
else
    using EDKit
end

Random.seed!(7)

L = 10
delta = 1.2
disorder = 0.8
fields = disorder .* randn(L)

sector = basis(L = L, N = L ÷ 2)
exchange = spin((1.0, "xx"), (1.0, "yy"), (delta, "zz"))

H = trans_inv_operator(exchange, 2, sector)
H += operator([h * spin("Z") for h in fields], collect(1:L), sector)

vals = eigvals(Hermitian(H)) |> sort
mid = length(vals) ÷ 2

println("Random-field XXZ chain in the zero-magnetization sector")
println("  system size                  = $L")
println("  sector dimension             = $(size(sector, 1))")
println("  mean gap ratio               = $(meangapratio(vals))")
println("  central five eigenvalues     = $(vals[mid-2:mid+2])")
