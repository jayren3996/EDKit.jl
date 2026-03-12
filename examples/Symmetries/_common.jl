using LinearAlgebra

const DEV = true

if DEV
    include("../../src/EDKit.jl")
    using .EDKit
else
    using EDKit
end

const L = 8
const HALF_FILLING = L ÷ 2
const XXZ_BOND = spin((1.0, "xx"), (1.0, "yy"), (0.6, "zz"))

function build_sector_hamiltonian(; kwargs...)
    B = basis(L = L; kwargs...)
    H = trans_inv_operator(XXZ_BOND, 2, B)
    B, H
end

function print_sector_summary(title; kwargs...)
    B, H = build_sector_hamiltonian(; kwargs...)
    vals = eigvals(Hermitian(H)) |> sort
    println(title)
    println("  sector kwargs   = ", (; kwargs...))
    println("  Hilbert dim     = ", size(B, 1))
    println("  matrix size     = ", size(H))
    println("  ground energy   = ", vals[1])
end
