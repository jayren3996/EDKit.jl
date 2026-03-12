using LinearAlgebra

const DEV = true

if DEV
    include("../../src/EDKit.jl")
    using .EDKit
else
    using EDKit
end

L = 8
half_filling = L ÷ 2
special_k = (0, L ÷ 2)

# This nearest-neighbor XXZ chain preserves:
# - translation
# - reflection parity
# - spin-flip
# - total Sz / particle number
bond = spin((1.0, "xx"), (1.0, "yy"), (0.6, "zz"))

sector_spectrum(; kwargs...) = begin
    B = basis(L = L; kwargs...)
    vals = iszero(size(B, 1)) ? Float64[] : eigvals(Hermitian(trans_inv_operator(bond, 2, B))) |> sort
    B, vals
end

function print_sector(label; kwargs...)
    B, vals = sector_spectrum(; kwargs...)
    energy = isempty(vals) ? "empty" : string(vals[1])
    println(rpad(label, 18), " dim = ", lpad(size(B, 1), 4), "   ground = ", energy)
end

function recombination_error(parts)
    ref_kwargs, child_kwargs = parts
    _, ref_vals = sector_spectrum(; ref_kwargs...)
    vals = Float64[]
    for kw in child_kwargs
        _, child_vals = sector_spectrum(; kw...)
        append!(vals, child_vals)
    end
    isempty(ref_vals) ? 0.0 : norm(sort(vals) - ref_vals)
end

println("Symmetry reduction examples with basis(...)\n")
println("System size L = $L")
println("Model: H = sum_i (XX + YY + 0.6 ZZ)")
println()

println("Single symmetries")
print_sector("none")
print_sector("N"; N = half_filling)
print_sector("p = +1"; p = 1)
print_sector("p = -1"; p = -1)
print_sector("z = +1"; z = 1)
print_sector("z = -1"; z = -1)
for k in 0:L-1
    print_sector("k = $k"; k = k)
end

println()
println("Two-symmetry combinations")
print_sector("N + p = +1"; N = half_filling, p = 1)
print_sector("N + p = -1"; N = half_filling, p = -1)
print_sector("N + z = +1"; N = half_filling, z = 1)
print_sector("N + z = -1"; N = half_filling, z = -1)
print_sector("p=+1, z=+1"; p = 1, z = 1)
print_sector("p=+1, z=-1"; p = 1, z = -1)
print_sector("p=-1, z=+1"; p = -1, z = 1)
print_sector("p=-1, z=-1"; p = -1, z = -1)
for k in 0:L-1
    print_sector("N + k = $k"; N = half_filling, k = k)
    print_sector("k = $k, z = +1"; k = k, z = 1)
    print_sector("k = $k, z = -1"; k = k, z = -1)
end
for k in special_k
    print_sector("k = $k, p = +1"; k = k, p = 1)
    print_sector("k = $k, p = -1"; k = k, p = -1)
end

println()
println("Three-symmetry combinations")
print_sector("N, p=+1, z=+1"; N = half_filling, p = 1, z = 1)
print_sector("N, p=+1, z=-1"; N = half_filling, p = 1, z = -1)
print_sector("N, p=-1, z=+1"; N = half_filling, p = -1, z = 1)
print_sector("N, p=-1, z=-1"; N = half_filling, p = -1, z = -1)
for k in 0:L-1
    print_sector("N, k=$k, z=+1"; N = half_filling, k = k, z = 1)
    print_sector("N, k=$k, z=-1"; N = half_filling, k = k, z = -1)
end
for k in special_k
    print_sector("N, k=$k, p=+1"; N = half_filling, k = k, p = 1)
    print_sector("N, k=$k, p=-1"; N = half_filling, k = k, p = -1)
    print_sector("k=$k, p=+1, z=+1"; k = k, p = 1, z = 1)
    print_sector("k=$k, p=+1, z=-1"; k = k, p = 1, z = -1)
    print_sector("k=$k, p=-1, z=+1"; k = k, p = -1, z = 1)
    print_sector("k=$k, p=-1, z=-1"; k = k, p = -1, z = -1)
end

println()
println("Four-symmetry combinations")
for k in special_k, p in (1, -1), z in (1, -1)
    print_sector("N, k=$k, p=$p, z=$z"; N = half_filling, k = k, p = p, z = z)
end

println()
println("Compatibility rules")
println("  N with z requires half filling: N = L/2 for base-2 systems.")
println("  k with p requires k = 0 or k = L/2.")
println("  k with p and z also requires k = 0 or k = L/2.")
println("  N + k + p + z therefore requires both half filling and k in {0, L/2}.")

println()
println("A few recombination checks")
println("  none <- N sectors error         = ", recombination_error(((;), [(; N = n) for n in 0:L])))
println("  none <- p sectors error         = ", recombination_error(((;), [(; p = p) for p in (1, -1)])))
println("  none <- z sectors error         = ", recombination_error(((;), [(; z = z) for z in (1, -1)])))
println("  none <- k sectors error         = ", recombination_error(((;), [(; k = k) for k in 0:L-1])))
println("  N <- (N,p) sectors error        = ", recombination_error(((; N = half_filling), [(; N = half_filling, p = p) for p in (1, -1)])))
println("  N <- (N,z) sectors error        = ", recombination_error(((; N = half_filling), [(; N = half_filling, z = z) for z in (1, -1)])))
println("  (N,k=0) <- parity sectors error = ", recombination_error(((; N = half_filling, k = 0), [(; N = half_filling, k = 0, p = p) for p in (1, -1)])))
println("  (N,k=0) <- flip sectors error   = ", recombination_error(((; N = half_filling, k = 0), [(; N = half_filling, k = 0, z = z) for z in (1, -1)])))
println("  (N,k=0) <- p,z sectors error    = ", recombination_error(((; N = half_filling, k = 0), [(; N = half_filling, k = 0, p = p, z = z) for p in (1, -1) for z in (1, -1)])))
