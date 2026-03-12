using LinearAlgebra
using ITensors
using ITensorMPS

const DEV = true

if DEV
    include("../src/EDKit.jl")
    using .EDKit
else
    using EDKit
end

L = 5
s = siteinds(2, L)
ps = siteinds("Pauli", L)

psi_vec = randn(ComplexF64, 2^L) |> normalize
psi_mps = vec2mps(psi_vec, s)
psi_back = mps2vec(psi_mps)

pmps = mps2pmps(psi_mps, ps)
rho_mpo = pmps2mpo(pmps, s)
pmpo = mpo2pmpo(rho_mpo, ps)

bond_dims = [linkdim(psi_mps, b) for b in 1:length(psi_mps)-1]

println("MPS and Pauli-space utilities")
println("  system size                 = $L")
println("  vec2mps round-trip error    = $(norm(psi_vec - psi_back))")
println("  max MPS bond dimension      = $(maximum(bond_dims; init = 1))")
println("  Pauli MPS length            = $(length(pmps))")
println("  density MPO length          = $(length(rho_mpo))")
println("  Pauli MPO length            = $(length(pmpo))")
