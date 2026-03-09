using LinearAlgebra

const DEV = true

if DEV
    include("../src/EDKit.jl")
    using .EDKit
else
    using EDKit
end

function pxp_obc_constraint(v::AbstractVector{<:Integer})
    all(v[i] == 0 || v[i + 1] == 0 for i in 1:length(v)-1)
end

L = 10

local_term = let projector = Diagonal([1, 1, 1, 0, 1, 1, 0, 0])
    projector * kron(I(2), spin("X"), I(2)) * projector
end

basis = ProjectedBasis(L = L, f = pxp_obc_constraint)
inds = [[i, i + 1, i + 2] for i in 1:L-2]
H = operator(fill(local_term, L - 2), inds, basis)

neel_state = [isodd(i) ? 1 : 0 for i in 1:L]
psi0 = productstate(neel_state, basis)

vals, vecs = eigen(Hermitian(H))
overlaps = abs2.(vecs' * psi0)
scar_index = argmax(overlaps)

println("PXP chain with open boundary conditions")
println("  system size                = $L")
println("  constrained Hilbert space  = $(size(basis, 1))")
println("  ground-state energy        = $(vals[1])")
println("  strongest Neel overlap     = $(overlaps[scar_index])")
println("  entropy of that eigenstate = $(ent_S(vecs[:, scar_index], 1:L÷2, basis))")
