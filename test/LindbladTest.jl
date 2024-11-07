include("../src/EDKit.jl")
using LinearAlgebra
using .EDKit
using Test

#---------------------------------------------------------------------------------------------------
# Tested model:
#     H  = ∑ᵢ rᵢ(XᵢXᵢ₊₁+YᵢYᵢ₊₁)  ⟷  -i/4 ∑ᵢ 2rᵢ(ωᵢωᵢ₊ₙ₊₁ - ωᵢ₊ₙωᵢ₊₁ - ωᵢ₊ₙ₊₁ωᵢ + ωᵢ₊₁ωᵢ₊ₙ)
#     Lᵢ = aᵢKᵢXᵢ + bᵢKᵢYᵢ       ⟷  aᵢωᵢ + bᵢωᵢ₊ₙ
#     Mᵢ = cᵢZᵢ                  ⟷  -i cᵢ/2 (ωᵢωᵢ₊ₙ - ωᵢ₊ₙωᵢ)
#
# Initial state:
#     |ψ₀⟩ = |↑↓⋯↑↓⟩             ⟷  |10⋯10⟩
#---------------------------------------------------------------------------------------------------
function singlebody(L, r, a, b, c)
    H1 = begin
        mat = zeros(2L, 2L)
        for i = 1:L-1
            mat[i, i+1+L] = +2r[i]
            mat[i+1+L, i] = -2r[i]
            mat[i+1, i+L] = +2r[i]
            mat[i+L, i+1] = -2r[i]
        end
        mat
    end
    L1 = begin
        mat = zeros(ComplexF64, 2L, L)
        for i = 1:L 
            mat[i, i] = a[i]
            mat[i+L, i] = b[i]
        end
        mat
    end
    M1 = begin
        mats = Vector{Matrix{Float64}}(undef, L)
        for i = 1:L 
            mat = zeros(2L,2L)
            mat[i, i+L] = +c[i] / 2 
            mat[i+L, i] = -c[i] / 2
            mats[i] = mat 
        end
        mats
    end
    QL = quadraticlindblad(H1, L1, M1)
    Γ = begin
        Z2 = [mod(i,2) for i=1:L]
        covariancematrix(Z2)
    end
    QL, Γ
end
#---------------------------------------------------------------------------------------------------
function manybody(L, r, a, b, c)
    H2 = begin
        mats = [spin((r[i],"XX"), (r[i],"YY"), D=2) for i=1:L-1]
        inds = [[i,i+1] for i=1:L-1]
        operator(mats, inds, L) |> Array
    end
    L2 = begin
        mats = Vector{Matrix{ComplexF64}}(undef, 2L)
        for i = 1:L 
            Ki = if i == 1
                I(1)
            elseif i == 2
                [-1 0; 0 1] 
            else
                kron(fill([-1 0; 0 1], i-1)...)
            end
            mat = spin((a[i], "X"), (b[i], "Y"), D=2)
            Ir = I(2^(L-i))
            mats[i] = kron(Ki, mat, Ir)
        end
        for i=1:L 
            mats[L+i] = operator(c[i] * spin("Z"), [i], L) |> Array
        end
        mats
    end
    EL = lindblad(H2, L2)
    ρ = begin
        Z2 = [mod(i-1,2) for i=1:L]
        ind = index(Z2)
        densitymatrix(ind, L)
    end
    EL, ρ
end

#---------------------------------------------------------------------------------------------------
# Helper
#---------------------------------------------------------------------------------------------------
function density(dm)
    L = round(Int, log(2, size(dm.ρ, 1)))
    n = zeros(L)
    Base.Threads.@threads for i=1:L 
        ni = operator([1 0;0 0], [i], L) |> Hermitian
        n[i] = expectation(ni, dm)
    end
    n
end

#---------------------------------------------------------------------------------------------------
# main
#---------------------------------------------------------------------------------------------------
function test(steps=30)
    L = 10
    r, a, b, c = rand(L-1), rand(ComplexF64, L), rand(ComplexF64, L), rand(L)
    # r, a, b, c = rand(L-1), rand(L), rand(L), rand(L)
    QL, Γ = singlebody(L, r, a, b, c)
    EL, ρ = manybody(L, r, a, b, c)

    n1 = zeros(steps, L)
    for i=1:steps 
        Γ = QL(Γ)
        n1[i, :] = fermioncorrelation(Γ, 1) |> diag |> real
    end

    n2 = zeros(steps, L)
    for i=1:steps 
        ρ = EL(ρ)
        n2[i,:] = density(ρ)
        println("Finish $i/$steps, error = $(norm(n2[i,:]-n1[i,:])).")
    end
    println("Total error = $(norm(n2-n1)).")
    n1, n2
end

n1, n2 = test()
