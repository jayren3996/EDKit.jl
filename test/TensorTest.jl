include("../src/EDKit.jl")
using Main.EDKit
using LinearAlgebra, ITensors, ITensorMPS
using Test

#----------------------------------------------------------------------------------------------------
# Basics
#----------------------------------------------------------------------------------------------------
@testset "check tensor" begin
    for L in 2:5, d = 2:4
        # random many-body state
        v = randn(ComplexF64, d^L) |> normalize
        T = EDKit.vec2tensor(v, base=d)
        for i in eachindex(v)
            ind = digits(i-1, base=d, pad=L) .+ 1 |> reverse
            @test v[i] ≈ T[ind...]
        end
    end
end
#----------------------------------------------------------------------------------------------------
@testset "MPS <--> vec" begin
    for L in 2:5, d = 2:4
        v = randn(ComplexF64, d^L) |> normalize
        T = EDKit.vec2tensor(v, base=d)

        s = siteinds(d, L)
        ψ = MPS(T, s)

        vec = mps2vec(ψ)
        for i in eachindex(v)
            @test vec[i] ≈ v[i]
        end
        psi = vec2mps(v, s)
        @test norm(psi) ≈ norm(ψ)
        @test abs(inner(psi, ψ)) ≈ norm(ψ)^2
    end
end
#----------------------------------------------------------------------------------------------------
@testset "op <--> mat" begin
    for L in 2:5, d = 2:4
        v = randn(ComplexF64, d^L) |> normalize
        T = EDKit.vec2tensor(v, base=d)

        s = siteinds(d, L)
        ψ = MPS(T, s)
        mat = Matrix(qr(randn(ComplexF64, d^2, d^2)).Q)
        v2 = kron(mat, I(d^(L-2))) * v 
        ψ2 = apply(mat2op(mat, s[1], s[2]), ψ) |> mps2vec 
        @test v2 ≈ ψ2 
        @test op2mat(EDKit.mat2op(mat, s[1], s[2]), s[1], s[2]) ≈ mat
    end
end
#----------------------------------------------------------------------------------------------------
@testset "Product MPS" begin
    for L = 2:4, d=2:4
        sites = siteinds(d, L)
        states = [normalize(rand(d)) for i in 1:L]

        v = kron(states...)
        ψ = productstate(sites, states)

        vec = mps2vec(ψ)
        for i in eachindex(v)
            @test vec[i] ≈ v[i]
        end
    end
end

#----------------------------------------------------------------------------------------------------
# AKLT Test
#----------------------------------------------------------------------------------------------------
λ = EDKit.λ
ITensors.op(::OpName"λ0",::SiteType"S=1") = λ(0)
ITensors.op(::OpName"λ1",::SiteType"S=1") = λ(1)
ITensors.op(::OpName"λ2",::SiteType"S=1") = λ(2)
ITensors.op(::OpName"λ3",::SiteType"S=1") = λ(3)
ITensors.op(::OpName"λ4",::SiteType"S=1") = λ(4)
ITensors.op(::OpName"λ5",::SiteType"S=1") = λ(5)
ITensors.op(::OpName"λ6",::SiteType"S=1") = λ(6)
ITensors.op(::OpName"λ7",::SiteType"S=1") = λ(7)
ITensors.op(::OpName"λ8",::SiteType"S=1") = λ(8)
ITensors.op(::OpName"++",::SiteType"S=1") = [0 0 1; 0 0 0; 0 0 0]
ITensors.op(::OpName"--",::SiteType"S=1") = [0 0 0; 0 0 0; 1 0 0]
#----------------------------------------------------------------------------------------------------
const AKLT_H2 = let h = zeros(9, 9)
    h[1, 1] = h[9, 9] = 1
    h[[2, 4], [2, 4]] .= 0.5
    h[[6, 8], [6, 8]] .= 0.5
    h[[3,5,7], [3,5,7]] .= [1; 2; 1] * [1 2 1] / 6
    h
end
#----------------------------------------------------------------------------------------------------
const AKLT_TENSOR = let A = zeros(2, 2, 3)
    A[1, 1, 2] = -1 / sqrt(2)
    A[2, 2, 2] = 1 / sqrt(2)
    A[1, 2, 1] = 1
    A[2, 1, 3] = -1
    A
end

#----------------------------------------------------------------------------------------------------
function aklt_mps(sites::AbstractVector; l::Int=0, r::Int=0)
    tensors = if iszero(l) && iszero(r)
        fill(AKLT_TENSOR, length(sites))
    else
        Tl = AKLT_TENSOR[[l], :, :]
        Tr = AKLT_TENSOR[:, [r], :]
        [Tl; fill(AKLT_TENSOR, length(sites)-2); Tr]
    end
    EDKit.pbcmps(sites, tensors)
end
#----------------------------------------------------------------------------------------------------
function aklt_tw(sites, k=length(sites)÷2)
    L = length(sites)
    os = OpSum()
    for i in 1:L 
        os += exp(-1im * k * 2π/L * (i-1)), "++", i
    end
    MPO(os, sites)
end
#----------------------------------------------------------------------------------------------------
@testset "AKLT Ground State" begin
    for L in 2:10
        B = TranslationalBasis(L=L, N=L, k=0, base=3)
        s = siteinds("S=1", L)
        psi = aklt_mps(s)
        ψ = mps2vec(psi, B)
        @test abs(norm(ψ)-1) < 1e-5

        H = trans_inv_operator(AKLT_H2, 2, B)
        @test norm(H * ψ) < 1e-5 
    end
end
#----------------------------------------------------------------------------------------------------
@testset "AKLT Scar Tower" begin
    L = 8
    s = siteinds("S=1", L)
    psi = aklt_mps(s)
    tw = aklt_tw(s)
    for j in 1:L÷2
        psi = apply(tw, psi) |> normalize!

        B = TranslationalBasis(L=L, N=L+2j, k=mod(j*L÷2, L), base=3)
        H = trans_inv_operator(AKLT_H2, 2, B)
        
        ψ = mps2vec(psi, B)
        @test abs(norm(ψ)-1) < 1e-5

        H = trans_inv_operator(AKLT_H2, 2, B)
        Hψ = H * ψ
        @test abs(norm(Hψ) - dot(ψ, Hψ) ) < 1e-5 
    end
end

#----------------------------------------------------------------------------------------------------
# Pauli Basis Test
#----------------------------------------------------------------------------------------------------
σ = EDKit.σ
@testset "Pauli Matrices" begin
    name = ["1", "X", "Y", "Z"]
    for L in 1:5
        for i in 0:4^L-1 
            inds = digits(i, base=4, pad=L)
            @test σ(inds) ≈ spin(prod(name[j+1] for j in inds))
        end
    end
end
#----------------------------------------------------------------------------------------------------
@testset "Pauli Coefficients" begin
    name = ["1", "X", "Y", "Z"]
    for L in 1:5
        c = randn(4^L)
        mat = sum(c[i] * pauli(i, L) for i in eachindex(c))
        plist = pauli_list(mat)
        @test plist ≈ c
    end
end
#----------------------------------------------------------------------------------------------------
@testset "Lindbladian" begin
    for i in 1:10
        A = randn(ComplexF64, 8, 8) |> Hermitian 
        Ac = commutation_mat(A)

        B = randn(ComplexF64, 8, 8) |> Hermitian 
        Bl = pauli_list(B)

        C = randn(ComplexF64, 8, 8)
        Cd = dissipation_mat(C)

        @test pauli(Ac * Bl) ≈ -1im * (A * B - B* A)
        @test pauli(Cd * Bl) ≈ C * B * C' - (C' * C * B + B * C' * C) / 2
    end
end

#----------------------------------------------------------------------------------------------------
# Pauli MPS
#----------------------------------------------------------------------------------------------------
@testset "Fidelity" begin
    for L in 2:7
        s = siteinds("S=1/2", L)
        ps = siteinds("Pauli", L)
        vec = rand(ComplexF64, 2^L) |> normalize!
        ρ = vec * vec' 
        pmps = vec2mps(pauli_list(ρ), ps) 
        mpo = pmps2mpo(pmps, s)
        ψ = vec2mps(vec, s)
        @test inner(ψ', mpo, ψ) ≈ 1.0
    end
end
#----------------------------------------------------------------------------------------------------
@testset "MPS -> PMPS" begin
    for L in 2:7
        s = siteinds("S=1/2", L)
        ps = siteinds("Pauli", L)
        vec = rand(ComplexF64, 2^L) |> normalize!
        ψ = vec2mps(vec, s)
        pmps = mps2pmps(ψ, ps)
        mpo = pmps2mpo(pmps, s)
        @test inner(ψ', mpo, ψ) ≈ 1.0
    end
end
