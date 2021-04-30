using LinearAlgebra
include("../src/EDKit.jl")
using .EDKit
using Test
const printstate = false
#-------------------------------------------------------------------------------------------------------------------------
# Spin-1/2 XY Model
#-------------------------------------------------------------------------------------------------------------------------
@testset "Spin-1/2 XY" begin
    L = 10
    if printstate
        println("------------------------------------------------------------")
        println("Spin-1/2 XY Model with L = $L")
        println("------------------------------------------------------------")
    end
    θ = 0.34
    expθ = exp(-1im*θ)
    mat = spin((expθ, "+-"), (1/expθ, "-+"), (1, "z1"), (1, "1z"))
    E = zeros(2^L)
    P = 0
    for n = 0:L, k = 0:L-1
        basis = translationalbasis(x->sum(x)==n, k, L)
        if printstate 
            println("N = $n, k = $k: $(size(basis,1)) states.") 
        end
        if (l = size(basis,1)) > 0
            vals = trans_inv_operator(mat, 2, basis) |> Array |> Hermitian |> eigvals
            E[P+1:P+l] = vals
            P += l
        end
    end
    @test P == 2^L
    vals = trans_inv_operator(mat, 2, L) |> Array |> Hermitian |> eigvals
    @test vals ≈ sort!(E)
end
#-------------------------------------------------------------------------------------------------------------------------
# Spin-1/2 Random Model
#-------------------------------------------------------------------------------------------------------------------------
@testset "Spin-1/2 Random" begin
    L = 10
    if printstate 
        println("------------------------------------------------------------")
        println("Spin-1/2 Random Model with L = $L")
        println("------------------------------------------------------------")
    end
    mat = begin
        M = zeros(ComplexF64, 8, 8)
        M[1,1] = rand()
        M[8,8] = rand()
        M[[2,3,5], [2,3,5]] = rand(ComplexF64, 3, 3) |> Hermitian
        M[[4,6,7], [4,6,7]] = rand(ComplexF64, 3, 3) |> Hermitian
        M
    end
    E = zeros(2^L)
    P = 0
    for n = 0:L, k = 0:L-1
        basis = translationalbasis(x->sum(x)==n, k, L)
        if (l = size(basis, 1)) > 0
            if printstate 
                println("N = $n, k = $k: $l states.")
            end
            vals = trans_inv_operator(mat, 3, basis) |> Array |> Hermitian |> eigvals
            E[P+1:P+l] = vals
            P += l
        end
    end
    @test P == 2^L
    vals = trans_inv_operator(mat, 3, L) |> Array |> Hermitian |> eigvals
    @test norm(vals - sort(E)) ≈ 0.0 atol = 1e-10

end

#-------------------------------------------------------------------------------------------------------------------------
# Spin-1 XY Model
#-------------------------------------------------------------------------------------------------------------------------
@testset "Spin-1 XY" begin
    L = 6
    if printstate 
        println("------------------------------------------------------------")
        println("Spin-1 XY Model with L = $L")
        println("------------------------------------------------------------")
    end
    h = 1
    mat = spin((1, "xx"), (1, "yy"), (h/2, "1z"), (h/2, "z1"), D=3)
    E = zeros(3^L)
    P = 0
    for n = 0:2L, k = 0:L-1
        basis = translationalbasis(x->sum(x)==n, k, L, base=3)
        if printstate 
            println("N = $n, k = $k: $(size(basis, 1)) states.")
        end
        if (l = size(basis, 1)) > 0
            vals = trans_inv_operator(mat, 2, basis) |> Array |> Hermitian |> eigvals
            E[P+1:P+l] = vals
            P += l
        end
    end
    @test P == 3^L
    vals = trans_inv_operator(mat, 2, L) |> Array |> Hermitian |> eigvals
    @test vals ≈ sort!(E)
end

#-------------------------------------------------------------------------------------------------------------------------
# AKLT Model
#-------------------------------------------------------------------------------------------------------------------------
@testset "AKLT" begin
    L = 6
    if printstate 
        println("------------------------------------------------------------")
        println("AKLT Model with L = $L")
        println("------------------------------------------------------------")
    end
    mat = begin
        ss = spin((1, "xx"), (1, "yy"), (1, "zz"), D=3) 
        1/2 * ss + 1/6 * ss^2 + 1/3 * I
    end
    E = zeros(3^L)
    P = 0
    for n = 0:2L, k = 0:L-1
        basis = translationalbasis(x->sum(x)==n, k, L, base=3)
        if printstate 
            println("N = $n, k = $k: $(size(basis, 1)) states.")
        end
        if (l = size(basis, 1)) > 0
            vals = trans_inv_operator(mat, 2, basis) |> Array |> Hermitian |> eigvals
            E[P+1:P+l] = vals
            P += l
        end
    end
    vals = trans_inv_operator(mat, 2, L) |> Array |> Hermitian |> eigvals
    @test vals ≈ sort!(E)
end

#-------------------------------------------------------------------------------------------------------------------------
# Spin-1 Random
#-------------------------------------------------------------------------------------------------------------------------
@testset "Spin-1 Random" begin
    L = 6
    if printstate 
        println("------------------------------------------------------------")
        println("Spin-1 Random Model with L = $L")
        println("------------------------------------------------------------")
    end
    mat = begin    
        m2 = [2, 4, 10]
        m1 = [3, 5, 7, 11, 13, 19]
        m0 = [6, 8, 12, 14, 16, 20, 22]
        p1 = [9, 15, 17, 21, 23, 25]
        p2 = [18, 24, 26]
        M = zeros(ComplexF64, 27, 27)
        M[1, 1]    = rand()
        M[27, 27]  = rand()
        M[m2, m2] .= rand(ComplexF64, 3, 3) |> Hermitian
        M[m1, m1] .= rand(ComplexF64, 6, 6) |> Hermitian
        M[m0, m0] .= rand(ComplexF64, 7, 7) |> Hermitian
        M[p1, p1] .= rand(ComplexF64, 6, 6) |> Hermitian
        M[p2, p2] .= rand(ComplexF64, 3, 3) |> Hermitian
        M
    end
    E = zeros(3^L)
    P = 0
    for n = 0:2L, k = 0:L-1
        basis = translationalbasis(x->sum(x)==n, k, L, base=3)
        if printstate 
            println("N = $n, k = $k: $(size(basis, 1)) states.")
        end
        if (l = size(basis, 1)) > 0
            vals = trans_inv_operator(mat, 3, basis) |> Array |> Hermitian |> eigvals
            E[P+1:P+l] = vals
            P += l
        end
    end
    vals = trans_inv_operator(mat, 3, L) |> Array |> Hermitian |> eigvals
    @test vals ≈ sort!(E)
end
