include("../src/EDKit.jl")
using .EDKit
using LinearAlgebra
using Test

#-------------------------------------------------------------------------------------------------------------------------
# Tesnsor Basis
#-------------------------------------------------------------------------------------------------------------------------
println("--------------------------------------------------")
println("                   Tensor Basis                   ")
println("--------------------------------------------------")

const X = [0 1; 1 0]
const Y = [0 -1; 1 0] * 1im
const Z = [1 0; 0 -1]

directkron(mat, i, j, L, base=2) = kron(I(base^(i-1)), mat, I(base^(L-j)))
boundarykron(m1, m2, i, j, base=2) = kron(m1, I(base^(j-i+1)), m2)

@testset "Base-2 Random Kron" begin
    L = 10
    for i=1:10
        mats = [rand(2^i, 2^i) for j=1:L-i+1]
        inds = [j:j+i-1 for j=1:L-i+1]
        M_test = operator(mats, inds, L) |> Array 
        M_targ = sum(directkron(mats[j], j, j+i-1, L) for j = 1:L-i+1)
        @test M_test ≈ M_targ
    end
end

@testset "Base-3 Random Kron" begin
    L = 6
    for i=1:6
        mats = [rand(3^i, 3^i) for j=1:L-i+1]
        inds = [j:j+i-1 for j=1:L-i+1]
        M_test = operator(mats, inds, L) |> Array 
        M_targ = sum(directkron(mats[j], j, j+i-1, L, 3) for j = 1:L-i+1)
        @test M_test ≈ M_targ
    end
end

@testset "Base-3 Boundary Kron" begin
    L = 6
    for i = 1:3, j = 1:3
        m1, m2 = rand(3^i, 3^i), rand(3^j, 3^j)
        mats = kron(m1, m2)
        inds = vcat(Vector(1:i), Vector(L-j+1:L))
        M_test = operator(mats, inds, L) |> Array 
        M_targ = boundarykron(m1, m2, i+1, L-j, 3)
        @test M_test ≈ M_targ
    end
end

@testset "Operator Multiplication" begin
    L = 10
    mat = kron(I(2), X, Y, Z)
    v = rand(2^L)
    opt = trans_inv_operator(mat, 4, L)
    v2 = opt * v
    optm = Array(opt)
    v3 = optm * v
    @test v2 ≈ v3
end

#-------------------------------------------------------------------------------------------------------------------------
# Projected Basis
#-------------------------------------------------------------------------------------------------------------------------
println("--------------------------------------------------")
println("                 Projected Basis                  ")
println("--------------------------------------------------")

@testset "Spin-1/2 XY" begin
    L = 10
    mat = spin((1, "xx"), (1, "yy"))
    E = zeros(2^L)
    P = 0
    for n = 0:L
        basis = ProjectedBasis(x->sum(x)==n, L)
        if (l = size(basis, 1)) > 0
            vals = trans_inv_operator(mat, 2, basis) |> Array |> Hermitian |> eigvals
            E[P+1:P+l] = vals
            P += l
        end
    end
    @test P == 2^L
    vals = trans_inv_operator(mat, 2, L) |> Array |> Hermitian |> eigvals
    @test vals ≈ sort!(E)
end

@testset "Spin-1/2 Random" begin
    L = 10
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
    for n = 0:L
        basis = ProjectedBasis(x->sum(x)==n, L)
        if (l = size(basis, 1)) > 0
            vals = trans_inv_operator(mat, 3, basis) |> Array |> Hermitian |> eigvals
            E[P+1:P+l] = vals
            P += l
        end
    end
    @test P == 2^L
    vals = trans_inv_operator(mat, 3, L) |> Array |> Hermitian |> eigvals
    @test vals ≈ sort!(E)
end

@testset "PXP_OBC" begin
    function fib(n::Integer)
        if n == 0
            return 1
        elseif n == 1
            return 2
        end
        a, b = 1, 2
        for m = 2:n
            a, b = b, a+b
        end
        b
    end
    mat = begin
        P = Diagonal([1, 1, 1, 0, 1, 1, 0, 0])
        P * kron(I(2), X, I(2)) * P
    end
    pxpf(v::Vector{<:Integer}) = all(v[i]==0 || v[i+1]==0 for i=1:length(v)-1)
    for L = 3:20
        basis = ProjectedBasis(pxpf, L)
        @test size(basis, 1) == fib(L)
    end
    for L = 3:10
        basis = ProjectedBasis(pxpf, L)
        mats = fill(mat, L-2)
        inds = [[i,i+1, i+2] for i=1:L-2]
        E = operator(mats, inds, basis) |> Array |> Hermitian |> eigvals
        H = operator(mats, inds, L) |> Array
        vals = H[basis.I, basis.I] |> Hermitian |> eigvals
        @test norm(vals - E) ≈ 0.0
    end
end

@testset "PXP_PBC Projected Basis" begin
    mat = begin
        P = Diagonal([1, 1, 1, 0, 1, 1, 0, 0])
        P * kron(I(2), X, I(2)) * P
    end
    pxpf(v::Vector{<:Integer}) = all(v[i]==0 || v[mod(i, length(v))+1]==0 for i=1:length(v))
    for L = 3:14
        basis = ProjectedBasis(pxpf, L)
        E = trans_inv_operator(mat, 3, basis) |> Array |> Hermitian |> eigvals
        H = trans_inv_operator(mat, 3, L) |> Array
        vals = H[basis.I, basis.I] |> Hermitian |> eigvals
        @test norm(vals - E) ≈ 0.0
    end
end

@testset "Spin-1 Random" begin
    L = 6
    R3 = begin    
        m2 = [2, 4, 10]
        m1 = [3, 5, 7, 11, 13, 19]
        m0 = [6, 8, 12, 14, 16, 20, 22]
        p1 = [9, 15, 17, 21, 23, 25]
        p2 = [18, 24, 26]
        mat = zeros(ComplexF64, 27, 27)
        mat[1, 1]    = rand()
        mat[27, 27]  = rand()
        mat[m2, m2] .= rand(ComplexF64, 3, 3) |> Hermitian
        mat[m1, m1] .= rand(ComplexF64, 6, 6) |> Hermitian
        mat[m0, m0] .= rand(ComplexF64, 7, 7) |> Hermitian
        mat[p1, p1] .= rand(ComplexF64, 6, 6) |> Hermitian
        mat[p2, p2] .= rand(ComplexF64, 3, 3) |> Hermitian
        mat
    end
    E = zeros(3^L)
    P = 0
    for n = 0:2L
        basis = ProjectedBasis(x->sum(x)==n, L, base=3)
        if (l = size(basis, 1)) > 0
            vals = trans_inv_operator(mat, 3, basis) |> Array |> Hermitian |> eigvals
            E[P+1:P+l] = vals
            P += l
        end
    end
    @test P == 3^L
    vals = trans_inv_operator(mat, 3, L) |> Array |> Hermitian |> eigvals
    @test vals ≈ sort!(E)
end

#-------------------------------------------------------------------------------------------------------------------------
# Translation Basis
#-------------------------------------------------------------------------------------------------------------------------
println("--------------------------------------------------")
println("               Translational Basis                ")
println("--------------------------------------------------")

@testset "Spin-1/2 XY" begin
    L = 10
    θ = 0.34
    expθ = exp(-1im*θ)
    mat = spin((expθ, "+-"), (1/expθ, "-+"), (1, "z1"), (1, "1z"))
    E = zeros(2^L)
    P = 0
    for n = 0:L, k = 0:L-1
        basis = TranslationalBasis(x->sum(x)==n, k, L)
        #println(size(basis,1))
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

@testset "Spin-1/2 Random" begin
    L = 10
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
        basis = TranslationalBasis(x->sum(x)==n, k, L)
        if (l = size(basis, 1)) > 0
            vals = trans_inv_operator(mat, 3, basis) |> Array |> Hermitian |> eigvals
            E[P+1:P+l] = vals
            P += l
        end
    end
    @test P == 2^L
    vals = trans_inv_operator(mat, 3, L) |> Array |> Hermitian |> eigvals
    @test norm(vals - sort(E)) ≈ 0.0 atol = 1e-10
end

@testset "Spin-1 XY" begin
    L = 6
    h = 1
    mat = spin((1, "xx"), (1, "yy"), (h/2, "1z"), (h/2, "z1"), D=3)
    E = zeros(3^L)
    P = 0
    for n = 0:2L, k = 0:L-1
        basis = TranslationalBasis(x->sum(x)==n, k, L, base=3)
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

@testset "AKLT" begin
    L = 6
    mat = begin
        ss = spin((1, "xx"), (1, "yy"), (1, "zz"), D=3) 
        1/2 * ss + 1/6 * ss^2 + 1/3 * I
    end
    E = zeros(3^L)
    P = 0
    for n = 0:2L, k = 0:L-1
        basis = TranslationalBasis(x->sum(x)==n, k, L, base=3)
        if (l = size(basis, 1)) > 0
            vals = trans_inv_operator(mat, 2, basis) |> Array |> Hermitian |> eigvals
            E[P+1:P+l] = vals
            P += l
        end
    end
    vals = trans_inv_operator(mat, 2, L) |> Array |> Hermitian |> eigvals
    @test vals ≈ sort!(E)
end

@testset "Spin-1 Random" begin
    L = 6
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
        basis = TranslationalBasis(x->sum(x)==n, k, L, base=3)
        if (l = size(basis, 1)) > 0
            vals = trans_inv_operator(mat, 3, basis) |> Array |> Hermitian |> eigvals
            E[P+1:P+l] = vals
            P += l
        end
    end
    vals = trans_inv_operator(mat, 3, L) |> Array |> Hermitian |> eigvals
    @test vals ≈ sort!(E)
end

#-------------------------------------------------------------------------------------------------------------------------
# Translation + Parity
#-------------------------------------------------------------------------------------------------------------------------
println("--------------------------------------------------")
println("             Translation Parity Basis             ")
println("--------------------------------------------------")

@testset "Spin-1/2 XY" begin
    mat = spin((1, "+-"), (1, "-+"), (1, "z1"), (1, "1z"))
    for L = 2:2:10
        for n = 0:L, k in [0, L ÷ 2]
            be = TranslationParityBasis(x->sum(x)==n, k, +1, L)
            bo = TranslationParityBasis(x->sum(x)==n, k, -1, L)
            ba = TranslationalBasis(x->sum(x)==n, k, L)
            ve = trans_inv_operator(mat, 2, be) |> Array |> Hermitian |> eigvals
            vo = trans_inv_operator(mat, 2, bo) |> Array |> Hermitian |> eigvals
            va = trans_inv_operator(mat, 2, ba) |> Array |> Hermitian |> eigvals
            E = sort(vcat(ve, vo))
            @test norm(E-va) ≈ 0.0 atol = 1e-12
        end
    end
end

@testset "Spin-1 XY" begin
    mat = spin((1, "+-"), (1, "-+"), (1, "z1"), (1, "1z"), D=3)
    for L = 2:2:6
        for n = 0:2L, k in [0, L ÷ 2]
            be = TranslationParityBasis(x->sum(x)==n, k, +1, L, base=3)
            bo = TranslationParityBasis(x->sum(x)==n, k, -1, L, base=3)
            ba = TranslationalBasis(x->sum(x)==n, k, L, base=3)
            ve = trans_inv_operator(mat, 2, be) |> Array |> Hermitian |> eigvals
            vo = trans_inv_operator(mat, 2, bo) |> Array |> Hermitian |> eigvals
            va = trans_inv_operator(mat, 2, ba) |> Array |> Hermitian |> eigvals
            E = sort(vcat(ve, vo))
            @test norm(E-va) ≈ 0.0 atol = 1e-12
        end
    end
end

@testset "PXP" begin
    mat = begin
        P = Diagonal([1, 1, 1, 0, 1, 1, 0, 0])
        P * kron(I(2), X, I(2)) * P
    end
    pxpf(v::Vector{<:Integer}) = all(v[i]==0 || v[mod(i, length(v))+1]==0 for i=1:length(v))
    for L = 4:2:20
        for k in [0, L ÷ 2]
            be = TranslationParityBasis(pxpf, k, +1, L)
            bo = TranslationParityBasis(pxpf, k, -1, L)
            ba = TranslationalBasis(pxpf, k, L)
            ve = trans_inv_operator(mat, 2, be) |> Array |> Hermitian |> eigvals
            vo = trans_inv_operator(mat, 2, bo) |> Array |> Hermitian |> eigvals
            va = trans_inv_operator(mat, 2, ba) |> Array |> Hermitian |> eigvals
            E = sort(vcat(ve, vo))
            @test norm(E-va) ≈ 0.0 atol = 1e-12
        end
    end
end

@testset "Ising Spin-flip" begin
    zz = spin((1, "zz"))
    x = spin((0.3,"x"))
    for L = 4:2:12
        E = zeros(2^L)
        for k in 0:L-1
            ba = TranslationalBasis(k, L)
            be = TranslationFlipBasis(k, 1, L)
            bo = TranslationFlipBasis(k, -1, L)
            va = trans_inv_operator(zz, 2, ba) + trans_inv_operator(x, 1, ba) |> Array |> Hermitian |> eigvals
            ve = trans_inv_operator(zz, 2, be) + trans_inv_operator(x, 1, be) |> Array |> Hermitian |> eigvals
            vo = trans_inv_operator(zz, 2, bo) + trans_inv_operator(x, 1, bo) |> Array |> Hermitian |> eigvals
            veo = sort(vcat(ve, vo))
            @test norm(va-veo) ≈ 0.0 atol = 1e-12
        end
    end
end
