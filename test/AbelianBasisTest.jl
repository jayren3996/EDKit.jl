include("../src/EDKit.jl")
using Main.EDKit, LinearAlgebra, Test, Profile

a = 1.0 
b = randn()
c = randn()
d = randn()
e = randn()
f = randn()

mat = spin(
    (a, "xx11"), (a, "yy11"), 
    (b, "x1x1"), (b, "y1y1"), 
    (c, "x11x"), (c, "y11y"), 
    (d, "zz11"),
    (e, "z1z1"), 
    (f, "z11z")
)

function add_value!(E, B)
    iszero(size(B, 1)) && return 
    vals = trans_inv_operator(mat, 4, B) |> Hermitian |> eigvals
    append!(E, vals)
end

@testset "1-sym-vals" begin
    for L in 4:10 
        B = basis(;L)
        E = trans_inv_operator(mat, 4, B) |> Hermitian |> eigvals
        
        # N total 
        Es = []
        for N in 0:L 
            B = basis(;L, N)
            add_value!(Es, B)
        end
        @test sort(Es) ≈ E

        # P
        Es = []
        for p in [1, -1] 
            B = basis(;L, p)
            add_value!(Es, B)
        end
        @test sort(Es) ≈ E 

        # Z
        Es = []
        for z in [1, -1] 
            B = basis(;L, z)
            add_value!(Es, B)
        end
        @test sort(Es) ≈ E 

        # k 
        Es = []
        for k in 0:L-1
            B = basis(;L, k)
            add_value!(Es, B)
        end
        @test sort(Es) ≈ E
    end
end

@testset "2-sym-vals" begin
    for L in 4:10 
        E = trans_inv_operator(mat, 4, L) |> Hermitian |> eigvals
        
        # N + P 
        Es = []
        for N in 0:L, p in [1, -1]
            B = basis(;L, N, p)
            add_value!(Es, B)
        end
        @test sort(Es) ≈ E 

        # N + k
        Es = []
        for N in 0:L, k in 0:L-1
            B = basis(;L, N, k)
            add_value!(Es, B)
        end
        @test sort(Es) ≈ E 

        # N + Z 

        # P + Z
        Es = []
        for p in [1, -1], z in [1, -1]
            B = basis(;L, p, z)
            add_value!(Es, B)
        end
        @test sort(Es) ≈ E 

        # k + Z
        Es = []
        for k in 0:L-1, z in [1, -1]
            B = basis(;L, k, z)
            add_value!(Es, B)
        end
        @test sort(Es) ≈ E 

        # k + P
        
    end
end


@testset "Half-k" begin
    for L in 4:2:10, k in [0, L÷2]
        B = basis(;L, k); iszero(size(B, 1)) && continue  
        E = trans_inv_operator(mat, 4, B) |> Hermitian |> eigvals

        # k + P 
        Es = []
        for p in [1, -1]
            B = basis(;L, k, p)
            add_value!(Es, B)
        end
        @test sort(Es) ≈ E 

        # k + P + Z
        Es = []
        for p in [1, -1], z in [1, -1]
            B = basis(;L, k, p, z)
            add_value!(Es, B)
        end
        @test sort(Es) ≈ E 

        # k + P + N
        Es = []
        for N in 0:L, p in [1, -1]
            B = basis(;L, N, k, p)
            add_value!(Es, B)
        end
        @test sort(Es) ≈ E 
    end
end

@testset "Half-filling" begin
    for L in 4:2:10
        N = L÷2
        B = basis(;L, N); iszero(size(B, 1)) && continue  
        E = trans_inv_operator(mat, 4, B) |> Hermitian |> eigvals

        # N + Z 
        Es = []
        for z in [1, -1]
            B = basis(;L, N, z)
            add_value!(Es, B)
        end
        @test sort(Es) ≈ E 

        # N + Z + P
        Es = []
        for z in [1, -1], p in [1, -1]
            B = basis(;L, N, p, z)
            add_value!(Es, B)
        end
        @test sort(Es) ≈ E 

        # N + Z + k
        Es = []
        for k in 0:L-1, z in [1, -1]
            B = basis(;L, N, k, z)
            add_value!(Es, B)
        end
        @test sort(Es) ≈ E 
    end
end

@testset "Symmetric" begin
    for L in 4:2:12, k in [0, L÷2]
        N = L÷2
        
        B = basis(;L, N, k); iszero(size(B, 1)) && continue  
        E = trans_inv_operator(mat, 4, B) |> Hermitian |> eigvals

        # N + k + P + Z
        Es = []
        for p in [1, -1], z in [1, -1]
            B = basis(;L, N, k, p, z)
            add_value!(Es, B)
        end
        @test sort(Es) ≈ E 
    end
end

function add_value!(M::AbstractMatrix, E, S, B)
    iszero(size(B, 1)) && return 
    e, v = trans_inv_operator(M, round(Int, log(2, size(M, 1))), B) |> Hermitian |> eigen
    append!(E, e)
    inds = 1:length(B)÷2 
    s = [ent_S(v[:, i], inds, B) for i in axes(v, 2)]
    append!(S, s)
end

function add_value!(Ms::Vector, E, S, B)
    iszero(size(B, 1)) && return 
    n = round(Int, log(2, size(Ms[1], 1)))
    L = length(B)
    inds = [mod.(i-1:i+n-2, L) .+ 1 for i in 1:L]
    e, v = operator(Ms, inds, B) |> Hermitian |> eigen
    append!(E, e)
    s = [ent_S(v[:, i], 1:L÷2, B) for i in axes(v, 2)]
    append!(S, s)
end

@testset "Relax k" begin
    for L in 4:8
        M = randn(ComplexF64, 8, 8) |> Hermitian |> Array  

        B = basis(;L); iszero(size(B, 1)) && continue  
        E, V = trans_inv_operator(M, 3, B) |> Hermitian |> eigen
        S = [ent_S(V[:, i], 1:L÷2, B) for i in axes(V, 2)] |> sort
        
        Es, Ss = [], []
        for k in 0:L-1
            B = basis(;L, k)
            add_value!(M, Es, Ss, B)
        end
        @test sort(Es) ≈ E
        @test norm(sort(Ss) - S) < 1e-7
    end
end



@testset "Relax N" begin
    a = 1.0 
    b = randn()
    c = randn()
    d = randn()
    e = randn()
    f = 5

    mat = spin(
        (a, "xx11"), (a, "yy11"), 
        (b, "x1x1"), (b, "y1y1"), 
        (c, "x11x"), (c, "y11y"), 
        (d, "zz11"),
        (e, "z1z1"), 
        (f, "z111")
    )
    for L in 4:8
        mats = [randn()*mat for i in 1:L]
        inds = [mod.(i-1:i+2, L) .+ 1 for i in 1:L]

        B = basis(;L); iszero(size(B, 1)) && continue  
        E, V = operator(mats, inds, B) |> Hermitian |> eigen
        S = [ent_S(V[:, i], 1:L÷2, B) for i in axes(V, 2)] |> sort
        
        Es, Ss = [], []
        for N in 0:L
            B = basis(;L, N)
            add_value!(mats, Es, Ss, B)
        end
        @test sort(Es) ≈ E
        @test norm(sort(Ss) - S) < 1e-7
        
    end
end

@testset "k -> N" begin
    a = 1.0 
    b = randn()
    c = randn()
    d = randn()
    e = randn()
    f = 5

    mat = spin(
        (a, "xx11"), (a, "yy11"), 
        (b, "x1x1"), (b, "y1y1"), 
        (c, "x11x"), (c, "y11y"), 
        (d, "zz11"),
        (e, "z1z1"), 
        (f, "z111")
    )

    for L in 4:8, k in 0:L-1
        B = basis(;L, k); iszero(size(B, 1)) && continue  
        E, V = trans_inv_operator(mat, 4, B) |> Hermitian |> eigen
        S = [ent_S(V[:, i], 1:L÷2, B) for i in axes(V, 2)] |> sort
        
        Es, Ss = [], []
        for N in 0:L
            B = basis(;L, N, k)
            add_value!(mat, Es, Ss, B)
        end
        @test sort(Es) ≈ E
        @test norm(sort(Ss) - S) < 1e-7
        
    end
end

@testset "N k -> P Z" begin
    for L in 4:2:12, k in [0, L÷2]
        N = L÷2

        # N + k + P + Z
        B = basis(;L, N, k); iszero(size(B, 1)) && continue  
        E, V = trans_inv_operator(mat, 4, B) |> Hermitian |> eigen
        S = [ent_S(V[:, i], 1:L÷2, B) for i in axes(V, 2)] |> sort
        
        Es, Ss = [], []
        for p in [1, -1], z in [1, -1]
            B = basis(;L, N, k, p, z)
            add_value!(mat, Es, Ss, B)
        end
        @test sort(Es) ≈ E
        @test sort(Ss) ≈ S
    end
end
