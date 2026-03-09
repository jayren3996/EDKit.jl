@testset "Core Operators And Helpers" begin
    @test Array(spin("X")) ≈ [0 1; 1 0]
    @test Array(spin("Y")) ≈ [0 -1im; 1im 0]
    @test Array(spin("xx")) ≈ kron(Array(spin("x")), Array(spin("x")))
    @test Array(spin((1.0, "xx"), (1.0, "yy"), (0.5, "zz"))) ≈ Array(spin("xx") + spin("yy") + 0.5 * spin("zz"))

    L = 4
    mat = spin((1.0, "xx"), (0.3, "zz"))
    inds = [2, 3]
    H = operator(mat, inds, L)
    @test Array(H) ≈ kron_embed(mat, inds, L)
    @test sparse(H) ≈ sparse(Array(H))

    bond = spin((1.0, "xx"), (1.0, "yy"))
    H1 = trans_inv_operator(bond, 2, L)
    H2 = operator(fill(bond, L), [[i, mod(i, L) + 1] for i in 1:L], L)
    @test Array(H1) ≈ Array(H2)

    v = randn(ComplexF64, 2^L)
    @test H1 * v ≈ Array(H1) * v
    @test EDKit.mul(H1, v) ≈ Array(H1) * v

    B = TensorBasis(L = L, base = 2)
    state = productstate([0, 1, 0, 1], B)
    @test count(!iszero, state) == 1
    @test state[index([0, 1, 0, 1], base = 2)] == 1

    bell = normalize([1.0, 0.0, 0.0, 1.0])
    @test ent_S(bell, [1], TensorBasis(L = 2, base = 2)) ≈ log(2)

    E = [0.0, 1.0, 3.0, 6.0]
    @test gapratio(E) ≈ [0.5, 2 / 3]
    @test meangapratio(E) ≈ 7 / 12
end
