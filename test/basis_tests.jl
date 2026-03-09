@testset "Basis Constructors And Symmetry Decomposition" begin
    L = 6
    Nhalf = L ÷ 2

    pxp(v) = all(v[i] == 0 || v[i + 1] == 0 for i in 1:length(v)-1)
    BP = ProjectedBasis(L = L, f = pxp)
    @test size(BP, 1) == 21

    BN = ProjectedBasis(L = L, N = Nhalf)
    @test size(BN, 1) == binomial(L, Nhalf)

    XXZ = spin((1.0, "xx"), (1.0, "yy"), (0.4, "zz"))
    full_vals = eigvals(Hermitian(trans_inv_operator(XXZ, 2, L))) |> sort

    vals_N = Float64[]
    for N in 0:L
        B = basis(L = L, N = N)
        append!(vals_N, eigvals(Hermitian(trans_inv_operator(XXZ, 2, B))))
    end
    @test sort(vals_N) ≈ full_vals

    vals_k = Float64[]
    for k in 0:L-1
        B = basis(L = L, k = k)
        append!(vals_k, eigvals(Hermitian(trans_inv_operator(XXZ, 2, B))))
    end
    @test sort(vals_k) ≈ full_vals

    Bk = basis(L = L, N = Nhalf, k = 0)
    vals_kp = Float64[]
    for p in (1, -1)
        B = basis(L = L, N = Nhalf, k = 0, p = p)
        append!(vals_kp, eigvals(Hermitian(trans_inv_operator(XXZ, 2, B))))
    end
    @test sort(vals_kp) ≈ eigvals(Hermitian(trans_inv_operator(XXZ, 2, Bk))) |> sort

    vals_kz = Float64[]
    for z in (1, -1)
        B = basis(L = L, N = Nhalf, k = 1, z = z)
        append!(vals_kz, eigvals(Hermitian(trans_inv_operator(XXZ, 2, B))))
    end
    @test sort(vals_kz) ≈ eigvals(Hermitian(trans_inv_operator(XXZ, 2, basis(L = L, N = Nhalf, k = 1)))) |> sort

    vals_a2 = Float64[]
    for k in 0:(L ÷ 2 - 1)
        B = TranslationalBasis(L = L, k = k, a = 2)
        append!(vals_a2, eigvals(Hermitian(trans_inv_operator(XXZ, 2, B))))
    end
    @test sort(vals_a2) ≈ full_vals
end
