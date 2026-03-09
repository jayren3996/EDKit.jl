@testset "Entanglement Across Symmetry Sectors" begin
    function entropy_spectrum(B::AbstractBasis, mat::AbstractMatrix; cut = 1:length(B) ÷ 2)
        vals, vecs = eigen(Hermitian(trans_inv_operator(mat, 2, B)))
        sort([ent_S(vecs[:, i], collect(cut), B) for i in axes(vecs, 2)])
    end

    @testset "Momentum Resolves Tensor-Basis Entropies" begin
        L = 6
        mat = randn(ComplexF64, 4, 4) |> Hermitian |> Array

        B0 = TensorBasis(L = L, base = 2)
        S0 = entropy_spectrum(B0, mat)

        Ss = Float64[]
        for k in 0:L-1
            B = TranslationalBasis(L = L, k = k)
            append!(Ss, entropy_spectrum(B, mat))
        end
        @test sort(Ss) ≈ S0
    end

    @testset "Momentum Resolves Tensor-Basis Entropies With a=2" begin
        L = 6
        mat = randn(ComplexF64, 4, 4) |> Hermitian |> Array

        B0 = TensorBasis(L = L, base = 2)
        S0 = entropy_spectrum(B0, mat)

        Ss = Float64[]
        for k in 0:(L ÷ 2 - 1)
            B = TranslationalBasis(L = L, k = k, a = 2)
            append!(Ss, entropy_spectrum(B, mat))
        end
        @test sort(Ss) ≈ S0
    end

    @testset "Parity Splits k=0 And k=pi Sectors" begin
        L = 6
        h = 0.37
        mat = spin((1.0, "xx"), (0.8, "yy"), (0.4, "zz"))

        function solve(B)
            H = trans_inv_operator(mat, 2, B) + trans_inv_operator(spin((h, "z")), 1, B)
            eigen(Hermitian(H)).vectors
        end

        for k in (0, L ÷ 2)
            Bt = TranslationalBasis(L = L, k = k)
            Bp = TranslationParityBasis(L = L, k = k, p = 1)
            Bm = TranslationParityBasis(L = L, k = k, p = -1)
            vt, vp, vm = solve(Bt), solve(Bp), solve(Bm)
            for inds in (collect(1:L÷2), [1, 3], [2, 4, 6])
                St = sort([ent_S(vt[:, i], inds, Bt) for i in axes(vt, 2)])
                Sp = [ent_S(vp[:, i], inds, Bp) for i in axes(vp, 2)]
                Sm = [ent_S(vm[:, i], inds, Bm) for i in axes(vm, 2)]
                @test St ≈ sort(vcat(Sp, Sm))
            end
        end
    end

    @testset "Spin-Flip Splits Momentum Sectors" begin
        L = 6
        mat = spin((1.0, "xx"), (0.8, "yy"), (0.3, "zz"))

        for k in 0:L-1
            Bt = TranslationalBasis(L = L, k = k)
            Bp = TranslationFlipBasis(L = L, k = k, p = 1)
            Bm = TranslationFlipBasis(L = L, k = k, p = -1)
            vt = eigen(Hermitian(trans_inv_operator(mat, 2, Bt))).vectors
            vp = eigen(Hermitian(trans_inv_operator(mat, 2, Bp))).vectors
            vm = eigen(Hermitian(trans_inv_operator(mat, 2, Bm))).vectors
            for inds in (collect(1:L÷2), [1, 2], [2, 4, 6])
                St = sort([ent_S(vt[:, i], inds, Bt) for i in axes(vt, 2)])
                Sp = [ent_S(vp[:, i], inds, Bp) for i in axes(vp, 2)]
                Sm = [ent_S(vm[:, i], inds, Bm) for i in axes(vm, 2)]
                @test St ≈ sort(vcat(Sp, Sm))
            end
        end
    end

    @testset "Abelian basis matches dedicated symmetry constructors" begin
        L = 6
        N = L ÷ 2
        mat = spin((1.0, "xx"), (0.7, "yy"), (0.2, "zz"))

        cases = [
            (basis(L = L, N = N), ProjectedBasis(L = L, N = N)),
            (basis(L = L, k = 1), TranslationalBasis(L = L, k = 1)),
            (basis(L = L, p = 1), ParityBasis(L = L, p = 1)),
            (basis(L = L, z = 1), FlipBasis(L = L, p = 1)),
            (basis(L = L, N = N, k = 1), TranslationalBasis(L = L, N = N, k = 1)),
            (basis(L = L, N = N, p = 1), ParityBasis(L = L, N = N, p = 1)),
            (basis(L = L, N = N, z = 1), FlipBasis(L = L, N = N, p = 1)),
            (basis(L = L, k = 0, p = 1), TranslationParityBasis(L = L, k = 0, p = 1)),
            (basis(L = L, k = 1, z = 1), TranslationFlipBasis(L = L, k = 1, p = 1)),
            (basis(L = L, p = 1, z = 1), ParityFlipBasis(L = L, p = 1, z = 1)),
            (basis(L = L, N = N, k = 0, p = 1), basis(L = L, N = N, k = 0, p = 1)),
            (basis(L = L, N = N, k = 1, z = 1), basis(L = L, N = N, k = 1, z = 1)),
            (basis(L = L, N = N, p = 1, z = 1), basis(L = L, N = N, p = 1, z = 1)),
            (basis(L = L, k = 0, p = 1, z = 1), basis(L = L, k = 0, p = 1, z = 1)),
            (basis(L = L, N = N, k = 0, p = 1, z = 1), basis(L = L, N = N, k = 0, p = 1, z = 1)),
        ]

        for (Ba, Bb) in cases
            iszero(size(Ba, 1)) && continue
            va = eigen(Hermitian(trans_inv_operator(mat, 2, Ba))).vectors
            vb = eigen(Hermitian(trans_inv_operator(mat, 2, Bb))).vectors
            Sa = sort([ent_S(va[:, i], collect(1:L÷2), Ba) for i in axes(va, 2)])
            Sb = sort([ent_S(vb[:, i], collect(1:L÷2), Bb) for i in axes(vb, 2)])
            @test Sa ≈ Sb
        end
    end
end
