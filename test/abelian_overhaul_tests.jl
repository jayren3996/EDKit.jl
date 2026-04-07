using EDKit
using LinearAlgebra
using Test

# Access internal functions for testing
import EDKit: compile_benes, apply_benes, apply_perm_int, BenesNetwork,
              AbelianOperator, check_min, shift_canonical!, check_min_int,
              shift_canonical_int, _gosper_next, _gosper_enumerate, _has_all_benes,
              apply_perm!, apply_inv!, init!, phase, order

@testset "Abelian Overhaul" begin

    @testset "Benes Network" begin
        @testset "Identity permutation" begin
            for L in [4, 8, 12, 16]
                perm = collect(1:L)
                bn = compile_benes(perm, L)
                for s in UInt64(0):UInt64(min(2^L - 1, 255))
                    @test apply_benes(bn, s, L) == s
                end
            end
        end

        @testset "Reverse permutation" begin
            for L in [4, 6, 8]
                perm = collect(L:-1:1)
                bn = compile_benes(perm, L)
                for s in UInt64(0):UInt64(2^L - 1)
                    expected = apply_perm_int(s, perm, L)
                    @test apply_benes(bn, s, L) == expected
                end
            end
        end

        @testset "Cyclic shift" begin
            for L in [4, 6, 8]
                # Shift by 1: site i -> site mod1(i-1, L)
                perm = Vector{Int}(undef, L)
                for i in 1:L
                    perm[i] = mod1(i - 1, L)
                end
                bn = compile_benes(perm, L)
                for s in UInt64(0):UInt64(2^L - 1)
                    expected = apply_perm_int(s, perm, L)
                    @test apply_benes(bn, s, L) == expected
                end
            end
        end

        @testset "Random permutations L=$L" for L in [4, 8, 12, 16]
            using Random
            rng = MersenneTwister(42 + L)
            for _ in 1:5
                perm = randperm(rng, L)
                bn = compile_benes(perm, L)
                # Test a subset of states
                ntest = min(2^L, 256)
                states = L <= 8 ? collect(UInt64(0):UInt64(2^L - 1)) :
                         UInt64.(rand(rng, 0:2^L-1, ntest))
                for s in states
                    expected = apply_perm_int(s, perm, L)
                    result = apply_benes(bn, s, L)
                    @test result == expected
                end
            end
        end
    end

    @testset "AbelianOperator Constructor" begin
        @testset "Translation generator" begin
            L = 6
            perm = Vector{Int}(undef, L)
            for i in 1:L
                perm[i] = mod1(i - 1, L)
            end
            ag = AbelianOperator(L, 0, perm)
            @test order(ag) == L
            @test length(ag.perm) == 1
            @test ag.perm[1] == perm
            @test !any(ag.inv[1])
        end

        @testset "Parity generator" begin
            L = 6
            perm = collect(L:-1:1)
            ag = AbelianOperator(2, 1, perm)
            @test order(ag) == 2
        end

        @testset "Spin-flip generator" begin
            L = 6
            perm = collect(1:L)
            inv_mask = trues(L)
            ag = AbelianOperator(2, 0, perm; inv=inv_mask)
            @test order(ag) == 2
            @test all(ag.inv[1])
        end

        @testset "Combination (direct product)" begin
            L = 6
            perm_t = Vector{Int}(undef, L)
            for i in 1:L; perm_t[i] = mod1(i - 1, L); end
            perm_p = collect(L:-1:1)
            ag1 = AbelianOperator(L, 0, perm_t)
            ag2 = AbelianOperator(2, 0, perm_p)
            ag = ag1 + ag2
            @test order(ag) == 2L
            @test length(ag.perm) == 2
        end

        @testset "Apply permutation to digit buffer" begin
            L = 4
            perm = [2, 3, 4, 1]  # perm[i]=j: value at site i goes to site j
            dgt = [1, 0, 0, 1]
            apply_perm!(dgt, perm)
            # dgt_new[2]=dgt_old[1]=1, dgt_new[3]=dgt_old[2]=0,
            # dgt_new[4]=dgt_old[3]=0, dgt_new[1]=dgt_old[4]=1
            @test dgt == [1, 1, 0, 0]
        end

        @testset "Odometer advances correctly" begin
            L = 4
            perm_t = Vector{Int}(undef, L)
            for i in 1:L; perm_t[i] = mod1(i - 1, L); end
            ag = AbelianOperator(L, 1, perm_t)
            init!(ag)
            dgt = [1, 0, 1, 0]
            original = copy(dgt)
            # Apply L times should return to original
            for _ in 1:L
                ag(dgt, 2)
            end
            @test dgt == original
        end
    end

    @testset "Integer Orbit Search" begin
        @testset "check_min_int matches check_min" begin
            L = 6
            perm_t = Vector{Int}(undef, L)
            for i in 1:L; perm_t[i] = mod1(i - 1, L); end
            ag = AbelianOperator(L, 0, perm_t)

            dgt = zeros(Int, L)
            for idx in 1:2^L
                change!(dgt, idx; base=Int64(2))
                ag_copy1 = deepcopy(ag)
                ag_copy2 = deepcopy(ag)
                q1, n1 = check_min(dgt, ag_copy1; base=2)
                q2, n2 = check_min_int(UInt64(idx - 1), ag_copy2, L)
                @test q1 == q2
                @test n1 == n2
            end
        end

        @testset "shift_canonical_int matches shift_canonical!" begin
            L = 6
            perm_t = Vector{Int}(undef, L)
            for i in 1:L; perm_t[i] = mod1(i - 1, L); end
            ag = AbelianOperator(L, 1, perm_t)

            dgt = zeros(Int, L)
            for idx in 1:2^L
                change!(dgt, idx; base=Int64(2))
                dgt_copy = copy(dgt)
                ag_copy1 = deepcopy(ag)
                ag_copy2 = deepcopy(ag)
                Im1, g1 = shift_canonical!(dgt_copy, ag_copy1; base=2)
                Im2, g2 = shift_canonical_int(UInt64(idx - 1), ag_copy2, L)
                @test Im1 == Im2
            end
        end

        @testset "Combined symmetry integer path" begin
            L = 6
            perm_t = Vector{Int}(undef, L)
            for i in 1:L; perm_t[i] = mod1(i - 1, L); end
            perm_p = collect(L:-1:1)
            ag = AbelianOperator(L, 0, perm_t) + AbelianOperator(2, 0, perm_p)

            dgt = zeros(Int, L)
            for idx in 1:2^L
                change!(dgt, idx; base=Int64(2))
                ag1 = deepcopy(ag)
                ag2 = deepcopy(ag)
                q1, n1 = check_min(dgt, ag1; base=2)
                q2, n2 = check_min_int(UInt64(idx - 1), ag2, L)
                @test q1 == q2
                @test n1 == n2
            end
        end
    end

    @testset "Gosper's Hack" begin
        @testset "_gosper_next" begin
            # Starting from 0b0111 (7), next with 3 bits in 4 positions
            s = UInt64(0b0111)
            s = _gosper_next(s)
            @test s == UInt64(0b1011)
            s = _gosper_next(s)
            @test s == UInt64(0b1101)
            s = _gosper_next(s)
            @test s == UInt64(0b1110)
        end

        @testset "_gosper_enumerate" begin
            # L=4, N=2: should give C(4,2) = 6 states
            result = _gosper_enumerate(4, 2)
            @test length(result) == 6
            # Check they all have exactly 2 bits set
            for idx in result
                @test count_ones(idx - 1) == 2
            end
            # Sorted
            @test issorted(result)

            # Edge cases
            @test _gosper_enumerate(4, 0) == [1]
            @test _gosper_enumerate(4, 4) == [16]

            # L=8, N=4: C(8,4) = 70
            result8 = _gosper_enumerate(8, 4)
            @test length(result8) == 70
        end
    end

    @testset "AbelianBasis Construction and Eigenvalues" begin
        L = 8
        XXZ = spin((1.0, "xx"), (1.0, "yy"), (0.4, "zz"))
        full_vals = eigvals(Hermitian(trans_inv_operator(XXZ, 2, L))) |> sort

        @testset "k-sectors reconstruct full spectrum" begin
            vals_k = Float64[]
            for k in 0:L-1
                B = basis(L=L, k=k)
                append!(vals_k, eigvals(Hermitian(trans_inv_operator(XXZ, 2, B))))
            end
            @test sort(vals_k) ≈ full_vals
        end

        @testset "k+p sectors reconstruct k=0 spectrum" begin
            Bk0 = basis(L=L, k=0)
            vals_k0 = eigvals(Hermitian(trans_inv_operator(XXZ, 2, Bk0))) |> sort

            vals_kp = Float64[]
            for p in (1, -1)
                B = basis(L=L, k=0, p=p)
                append!(vals_kp, eigvals(Hermitian(trans_inv_operator(XXZ, 2, B))))
            end
            @test sort(vals_kp) ≈ vals_k0
        end

        @testset "N+k sectors reconstruct N sector" begin
            Nhalf = L ÷ 2
            BN = basis(L=L, N=Nhalf)
            vals_N = eigvals(Hermitian(trans_inv_operator(XXZ, 2, BN))) |> sort

            vals_Nk = Float64[]
            for k in 0:L-1
                B = basis(L=L, N=Nhalf, k=k)
                sz = size(B, 1)
                sz == 0 && continue
                append!(vals_Nk, eigvals(Hermitian(trans_inv_operator(XXZ, 2, B))))
            end
            @test sort(vals_Nk) ≈ vals_N
        end

        @testset "N+k+z sectors at half-filling" begin
            Nhalf = L ÷ 2
            BNk = basis(L=L, N=Nhalf, k=1)
            vals_Nk = eigvals(Hermitian(trans_inv_operator(XXZ, 2, BNk))) |> sort

            vals_Nkz = Float64[]
            for z in (1, -1)
                B = basis(L=L, N=Nhalf, k=1, z=z)
                sz = size(B, 1)
                sz == 0 && continue
                append!(vals_Nkz, eigvals(Hermitian(trans_inv_operator(XXZ, 2, B))))
            end
            @test sort(vals_Nkz) ≈ vals_Nk
        end
    end

    @testset "Full Symmetric Basis (N+k+p+z)" begin
        L = 8
        Nhalf = L ÷ 2
        XXZ = spin((1.0, "xx"), (1.0, "yy"), (0.4, "zz"))

        # k=0 sector with N, p, z
        BNk0 = basis(L=L, N=Nhalf, k=0)
        vals_Nk0 = eigvals(Hermitian(trans_inv_operator(XXZ, 2, BNk0))) |> sort

        vals_full_sym = Float64[]
        for p in (1, -1), z in (1, -1)
            B = basis(L=L, N=Nhalf, k=0, p=p, z=z)
            sz = size(B, 1)
            sz == 0 && continue
            append!(vals_full_sym, eigvals(Hermitian(trans_inv_operator(XXZ, 2, B))))
        end
        @test sort(vals_full_sym) ≈ vals_Nk0
    end

    @testset "Symmetries keyword" begin
        L = 6
        # Translation as a custom symmetry should give same result as k keyword
        perm_t = Vector{Int}(undef, L)
        for i in 1:L; perm_t[i] = mod1(i - 1, L); end

        B1 = basis(L=L, k=0)
        B2 = basis(L=L, symmetries=[(L, 0, perm_t)])
        @test size(B1) == size(B2)
        @test B1.I == B2.I
    end

end
