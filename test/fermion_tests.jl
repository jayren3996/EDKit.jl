@testset "Spinless Fermion Basis" begin
    @testset "Basis construction" begin
        B = SpinlessFermionBasis(L = 4)
        @test size(B, 1) == 2^4
        @test size(B, 2) == 2^4
        @test B.B == 2
        @test length(B.I) == 2^4

        Bn = SpinlessFermionBasis(L = 6, N = 3)
        @test size(Bn, 1) == binomial(6, 3)

        # N=1 must mean exactly one occupied site (dgt=1), not one empty site.
        # binomial(6,1) is symmetric so size alone wouldn't catch a swapped
        # convention — check digit sums directly.
        Bone = SpinlessFermionBasis(L = 6, N = 1)
        @test size(Bone, 1) == 6
        for i in 1:size(Bone, 1)
            change!(Bone, i)
            @test sum(Bone.dgt) == 1
        end

        # copy gives an independent digit buffer
        B = SpinlessFermionBasis(L = 4)
        Bc = copy(B)
        @test Bc isa SpinlessFermionBasis
        @test Bc.I === B.I    # share the representative list
        @test Bc.dgt !== B.dgt # but own dgt

        # index round-trip
        change!(B, 5)
        coeff, pos = index(B)
        @test coeff == 1
        @test pos == 5
    end

    @testset "Number and nearest-neighbor hop" begin
        @test fermion("n") ≈ [0.0 0.0; 0.0 1.0]

        expected = kron(Array(EDKit.spin_Sm(2)), Array(EDKit.spin_Sp(2)))
        @test fermion("+-", 2) ≈ expected

        B = SpinlessFermionBasis(L = 4)
        Bt = TensorBasis(L = 4, base = 2)
        n_op_fermion = Array(fermion_operator("n", [2], B))
        n_op_spin = Array(operator(diagm([0.0, 1.0]), [2], Bt))
        @test n_op_fermion ≈ n_op_spin

        # The bare local matrices are accessible via fermion() (no JW prefix);
        # fermion_operator now embeds them with the proper JW string — see
        # the "JW-embedded single-site c† and c" testset below for that path.
        @test fermion("+") ≈ Array(EDKit.spin_Sm(2))
        @test fermion("-") ≈ Array(EDKit.spin_Sp(2))
    end

    @testset "Non-adjacent hop carries the JW string" begin
        # c†_1 c_3 on a 3-site block: σ⁻_1 σ_z^2 σ⁺_3
        σ⁺ = Array(EDKit.spin_Sp(2))
        σ⁻ = Array(EDKit.spin_Sm(2))
        σ_z = [1.0 0.0; 0.0 -1.0]
        expected = kron(kron(σ⁻, σ_z), σ⁺)
        @test fermion("+-", 3) ≈ expected

        # c†_1 c_4 on a 4-site block: σ⁻_1 σ_z^2 σ_z^3 σ⁺_4
        expected4 = kron(kron(kron(σ⁻, σ_z), σ_z), σ⁺)
        @test fermion("+-", 4) ≈ expected4

        # Compare a many-body matrix built via fermion_operator against the
        # explicit JW spin construction on a TensorBasis.
        L = 5
        Bf = SpinlessFermionBasis(L = L)
        Bt = TensorBasis(L = L, base = 2)
        # c†_1 c_4
        H_fermion = Array(fermion_operator("+-", [1, 4], Bf))
        H_spin = Array(operator(kron(σ⁻, σ_z, σ_z, σ⁺), [1, 2, 3, 4], Bt))
        @test H_fermion ≈ H_spin
        # c†_2 c_5
        H_fermion2 = Array(fermion_operator("+-", [2, 5], Bf))
        H_spin2 = Array(operator(kron(σ⁻, σ_z, σ_z, σ⁺), [2, 3, 4, 5], Bt))
        @test H_fermion2 ≈ H_spin2

        # Pair operators: "--" has a leading − from σ⁺σ_z = −σ⁺, while "++" has
        # no overall sign because σ⁻σ_z = +σ⁻.
        @test fermion("--", 2) ≈ -kron(σ⁺, σ⁺)
        @test fermion("++", 2) ≈ kron(σ⁻, σ⁻)
        @test fermion("--", 4) ≈ -kron(kron(kron(σ⁺, σ_z), σ_z), σ⁺)
        @test fermion("++", 4) ≈ kron(kron(kron(σ⁻, σ_z), σ_z), σ⁻)

        # Swapped sites: c†_3 c_1 = (c†_1 c_3)† = +σ⁺_1 σ_z^2 σ⁻_3.
        # Verify both against the hand-built JW spin operator.
        Bf3 = SpinlessFermionBasis(L = L)
        Bt3 = TensorBasis(L = L, base = 2)
        forward = Array(fermion_operator("+-", [1, 3], Bf3))
        forward_spin = Array(operator(kron(σ⁻, σ_z, σ⁺), [1, 2, 3], Bt3))
        @test forward ≈ forward_spin
        reverse = Array(fermion_operator("+-", [3, 1], Bf3))
        reverse_spin = Array(operator(kron(σ⁺, σ_z, σ⁻), [1, 2, 3], Bt3))
        @test reverse ≈ reverse_spin

        # Fermionic anticommutation: {c†_1, c†_3} = 0 (cross-check on the signs
        # produced by both the natural and the swapped branch of "++").
        c13 = Array(fermion_operator("++", [1, 3], Bf3))
        c31 = Array(fermion_operator("++", [3, 1], Bf3))
        @test c13 + c31 ≈ zeros(2^L, 2^L) atol = 1e-12

        # "-+" operator: c_1 c_span†. Leading − sign from σ⁺ σ_z = -σ⁺ at site 1.
        @test fermion("-+", 2) ≈ -kron(σ⁺, σ⁻)
        @test fermion("-+", 3) ≈ -kron(kron(σ⁺, σ_z), σ⁻)
        # Many-body matrix check across a non-trivial span.
        Lpm = 5
        Bpm = SpinlessFermionBasis(L = Lpm)
        Btpm = TensorBasis(L = Lpm, base = 2)
        @test Array(fermion_operator("-+", [1, 4], Bpm)) ≈
              Array(operator(-kron(kron(kron(σ⁺, σ_z), σ_z), σ⁻), [1, 2, 3, 4], Btpm))
        # Swap branch for "-+": (c_lo c†_hi)† = c†_hi c_lo... but for "-+" the swap
        # rule is the Hermitian conjugate (same as "+-"), so [3,1] gives the
        # conjugate of [1,3].
        @test Array(fermion_operator("-+", [3, 1], Bpm)) ≈
              Array(fermion_operator("-+", [1, 3], Bpm))'
    end

    @testset "Hermiticity and free-fermion spectrum" begin
        L = 6
        Bf = SpinlessFermionBasis(L = L, N = 1)  # single-particle sector

        # H = -Σ_i (c†_i c_{i+1} + h.c.) with PBC
        t = 1.0
        H = sum(
            -t * (
                fermion_operator("+-", [i, mod1(i + 1, L)], Bf) +
                fermion_operator("+-", [mod1(i + 1, L), i], Bf)
            )
            for i in 1:L
        )
        Hm = Array(H)
        @test Hm ≈ Hm'  # Hermitian

        # Single-particle spectrum: -2t cos(2π k / L) for k = 0, 1, …, L-1
        eigvals_h = eigvals(Hermitian(Hm))
        analytic = sort([-2 * t * cos(2π * k / L) for k in 0:L-1])
        @test eigvals_h ≈ analytic atol = 1e-10
    end

    @testset "Single-site z and I operators" begin
        @test fermion("z") ≈ [-0.5 0.0; 0.0 0.5]
        @test fermion("I") ≈ [1.0 0.0; 0.0 1.0]

        L = 5
        Bf = SpinlessFermionBasis(L = L)
        Bt = TensorBasis(L = L, base = 2)
        @test Array(fermion_operator("z", [3], Bf)) ≈
              Array(operator([-0.5 0.0; 0.0 0.5], [3], Bt))
        @test Array(fermion_operator("I", [3], Bf)) ≈
              Array(operator(Matrix{Float64}(I, 2, 2), [3], Bt))
    end

    @testset "JW-embedded single-site c† and c" begin
        σ⁺ = Array(EDKit.spin_Sp(2))
        σ⁻ = Array(EDKit.spin_Sm(2))
        σ_z = [1.0 0.0; 0.0 -1.0]
        L = 5
        Bf = SpinlessFermionBasis(L = L)
        Bt = TensorBasis(L = L, base = 2)

        # c†_1 has no JW prefix
        @test Array(fermion_operator("+", [1], Bf)) ≈ Array(operator(σ⁻, [1], Bt))
        # c†_3 = σ_z^1 σ_z^2 σ⁻_3
        @test Array(fermion_operator("+", [3], Bf)) ≈
              Array(operator(kron(σ_z, σ_z, σ⁻), [1, 2, 3], Bt))
        # c_4 = σ_z^1 σ_z^2 σ_z^3 σ⁺_4
        @test Array(fermion_operator("-", [4], Bf)) ≈
              Array(operator(kron(σ_z, σ_z, σ_z, σ⁺), [1, 2, 3, 4], Bt))

        # Canonical anticommutation {c_i, c†_j} = δ_ij × I
        for i in 1:L, j in 1:L
            a = Array(fermion_operator("+", [i], Bf))
            b = Array(fermion_operator("-", [j], Bf))
            expected = i == j ? Matrix{Float64}(I, 2^L, 2^L) : zeros(2^L, 2^L)
            @test a * b + b * a ≈ expected atol = 1e-12
        end
    end

    @testset "Density-density operator nn" begin
        n_mat = [0.0 0.0; 0.0 1.0]
        Id2 = Matrix{Float64}(I, 2, 2)
        @test fermion("nn", 2) ≈ kron(n_mat, n_mat)
        @test fermion("nn", 4) ≈ kron(kron(kron(n_mat, Id2), Id2), n_mat)

        L = 5
        Bf = SpinlessFermionBasis(L = L)
        Bt = TensorBasis(L = L, base = 2)
        @test Array(fermion_operator("nn", [2, 5], Bf)) ≈
              Array(operator(kron(kron(kron(n_mat, Id2), Id2), n_mat), [2, 3, 4, 5], Bt))
        # nn is symmetric in its arguments
        @test Array(fermion_operator("nn", [2, 5], Bf)) ≈
              Array(fermion_operator("nn", [5, 2], Bf))
        # nn equals product of two single-site n operators
        n2 = Array(fermion_operator("n", [2], Bf))
        n5 = Array(fermion_operator("n", [5], Bf))
        @test Array(fermion_operator("nn", [2, 5], Bf)) ≈ n2 * n5
    end

    @testset "Density-fraction nf and multi-sector N" begin
        # nf shorthand: nf = N / L
        @test SpinlessFermionBasis(L = 6, nf = 1 / 2).I == SpinlessFermionBasis(L = 6, N = 3).I
        @test SpinlessFermionBasis(L = 8, nf = 1 / 4).I == SpinlessFermionBasis(L = 8, N = 2).I

        # Multi-sector N as a vector: union of fixed-N sectors, sorted reps
        B12 = SpinlessFermionBasis(L = 6, N = [1, 2])
        @test size(B12, 1) == binomial(6, 1) + binomial(6, 2)
        @test issorted(B12.I)
        for i in 1:size(B12, 1)
            change!(B12, i)
            @test sum(B12.dgt) in (1, 2)
        end

        # Error on conflicting or invalid spec
        @test_throws ErrorException SpinlessFermionBasis(L = 6, N = 3, nf = 0.5)
        @test_throws ErrorException SpinlessFermionBasis(L = 6, nf = 0.3)
    end

    @testset "Tight-binding ring PBC at N ≥ 2 (long-way JW on wrap)" begin
        # Free-fermion eigenvalues at fixed N on the L-site ring: pick any N
        # single-particle modes from ε_k = -2t cos(2π k / L) and sum them.
        # The N = 1 sector is already covered above; this set targets N ≥ 2,
        # where a missing JW string on the L→1 bond gives wrong eigenvalues.
        function tight_binding_spectrum(L::Integer, N::Integer; t::Real=1.0)
            modes = [-2 * t * cos(2π * k / L) for k in 0:L-1]
            spectrum = Float64[]
            for bits in 0:(2^L - 1)
                count_ones(bits) == N || continue
                E = sum(modes[i+1] for i in 0:L-1 if (bits >> i) & 1 == 1; init = 0.0)
                push!(spectrum, E)
            end
            sort!(spectrum)
        end

        for (L, N) in [(4, 2), (6, 3)]
            Bf = SpinlessFermionBasis(L = L, N = N)

            # Explicit per-bond construction (mod1 wraps include the L→1 bond)
            H_explicit = sum(
                fermion_operator("+-", [i, mod1(i + 1, L)], Bf) +
                fermion_operator("+-", [mod1(i + 1, L), i], Bf)
                for i in 1:L
            )
            Hm_explicit = -Array(H_explicit)
            @test Hm_explicit ≈ Hm_explicit'

            # Same Hamiltonian via trans_inv_fermion_operator
            H_hop = trans_inv_fermion_operator("+-", [1, 2], Bf)
            Hm_trans = -Array(H_hop + adjoint(H_hop))
            @test Hm_trans ≈ Hm_trans'

            @test Hm_trans ≈ Hm_explicit atol = 1e-10
            @test sort(eigvals(Hermitian(Hm_explicit))) ≈ tight_binding_spectrum(L, N) atol = 1e-10
        end

        # Range form sugar: span Integer agrees with explicit support vector.
        Bf = SpinlessFermionBasis(L = 5, N = 2)
        @test Array(trans_inv_fermion_operator("+-", 2, Bf)) ≈
              Array(trans_inv_fermion_operator("+-", [1, 2], Bf))
    end

    @testset "fermion_operator errors on symmetry-reduced bases" begin
        L = 6
        Bt = TranslationalBasis(L = L, k = 0, base = 2)
        Bp = ParityBasis(L = L, p = 1, base = 2)

        # JW-bearing operators must error — silently embedding them on a
        # permutation-symmetric basis would give wrong eigenvalues.
        for op in ("+", "-")
            @test_throws ErrorException fermion_operator(op, [1], Bt)
            @test_throws ErrorException fermion_operator(op, [1], Bp)
        end
        for op in ("+-", "-+", "++", "--")
            @test_throws ErrorException fermion_operator(op, [1, 2], Bt)
            @test_throws ErrorException fermion_operator(op, [1, 2], Bp)
        end

        # Diagonal / identity operators carry no JW string and commute with
        # any permutation, so they remain valid.
        @test fermion_operator("n", [1], Bt) isa EDKit.Operator
        @test fermion_operator("z", [2], Bt) isa EDKit.Operator
        @test fermion_operator("I", [3], Bp) isa EDKit.Operator
        @test fermion_operator("nn", [1, 3], Bt) isa EDKit.Operator

        # And the trans_inv helper itself does not accept those bases (the
        # method only dispatches on SpinlessFermionBasis).
        @test_throws MethodError trans_inv_fermion_operator("+-", [1, 2], Bt)
    end
end
