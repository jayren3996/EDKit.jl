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

        # Reject bare "+" / "-" in fermion_operator (would silently miss the JW string).
        @test_throws ErrorException fermion_operator("+", [2], B)
        @test_throws ErrorException fermion_operator("-", [3], B)
        # The bare matrices are still accessible via fermion() itself.
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
    end
end
