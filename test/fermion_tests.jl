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
end
