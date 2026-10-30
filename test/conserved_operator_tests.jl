@testset "Conserved and slow operators" begin
    B = TensorBasis(L=4, base=2)
    H = operator(spin("z"), [1], B)
    candidates = [operator(spin(s), [1], B) for s in ("z", "x", "y")]
    result = conserved_operators(H, candidates)
    @test result.values ≈ [0, 1, 1] atol=1e-12
    @test result.vectors' * result.overlap * result.vectors ≈ I
    @test norm(Matrix(H) * Matrix(candidates[1]) - Matrix(candidates[1]) * Matrix(H)) < 1e-12

    # Linear dependence and scaling must not introduce fake conserved charges.
    redundant = conserved_operators(H, [candidates[1], 2 * candidates[1], candidates[2]])
    @test length(redundant.values) == 2
    @test redundant.values ≈ [0, 1] atol=1e-12
    @test conserved_operators(H, [100 * candidates[1], candidates[2]]).values ≈ [0, 1] atol=1e-12

    # Complex coefficients can combine Hermitian candidates into a non-Hermitian eigenoperator.
    plus = [0.0 1; 0 0]
    minus = plus'
    complex_result = conserved_operators([1.0 0; 0 -1], [plus, minus]; coefficients=:complex)
    @test complex_result.values ≈ [4, 4]
    @test_throws ArgumentError conserved_operators([1.0 0; 0 -1], [plus])
    @test_throws ArgumentError conserved_operators(H, [0 * candidates[1]])
    @test_throws DimensionMismatch conserved_operators(H, [zeros(2, 2)])

    # Projected operators are multiplied in the constrained space itself.
    constraint(dgt) = all(i -> !(dgt[i] == 1 && dgt[mod1(i + 1, length(dgt))] == 1), eachindex(dgt))
    Bp = ProjectedBasis(L=5, f=constraint, threaded=false)
    Hp = trans_inv_operator(spin("x"), 1, Bp)
    cp = [trans_inv_operator(spin(s), length(s), Bp) for s in ("x", "z", "zz")]
    projected = conserved_operators(Hp, cp)
    @test minimum(abs, projected.values) < 1e-11
    @test projected.commutator[1, 1] ≈ 0 atol=1e-12

    # XXZ has a conserved range-three energy current, independent of H.
    L = 6
    Bxxz = TensorBasis(L=L, base=2)
    h = Matrix(spin((1.0, "xx"), (1.0, "yy"), (0.7, "zz")))
    Hxxz = trans_inv_operator(h, 2, Bxxz)
    id = Matrix{ComplexF64}(I, 2, 2)
    h12, h23 = kron(h, id), kron(id, h)
    j = im * (h12 * h23 - h23 * h12)
    current = trans_inv_operator(j, 3, Bxxz)
    ops = [current, trans_inv_operator(spin("xx"), 2, Bxxz),
           trans_inv_operator(spin("zz"), 2, Bxxz), trans_inv_operator(spin("x"), 1, Bxxz)]
    xxz = conserved_operators(Hxxz, ops)
    @test norm(Matrix(current)) > 1
    @test norm(Matrix(Hxxz) * Matrix(current) - Matrix(current) * Matrix(Hxxz)) < 1e-12
    @test xxz.values ≈ [0, 0.045, 0.5, 1.245] atol=1e-11
    @test abs(dot(xxz.vectors[:, 1], xxz.overlap[:, 1])) / sqrt(real(xxz.overlap[1, 1])) ≈ 1

    # Recover the current from individual strings, without supplying the current.
    strings = ("xyz", "xzy", "yxz", "yzx", "zxy", "zyx")
    blind_ops = [trans_inv_operator(spin(s), 3, Bxxz) for s in strings]
    blind = conserved_operators(Hxxz, blind_ops)
    null_vectors = blind.vectors[:, findall(v -> abs(v) < 1e-11, blind.values)]
    overlaps = [real(dot(Matrix(o), Matrix(current))) / size(Hxxz, 1) for o in blind_ops]
    coeffs = null_vectors * (null_vectors' * overlaps)
    recovered = sum(coeffs[i] * Matrix(blind_ops[i]) for i in eachindex(blind_ops))
    @test norm(recovered - Matrix(current)) / norm(Matrix(current)) < 1e-10
end

@testset "Shared operator-search optimization" begin
    X = [0.0 1; 1 0]
    Y = ComplexF64[0 -im; im 0]
    Z = [1.0 0; 0 -1]
    Id = Matrix{Float64}(I, 2, 2)

    # Rank is determined after scaling individual physical operators.
    scaled = conserved_operators(Z, [1e6 * Z, X])
    @test scaled.values ≈ [0, 4] atol=1e-12
    @test scaled.rank == 2
    @test scaled.residuals .^ 2 ≈ scaled.values atol=1e-12
    @test scaled.vectors' * scaled.overlap * scaled.vectors ≈ I

    # Include zero and dependent directions, including more candidates than
    # the dimension of the complex operator space.
    redundant = conserved_operators(Z, [Z, 2Z, X, Y, Id, zeros(2, 2)]; coefficients=:complex)
    @test redundant.rank == 4
    @test size(redundant.discarded) == (6, 2)
    @test norm(hcat(vec(Z), vec(2Z), vec(X), vec(Y), vec(Id), zeros(4)) * redundant.discarded) < 1e-12

    # Exclusions use the projection of even an out-of-span operator.
    excluded = conserved_operators(Z, [Z, X]; exclude=[Z + Y])
    @test excluded.excluded_rank == 1
    @test excluded.values ≈ [4] atol=1e-12
    @test abs(dot(vec(Z + Y), vec(excluded.vectors[1, 1] * Z + excluded.vectors[2, 1] * X))) < 1e-12
    @test conserved_operators(Z, [Z, X]; exclude=[Y]).excluded_rank == 0
    @test_throws ArgumentError conserved_operators(Z, [Z]; exclude=[Z])
    @test_throws ArgumentError conserved_operators(Z, [Z]; metric=:coefficients, exclude=[Z])

    # Different bases of equal dimension cannot safely be mixed.
    Ba = ProjectedBasis(L=2, f=d -> d[1] == 0, threaded=false)
    Bb = ProjectedBasis(L=2, f=d -> d[1] == 1, threaded=false)
    @test size(Ba) == size(Bb)
    @test_throws ArgumentError conserved_operators(operator(spin("z"), [1], Ba), [operator(spin("z"), [1], Bb)])
    Ba2 = ProjectedBasis(L=2, f=d -> d[1] == 0, threaded=false)
    @test conserved_operators(operator(spin("z"), [1], Ba), [operator(spin("z"), [1], Ba2)]).rank == 1

    # Each state is an eigenstate of Z, but the mixture has nonzero variance.
    states = [2.0 0; 0 3.0]
    each = eigenstate_operators([Z], states)
    mixture = eigenstate_operators([Z], states; objective=:mixture)
    @test each.values ≈ [0] atol=1e-12
    @test mixture.values ≈ [1] atol=1e-12
    @test Matrix(covmat([Z], Matrix{Float64}(I, 2, 2))) ≈ Matrix(mixture.covariance)
    @test eigenstate_operators([Z], states; objective=:mixture, weights=[1, 3]).values ≈ [0.75]
    @test_throws ArgumentError eigenstate_operators([Z], zeros(2))
    @test_throws ArgumentError eigenstate_operators([Z], states; weights=[0, 0])

    # A metric built from physical operators removes dependent directions.
    qim = eigenstate_operators([Id, Z, 2Z, X], states[:, 1]; metric=:operator, exclude=[Id])
    @test qim.rank == 3
    @test qim.excluded_rank == 1
    @test qim.values ≈ [0, 1] atol=1e-12
    @test qim.vectors' * qim.overlap * qim.vectors ≈ I
    @test qim.residuals .^ 2 ≈ qim.values atol=1e-12

    # A wide residual map has an exact null space beyond the thin SVD.
    wide = eigenstate_operators([Id, X, Y, Z], [1.0, 0]; coefficients=:complex)
    @test length(wide.values) == 4
    @test count(v -> v < 1e-12, wide.values) == 3
    @test all(wide.residuals .^ 2 .≈ wide.values)
end

@testset "Eigenstate search with EDKit operators" begin
    B = TensorBasis(L=3, base=2)
    z = operator(spin("z"), [1], B)
    x = operator(spin("x"), [1], B)
    psi = zeros(ComplexF64, size(B, 1))
    psi[1] = 2
    fit = eigenstate_operators([z, x], psi)
    @test fit.values ≈ [0, 0.25] atol=1e-12
    @test fit.residuals .^ 2 ≈ fit.values atol=1e-12
end
