@testset "Lindblad Dynamics" begin
    @testset "Single-Site Amplitude Damping" begin
        γ = 0.4
        jump = sqrt(γ) * ComplexF64[0 1; 0 0]
        lb = lindblad(zeros(ComplexF64, 2, 2), [jump])
        dm = densitymatrix(ComplexF64[0.0, 1.0])

        dt = 0.02
        steps = 50
        projector_excited = ComplexF64[0 0; 0 1]
        projector_ground = ComplexF64[1 0; 0 0]

        for _ in 1:steps
            dm = lb(dm, dt, order = 8)
        end

        t = dt * steps
        p_excited = real(expectation(projector_excited, dm))
        p_ground = real(expectation(projector_ground, dm))

        @test tr(dm.ρ) ≈ 1.0 atol = 1e-10
        @test minimum(eigvals(Hermitian(dm.ρ))) ≥ -1e-10
        @test p_excited ≈ exp(-γ * t) atol = 5e-5
        @test p_ground ≈ 1 - exp(-γ * t) atol = 5e-5
    end

    @testset "Quadratic and Many-Body XX Loss Match" begin
        ql_empty = quadraticlindblad(zeros(4, 4), zeros(ComplexF64, 4, 2), Matrix{Float64}[])
        @test ql_empty isa EDKit.QuardraticLindblad
        @test ql_empty(covariancematrix([1, 0]), 0.1, order = 4) isa EDKit.CovarianceMatrix

        function quadratic_chain(L, r, a, b, c)
            H1 = zeros(2L, 2L)
            for i in 1:L-1
                H1[i, i + 1 + L] = 2r[i]
                H1[i + 1 + L, i] = -2r[i]
                H1[i + 1, i + L] = 2r[i]
                H1[i + L, i + 1] = -2r[i]
            end

            L1 = zeros(ComplexF64, 2L, L)
            for i in 1:L
                L1[i, i] = a[i]
                L1[i + L, i] = b[i]
            end

            M1 = Matrix{Float64}[]
            for i in 1:L
                mat = zeros(2L, 2L)
                mat[i, i + L] = c[i] / 2
                mat[i + L, i] = -c[i] / 2
                push!(M1, mat)
            end

            ql = quadraticlindblad(H1, L1, M1)
            occ = [mod(i, 2) for i in 1:L]
            cm = covariancematrix(occ)
            ql, cm
        end

        function manybody_chain(L, r, a, b, c)
            mats = [spin((r[i], "XX"), (r[i], "YY"), D = 2) for i in 1:L-1]
            inds = [[i, i + 1] for i in 1:L-1]
            H = Array(operator(mats, inds, L))

            jumps = Matrix{ComplexF64}[]
            for i in 1:L
                string_z = if i == 1
                    Matrix{Float64}(I, 1, 1)
                elseif i == 2
                    [-1.0 0.0; 0.0 1.0]
                else
                    reduce(kron, fill([-1.0 0.0; 0.0 1.0], i - 1))
                end
                local_jump = spin((a[i], "X"), (b[i], "Y"), D = 2)
                right_eye = Matrix{Float64}(I, 2^(L - i), 2^(L - i))
                push!(jumps, kron(string_z, local_jump, right_eye))
            end
            for i in 1:L
                push!(jumps, Array(operator(c[i] * spin("Z"), [i], L)))
            end

            lb = lindblad(H, jumps)
            pattern = [mod(i - 1, 2) for i in 1:L]
            dm = densitymatrix(index(pattern), L)
            lb, dm
        end

        function densities(dm, L)
            vals = zeros(L)
            for i in 1:L
                n̂ = Hermitian(Array(operator([1.0 0.0; 0.0 0.0], [i], L)))
                vals[i] = expectation(n̂, dm)
            end
            vals
        end

        L = 4
        r = [0.7, 0.4, 0.6]
        a = ComplexF64[0.2, 0.1, 0.25, 0.15]
        b = ComplexF64[0.0, 0.05im, 0.0, -0.03im]
        c = [0.1, 0.0, 0.08, 0.0]

        ql, cm = quadratic_chain(L, r, a, b, c)
        lb, dm = manybody_chain(L, r, a, b, c)

        dt = 0.01
        steps = 12
        for _ in 1:steps
            cm = ql(cm, dt, order = 8)
            dm = lb(dm, dt, order = 8)
        end

        n_quadratic = real.(diag(fermioncorrelation(cm, 1)))
        n_manybody = densities(dm, L)

        @test tr(dm.ρ) ≈ 1.0 atol = 1e-8
        @test norm(n_manybody - n_quadratic) ≤ 1e-5
    end
end
