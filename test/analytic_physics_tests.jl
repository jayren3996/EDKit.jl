@testset "Analytic Physical Spectra" begin

    # ---------------------------------------------------------------------------
    # PBC Transverse-Field Ising Model: H = -J Σ σ^x_i σ^x_{i+1} - h Σ σ^z_i
    #
    # Diagonalised exactly via Jordan-Wigner + Bogoliubov (Pfeuty, Annals of
    # Physics 57, 79 (1970)). The L-site spectrum splits by Z₂ spin parity
    # P = Π σ^z:
    #   - P = +1 (anti-periodic / Neveu-Schwarz) momenta k = π(2m+1)/L
    #   - P = -1 (periodic     / Ramond)        momenta k = 2π m / L
    # In each sector the quasiparticle energies are ε_k = 2 √(J² + h² - 2Jh cos k)
    # and the many-body energy is Σ_k ε_k (n_k - 1/2) for n_k ∈ {0,1}, filtered
    # to states with the η-fermion parity dictated by the sector.
    #
    # The R sector has the unpaired k = 0 mode whose bare single-particle
    # energy 2(h - J) changes sign at h = J. In the ferromagnetic phase (h ≤ J)
    # that mode is filled in the BdG vacuum, flipping the required η-parity in
    # the R sector. At h = J the spectrum is degenerate and either convention
    # gives the same multiset. Empirically verified across L=4,6,8 and 30
    # (J, h) points — see commit message for the parity-convention search.
    # ---------------------------------------------------------------------------
    @testset "TFIM exact spectrum (Pfeuty)" begin
        function tfim_hamiltonian(L::Integer, J::Real, h::Real)
            Bt = TensorBasis(L = L, base = 2)
            H_xx = trans_inv_operator(Array(EDKit.spin("XX")), 2, Bt)
            H_z  = trans_inv_operator(Array(EDKit.spin("Z")),  1, Bt)
            Array(-J * H_xx - h * H_z)
        end

        function tfim_analytic_spectrum(L::Integer, J::Real, h::Real)
            iseven(L) || error("Pfeuty formula coded for even L only (got $L).")
            ε(k) = 2 * sqrt(J^2 + h^2 - 2 * J * h * cos(k))
            k_NS = [π * (2m + 1) / L for m in 0:L-1]
            k_R  = [2π * m / L      for m in 0:L-1]
            function sector(ks::Vector{<:Real}, want_even_η::Bool)
                eps = [ε(k) for k in ks]
                out = Float64[]
                for bits in 0:(2^L - 1)
                    iseven(count_ones(bits)) == want_even_η || continue
                    E = sum(eps[m + 1] * (((bits >> m) & 1) - 0.5) for m in 0:L-1)
                    push!(out, E)
                end
                out
            end
            R_even = h <= J  # ferromagnetic phase flips R-sector η-parity
            sort(vcat(sector(k_NS, true), sector(k_R, R_even)))
        end

        # Cover both phases, the QPT, and the trivial-limit reductions.
        cases = [
            (1.0, 0.0),    # pure Ising
            (1.0, 0.3),    # deep ferromagnetic
            (1.0, 0.7),    # ferromagnetic
            (1.0, 1.0),    # quantum critical point
            (1.0, 1.5),    # paramagnetic
            (1.0, 2.0),    # deep paramagnetic
            (0.0, 1.0),    # pure transverse field
            (2.0, 0.5),    # rescaled J
        ]
        for L in (4, 6, 8), (J, h) in cases
            H = tfim_hamiltonian(L, J, h)
            @test H ≈ H'
            spec_num = sort(eigvals(Hermitian(H)))
            spec_an  = tfim_analytic_spectrum(L, J, h)
            @test spec_num ≈ spec_an atol = 1e-9
        end

        # Spot-check the ground-state energy against a closed form: at h = 0
        # the ground state is the doubly-degenerate ferromagnet, E₀ = -J L.
        for L in (4, 6, 8)
            H = tfim_hamiltonian(L, 1.0, 0.0)
            @test minimum(eigvals(Hermitian(H))) ≈ -1.0 * L atol = 1e-10
        end

        # And at J = 0 the ground state is fully polarised along σ^z,
        # E₀ = -h L.
        for L in (4, 6, 8), h in (0.7, 1.5)
            H = tfim_hamiltonian(L, 0.0, h)
            @test minimum(eigvals(Hermitian(H))) ≈ -h * L atol = 1e-10
        end
    end

    # ---------------------------------------------------------------------------
    # Pure classical Ising at h = 0: spectrum is determined by counting
    # antialigned bonds. For an L-site PBC ring with k antialigned bonds
    # (k must be even for closure), the energy is -J(L - 2k) and the number
    # of configurations is 2 × C(L, k).
    # ---------------------------------------------------------------------------
    @testset "Pure Ising domain-wall enumeration (h = 0)" begin
        for L in (4, 6, 8), J in (0.7, 1.5)
            Bt = TensorBasis(L = L, base = 2)
            H = -J * Array(trans_inv_operator(Array(EDKit.spin("XX")), 2, Bt))
            spec = sort(eigvals(Hermitian(H)))
            analytic = Float64[]
            for k in 0:2:L
                E = -J * (L - 2k)
                degeneracy = 2 * binomial(L, k)
                append!(analytic, fill(E, degeneracy))
            end
            sort!(analytic)
            @test spec ≈ analytic atol = 1e-10
        end
    end

    # ---------------------------------------------------------------------------
    # Pure transverse field at J = 0: each spin contributes ±h independently.
    # For m occupied (dgt=1) sites out of L the energy is -h(L - 2m) with
    # multiplicity C(L, m).
    # ---------------------------------------------------------------------------
    @testset "Pure transverse field (J = 0)" begin
        for L in (4, 6, 8), h in (0.5, 1.7)
            Bt = TensorBasis(L = L, base = 2)
            H = -h * Array(trans_inv_operator(Array(EDKit.spin("Z")), 1, Bt))
            spec = sort(eigvals(Hermitian(H)))
            analytic = Float64[]
            for m in 0:L
                append!(analytic, fill(-h * (L - 2m), binomial(L, m)))
            end
            sort!(analytic)
            @test spec ≈ analytic atol = 1e-10
        end
    end
end
