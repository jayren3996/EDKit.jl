function completebasis(f, L::Integer)
    [translationalbasis(x -> f(x) && sum(x)==L-n, k, L) for n in 0:L, k in 0:L-1]
end

completebasis(L::Integer) = completebasis(x -> true, L)

export block_diagonalize
function block_diagonalize(H::Operator; path::String="temp", basis = completebasis(length(H.B)))
    mkpath(path)
    M, I = H.M, H.I
    for n in 0:L, k in 0:L-1
        H = Hamiltonian(M, I, basis[n+1, k+1]) |> Array |> Hermitian
        E, V = if size(H, 1) == 0
            zeros(Float64, 0), zeros(ComplexF64, 0, 0)
        else
            eigen(H)
        end
        @save "$path/$n-$k" E V
    end
end

