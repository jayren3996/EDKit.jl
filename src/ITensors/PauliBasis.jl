const PAULI_CONVERSION = let T = zeros(ComplexF64, 4, 2, 2)
    for i in 1:4 
        T[i, :, :] += conj.(PAULI[i]) / 2
    end
    T
end
const PAULI_INV_CONV = let T = zeros(ComplexF64, 2, 2, 4)
    for i in 1:4 
        T[:, :, i] += PAULI[i]
    end
    T
end
#---------------------------------------------------------------------------------------------------
export pauli
"""
    pauli(i, L=1) 

Return ith (size-L) Pauli matrices.
"""
function pauli(i::Integer, L::Integer=1) 
    inds = digits(i-1, base=4, pad=L) |> reverse!
    σ(inds)
end
#---------------------------------------------------------------------------------------------------
function pauli(Is::AbstractVector{<:Real})
    L = round(Integer, log(4, length(Is)))
    @assert 4^L == length(Is)
    out = zeros(ComplexF64, 2^L, 2^L)
    for i in eachindex(Is)
        out += Is[i] * pauli(i, L)
    end
    out
end
#---------------------------------------------------------------------------------------------------
"""
Compute tr(A * B), where A is a Pauli matrix.
"""
function pauli_mtr(A::SparseMatrixCSC, B::AbstractMatrix)
    rows, vals = rowvals(A), nonzeros(A)
    res = zero(promote_type(eltype(A), eltype(B)))
    for j in axes(A, 2)
        val, row = vals[j], rows[j]
        res += val * B[j, row]
    end
    #res = tr(A * B)
    res
end
#---------------------------------------------------------------------------------------------------
export pauli_list
"""
    pauli_list(A)

Return pauli components of a 2×2 matrix `A`.
"""
function pauli_list(A::AbstractMatrix, T::DataType=Float64)
    n = round(Integer, log(2, size(A, 1)))
    N = 2^n
    @assert N == size(A, 1) "dimension D = $(size(A, 1)) not right."
    list = Vector{T}(undef, N^2)
    for i in eachindex(list)
        c = pauli_mtr(pauli(i, n), A) / N 
        list[i] = T <: Real ? real(c) : c 
    end
    return list
end
#---------------------------------------------------------------------------------------------------
"""
    pauli_op_mat(f)

Matrix representation for f(ρ).
"""
function pauli_op_mat(f, L::Integer=1)
    N = 4^L
    mat = Matrix{Float64}(undef, N, N)
    for i in axes(mat, 2)
        mat[:, i] = f(pauli(i, L)) |> pauli_list
    end
    return mat
end

#---------------------------------------------------------------------------------------------------
"""
    dissipation(L,ρ)

Compute dissipation
    D[L]ρ = L⋅ρ⋅L⁺ - 1/2{L⁺L,ρ}
"""
function dissipation(L::AbstractMatrix, ρ::AbstractMatrix)
    Lρ = L * ρ
    Ld = L'
    LLρ = Ld * Lρ
    return Lρ * Ld - (LLρ + LLρ') / 2
end
#---------------------------------------------------------------------------------------------------
"""
    commutation(L,ρ)

Compute commutation
    -i[H, ρ]
"""
function commutation(H::AbstractMatrix, ρ::AbstractMatrix)
    Hρ = -1im * H * ρ 
    Hρ + Hρ'
end
#---------------------------------------------------------------------------------------------------
export dissipation_mat
"""
    dissipation_mat(L)

Matrix representation for D[L].
"""
function dissipation_mat(L::AbstractMatrix)
    n = round(Integer, log(2, size(L, 1)))
    @assert 2^n == size(L, 1) "dimension D = $(size(L, 1)) not right."
    f = ρ -> dissipation(L,ρ)
    pauli_op_mat(f, n)
end
#---------------------------------------------------------------------------------------------------
export commutation_mat
"""
    commutation_mat(L)

Matrix representation for -i[H,⋅].
"""
function commutation_mat(H::AbstractMatrix)
    n = round(Integer, log(2, size(H, 1)))
    @assert 2^n == size(H, 1) "dimension D = $(size(H, 1)) not right."
    f = ρ -> commutation(H,ρ)
    pauli_op_mat(f, n)
end


#---------------------------------------------------------------------------------------------------
# Define Pauli basis
#---------------------------------------------------------------------------------------------------
ITensors.space(::SiteType"Pauli") = 4

ITensors.state(::StateName"I", ::SiteType"Pauli") = [1, 0, 0, 0]
ITensors.state(::StateName"X", ::SiteType"Pauli") = [0, 1, 0, 0]
ITensors.state(::StateName"Y", ::SiteType"Pauli") = [0, 0, 1, 0]
ITensors.state(::StateName"Z", ::SiteType"Pauli") = [0, 0, 0, 1]
ITensors.state(::StateName"Up", ::SiteType"Pauli") = [1/2, 0, 0, 1/2]
ITensors.state(::StateName"Dn", ::SiteType"Pauli") = [1/2, 0, 0, -1/2]
#---------------------------------------------------------------------------------------------------
# MPS/MPO
#---------------------------------------------------------------------------------------------------
export mps2pmps
function mps2pmps(ψ::MPS, S::AbstractVector)
    s = siteinds(ψ)
    L = length(s)
    psi = MPS(L)
    
    psi[1] = begin
        l1 = linkind(ψ, 1)
        Cl = combiner(l1, l1', tags="Link,l=1")
        C = ITensor(PAULI_CONVERSION, S[1], s[1]', s[1])
        ψ[1]' * conj(ψ[1]) * C * Cl |> real
    end
    
    for i in 2:L-1 
        li = linkind(ψ, i)
        Cl2 = combiner(li, li', tags="Link,l=$i")
        C = ITensor(PAULI_CONVERSION, S[i], s[i]', s[i])
        psi[i] = ψ[i]' * conj(ψ[i]) * C * Cl * Cl2 |> real
        Cl = Cl2
    end

    psi[L] = begin
        C = ITensor(PAULI_CONVERSION, S[L], s[L]', s[L])
        ψ[L]' * conj(ψ[L]) * C * Cl |> real
    end
    orthogonalize!(psi, 1)
end
#---------------------------------------------------------------------------------------------------
export pmps2mpo
"""
    pmps2mpo(ψ, s)

Convert Pauli MPS to MPO.
"""
function pmps2mpo(ψ::MPS, s::AbstractVector)
    L = length(s)
    S = siteinds(ψ)
    @assert length(ψ) == L 
    O = MPO(L) 
    for i in eachindex(s)
        si = s[i]
        C = ITensor(PAULI_INV_CONV, si', si, S[i])
        O[i] = C * ψ[i]
    end
    O
end
#---------------------------------------------------------------------------------------------------
export mpo2pmps
"""
    mpo2pmps(O, S)

Convert MPO to Pauli MPS.
"""
function mpo2pmps(O::MPO, S::AbstractVector)
    s = siteinds(O)
    L = length(s)
    ψ = MPS(L) 
    for i in eachindex(s)
        C = ITensor(PAULI_CONVERSION, S[i], s[i]...)
        ψ[i] = C * O[i] |> real
    end
    orthogonalize!(ψ, 1)
end

#---------------------------------------------------------------------------------------------------
# Function on Pauli basis
#---------------------------------------------------------------------------------------------------
"""
Return the operator norm (N = |⟨I|ρ⟩|).

The renormed (ρ/N = I + ⋯) is the correct density matrix. 
"""
function density_norm(ψ::MPS)
    s = siteinds(ψ)
    V = ITensor(1.0)
    for j in eachindex(s)
        V *= ψ[j] * state(s[j], 1)
    end
    scalar(V) |> abs
end
#---------------------------------------------------------------------------------------------------
function density_expect(ψ::MPS, o::Integer)
    s = siteinds(ψ)
    L = let V = ITensor(1.0)
        l = Vector{ITensor}(undef, length(s)); l[1] = V
        for j in 2:length(s) 
            V *= ψ[j-1] * state(s[j-1], 1)
            l[j] = V
        end
        l
    end
    R = let V = ITensor(1.0)
        l = Vector{ITensor}(undef, length(s)); l[end] = V
        for j in length(s)-1:-1:1
            V *= ψ[j+1] * state(s[j+1], 1)
            l[j] = V
        end
        l
    end
    out = Vector{Float64}(undef, length(s))
    for j in eachindex(s)
        V = ψ[j] * state(s[j], o)
        out[j] = L[j] * V * R[j] |> scalar |> real 
    end
    N = L[2]*R[1] |> scalar |> real 
    out ./ N
end
