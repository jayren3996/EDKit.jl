#---------------------------------------------------------------------------------------------------
# Many-body Lindblad Form
#---------------------------------------------------------------------------------------------------
export lindblad, densitymatrix, expectation
"""
Lindblad Equation:
    ∂ₜρ = -i[H, ρ] + ∑ᵢ LᵢρLᵢ⁺ - 1/2 ∑ᵢ {Lᵢ⁺Lᵢ, ρ} 
The object `Lindblad` store the information of the supper operator. 
"""
struct Lindblad{T1, T2}
    H::Matrix{T1}
    L::Vector{Matrix{T2}}
end

"""
Construction for `Lindblad`.
"""
lindblad(H, L) = Lindblad(Array(H), [Array(l) for l in L])
eltype(::Lindblad{T1, T2}) where {T1, T2} = promote_type(T1, T2)

#---------------------------------------------------------------------------------------------------
# Density Matrix
#---------------------------------------------------------------------------------------------------
"""
Density Matrix under Many-body Basis
"""
struct DensityMatrix{T}
    ρ::Matrix{T}
end
#---------------------------------------------------------------------------------------------------
densitymatrix(ρ::AbstractArray) = DensityMatrix(ρ)
densitymatrix(ψ::AbstractVector) = DensityMatrix(ψ * ψ')
function densitymatrix(i::Integer, L::Integer; base::Integer=2)
    N = base ^ L
    ρ = zeros(N, N)
    ρ[i, i] = 1
    DensityMatrix(ρ)
end
#---------------------------------------------------------------------------------------------------
Array(dm::DensityMatrix) = Hermitian(dm.ρ)
"""
Normalize the density operator so that tr[ρ]=1.
"""
function LinearAlgebra.normalize!(dm::DenseMatrix)
    dm.ρ ./= tr(dm.ρ)
end
#---------------------------------------------------------------------------------------------------
expectation(O::AbstractMatrix, dm::DensityMatrix) = tr(O * dm.ρ)
expectation(O::Hermitian, dm::DensityMatrix) = real(tr(O * dm.ρ))
#---------------------------------------------------------------------------------------------------
function entropy(dm::DensityMatrix; α::Real=1, cutoff::Real=1e-20)
    λ = eigvals(Hermitian(dm.ρ))
    entropy(λ, α=α, cutoff=cutoff)
end

#---------------------------------------------------------------------------------------------------
# Lindblad Evolution
#---------------------------------------------------------------------------------------------------
function *(lb::Lindblad, ρ::Matrix)
    H, L = lb.H, lb.L
    out = -1im * (H * ρ - ρ * H)
    LdL = zeros(eltype(L[1]), size(L[1]))
    for l in L 
        ld = l'
        out += l * ρ * ld
        LdL += ld * l
    end
    LdL ./= -2
    out += LdL * ρ
    out += ρ * LdL 
    out
end
#---------------------------------------------------------------------------------------------------
function (lb::Lindblad)(dm::DensityMatrix, dt::Real=0.05; order::Integer=5)
    ρ = dm.ρ
    mat = Array(ρ)
    dρ = dt * (lb * ρ)
    mat += dρ
    for i = 2:order 
        dρ = (dt/i) * (lb * dρ)
        mat += dρ
    end
    DensityMatrix(mat)
end


#---------------------------------------------------------------------------------------------------
# Single-body Lindblad Form
#---------------------------------------------------------------------------------------------------
export quadraticlindblad, covariancematrix,majoranaform, fermioncorrelation
struct QuardraticLindblad{T1 <: Real, T2 <: Real, T3 <: Real} 
    X::Matrix{T1}
    Y::Matrix{T2}
    Z::Vector{Matrix{T3}}
end
#---------------------------------------------------------------------------------------------------
function quadraticlindblad(
    H::AbstractMatrix, 
    L::AbstractMatrix, 
    M::AbstractVector{<:AbstractMatrix}
)
    B = L * L' 
    X = H - 2 * real(B) + 8 * sum(Ms^2 for Ms in M)
    Y = 4 * imag(B)
    Z = [4 * Ms for Ms in M]
    QuardraticLindblad(X, Y, Z)
end

function quadraticlindblad(
    H::AbstractMatrix, 
    L::AbstractMatrix
)
    B = L * L' 
    X = H - 2 * real(B)
    Y = 4 * imag(B)
    Z = Matrix{Float64}[]
    QuardraticLindblad(X, Y, Z)
end

#---------------------------------------------------------------------------------------------------
# Majorana Covariant Matrix
#---------------------------------------------------------------------------------------------------
struct CovarianceMatrix{T <: Real}
    Γ::Matrix{T}
    N::Integer
end
#---------------------------------------------------------------------------------------------------
function covariancematrix(Γ::AbstractMatrix{<:Real})
    N = size(Γ, 1) ÷ 2 
    CovarianceMatrix(Γ, N)
end

function covariancematrix(n::AbstractVector{<:Integer})
    N = length(n)
    D = diagm(-2 * n .+ 1)
    Z = zeros(N, N)
    CovarianceMatrix([Z D; -D Z], N)
end
#---------------------------------------------------------------------------------------------------
function fermioncorrelation(cm::CovarianceMatrix)
    Γ, n = cm.Γ, cm.N
    Γ11, Γ12, Γ21, Γ22 = Γ[1:n,1:n], Γ[1:n,n+1:2n], Γ[n+1:2n,1:n], Γ[n+1:2n,n+1:2n]
    A = (Γ21 - Γ12 + 1im * Γ11 + 1im * Γ22) / 4 + I / 2
    B = (Γ21 + Γ12 + 1im * Γ11 - 1im * Γ22) / 4
    A, B
end

function fermioncorrelation(cm::CovarianceMatrix, i::Integer)
    Γ, n = cm.Γ, cm.N
    Γ11, Γ12, Γ21, Γ22 = Γ[1:n,1:n], Γ[1:n,n+1:2n], Γ[n+1:2n,1:n], Γ[n+1:2n,n+1:2n]
    if i == 1
        (Γ21 - Γ12 + 1im * Γ11 + 1im * Γ22) / 4 + I / 2
    elseif i == 2
        (Γ21 + Γ12 + 1im * Γ11 - 1im * Γ22) / 4
    elseif i == 3
        (-Γ21 - Γ12 + 1im * Γ11 - 1im * Γ22) / 4
    else
        error("i should equals to 1, 2, or 3.")
    end
end


#---------------------------------------------------------------------------------------------------
# Evolution of Covariance Matrix
#---------------------------------------------------------------------------------------------------
function *(ql::QuardraticLindblad, Γ::AbstractMatrix{<:Real})
    transpose(ql.X) * Γ + Γ * ql.X + sum(transpose(Zs) * Γ * Zs for Zs in ql.Z)
end
#---------------------------------------------------------------------------------------------------
function (ql::QuardraticLindblad)(cm::CovarianceMatrix, dt::Real=0.05; order::Integer=5)
    Γ = cm.Γ
    mat = Array(Γ)
    dΓ = dt * (ql * Γ + ql.Y)
    mat += dΓ
    for i=2:order 
        dΓ = (dt/i) * (ql * dΓ)
        mat += dΓ
    end
    CovarianceMatrix(mat, cm.N)
end

#---------------------------------------------------------------------------------------------------
"""
    majoranaform(A::AbstractMatrix, B::AbstractMatrix)

Return the Majorana quadratic form
    Ĥ = -i/4 ∑ Hᵢⱼ ωᵢωⱼ 
from the fermion quadratic form
    Ĥ = 1/2 ∑(Aᵢⱼ cᵢ⁺cⱼ + Bᵢⱼcᵢ⁺cⱼ⁺ + h.c.).
"""
function majoranaform(A::AbstractMatrix, B::AbstractMatrix)
    AR, AI, BR, BI = real(A), imag(A), real(B), imag(B)
    [-AI-BI AR-BR; -AR-BR -AI+BI]
end

