#---------------------------------------------------------------------------------------------------
# Many-body Lindblad Form
#---------------------------------------------------------------------------------------------------
export lindblad, evo
"""
Lindblad Equation:
    ∂ₜρ = -i[H, ρ] + ∑ᵢ LᵢρLᵢ⁺ - 1/2 ∑ᵢ {Lᵢ⁺Lᵢ, ρ} 
The object `Lindblad` store the information of the supper operator. 
"""
struct Lindblad{T1, T2}
    H::Matrix{T1}
    L::Vector{Matrix{T2}}
end

lindblad(H, L) = Lindblad(Array(H), [Array(l) for l in L])
Base.eltype(::Lindblad{T1, T2}) where {T1, T2} = promote_type(T1, T2)
#---------------------------------------------------------------------------------------------------
function *(lb::Lindblad, ρ::AbstractMatrix)
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
function evo(
    lb::Lindblad, ρ::AbstractMatrix,
    dt::Real=0.05; order::Integer=5
)
    mat = Array(ρ)
    dρ = dt * (lb * ρ)
    mat += dρ
    for i = 2:order 
        dρ = (dt/i) * (lb * dρ)
        mat += dρ
    end
    mat
end

#---------------------------------------------------------------------------------------------------
# Single-body Lindblad Form
#---------------------------------------------------------------------------------------------------
export quadraticlindblad
struct QuardraticLindblad{T1 <: Real, T2 <: Real, T3 <: Real} 
    X::Matrix{T1}
    Y::Matrix{T2}
    Z::Vector{Matrix{T3}}
end

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
function *(ql::QuardraticLindblad, Γ::AbstractMatrix{<:Real})
    transpose(ql.X) * Γ + Γ * ql.X + sum(transpose(Zs) * Γ * Zs for Zs in ql.Z)
end
#---------------------------------------------------------------------------------------------------
function evo(
    ql::QuardraticLindblad, Γ::AbstractMatrix{<:Real},
    dt::Real=0.05; order::Integer=5
)
    mat = Array(Γ)
    dΓ = dt * (ql * Γ + ql.Y)
    mat += dΓ
    for i=2:order 
        dΓ = (dt/i) * (ql * dΓ)
        mat += dΓ
    end
    mat
end

#---------------------------------------------------------------------------------------------------
# Majorana Form
#---------------------------------------------------------------------------------------------------
export majoranaform, fermioncorrelation
"""
Majorana form:
    Ĥ = -i/4 ∑ Hᵢⱼ ωᵢωⱼ 
"""
function majoranaform(A::AbstractMatrix, B::AbstractMatrix)
    AR, AI, BR, BI = real(A), imag(A), real(B), imag(B)
    [-AI-BI AR-BR; -AR-BR -AI+BI]
end
#---------------------------------------------------------------------------------------------------
function fermioncorrelation(Γ::AbstractMatrix)
    n = size(Γ, 1) ÷ 2 
    Γ11, Γ12, Γ21, Γ22 = Γ[1:n,1:n], Γ[1:n,n+1:2n], Γ[n+1:2n,1:n], Γ[n+1:2n,n+1:2n]
    A = (Γ21 - Γ12 + 1im * Γ11 + 1im * Γ22) / 4 + I / 2
    B = (Γ21 + Γ12 + 1im * Γ11 - 1im * Γ22) / 4
    A, B
end
