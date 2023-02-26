#--------------------------------------------------------------------------------
# Gaussian States
#--------------------------------------------------------------------------------
export GaussianState
"""
Particle-number-preserving fermionic Gaussian state

A Gaussian state `s` is represented by a matrix `B`,
    |ψ⟩ = ∏ bₖ⁺|0⟩ = ∑ cᵢ⁺ Bᵢₖ.
A `GaussianState` object can be access as a matrix.
The element `s[i,j]` is the two-point function ⟨cᵢ⁺cⱼ⟩.

The boolean value `N` keeps track of the orthogonality of the representation.
"""
struct GaussianState{T<:Number}
    B::Matrix{T}    # Matrix representation of free modes
    N::Bool         # Whether the Gaussian state is normalized
end

#--------------------------------------------------------------------------------
# Initialization
"""
    GaussianState(dtype::DataType, pos::AbstractVector{<:Integer}, L::Integer)

Initialize GaussianState with given data type, particle positions and system length.
"""
function GaussianState(dtype::DataType, pos::AbstractVector{<:Integer}, L::Integer)
    B = zeros(dtype, L, length(pos))
    for i = 1:length(pos)
        B[pos[i], i] = 1
    end
    GaussianState(B, true)
end

"""
    GaussianState(pos::AbstractVector{<:Integer}, L::Integer) 

Initialize GaussianState with given particle positions and system length.
"""
function GaussianState(pos::AbstractVector{<:Integer}, L::Integer) 
    GaussianState(ComplexF64, pos, L)
end

"""
    GaussianState(dtype::DataType, vec::AbstractVector{Bool})

Initialize GaussianState with vector indicating particle positions.
"""
function GaussianState(dtype::DataType, vec::AbstractVector{Bool})
    pos = Int64[]
    for i = 1:length(vec) 
        if vec[i] 
            push!(pos, i)
        end
    end
    GaussianState(dtype, pos, length(vec))
end

function GaussianState(vec::AbstractVector{Bool}) 
    GaussianState(ComplexF64, vec)
end

"""
    GaussianState(;L::Integer, N::Integer, config::String="Z2")

More initial configuration.
"""
function GaussianState(;L::Integer, N::Integer, config::String="Z2")
    pos = if config == "center"
        l = (L-N)÷2
        l+1:l+N
    elseif config == "left"
        1:N
    elseif config == "right"
        L-N+1:L
    elseif config == "random"
        shuffle(1:L)[1:N]
    elseif config[1] == 'Z'
        n = parse(Int, config[2])
        1:n:1+n*(N-1)
    else
        error("Invalid configureation: $config.")
    end
    GaussianState(pos, L)
end

#--------------------------------------------------------------------------------
# Properties
Base.eltype(::GaussianState{T}) where T = T
Base.length(s::GaussianState) = size(s.B, 1)

"""
Compute correlation ⟨cᵢ⁺cⱼ⟩
"""
Base.getindex(s::GaussianState, i::Integer, j::Integer) = dot(s.B[i,:], s.B[j,:])

"""
Compute correlation matrix
"""
Base.getindex(s::GaussianState, ::Colon) = conj(s.B) * transpose(s.B)
Base.getindex(s::GaussianState, ::Colon, ::Colon) = conj(s.B) * transpose(s.B)

"""
Compute particle density on each site.
"""
function LinearAlgebra.diag(s::GaussianState)
    n = Vector{Float64}(undef,length(s))
    for i in 1:length(s)
        n[i] = real(s[i,i])
    end
    n
end

"""
Return particle number.
"""
LinearAlgebra.rank(s::GaussianState) = size(s.B, 2)

"""
Return the normalized GaussianState.
"""
function LinearAlgebra.normalize(s::GaussianState) 
    if s.N 
        s
    else
        GaussianState(orthogonalize(s.B), true)
    end
end

orthogonalize(B::AbstractMatrix) = Matrix(qr(B).Q)
#--------------------------------------------------------------------------------
"""
    entropy(s::GaussianState, i::AbstractVector{<:Integer})

Entaglement entropy of gaussian state with system A chosen to be the sites {i}.
"""
function entropy(s::GaussianState, i::AbstractVector{<:Integer})
    B = eltype(s) <: Real ? Float64.(s.B[i,:]) : ComplexF64.(s.B[i,:])
    vals = svdvals(B).^2
    EE = 0.0
    for x in vals
        if x < 1.0
            x < 1e-14 && continue
            EE -= x*log(x)+(1-x)*log(1-x)
        else
            @assert x-1 < 1e-6 "Got a Schmidt value λ = $x."
        end
    end
    EE
end

"""
    entropy(s::GaussianState, i::AbstractVector{<:Integer}, j::AbstractVector{<:Integer})

Mutual information for GaussianState.
"""
function entropy(
    s::GaussianState, 
    i::AbstractVector{<:Integer}, 
    j::AbstractVector{<:Integer}
)
    SA = entropy(s, i)
    SB = entropy(s, j)
    SAB = entropy(s, vcat(i,j))
    SA + SB - SAB
end


#--------------------------------------------------------------------------------
# Evolution under unitary matrices
"""
Data type for Unitary matrices
"""
struct Unitary{T <: AbstractMatrix} 
    M::T
end
#--------------------------------------------------------------------------------
export evo_operator
"""
    evo_operator(H::Hermitian, dt::Real)

Exponential for Hermitian matrix.
"""
function evo_operator(H::Hermitian, dt::Real)
    vals, vecs = eigen(H)
    D = exp.(-dt * im * vals) |> Diagonal
    Unitary(vecs * D * vecs')
end

"""
    evo_operator(H::AbstractMatrix, dt::Real; order::Integer=10)

Exponential for general matrix using series expansion.
"""
function evo_operator(H::AbstractMatrix, dt::Real; order::Integer=10)
    N = round(Int, 10 * dt * maximum(abs, H))
    iszero(N) && (N = 1)
    expm(-dt*im*H/N, order)^N
end

"""
    expm(A::AbstractMatrix, order::Integer=10)

Matrix exponential using Taylor expansion.
"""
function expm(A::AbstractMatrix, order::Integer=10)
    mat = I + A / order
    order -= 1
    while order > 0
        mat = A * mat
        mat ./= order
        mat += I
        order -= 1
    end
    mat
end

"""
    expv(A::AbstractMatrix, v::AbstractVecOrMat, order::Integer=10)

Compute exp(A)*v using Taylor expansion
"""
function expv(A::AbstractMatrix, v::AbstractVecOrMat, order::Integer=10)
    vec = v + A * v / order
    order -= 1
    while order > 0
        vec = A * vec
        vec ./= order
        vec += v
        order -= 1
    end
    vec
end

export rand_unitary
"""
    rand_unitary(N::Integer)

Generate random unitary matrix. 
"""
function rand_unitary(N::Integer)
    F = qr(randn(ComplexF64, N,N)) 
    Unitary(Matrix(F.Q))
end
#--------------------------------------------------------------------------------
"""
Multiply Unitary matrix to GaussianState.
"""
function *(U::Unitary, s::GaussianState) 
    GaussianState(U.M * s.B, true)
end

"""
Multiply general matrix to GaussianState, and by default normalize the output.
"""
function *(M::AbstractMatrix, s::GaussianState; normalize::Bool=true)
    B_new = normalize ? orthogonalize(M * s.B) : M * s.B
    GaussianState(B_new, true)
end

export apply!
"""
    apply!(U::Unitary, s::GaussianState, ind::AbstractVector{<:Integer})

Apply local unitary gate to GaussianState `s` on sites `inds`.
"""
function apply!(U::Unitary, s::GaussianState, inds::AbstractVector{<:Integer})
    s.B[inds, :] = U.M * s.B[inds, :]
    s
end

#--------------------------------------------------------------------------------
# Quantum Jumps
#--------------------------------------------------------------------------------
"""
Quasi Mode d⁺ = ∑ Vⱼ cⱼ⁺
"""
struct QuasiMode{T<:Number}
    I::Vector{Int64}
    V::Vector{T}
    L::Int64
end
#--------------------------------------------------------------------------------
vector(qm::QuasiMode) = sparsevec(qm.I, qm.V, qm.L)
inner(qm::QuasiMode, s::GaussianState) = vec(qm.V' * s.B[qm.I, :])
#--------------------------------------------------------------------------------
function particle_exp(qm, γdt)
    v = hcat(qm.V, nullspace(qm.V'))
    d = ones(length(qm.I))
    d[1] = exp(-γdt/2)
    v * Diagonal(d) * v'
end
#--------------------------------------------------------------------------------
function hole_exp(qm, γdt)
    v = hcat(qm.V, nullspace(qm.V'))
    d = fill(exp(-γdt/2), length(qm.I))
    d[1] = 1
    v * Diagonal(d) * v'
end
#--------------------------------------------------------------------------------
abstract type QuantumJump end
export QJParticle, QJHole, QJDrain, QJSource
struct QJParticle{T1, T2, T3<:Real} <: QuantumJump
    M::QuasiMode{T1}
    P::Matrix{T2}
    γdt::T3
end
struct QJHole{T1, T2, T3<:Real} <: QuantumJump
    M::QuasiMode{T1}
    P::Matrix{T2}
    γdt::T3
end
struct QJDrain{T1, T2, T3<:Real} <: QuantumJump
    M::QuasiMode{T1}
    P::Matrix{T2}
    γdt::T3
end
struct QJSource{T1, T2, T3<:Real} <: QuantumJump
    M::QuasiMode{T1}
    P::Matrix{T2}
    γdt::T3
end
#--------------------------------------------------------------------------------
function QJParticle(I::AbstractVector{<:Integer}, V::AbstractVector, L::Integer, γdt::Real)
    M = QuasiMode(I, V, L) 
    P = particle_exp(M, γdt)
    QJParticle(M, P, γdt)
end
function QJHole(I::AbstractVector{<:Integer}, V::AbstractVector, L::Integer, γdt::Real)
    M = QuasiMode(I, V, L) 
    P = hole_exp(M, γdt)
    QJHole(M, P, γdt)
end
function QJDrain(I::AbstractVector{<:Integer}, V::AbstractVector, L::Integer, γdt::Real)
    M = QuasiMode(I, V, L) 
    P = particle_exp(M, γdt)
    QJDrain(M, P, γdt)
end
function QJSource(I::AbstractVector{<:Integer}, V::AbstractVector, L::Integer, γdt::Real)
    M = QuasiMode(I, V, L) 
    P = hole_exp(M, γdt)
    QJSource(M, P, γdt)
end
#--------------------------------------------------------------------------------
function *(qj::QuantumJump, s::GaussianState)
    Q, s2 = jump(qj, s)
    Q && return s2 
    B = s2.B
    B[qj.M.I, :] .= qj.P * B[qj.M.I, :]
    GaussianState(orthogonalize(B), true)
end
#--------------------------------------------------------------------------------
function *(qjs::AbstractVector{<:QuantumJump}, s::GaussianState)
    N = length(qjs)
    Qs = Vector{Bool}(undef, N)
    for i in 1:N 
        Qs[i], s = jump(qjs[i], s)
    end
    B = s.B
    normQ = true
    for i in 1:N 
        qj = qjs[i]
        Qs[i] && continue
        B[qj.M.I, :] .= qj.P * B[qj.M.I, :]
        normQ = false
    end
    normQ ? S : GaussianState(orthogonalize(B), true)
end
#--------------------------------------------------------------------------------
function jump(qj::QJParticle, s::GaussianState)
    p = inner(qj.M, s)
    rand() > real(dot(p, p)) * qj.γdt && return (false, s)
    B = replace_vector(s.B, vector(qj.M), p)
    true, GaussianState(B, true)
end
function jump(qj::QJHole, s::GaussianState)
    p = inner(qj.M, s)
    rand() > (1-real(dot(p, p))) * qj.γdt && return (false, s)
    B = avoid_vector(s.B, vector(qj.M), p)
    true, GaussianState(B, true)
end
function jump(qj::QJDrain, s::GaussianState)
    p = inner(qj.M, s)
    rand() > real(dot(p, p)) * qj.γdt && return (false, s)
    B = delete_vector(s.B, vector(qj.M), p)
    true, GaussianState(B, true)
end
function jump(qj::QJSource, s::GaussianState)
    p = inner(qj.M, s)
    rand() > (1-real(dot(p, p))) * qj.γdt && return (false, s)
    B = insert_vector(s.B, vector(qj.M))
    true, GaussianState(B, true)
end
#--------------------------------------------------------------------------------
export measure
"""
Projective Measure
"""
function measure(qm::QuasiMode, s::GaussianState)
    v = vector(qm)
    p = (v' * s.B)[:]
    if rand() < real(dot(p, p))
        B = replace_vector(s.B, v, p)
        true, GaussianState(B, true)
    else
        B = avoid_vector(s.B, v, p)
        false, GaussianState(B, true)
    end
end

#--------------------------------------------------------------------------------
# Conditional Gates
#--------------------------------------------------------------------------------
export ConditionalJump
struct ConditionalJump{T<:QuantumJump, Tm<:AbstractMatrix}
    J::T
    U::Tm
end
#--------------------------------------------------------------------------------
function *(cj::ConditionalJump, s::GaussianState)
    q, s2 = jump(cj.J, s)
    ind = cj.J.M.I
    if q 
        s2.B[ind, :] .= cj.U * s2.B[ind, :] 
        s2
    else
        s2.B[ind, :] .= cj.J.P * s2.B[ind, :]
        GaussianState(orthogonalize(s2.B), true)
    end
end
#--------------------------------------------------------------------------------
function *(cjs::AbstractVector{<:ConditionalJump}, s::GaussianState)
    N = length(cjs)
    Qs = Vector{Bool}(undef, N)
    for i in 1:N
        Qs[i], s = jump(cjs[i].J, s)
    end
    normQ = true
    B = s.B
    for i in 1:N 
        ind = cjs[i].J.M.I
        if Qs[i]
            B[ind, :] .= cjs[i].U * B[ind, :] 
        else
            B[ind, :] .= cjs[i].J.P * B[ind, :]
            normQ = false
        end
    end
    normQ ? s : GaussianState(orthogonalize(B), true)
end
#--------------------------------------------------------------------------------
export ConditionalMeasure
struct ConditionalMeasure{T<:QuasiMode, Tm1<:AbstractMatrix, Tm2<:AbstractMatrix}
    M::T
    Ut::Tm1
    Uf::Tm2
end
#--------------------------------------------------------------------------------
function *(cm::ConditionalMeasure, s::GaussianState)
    q, s2 = measure(cm.M, s)
    ind = cj.M.I
    s2.B[ind, :] .= (q ? cm.Ut : cm.Uf) * s2.B[ind, :] 
    s2
end


#--------------------------------------------------------------------------------
# Helper
#--------------------------------------------------------------------------------
"""
Given a vector space spanned by column of `A`, and a vector `v`, find the nullspace 
of `v` inside A. The inpute `vA` is the pre-computed inner product ⟨v|Aᵢ⟩.

The algorithm first choose the pivot vector |Aₚ⟩ with largest overlap ‖⟨v|Aₚ⟩‖ with
vector `v`, then calculate the vectors 
    |Ãᵢ⟩ := |Aᵢ⟩ - (⟨v|Aᵢ⟩/⟨v|Aₚ⟩) * |Aₚ⟩,  i ≠ p.
The set of vectors {|Ãᵢ⟩} span the null space of `v`, though not orthogonal.
"""
function delete_vector(A::AbstractMatrix, v::AbstractVector, vA::AbstractVector)
    mat = Matrix{eltype(vA)}(undef, size(A, 1), size(A, 2)-1)
    p = argmax(abs.(vA))
    Ap = A[:, p] / vA[p]
    mat[:, 1] .= v 
    for i = 1:p-1
        mat[:, i] .= A[:, i] - vA[i] * Ap
    end
    for i = p+1:size(A, 2)
        mat[:, i-1] .= A[:, i] - vA[i] * Ap
    end
    orthogonalize(mat)
end
#--------------------------------------------------------------------------------
function replace_vector(A::AbstractMatrix, v::AbstractVector, vA::AbstractVector)
    mat = Matrix{eltype(vA)}(undef, size(A, 1), size(A, 2))
    mat[:, 1] .= v 
    mat[:, 2:end] .= delete_vector(A, v, vA)
    mat
end
#--------------------------------------------------------------------------------
"""
Compute null vectors:
    |A'ᵢ⟩ = |Aᵢ⟩ - ⟨v|Aᵢ⟩|v⟩
"""
function avoid_vector(A::AbstractMatrix, v::AbstractVector, vA::AbstractVector)
    mat = Matrix{eltype(vA)}(undef, size(A, 1), size(A, 2))
    for i = 1:size(A, 2)
        mat[:, i] .= A[:, i] - vA[i] * v
    end
    orthogonalize(mat)
end
#--------------------------------------------------------------------------------
function insert_vector(A::AbstractMatrix, v::AbstractVector)
    mat = hcat(v, A)
    orthogonalize(mat)
end

