## Lindblad Equation
The Lindblad master equation for Makovian open system is
```math
\frac{d}{dt} \hat\rho = -i[\hat H, \hat \rho] + \sum_{\mu=1}^{m} \hat L_\mu \hat\rho \hat L_\mu^\dagger -\frac{1}{2} \sum_{\mu=1}^{m} \{\hat L_\mu^\dagger \hat L_\mu, \hat \rho\}
```
The Lindblad super-operator is represented by object `Lindblad`:
```julia
struct Lindblad{T1, T2}
    H::Matrix{T1}
    L::Vector{Matrix{T2}}
end
```
To create the object, simply use the command `lindblad(H, L)` where `H` is the Hamiltonian and `L` is the list of jump operators. The Lindblad operator act on the `DensityMatrix` object:
```julia
struct DensityMatrix{T}
    ρ::Matrix{T}
end
```
which is created by the following constructors:
```julia
densitymatrix(ρ::AbstractArray) = DensityMatrix(ρ)
densitymatrix(ψ::AbstractVector) = DensityMatrix(ψ * ψ')
densitymatrix(i::Integer, L::Integer; base::Integer=2)
```
Given a Lindblad super-operator and a density matrix, the time evolution is simulated by
```julia
(lb::Lindblad)(dm::DensityMatrix, dt::Real=0.05; order::Integer=5)
```

## Free Fermion Case
When the jump operators contain only the linear Majorana operators, the Lindblad equation preserve Gaussianity. For jump operator contains up to quadratic Majorana terms, the evolution will break the Gaussian form, however, the $2n$-point correlation is still solvable for free fermion system. It is better to expressed the problem in Majorana operators rather than the ordinary fremion operators, the former is defined here as
```math
\left[\begin{array}{c} \omega_{i} \\ \omega_{i+N} \end{array}\right]
= \left[\begin{array}{cc} 
	1 & 1 \\ 
	i & -i 
\end{array}\right] 
\left[\begin{array}{c} 
	c_i \\ c_i^\dagger 
\end{array}\right].
```
A fermion bilinear of the form
```math
\hat H_{\mathrm{free}} = \sum_{i,j=1}^N A_{ij} c_i^\dagger c_j + \frac{1}{2}\sum_{i,j=1}^N B_{ij} c_i c_j + \frac{1}{2}\sum_{i,j=1}^N B_{ij}^* c_j^\dagger c_i^\dagger
```
can be brought to Majorana form
```math
\hat H_{\mathrm{free}} = -\frac{i}{4} \sum_{i,j=1}^{2N} H_{ij} \omega_i \omega_j,
```
where the single-body matrix $H$ is a `2N * 2N` real anti-symmetric matrix:
```math
	H = \left[\begin{array}{cc} 
		-A^I - B^I & A^R - B^R \\
    	-A^R - B^R &  -A^I + B^I 
	\end{array}\right].
```
This can be done using the function `majoranaform(A::AbstractMatrix, B::AbstractMatrix)`.


We assume that the jump operator has up to quadratic Majorana terms. In particular, we denote the linear terms and the Hermitian quadratic terms as
```math
\hat L_r = \sum_{j=1}^{2N} L^r_{j} \omega_j, \quad \hat L_s = \sum_{j,k=1}^{2N} M^s_{jk} \omega_j \omega_k.
```
Consider the *covariance matrix*
```math
\Gamma_{ij} = i\langle \omega_i \omega_j\rangle,
```
its equation of motion is:
```math
\partial_t \Gamma = X^T\cdot\Gamma + \Gamma \cdot X + \sum_s (Z^s)^T \cdot \Gamma\cdot Z^s + Y,
```
where
```math
X = H - 2B^R + 8 \sum_s (\mathrm{Im} M^s)^2, \quad Y = 4B^I, \quad Z = 4 \mathrm{Im} M^s.
```
The Lindblad super-operator for this free fermion system is stored in the object
```julia
struct QuardraticLindblad{T1 <: Real, T2 <: Real, T3 <: Real} 
    X::Matrix{T1}
    Y::Matrix{T2}
    Z::Vector{Matrix{T3}}
end
```
There are two constructors for the `QuardraticLindblad`:
```julia
quadraticlindblad(H::AbstractMatrix, L::AbstractMatrix, M::AbstractVector{<:AbstractMatrix})
quadraticlindblad(H::AbstractMatrix, L::AbstractMatrix)
```
The `QuardraticLindblad` object act on the `CovarianceMatrix` object:
```julia
struct CovarianceMatrix{T <: Real}
    Γ::Matrix{T}
    N::Integer
end
```
created by
```julia
covariancematrix(Γ::AbstractMatrix{<:Real})
covariancematrix(n::AbstractVector{<:Integer})
```
The evolution of the covariance matrix is simulated by the following command:
```julia
(ql::QuardraticLindblad)(cm::CovarianceMatrix, dt::Real=0.05; order::Integer=5)
```