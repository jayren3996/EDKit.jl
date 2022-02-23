# EDKit.jl

Julia package for general many-body exact diagonalization calculation. The package provide a general Hamiltonian constructing routine for specific symmetry sectors. The functionalities can be extended with user-defined bases.

## Installation

Run the following script in the ```Pkg REPL``` environment:

```julia
pkg> add EDKit
```

## Documentation

https://docs.juliahub.com/EDKit/JkYPS/0.3.11/

## Examples

### Radom XXZ Model with Random

Consider the Hamiltonian `H = ∑ᵢ (σᵢˣσᵢ₊₁ˣ + σᵢʸσᵢ₊₁ʸ + hᵢσᵢᶻσᵢ₊₁ᶻ)`. We choose the system size to be `L=10`. The Hamiltonian need 3 generic information: 

1. Local operators represented by matrices;
2. Site indices where each local operator acts on;
3. Basis, if use the default tensor-product basis, only need to provide the system size.

The following script generate the information we need to generate XXZ Hamiltonian:

```julia
L = 10
mats = [
    fill(spin("XX"), L);
    fill(spin("YY"), L);
    [randn() * spin("ZZ") for i=1:L]
]
inds = [
    [[i, mod(i, L)+1] for i=1:L];
    [[i, mod(i, L)+1] for i=1:L];
    [[i, mod(i, L)+1] for i=1:L]
]
H = operator(mats, inds, L)
```

Then we can use the constructor `operator` to create Hamiltonian:

```julia
julia> H = operator(mats, inds, L)
Operator of size (1024, 1024) with 10 terms.
```

The constructor return an `Operator` object, which is a linear operator that can act on vector/ matrix. For example, we can act `H` on the ferromagnetic state:

```julia
julia> ψ = zeros(2^L); ψ[1] = 1; H * random_state
1024-element Vector{Float64}:
 -1.5539463277491536
  5.969061189628827
  3.439873269795492
  1.6217619009059376
  0.6101231697221667
  6.663735992405236
  ⋮
  5.517409105968883
  0.9498121684380652
 -0.0004996659995972763
  2.6020967735388734
  4.99027405325114
 -1.4831032210847952
```

If we need a matrix representation of the Hamitonian, we can convert `H` to julia array by:

```julia
julia> Array(H)
1024×1024 Matrix{Float64}:
 -1.55617  0.0       0.0       0.0     …   0.0      0.0       0.0
  0.0      4.18381   2.0       0.0         0.0      0.0       0.0
  0.0      2.0      -1.42438   0.0         0.0      0.0       0.0
  0.0      0.0       0.0      -1.5901      0.0      0.0       0.0
  0.0      0.0       2.0       0.0         0.0      0.0       0.0
  0.0      0.0       0.0       2.0     …   0.0      0.0       0.0
  ⋮                                    ⋱                     
  0.0      0.0       0.0       0.0         0.0      0.0       0.0
  0.0      0.0       0.0       0.0         2.0      0.0       0.0
  0.0      0.0       0.0       0.0     …   0.0      0.0       0.0
  0.0      0.0       0.0       0.0        -1.42438  2.0       0.0
  0.0      0.0       0.0       0.0         2.0      4.18381   0.0
  0.0      0.0       0.0       0.0         0.0      0.0      -1.55617
```

Or use the function `sparse` to create the sparse matrix (requires the module `SparseArrays` being imported):

```julia
julia> sparse(H)
1024×1024 SparseMatrixCSC{Float64, Int64} with 6144 stored entries:
⠻⣦⣄⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠳⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⢹⡻⣮⡳⠄⢠⡀⠀⠀⠀⠀⠀⠀⠈⠳⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠙⠎⢿⣷⡀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠈⠳⣄⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠲⣄⠈⠻⣦⣄⠙⠀⠀⠀⢦⡀⠀⠀⠀⠈⠳⣄⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠈⠳⣄⠙⡻⣮⡳⡄⠀⠀⠙⢦⡀⠀⠀⠀⠈⠳⣄⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠙⠮⢻⣶⡄⠀⠀⠀⠙⢦⡀⠀⠀⠀⠈⠳⣄⠀
⢤⡀⠀⠀⠀⠀⠠⣄⠀⠀⠀⠉⠛⣤⣀⠀⠀⠀⠙⠂⠀⠀⠀⠀⠈⠓
⠀⠙⢦⡀⠀⠀⠀⠈⠳⣄⠀⠀⠀⠘⠿⣧⡲⣄⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠙⢦⡀⠀⠀⠀⠈⠳⣄⠀⠀⠘⢮⡻⣮⣄⠙⢦⡀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠈⠳⠀⠀⠀⣄⠙⠻⣦⡀⠙⠦⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠈⠳⣄⠈⢿⣷⡰⣄⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠈⠃⠐⢮⡻⣮⣇⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠙⠻⣦
```

### Translational-invariant Systems

Consider the AKLT model `H = ∑ᵢ S⃗ᵢ⋅S⃗ᵢ₊₁ + 1/3 ∑ᵢ (S⃗ᵢ⋅S⃗ᵢ₊₁)²`, with system size chosen to be `L=8`. The Hamiltonian operator for this translational-invariant Hamiltonian can be constructed using the `trans_inv_operator` function:

```julia
L = 8
SS = spin((1, "xx"), (1, "yy"), (1, "zz"), D=3)
mat = SS + 1/3 * SS^2
H = trans_inv_operator(mat, 1:2, L)
```

The second input specifies the indices the operators act on.

Because of the translational symmetry, we can simplify the problem by considering the symmetry. We construct a translational-symmetric basis by:

```julia
B = TranslationalBasis(0, 8, base=3)
```

Here the first argument labels the momentum `k = 0,...,L-1`, the second argument is the length of the system. The function `TranslationalBasis` return a basis object containing 834 states. We can obtain the Hamiltonian in this sector by:

```julia
julia> H = trans_inv_operator(mat, 1:2, B)
Operator of size (834, 834) with 8 terms.
```

In addition, we can take into account the total `Sz` conservation, by constructing the basis

```julia
B = TranslationalBasis(x -> sum(x) == 8, 0, 8, base=3)
```

where the first argument is the selection function. The function `(x -> sum(x) == 8)` means we select those states whose total `Sz` equalls 0 (note that we use 0,1,2 to label the `Sz=1,0,-1` states). This gives a further reduced Hamiltonian matrix:

```julia
julia> H = trans_inv_operator(mat, 1:2, B)
Operator of size (142, 142) with 8 terms.
```

We can go on step further by considering the spatial reflection symmetry.

```julia
B = TranslationParityBasis(x -> sum(x) == 8, 0, 1, L, base=3)
```

where the second argument is the momentum, the third argument is the parity `p = ±1`.

```julia
julia> H = trans_inv_operator(mat, 1:2, B)
Operator of size (84, 84) with 8 terms.
```



## Basis Object

In `EDKit.jl`, the fundamental objects are basis and operator. The `AbstractBasis` is the abstract type of basis. Currently there are 4 concrete basis:

1. `TensorBasis`: Ordinary basis without any symmetry.
2. `ProjectedBasis`: Basis for subspace that is spanned only by product states.
3. `TranslationalBasis`: Basis for translational symmetric Hamiltonian.
4. `TranslationalParityBasis` : Basis with translational as well as reflection symmetry. The momensum should be 0 or π.

The basis object can be extended. To construct linear operation, we need to define 3 functions for a new basis type:

1. `size(b::AbstractBasis)` : Size of matrix representations of the operators in this subspace.
2. `change!(b::AbstractBasis, i::Integer)` : Return the normalization of ith state and change the digits to ith states in this subspace.
3. `index(b::AbstractBasis)` : Return the coefficient and index of the digits.

Optionally, we can define `eltype` for a basis object (default is `ComplexF64`).

If the calculation is done on the entire Hilbert space, the basis object need not be explicitly constructed. The `Operator` will use `TensorBasis` by default. The construction of other basis with symmetry concern are discussed below.

In addition, if the entaglement entropy is needed, the user-defined basis should implement a function `schmidt!(target, v, Ainds, b::AbstractBasis)`.



## Operator Object

In `EDKit.jl` , a many-body operator is represented by the type `Operator`:

```julia
struct Operator{Tv<:Number, Tb<:AbstractBasis}
    M::Vector{SparseMatrixCSC{Tv, Int}}
    I::Vector{Vector{Int}}
    B::Tb
end
```

In this definition, `M` is the list of matrix representations of local operators, `I` is the list of indices of sites it acts on.

### Construction

To construct an `Operator` object, we need 3 inputs:

```julia
M <: AbstractVector{<:AbstractMatrix}
I <: AbstractVector{<:AbstractVector{<:Integer}}
B <: AbstractBasis
```

where `M` is the list of matrices representing the local operators, `I` is the list of vectors representing the sites it acts on. `B` is a basis object. If use `TensorBasis`, we can directly using the constructing method

```julia
operator(mats::AbstractVector{<:AbstractMatrix}, inds::AbstractVector{<:AbstractVector}, L::Integer)
```

### Convert to matrix

An `Operator` object is basically like a matrix, and it can be convert to dense matrix using the function

```julia
Array(opt::Operation)
```

It can also write to a given matrix with correct dimension using function

```julia
addto!(M::AbstractMatrix, opt::Operator)
```

Note that to get correct answer, `M` should de initialized as a zero matrix.

### Multiply to vector or matrix

We can directly using

```julia
O::Operator * M::AbstractVecOrMat
```

to do the multiplycation. Or, use the function

```julia
mul!(target::AbstractVecOrMat, opt::Operator, v::AbstractVecOrMat)
```

to modify `target` (similarly, `target` should be initialized as a zero vector/matrix).

### Compute entaglement entropy

After obtaining Hamiltonian in a symmetry sector. We can calculate the entaglement entropy of an eigenvector `v` (assume the system size is `2L`, and the entropy cut is at the middel of the chain) by

```julia
ent_S(v::AbstractVector, 1:L, b::AbstractBasis)
```



## Concrete Implementations of Basis

Here we introduce 3 concrete implementation of `AbstractBasis`.

### TensorBasis

The type `TensorBasis` has the fields:

```julia
struct TensorBasis <: AbstractBasis
    dgt::Vector{Int}
    B::Int
end
```

We can use the function

```julia
tensorbasis(L::Integer; base::Integer=2)
```

to construct a basis, though in most cases it is not necessary.

### ProjectedBasis

The type `ProjectedBasis` has the fields:

```julia
struct ProjectedBasis <: AbstractBasis
    dgt::Vector{Int}
    I::Vector{Int}
    B::Int
end
```

We can use the function

```julia
projectedbasis(f, L::Integer; base::Integer=2)
```

In the definition, `f` is a function acting on digits that tells whether a given digits is valid in this basis. For eaxample, consider a S=1/2 chain with L=6 (conserve megnetic quantum number). If we consider the subspace spaned by those states with 3 up-spins, we can create the basis for the subspace by

```julia
projectedbasis(x -> sum(x)==3, 6; base=2)
```

### TranslationalBasis

The type `TranslationalBasis` has the fields:

```julia
struct TranslationalBasis <: AbstractBasis
    dgt::Vector{Int}
    I::Vector{Int}
    R::Vector{Float64}
    C::ComplexF64
    B::Int
end
```

We can use the function

```julia
translationalbasis(k::Integer, L::Integer; base::Integer=2)
translationalbasis(f, k::Integer, L::Integer; base::Integer=2)
```

In the definition, `k` is the momentum sector we are interested, and `f` is a function acting on digits that tells whether a given digits is valid in this basis. Note that `k` is an integer, representing the momentum at `2πk/L`. For example, consider a `S=1/2` chain with `L=6` (with translational symmetry and conserve magnetic quantum number). If we consider the subspace spaned by those states with 3 up-spins with zero momentum, we can create the basis for the subspace by

```julia
translationalbasis(x -> sum(x)==3, 0, 6; base=2)
```
