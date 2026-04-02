export AbstractBasis, content, base, index, change!
"""
    AbstractBasis

Abstract supertype for all basis objects used by EDKit.

Every concrete basis stores a mutable digit buffer `dgt` together with enough
metadata to move between three related descriptions of a many-body state:

- a concrete digit string such as `[0, 1, 0, 1]`,
- a representative index stored in the basis contents,
- and, for symmetry-reduced bases, the normalization or phase needed to recover
  the corresponding physical state.

The two most important operations on any basis are:

- [`index`](@ref): interpret the current `dgt` buffer as coordinates in the
  basis and return a coefficient/index pair,
- [`change!`](@ref): load the `i`th stored basis state into `dgt` and return
  its normalization.

This abstraction is what allows operator construction, matrix assembly,
entanglement calculations, and inter-basis maps to work with the same code
across full and symmetry-reduced Hilbert spaces.
"""
abstract type AbstractBasis end

"""
    AbstractOnsiteBasis

Abstract supertype for bases whose stored coordinates are ordinary product-state
coordinates without nontrivial orbit structure.

These bases behave like direct tensor-product Hilbert spaces:

- `content(b, i)` is just the integer encoding of the `i`th product state,
- `norm(b, i) == 1`,
- `order(b) == 1`,
- and `eltype(b)` can remain real/integer because no symmetry phase factors are
  introduced by the basis itself.

`TensorBasis` and `ProjectedBasis` are the main concrete subtypes. Many generic
algorithms use this supertype to detect the simpler "no orbit normalization"
path.
"""
abstract type AbstractOnsiteBasis <: AbstractBasis end

"""
    AbstractPermuteBasis

Abstract supertype for symmetry-reduced bases built from permutation orbits of
product states.

Concrete subtypes typically keep one canonical representative per orbit under a
group action such as translation, parity, or spin flip. Their `index` methods
therefore return both:

- the stored representative index,
- and the coefficient/phase relating the current `dgt` configuration to the
  normalized symmetry eigenstate.

This category is used throughout EDKit to dispatch onto algorithms that must
account for orbit lengths, phase conventions, or canonical representative
selection.
"""
abstract type AbstractPermuteBasis <: AbstractBasis end

"""
    AbstractTranslationalParityBasis

Abstract supertype for bases that simultaneously track translation structure
and a reflection-like symmetry sector.

This covers the parity-resolved momentum bases, such as
`TranslationParityBasis` and `TranslationFlipBasis`, where the representative
search and normalization rules are more intricate than in pure translation or
pure parity sectors. Generic Schmidt and basis-construction helpers dispatch on
this supertype when they need the specialized translational-parity logic.
"""
abstract type AbstractTranslationalParityBasis <: AbstractPermuteBasis end

#-------------------------------------------------------------------------------------------------------------------------
# Default function definitions
#-------------------------------------------------------------------------------------------------------------------------
"""
    content(b::AbstractBasis, i::Integer)

Return the stored representative label for the `i`th basis state.

For most symmetry-reduced bases this is not the full state vector itself, but an
integer encoding of the canonical representative digit string associated with
that basis element. For `TensorBasis`, this is simply the ordinary tensor-basis
index.
"""
content(b::AbstractBasis, i::Integer) = b.I[i]
#-------------------------------------------------------------------------------------------------------------------------
"""
    norm(b::AbstractBasis, i::Integer)

Return the normalization factor associated with the `i`th basis state.

For onsite bases this is `1`. For symmetry-reduced bases this typically stores
the orbit-size-dependent normalization needed to convert between canonical
representatives and normalized symmetry eigenstates.
"""
norm(b::AbstractBasis, i::Integer) = b.R[i]
norm(::AbstractOnsiteBasis, ::Integer) = 1
#-------------------------------------------------------------------------------------------------------------------------
"""
    eltype(b::AbstractBasis)

Return the coefficient type naturally associated with the basis.

For onsite bases this defaults to `Int64` because no symmetry phases appear. For
most reduced bases the default is `ComplexF64`, reflecting that basis vectors
may carry complex momentum or symmetry phases.

This is used heavily in type promotion for operators, matrix assembly, and
state-space maps.
"""
eltype(::AbstractBasis) = ComplexF64
eltype(::AbstractOnsiteBasis) = Int64
#-------------------------------------------------------------------------------------------------------------------------
"""
    length(b::AbstractBasis)

Return the number of sites encoded by the basis digit buffer `b.dgt`.

This is the physical system size `L`, not the Hilbert-space dimension.
"""
length(b::AbstractBasis) = length(b.dgt)
#-------------------------------------------------------------------------------------------------------------------------
"""
    size(b::AbstractBasis[, dim])

Return the matrix shape induced by the basis.

For ordinary state bases this is the Hilbert-space dimension on both axes. Some
derived bases, such as `DoubleBasis`, override this behavior to describe
rectangular maps.
"""
function size(b::AbstractBasis, dim::Integer)
    isone(dim) && return length(b.I)
    isequal(dim, 2) ? length(b.I) : 1
end
size(b::AbstractBasis) = (size(b, 1), size(b, 2))
#-------------------------------------------------------------------------------------------------------------------------
"""
    change!(b::AbstractBasis, i::Integer)

Load the `i`th basis state into `b.dgt` and return its normalization factor.

This is the inverse companion to [`index`](@ref) at the basis-object level:
`change!` goes from stored basis coordinates to the current working digit
representation, whereas `index` interprets the current digits in basis
coordinates.

The method mutates `b.dgt`.
"""
function change!(b::AbstractBasis, i::Integer) 
    change!(b.dgt, content(b, i), base=b.B)
    norm(b, i)
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    int_type(b::AbstractBasis)

Return the integer type used to store representative indices for the basis.

This matters for large Hilbert spaces where the canonical representative labels
may exceed the safe range of default machine integers.
"""
int_type(b::AbstractBasis) = eltype(b.I)
#-------------------------------------------------------------------------------------------------------------------------
"""
    order(b::AbstractOnsiteBasis)

Return the orbit order associated with the basis.

For onsite bases there is no nontrivial symmetry orbit, so the order is `1`.
Reduced bases override this when they identify multiple product states into one
symmetry-adapted basis vector.
"""
order(b::AbstractOnsiteBasis) = 1
#-------------------------------------------------------------------------------------------------------------------------
# TensorBasis
#-------------------------------------------------------------------------------------------------------------------------
export TensorBasis
"""
    TensorBasis

Full tensor-product basis with no symmetry reduction.

The digit buffer `dgt` stores one local state label per site in base `B`. The
corresponding global index is computed by:

    I(i₁, i₂, ⋯, iₙ) = i₁ B^(n-1) + i₂ B^(n-2) + ⋯ + iₙ

`TensorBasis` is the reference basis against which more structured bases can be
understood. Many EDKit workflows do not need to build it explicitly, but it is
the conceptual full-space embedding used by maps, projections, and matrix
assembly.
"""
struct TensorBasis <: AbstractOnsiteBasis
    dgt::Vector{Int64}
    B::Int64
end

"""
    TensorBasis(; L, base=2)

Construct the full tensor-product basis on `L` sites with local base `base`.

Arguments:
- `L`: number of sites.
- `base`: local Hilbert-space dimension, for example `2` for spin-1/2 systems.

Returns:
- A fresh [`TensorBasis`](@ref) with a zero-initialized working digit buffer.
"""
function TensorBasis(;L::Integer, base::Integer=2)
    dgt = zeros(Int64, L)
    B = Int64(base)
    TensorBasis(dgt, B)
end
#-------------------------------------------------------------------------------------------------------------------------
content(::TensorBasis, i::Integer) = i
norm(b::TensorBasis, ::Integer) = 1
eltype(::TensorBasis) = Int64
size(b::TensorBasis, i::Integer) = isone(i) || isequal(i, 2) ? b.B^length(b.dgt) : 1
size(b::TensorBasis) = (l=b.B^length(b.dgt); (l, l))
index(b::TensorBasis) = 1, index(b.dgt, base=b.B)
int_type(::TensorBasis) = Int64
#-------------------------------------------------------------------------------------------------------------------------
"""
    copy(b::TensorBasis)

Return a copy of `b` with an independent digit buffer.

This matters because many algorithms mutate `dgt` while iterating over basis
states.
"""
function copy(b::TensorBasis) 
    TensorBasis(deepcopy(b.dgt), b.B)
end

#-------------------------------------------------------------------------------------------------------------------------
# Index functions
#-------------------------------------------------------------------------------------------------------------------------
"""
    index(dgt::AbstractVector{T}; base::Integer=2, dtype::DataType=T) where T <: Integer

Convert a full digit string into EDKit's 1-based integer index.

The mapping is:

    number = ∑ᵢ bits[i] * base^(L-i) + 1

Arguments:
- `dgt`: digit representation of a product state.
- `base`: local base used to interpret the digits.

Returns:
- The 1-based linear index corresponding to `dgt`.

Notes:
- This is the low-level indexing convention used throughout the package.
- The computation uses iterative polynomial evaluation rather than explicit
  powers for speed.
"""
@inline function index(dgt::AbstractVector{<:Integer}; base::T=2) where T <: Integer
    N = zero(T)
    if base == 2
        @inbounds for i = 1:length(dgt)
            N = (N << 1) | dgt[i]
        end
    else
        @inbounds for i = 1:length(dgt)
            N *= base
            N += dgt[i]
        end
    end
    N + one(T)
end

#-------------------------------------------------------------------------------------------------------------------------
"""
    index(dgt::AbstractVector{T}, sites::AbstractVector{<:Integer}; base::Integer=2) where T <: Integer

Convert a subset of sites from a digit buffer into a 1-based local index.

Arguments:
- `dgt`: full digit buffer.
- `sites`: site positions to read, in the order they should appear in the local
  tensor product.
- `base`: local base used for interpretation.

Returns:
- The 1-based index of the subsystem specified by `sites`.

This is a central helper in local operator application, where a sparse local
matrix acts only on a selected subset of sites.
"""
@inline function index(dgt::AbstractVector{T}, sites::AbstractVector{<:Integer}; base::T=2) where T <: Integer
    N = zero(T)
    if base == 2
        @inbounds for i in sites
            N = (N << 1) | dgt[i]
        end
    else
        @inbounds for i in sites
            N *= base
            N += dgt[i]
        end
    end
    N + one(T)
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    change!(dgt::AbstractVector{<:Integer}, ind::Integer; base::Integer=2)

Overwrite `dgt` with the product-state digits corresponding to `ind`.

Arguments:
- `dgt`: mutable digit buffer to overwrite.
- `ind`: 1-based product-state index.
- `base`: local base used for the decoding.

Returns:
- The mutated digit buffer `dgt`.

This is the inverse of `index(dgt; base=...)`.
"""
@inline function change!(dgt::AbstractVector{T}, ind::T; base::T=2) where T
    N = ind - one(T)
    if base == 2
        @inbounds for i = length(dgt):-1:1
            dgt[i] = N & one(T)
            N >>= 1
        end
    else
        @inbounds for i = length(dgt):-1:1
            N, dgt[i] = divrem(N, base)
        end
    end
end
#-------------------------------------------------------------------------------------------------------------------------
"""
    change!(dgt::AbstractVector{<:Integer}, sites::AbstractVector{<:Integer}, ind::Integer; base::Integer=2)

Overwrite only the selected entries of `dgt` with the subsystem digits encoded
by `ind`.

Arguments:
- `dgt`: mutable full-system digit buffer.
- `sites`: site positions to overwrite.
- `ind`: 1-based subsystem index.
- `base`: local base used for decoding.

Returns:
- The mutated digit buffer `dgt`.

This is the inverse of `index(dgt, sites; base=...)` and is heavily used
inside matrix-free operator application.
"""
@inline function change!(dgt::AbstractVector{T}, sites::AbstractVector{<:Integer}, ind::Integer; base::T=2) where T
    N = ind - one(T)
    if base == 2
        @inbounds for i = length(sites):-1:1
            dgt[sites[i]] = N & one(T)
            N >>= 1
        end
    else
        @inbounds for i = length(sites):-1:1
            N, dgt[sites[i]] = divrem(N, base)
        end
    end
end

#-------------------------------------------------------------------------------------------------------------------------
# Integer-state helpers for orbit search (avoid circshift!/reverse! on digit arrays)
#-------------------------------------------------------------------------------------------------------------------------
"""
    _int_translate(state, pow_shift, base_shift)

Integer equivalent of `circshift!(dgt, A)` followed by `index(dgt)`.

Operates on a 0-based integer state.  `pow_shift = base^(L-A)` and
`base_shift = base^A` must be precomputed by the caller.
"""
@inline function _int_translate(state::T, pow_shift::T, base_shift::T) where T <: Integer
    tail = state % base_shift
    tail * pow_shift + state ÷ base_shift
end

"""
    _int_reverse(state, L, base)

Integer equivalent of reversing the digit order (parity/reflection).

Operates on a 0-based integer state.
"""
@inline function _int_reverse(state::T, L::Int, base::T) where T <: Integer
    r = zero(T)
    v = state
    @inbounds for _ in 1:L
        r = r * base + v % base
        v ÷= base
    end
    r
end

"""
    _int_spinflip(state, maxstate)

Integer equivalent of spin-flip (complement each digit: d → base-1-d).

Operates on a 0-based integer state.  `maxstate = base^L - 1`.
"""
@inline _int_spinflip(state::T, maxstate::T) where T <: Integer = maxstate - state
