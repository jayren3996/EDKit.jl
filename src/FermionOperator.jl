export fermion, fermion_operator, trans_inv_fermion_operator, jw_string_required, hubbard

const _JW_BEARING_OPS  = ("+", "-", "+-", "-+", "++", "--")
const _JW_FREE_OPS     = ("n", "z", "I", "nn")

"""
    jw_string_required(op::AbstractString) -> Bool

Whether the fermion operator string `op` carries a Jordan-Wigner string when
embedded into a many-body basis.

JW-string-bearing operators (`"+"`, `"-"`, `"+-"`, `"-+"`, `"++"`, `"--"`) do
not commute with spatial permutation symmetries. The diagonal operators
`"n"`, `"z"`, `"I"`, and `"nn"` carry no JW string and remain safe on any
onsite basis. Raises an error for unrecognised op strings rather than
silently defaulting either way.
"""
function jw_string_required(op::AbstractString)
    op in _JW_BEARING_OPS && return true
    op in _JW_FREE_OPS    && return false
    error("Unknown fermion operator string \"$op\". Recognised: " *
          "$(_JW_BEARING_OPS) (JW-bearing) and $(_JW_FREE_OPS) (no JW).")
end

"""
    fermion(op::AbstractString[, span::Integer]; convention::Symbol=:left)

Return the local matrix representation of a spinless-fermion operator on a
contiguous block of `span` sites.

Single-site operators (`span = 1`):
- `"n"` — number operator `n = c†c` (`diag(0, 1)`)
- `"z"` — `n - 1/2` (`diag(-1/2, 1/2)`)
- `"I"` — identity
- `"+"` — bare creation `c†` (no JW string; embed via [`fermion_operator`](@ref))
- `"-"` — bare annihilation `c` (no JW string; embed via [`fermion_operator`](@ref))

Two-character operators (`span ≥ 2`):
- `"+-"` — hop `c†_1 c_span`
- `"-+"` — `c_1 c†_span`
- `"++"` — pair creation `c†_1 c†_span`
- `"--"` — pair annihilation `c_1 c_span`
- `"nn"` — density-density `n_1 n_span` (no JW string; both endpoints diagonal)

For the four `±` two-character operators, the Jordan-Wigner `σ_z` string is
inserted on sites `2, 3, …, span-1` whenever `span > 2`. An overall ±1 sign
is folded in to account for the on-site `σ⁺σ_z = -σ⁺` / `σ⁻σ_z = +σ⁻`
reduction at site 1.

`convention` selects the Jordan-Wigner direction. Only `:left` is implemented;
the kwarg is plumbed through `fermion`, `fermion_operator`, and
`trans_inv_fermion_operator` for forward compatibility with a future
`:right` implementation.

Embed the returned matrix into a many-body basis with [`operator`](@ref) or
[`fermion_operator`](@ref).
"""
function fermion(op::AbstractString, span::Integer=length(op); convention::Symbol=:left)
    convention === :left || error("Only the :left Jordan-Wigner convention is implemented (got $convention).")

    # Single-site diagonal / identity operators (no JW string).
    if op == "n"
        span == 1 || error("\"n\" is a single-site operator (got span=$span).")
        return [0.0 0.0; 0.0 1.0]
    elseif op == "z"
        span == 1 || error("\"z\" is a single-site operator (got span=$span).")
        return [-0.5 0.0; 0.0 0.5]
    elseif op == "I"
        span == 1 || error("\"I\" is a single-site operator (got span=$span).")
        return Matrix{Float64}(I, 2, 2)
    end

    if length(op) == 1
        c = op[1]
        span == 1 || error("Single-character operator \"$op\" is a one-site operator (got span=$span).")
        return c == '+' ? Array(spin_Sm(2)) :
               c == '-' ? Array(spin_Sp(2)) :
               error("Unsupported single-site fermion operator: \"$op\". Supported: \"n\", \"z\", \"I\", \"+\", \"-\".")
    end

    length(op) == 2 || error("Only 1- or 2-character fermion operator strings are supported in this MVP (got \"$op\").")
    span >= 2 || error("Two-fermion operators require span >= 2 (got $span).")

    # Density-density "nn" — diagonal at both endpoints, no JW string.
    # Use sparse arrays so the intermediate is O(2^(span-2)) nnz instead of a
    # dense 2^span × 2^span matrix.
    if op == "nn"
        n_sp  = sparse([0.0 0.0; 0.0 1.0])
        Id_sp = sparse(1.0 * I, 2, 2)
        mat = n_sp
        for k in 2:span-1
            mat = kron(mat, Id_sp)
        end
        mat = kron(mat, n_sp)
        return mat
    end

    # Two-character ±±/+-/-+ ops. The JW string is a diagonal of ±1 and each
    # endpoint is the rank-2 σ⁻ / σ⁺, so the full local matrix has only
    # 2^(span-1) nonzeros. Building it sparse avoids the 2^span × 2^span
    # dense intermediate that would OOM for span ≳ 20.
    σ⁺_sp = sparse(Array(spin_Sp(2)))
    σ⁻_sp = sparse(Array(spin_Sm(2)))
    σ_z_sp = sparse([1.0 0.0; 0.0 -1.0])  # spin_Sz gives Sz = σ_z/2; we want σ_z = diag(1, -1)

    site1 = op[1] == '+' ? σ⁻_sp : op[1] == '-' ? σ⁺_sp :
            error("Unsupported operator char: '$(op[1])' in \"$op\".")
    site2 = op[2] == '+' ? σ⁻_sp : op[2] == '-' ? σ⁺_sp :
            error("Unsupported operator char: '$(op[2])' in \"$op\".")

    # Sign from the on-site product at site 1:
    #   σ⁻ σ_z = +σ⁻  → "+? " → +1
    #   σ⁺ σ_z = −σ⁺  → "-? " → −1
    sign = op[1] == '+' ? 1.0 : -1.0

    mat = site1
    for k in 2:span-1
        mat = kron(mat, σ_z_sp)
    end
    mat = kron(mat, site2)
    sign * mat
end

"""
    fermion_operator(op::AbstractString, sites::AbstractVector{<:Integer},
                     B::AbstractBasis; convention::Symbol=:left)

Build the full-system [`Operator`](@ref) corresponding to a fermion operator
`op` placed on the supplied `sites`.

Supported `op` strings:
- `"n"`, `"z"`, `"I"` — single-site diagonal operators (no JW string).
  Require `length(sites) == 1`.
- `"+"`, `"-"` — single-site `c†_i` / `c_i` embedded with the full
  left-string Jordan-Wigner prefix `σ_z^1 ⋯ σ_z^{i-1}` on sites `1:i`.
  Require `length(sites) == 1`.
- `"+-"`, `"-+"`, `"++"`, `"--"` — two-site fermion products. Require
  `length(sites) == 2`. The operator's support is the contiguous block
  `min(sites):max(sites)`, with the JW `σ_z` string baked in on intermediate
  sites by [`fermion`](@ref).
- `"nn"` — two-site density-density `n_i n_j`. Require `length(sites) == 2`,
  `sites[1] != sites[2]`. Order-insensitive (`n_i` and `n_j` commute).

If `sites` is unsorted (e.g. `[3, 1]` for `c†_3 c_1`), the operator is
adjusted by the appropriate fermionic sign so that the returned matrix
represents the operator literally written by the caller.
"""
function fermion_operator(op::AbstractString, sites::AbstractVector{<:Integer},
        B::AbstractBasis; convention::Symbol=:left)
    convention === :left || error("Only the :left convention is implemented.")

    # Spinless fermion ops are defined on a 2-state local Hilbert space
    # (empty / occupied). A higher-base basis would silently produce a
    # dimension mismatch in `operator(...)` further down.
    B.B == 2 || error("fermion_operator requires a base=2 basis (got base=$(B.B)).")

    # Validate site indices against the basis size — internal kernels use
    # @inbounds and can silently miscompute (or corrupt the dgt buffer) on
    # out-of-range entries.
    L = length(B.dgt)
    isempty(sites) && error("`sites` must be non-empty (got $sites).")
    all(s -> 1 <= s <= L, sites) ||
        error("Each entry of `sites` must lie in 1:$L for the supplied basis (got $sites).")

    # JW-bearing operators (c†, c, c†c†, cc, c†c, cc†) do not commute with
    # spatial permutation symmetries. Embedding them on an AbstractPermuteBasis
    # would silently produce wrong eigenvalues. Diagonal/identity operators
    # carry no JW string and remain valid on any onsite permutation — gated
    # via `jw_string_required` so new diagonal ops are safe by default.
    if B isa AbstractPermuteBasis && !(B isa TranslationalFermionBasis) && jw_string_required(op)
        error("fermion_operator(\"$op\", $(sites), ::$(typeof(B))) is not supported: " *
              "the Jordan-Wigner string for c†/c does not commute with the symmetry of " *
              "$(typeof(B)). These AbstractPermuteBasis subtypes are rejected — " *
              "TranslationalBasis, ParityBasis, FlipBasis, ParityFlipBasis, " *
              "TranslationParityBasis, TranslationFlipBasis, and AbelianBasis. Use " *
              "SpinlessFermionBasis (with N=… or nf=…) for the full-space sector, or " *
              "TranslationalFermionBasis (momentum-resolved spinless fermions) with " *
              "trans_inv_fermion_operator. See the \"Symmetry caveats\" section of the " *
              "spinless fermions manual.")
    end

    # Single-site diagonal / identity operators — no JW string.
    if op == "n" || op == "z" || op == "I"
        length(sites) == 1 || error("\"$op\" takes exactly one site (got $(length(sites))).")
        return operator(fermion(op), [sites[1]], B)
    end

    # Single-site c† or c — embed with explicit left-string JW prefix.
    if length(op) == 1
        c = op[1]
        c == '+' || c == '-' || error("Unsupported single-site operator \"$op\". Supported: \"n\", \"z\", \"I\", \"+\", \"-\".")
        length(sites) == 1 || error("\"$op\" takes exactly one site (got $(length(sites))).")
        i = sites[1]
        if i == 1
            return operator(fermion(op), [1], B)
        end
        # Sparse JW prefix: i−1 σ_z factors (each diagonal) kron'd with the
        # endpoint σ⁻ / σ⁺. The full local matrix has 2^(i-1) nonzeros, so
        # a sparse build is O(2^(i-1)) memory vs the 2^i × 2^i dense
        # intermediate (which OOMs for i ≳ 20).
        σ_z_sp   = sparse([1.0 0.0; 0.0 -1.0])
        endpoint = sparse(fermion(op))
        mat = σ_z_sp
        for k in 2:i-1
            mat = kron(mat, σ_z_sp)
        end
        mat = kron(mat, endpoint)
        return operator(mat, collect(1:i), B)
    end

    length(op) == 2 || error("Only 1- or 2-character fermion operator strings are supported in this MVP (got \"$op\").")
    length(sites) == 2 || error("Two-fermion operators take exactly two sites.")
    i, j = sites[1], sites[2]

    # Density-density "nn" — n_i and n_j commute and are diagonal, so order
    # of `sites` is irrelevant and no fermion sign is incurred.
    if op == "nn"
        i == j && error("\"nn\" on the same site reduces to \"n\" (n^2 = n). Use fermion_operator(\"n\", [i], B) instead.")
        lo, hi = minmax(i, j)
        span = hi - lo + 1
        return operator(fermion("nn", span), collect(lo:hi), B)
    end

    if i == j
        if op == "+-"
            error("\"+-\" on the same site is the number operator: c†_$i c_$i = n_$i. " *
                  "Use fermion_operator(\"n\", [$i], B) instead.")
        elseif op == "-+"
            error("\"-+\" on the same site equals I − n: c_$i c†_$i = I − n_$i. " *
                  "Construct as I − fermion_operator(\"n\", [$i], B).")
        else
            # "++" and "--" only — both vanish by Pauli exclusion.
            error("\"$op\" on the same site is identically zero by Pauli exclusion (c†² = c² = 0).")
        end
    end
    swapped = i > j
    lo, hi = minmax(i, j)
    span = hi - lo + 1

    local_op = if swapped
        # Swapped-site reductions (lo < hi, sites[1] > sites[2]):
        #  "+-" with i>j: c†_hi c_lo = (c†_lo c_hi)†          (Hermitian conjugate, no sign)
        #  "-+" with i>j: c_hi c†_lo = (c_lo c†_hi)†          (Hermitian conjugate, no sign)
        #  "++" with i>j: c†_hi c†_lo = -c†_lo c†_hi          (fermionic anticommutation)
        #  "--" with i>j: c_hi c_lo  = -c_lo c_hi             (fermionic anticommutation)
        if op == "+-" || op == "-+"
            fermion(op, span)'
        elseif op == "++" || op == "--"
            -fermion(op, span)
        else
            error("Unsupported swap-mode fermion operator: \"$op\"")
        end
    else
        fermion(op, span)
    end
    operator(local_op, collect(lo:hi), B)
end

"""
    trans_inv_fermion_operator(op, support, B; convention=:left)
    trans_inv_fermion_operator(op, span::Integer, B; convention=:left)

Build a translationally-invariant fermion operator by summing
`fermion_operator(op, mod1.(support .+ t, L), B)` over `t = 0, 1, …, L-1`.

Unlike [`trans_inv_operator`](@ref) — which duplicates the same local matrix
on every translation — this helper rebuilds the Jordan-Wigner string for each
translated bond, so the wrap-around term carries the correct long-way JW
chain through sites `2, …, L-1`. Without that chain the boundary bond is
silently wrong at particle number `N ≥ 2` (the single-particle sector is
accidentally correct).

Accepts any [`AbstractOnsiteBasis`](@ref) — [`SpinlessFermionBasis`](@ref),
[`TensorBasis`](@ref), and `base=2` [`ProjectedBasis`](@ref) all work, since
the JW string is purely an on-site permutation in those representations.
Symmetry-reduced fermion bases (`AbstractPermuteBasis`) are not yet
supported; see the "Symmetry caveats" section of the spinless fermions
manual.

# Example: tight-binding ring with PBC
```julia
L = 6
B = SpinlessFermionBasis(L = L, N = L ÷ 2)
H_hop = trans_inv_fermion_operator("+-", [1, 2], B)
H = -(H_hop + adjoint(H_hop))      # = -Σ_i (c†_i c_{i+1} + h.c.)
```

!!! warning "L = 2 ring is a special case"
    On a 2-site ring there is only one bond, but `trans_inv_fermion_operator`
    still sums over `t = 0, 1`, producing `c†_1 c_2 + c†_2 c_1` (already
    Hermitian). Applying the `H_hop + adjoint(H_hop)` recipe above then
    double-counts. For `L = 2`, build the Hamiltonian directly as
    `H = -H_hop`.
"""
function trans_inv_fermion_operator(op::AbstractString,
        support::AbstractVector{<:Integer},
        B::Union{AbstractOnsiteBasis, TranslationalFermionBasis};
        convention::Symbol=:left)
    L = length(B.dgt)
    isempty(support) && error("`support` must be non-empty (got $support).")
    all(s -> 1 <= s <= L, support) ||
        error("Each entry of `support` must lie in 1:$L (got $support). " *
              "Out-of-range indices would be silently wrapped by mod1 into a wrong operator.")
    # Seed `total` with the t=0 term so its inferred type is `Operator` rather
    # than `Union{Nothing, Operator}` — keeps the return type stable for
    # downstream Array / eigen / sparse! calls.
    total = fermion_operator(op, mod1.(support, L), B; convention=convention)
    for t in 1:L-1
        sites = mod1.(support .+ t, L)
        total = total + fermion_operator(op, sites, B; convention=convention)
    end
    total
end

trans_inv_fermion_operator(op::AbstractString, span::Integer,
        B::Union{AbstractOnsiteBasis, TranslationalFermionBasis}; kwargs...) =
    trans_inv_fermion_operator(op, collect(1:span), B; kwargs...)

"""
    fermion_operator(op::AbstractString, site::Integer, B::AbstractBasis; kwargs...)

Single-site shorthand: equivalent to `fermion_operator(op, [site], B; kwargs...)`.
"""
fermion_operator(op::AbstractString, site::Integer, B::AbstractBasis; kwargs...) =
    fermion_operator(op, [site], B; kwargs...)

#-----------------------------------------------------------------------------------------------------
# Spinful / multi-species fermion operators
#-----------------------------------------------------------------------------------------------------
"""
    fermion_operator(op, sitespins::AbstractVector{<:Tuple}, B::SpinfulFermionBasis; convention=:left)
    fermion_operator(op, sitespin::Tuple, B::SpinfulFermionBasis; convention=:left)

Spin-aware embedding for a [`SpinfulFermionBasis`](@ref): each `(site, spin)`
tuple is mapped to its mode index via [`fermionmode`](@ref) and the call is
forwarded to the mode-indexed [`fermion_operator`](@ref). `spin` accepts an
`Integer` species or, for spin-½, the symbols `:↑`/`:↓` (`:up`/`:down`).

# Example
```julia
B = SpinfulFermionBasis(L=4, S=2, N=(2, 2))
hop = fermion_operator("+-", [(1, :↑), (2, :↑)], B)   # c†_{1↑} c_{2↑}
```
"""
function fermion_operator(op::AbstractString, sitespins::AbstractVector{<:Tuple},
        B::SpinfulFermionBasis; convention::Symbol=:left)
    modes = [fermionmode(B, sp[1], sp[2]) for sp in sitespins]
    fermion_operator(op, modes, B; convention=convention)
end

fermion_operator(op::AbstractString, sitespin::Tuple, B::SpinfulFermionBasis; kwargs...) =
    fermion_operator(op, [sitespin], B; kwargs...)

"""
    hubbard(B::SpinfulFermionBasis; t=1.0, U=0.0, μ=0.0, boundary=:periodic)

Build the single-band Fermi–Hubbard Hamiltonian on a spin-½
[`SpinfulFermionBasis`](@ref) (`S = 2`):

    H = −t Σ_{⟨ij⟩,σ} (c†_{iσ} c_{jσ} + h.c.)  +  U Σ_i n_{i↑} n_{i↓}  −  μ Σ_{iσ} n_{iσ}

Keyword arguments:
- `t`        : nearest-neighbour hopping amplitude.
- `U`        : on-site interaction.
- `μ`        : chemical potential.
- `boundary` : `:periodic` (ring) or `:open` (chain).

On a 2-site `:periodic` lattice the single bond is counted once (the generic
ring would otherwise traverse it twice). Returns an [`Operator`](@ref).
"""
function hubbard(B::SpinfulFermionBasis; t::Real=1.0, U::Real=0.0, μ::Real=0.0,
        boundary::Symbol=:periodic)
    B.S == 2 || error("hubbard requires a spin-½ basis (S=2); got S=$(B.S). Build multi-species Hamiltonians term by term with fermion_operator.")
    boundary in (:periodic, :open) || error("boundary must be :periodic or :open (got $boundary).")
    L = B.L

    bonds = if boundary === :open
        [(i, i + 1) for i in 1:L-1]
    elseif L == 2
        [(1, 2)]                                  # ring with one distinct bond
    else
        [(i, mod1(i + 1, L)) for i in 1:L]
    end

    total = nothing
    for (i, j) in bonds, σ in (:↑, :↓)
        hop = fermion_operator("+-", [(i, σ), (j, σ)], B)   # c†_{iσ} c_{jσ}
        total = total + (-t) * (hop + adjoint(hop))
    end
    if !iszero(U)
        for i in 1:L
            total = total + U * fermion_operator("nn", [(i, :↑), (i, :↓)], B)
        end
    end
    if !iszero(μ)
        for i in 1:L, σ in (:↑, :↓)
            total = total + (-μ) * fermion_operator("n", [(i, σ)], B)
        end
    end
    total
end
