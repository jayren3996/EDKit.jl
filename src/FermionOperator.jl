export fermion, fermion_operator

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
inserted on sites `2, 3, …, span-1` whenever `span > 2`, plus on site 1 of the
kron product whenever the second operator brings a `σ_z` chain through.

`convention` selects the Jordan-Wigner direction. Only `:left` is implemented
in this MVP.

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
    if op == "nn"
        n_mat = [0.0 0.0; 0.0 1.0]
        Id2 = Matrix{Float64}(I, 2, 2)
        factors = Vector{Matrix{Float64}}(undef, span)
        factors[1] = n_mat
        for k in 2:span-1
            factors[k] = Id2
        end
        factors[span] = n_mat
        mat = factors[1]
        for k in 2:span
            mat = kron(mat, factors[k])
        end
        return mat
    end

    σ⁺ = Array(spin_Sp(2))
    σ⁻ = Array(spin_Sm(2))
    σ_z = Array(spin_Sz(2)) * 2  # spin_Sz returns Sz = σ_z / 2; we want σ_z = diag(1, -1)

    site1 = op[1] == '+' ? σ⁻ : op[1] == '-' ? σ⁺ : error("Unsupported operator char: '$(op[1])' in \"$op\".")
    site2 = op[2] == '+' ? σ⁻ : op[2] == '-' ? σ⁺ : error("Unsupported operator char: '$(op[2])' in \"$op\".")

    factors = Vector{Matrix{Float64}}(undef, span)
    factors[1] = site1
    for k in 2:span-1
        factors[k] = σ_z
    end
    factors[span] = site2

    # Sign comes from the on-site product at site 1:
    #   σ⁻ σ_z = +σ⁻  → "+? " → +1
    #   σ⁺ σ_z = −σ⁺  → "-? " → −1
    sign = op[1] == '+' ? 1.0 : -1.0

    mat = factors[1]
    for k in 2:span
        mat = kron(mat, factors[k])
    end
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
        σ_z = [1.0 0.0; 0.0 -1.0]
        endpoint = fermion(op)
        factors = Vector{Matrix{Float64}}(undef, i)
        for k in 1:i-1
            factors[k] = σ_z
        end
        factors[i] = endpoint
        mat = factors[1]
        for k in 2:i
            mat = kron(mat, factors[k])
        end
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

    i == j && error("\"$op\" on the same site is identically zero (Pauli exclusion).")
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
