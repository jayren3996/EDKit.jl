export fermion, fermion_operator

"""
    fermion(op::AbstractString[, span::Integer]; convention::Symbol=:left)

Return the local matrix representation of a spinless-fermion operator on a
contiguous block of `span` sites.

The operator string uses one character per operator:
- `"n"` — single-site number operator `n = c†c`
- `"+"` — creation operator `c†`
- `"-"` — annihilation operator `c`

Multi-character strings build two-site fermion products:
- `"+-"` — hop `c†_1 c_span`
- `"-+"` — `c_1 c†_span`
- `"++"` — pair creation `c†_1 c†_span`
- `"--"` — pair annihilation `c_1 c_span`

The Jordan-Wigner `σ_z` string is inserted on sites `2, 3, …, span-1` whenever
`span > 2`, plus on site 1 of the kron product whenever the second operator
brings a `σ_z` chain through.

`convention` selects the Jordan-Wigner direction. Only `:left` is implemented
in this MVP. The default `span` for a one-character string is `1`; for a
two-character string it must be supplied explicitly and is `>= 2`.

Embed the returned matrix into a many-body basis with [`operator`](@ref) or
[`fermion_operator`](@ref).
"""
function fermion(op::AbstractString, span::Integer=length(op); convention::Symbol=:left)
    convention === :left || error("Only the :left Jordan-Wigner convention is implemented (got $convention).")

    if op == "n"
        span == 1 || error("\"n\" is a single-site operator (got span=$span).")
        return [0.0 0.0; 0.0 1.0]
    end

    if length(op) == 1
        c = op[1]
        span == 1 || error("Single-character operator \"$op\" is a one-site operator (got span=$span).")
        return c == '+' ? Array(spin_Sm(2)) :
               c == '-' ? Array(spin_Sp(2)) :
               error("Unsupported single-site fermion operator: \"$op\". Supported: \"n\", \"+\", \"-\".")
    end

    length(op) == 2 || error("Only 1- or 2-character fermion operator strings are supported in this MVP (got \"$op\").")
    span >= 2 || error("Two-fermion operators require span >= 2 (got $span).")

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
- `"n"` — single-site number operator `n_i` (requires `length(sites) == 1`).
- `"+-"` — hop `c†_i c_j`.
- `"-+"` — `c_i c†_j`.
- `"++"` — pair creation `c†_i c†_j`.
- `"--"` — pair annihilation `c_i c_j`.

Two-character strings require `length(sites) == 2`; the operator's support is
the contiguous block `min(sites):max(sites)`, with the Jordan-Wigner `σ_z`
string baked in on intermediate sites by [`fermion`](@ref).

If `sites` is unsorted (e.g. `[3, 1]` for `c†_3 c_1`), the operator is
adjusted by the appropriate fermionic sign so that the returned matrix
represents the operator literally written by the caller.

Bare single-site `c†_i` or `c_i` (i.e. `"+"` or `"-"` alone) are intentionally
not supported through this entry point — they would require an explicit
Jordan-Wigner string and are an uncommon construction. If you need the raw
2×2 matrix, call [`fermion`](@ref) with `"+"` or `"-"` and embed it yourself.
"""
function fermion_operator(op::AbstractString, sites::AbstractVector{<:Integer},
        B::AbstractBasis; convention::Symbol=:left)
    convention === :left || error("Only the :left convention is implemented.")

    if op == "n"
        length(sites) == 1 || error("\"n\" takes exactly one site (got $(length(sites))).")
        return operator(fermion("n"), [sites[1]], B)
    end

    if length(op) == 1
        error(
            "fermion_operator only supports \"n\" or two-character operators (\"+-\", \"-+\", \"++\", \"--\"). " *
            "Bare \"+\" / \"-\" would require an explicit Jordan-Wigner string and is not exposed in this MVP — " *
            "if you need the local matrix, call fermion(\"+\") or fermion(\"-\") and embed it yourself.",
        )
    end

    length(op) == 2 || error("Only 1- or 2-character fermion operator strings are supported in this MVP (got \"$op\").")
    length(sites) == 2 || error("Two-fermion operators take exactly two sites.")
    i, j = sites[1], sites[2]
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
