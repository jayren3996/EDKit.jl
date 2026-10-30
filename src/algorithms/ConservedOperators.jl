export conserved_operators, eigenstate_operators

# Embed complex vectors into a real vector space when coefficients are real.
_search_columns(xs, real_coefficients) = begin
    A = hcat((vec(x) for x in xs)...)
    real_coefficients ? vcat(real.(A), imag.(A)) : A
end

function _search_options(coefficients, metric, rtol)
    coefficients in (:real, :complex) || throw(ArgumentError("coefficients must be :real or :complex"))
    metric in (:operator, :coefficients) || throw(ArgumentError("metric must be :operator or :coefficients"))
    isfinite(rtol) && 0 <= rtol < 1 || throw(ArgumentError("rtol must be finite and satisfy 0 <= rtol < 1"))
end

# Check cases whose coordinate equivalence is unambiguous. Other basis types
# require a shared basis object; distinct objects cannot safely be compared by
# their mutable work buffers or just by dimension.
function _same_search_basis(a::AbstractBasis, b::AbstractBasis)
    a === b && return true
    a isa TensorBasis && b isa TensorBasis && return length(a.dgt) == length(b.dgt) && a.B == b.B
    a isa ProjectedBasis && b isa ProjectedBasis && return length(a.dgt) == length(b.dgt) && a.B == b.B && a.I == b.I
    false
end
function _check_search_bases(inputs)
    ops = [x for x in inputs if x isa Operator]
    for o in ops
        o.B isa DoubleBasis && throw(ArgumentError("DoubleBasis operators require an explicit matrix with a common domain and codomain"))
    end
    isempty(ops) && return
    all(o -> _same_search_basis(ops[1].B, o.B), ops) ||
        throw(ArgumentError("EDKit operators must use the same basis; reuse the basis object for symmetry bases"))
end

# B is the chosen metric map, R the residual map, and E contains metric-space
# vectors to exclude. All columns are already represented over the selected
# real or complex coefficient field.
function _solve_operator_search(B, R, E; rtol)
    n = size(B, 2)
    scales = [norm(view(B, :, i)) for i in 1:n]
    nz = findall(!iszero, scales)
    isempty(nz) && throw(ArgumentError("all candidate directions have zero norm in the chosen metric"))
    Bs = B[:, nz] * Diagonal(1 ./ scales[nz])
    f = svd(Bs; full=size(Bs, 1) < size(Bs, 2))
    rank = count(x -> x > rtol * first(f.S), f.S)
    rank > 0 || throw(ArgumentError("no candidate directions remain after rank truncation"))
    W = zeros(eltype(f.V), n, rank)
    W[nz, :] = Diagonal(1 ./ scales[nz]) * f.V[:, 1:rank] * Diagonal(1 ./ f.S[1:rank])

    discarded = zeros(eltype(W), n, n - rank)
    for (k, i) in enumerate(findall(iszero, scales))
        discarded[i, k] = 1
    end
    offset = n - length(nz)
    if rank < length(nz)
        discarded[nz, offset+1:end] = Diagonal(1 ./ scales[nz]) * f.V[:, rank+1:end]
    end

    excluded_rank = 0
    if size(E, 2) > 0
        # Normalize constraints before rank detection, without changing their null space.
        F = E' * B * W
        rownorms = [norm(view(F, i, :)) for i in axes(F, 1)]
        # Treat roundoff-level projection of an orthogonal exclusion as zero.
        enorms = [norm(view(E, :, i)) for i in axes(E, 2)]
        nonzero = [i for i in eachindex(rownorms) if rownorms[i] > max(rtol, eps(Float64) * 10) * enorms[i]]
        if !isempty(nonzero)
            F = Diagonal(1 ./ rownorms[nonzero]) * F[nonzero, :]
            ef = svd(F; full=size(F, 1) < size(F, 2))
            excluded_rank = count(x -> x > rtol * first(ef.S), ef.S)
            excluded_rank < rank || throw(ArgumentError("exclusions remove the entire candidate space"))
            W = W * ef.V[:, excluded_rank+1:end]
        end
    end
    A = R * W
    s = svd(A; full=size(A, 1) < size(A, 2))
    k = size(W, 2)
    sigma = vcat(s.S, zeros(eltype(s.S), k - length(s.S)))
    order = sortperm(sigma)
    vectors = W * s.V[:, order]
    residuals = [norm(R * vectors[:, i]) for i in axes(vectors, 2)]
    (; values=sigma[order].^2, vectors, residuals, rank, discarded,
       singular_values=f.S, rank_tolerance=rtol, excluded_rank)
end

"""
    conserved_operators(H, candidates; coefficients=:real, metric=:operator,
                        exclude=[], rtol=1e-10)

Find combinations `Q = sum(c[a] * candidates[a])` minimizing the squared
Hilbert–Schmidt norm of `[H,Q]` divided by the chosen norm of `Q`. Inputs may
be matrices or EDKit `Operator`s acting in the same nonempty basis. Operators
are converted to dense matrices; this is intended for small systems.

`metric=:operator` (default) normalizes the Hilbert–Schmidt norm of `Q`;
`:coefficients` uses the Euclidean norm of `c`, which depends on candidate
scaling and may admit zero physical operators. `exclude` is a vector of
operators to which the result must be Hilbert–Schmidt orthogonal; it requires
`metric=:operator`. `coefficients=:real` requires Hermitian candidates;
`:complex` permits general candidates. `rtol` is a relative threshold on
singular values after candidate normalization.

Returns a named tuple containing `values`, coefficient columns `vectors`,
`commutator` and `overlap` Gram matrices, `residuals`, retained `rank`,
`discarded` coefficient directions, `singular_values`, `rank_tolerance`, and
`excluded_rank`. `overlap` is the chosen metric. The squared residuals equal
`values` up to roundoff. Each Gram matrix uses the normalized trace `tr / D`.
For `metric=:coefficients`, the overlap is the identity instead.

A projected or symmetry-reduced basis measures conservation within that
subspace, not necessarily in the full Hilbert space. Matrix inputs carry no
basis information, so the caller must ensure their coordinates agree.
"""
function conserved_operators(H, candidates::AbstractVector; coefficients::Symbol=:real,
                             metric::Symbol=:operator, exclude::AbstractVector=[], rtol::Real=1e-10)
    _search_options(coefficients, metric, rtol)
    isempty(candidates) && throw(ArgumentError("at least one candidate is required"))
    !isempty(exclude) && metric !== :operator && throw(ArgumentError("exclude requires metric=:operator"))
    _check_search_bases(vcat(Any[H], candidates, exclude))
    h = Matrix(H)
    d = size(h, 1)
    d > 0 && size(h, 2) == d || throw(ArgumentError("H must be a nonempty square matrix"))
    ops = [Matrix(o) for o in candidates]
    excluded = [Matrix(o) for o in exclude]
    all(o -> size(o) == (d, d), [ops; excluded]) || throw(DimensionMismatch("all operators must have the same dimensions as H"))
    if coefficients === :real
        all(o -> isapprox(o, o'; rtol=1e-10, atol=0), ops) || throw(ArgumentError("real coefficients require Hermitian candidates"))
    end
    real_coefficients = coefficients === :real
    physical = _search_columns(ops, real_coefficients) / sqrt(d)
    R = _search_columns([h * o - o * h for o in ops], real_coefficients) / sqrt(d)
    B = metric === :operator ? physical : Matrix{eltype(physical)}(I, length(ops), length(ops))
    E = isempty(excluded) ? zeros(eltype(B), size(B, 1), 0) : _search_columns(excluded, real_coefficients) / sqrt(d)
    sol = _solve_operator_search(B, R, E; rtol)
    (; sol..., commutator=Hermitian(R' * R), overlap=Hermitian(B' * B))
end

"""
    eigenstate_operators(candidates, states; objective=:each_state,
                         metric=:coefficients, coefficients=:real, weights=nothing,
                         exclude=[], rtol=1e-10)

Find candidate operators for which the supplied state vector, or columns of a
state matrix, have small variance. States are normalized individually and
nonnegative `weights` are normalized to sum to one. `objective=:each_state`
minimizes the weighted average of individual variances, allowing different
eigenvalues for different states. `objective=:mixture` uses the variance in
the corresponding mixed state and therefore requires a common eigenvalue for
zero variance. Existing `covmat` and `qimsolve` retain their original behavior.

The default `metric=:coefficients` avoids materializing operator matrices but
depends on candidate scaling and can admit combinations representing the zero
operator. `metric=:operator` uses the normalized Hilbert–Schmidt norm and
removes dependent physical operators, but materializes dense matrices.
`exclude` requires `metric=:operator`. Real coefficients do not by themselves
ensure Hermiticity: supply Hermitian candidates when a Hamiltonian is desired.

Returns the same diagnostics as [`conserved_operators`](@ref), except that the
objective Gram matrix is named `covariance` rather than `commutator`. The
search acts in the basis represented by the inputs; projected or symmetry
bases yield only an in-subspace statement.
"""
function eigenstate_operators(candidates::AbstractVector, states::AbstractVecOrMat;
        objective::Symbol=:each_state, metric::Symbol=:coefficients, coefficients::Symbol=:real,
        weights=nothing, exclude::AbstractVector=[], rtol::Real=1e-10)
    _search_options(coefficients, metric, rtol)
    objective in (:each_state, :mixture) || throw(ArgumentError("objective must be :each_state or :mixture"))
    isempty(candidates) && throw(ArgumentError("at least one candidate is required"))
    !isempty(exclude) && metric !== :operator && throw(ArgumentError("exclude requires metric=:operator"))
    _check_search_bases([candidates; exclude])
    V = states isa AbstractVector ? reshape(states, :, 1) : states
    d, m = size(V)
    d > 0 && m > 0 || throw(ArgumentError("states must be nonempty"))
    all(o -> size(o) == (d, d), [candidates; exclude]) || throw(DimensionMismatch("operators and states must have matching dimensions"))
    norms = [norm(view(V, :, i)) for i in 1:m]
    all(x -> isfinite(x) && x > 0, norms) || throw(ArgumentError("states must have finite, nonzero norms"))
    P = V * Diagonal(1 ./ norms)
    w = isnothing(weights) ? fill(1 / m, m) : collect(weights)
    length(w) == m && all(x -> x isa Real && isfinite(x) && x >= 0, w) && sum(w) > 0 ||
        throw(ArgumentError("weights must be nonnegative, finite, match the states, and have positive sum"))
    w = w ./ sum(w)
    actions = [o * P for o in candidates]
    residuals = map(actions) do action
        means = [dot(view(P, :, i), view(action, :, i)) for i in 1:m]
        objective === :mixture && (means .= sum(w .* means))
        (action - P * Diagonal(means)) * Diagonal(sqrt.(w))
    end
    real_coefficients = coefficients === :real
    R = _search_columns(residuals, real_coefficients)
    if metric === :operator
        ops = [Matrix(o) for o in candidates]
        B = _search_columns(ops, real_coefficients) / sqrt(d)
        E = isempty(exclude) ? zeros(eltype(B), size(B, 1), 0) : _search_columns([Matrix(o) for o in exclude], real_coefficients) / sqrt(d)
    else
        B = Matrix{eltype(R)}(I, length(candidates), length(candidates))
        E = zeros(eltype(B), size(B, 1), 0)
    end
    sol = _solve_operator_search(B, R, E; rtol)
    (; sol..., covariance=Hermitian(R' * R), overlap=Hermitian(B' * B))
end
