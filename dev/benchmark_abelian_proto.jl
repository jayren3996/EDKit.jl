# Prototype: integer-based shift_canonical! for AbelianBasis
# Test correctness against current implementation, then benchmark.

using EDKit
using LinearAlgebra, SparseArrays
using BenchmarkTools

BenchmarkTools.DEFAULT_PARAMETERS.samples = 10
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 3

# ============================================================================
# Integer group actions (work on 0-based integer state, any base)
# ============================================================================

# Translation: circshift!(dgt, a) on integer
# circshift! shifts elements RIGHT by a positions:
#   [d1,d2,...,dL] -> [d(L-a+1),...,dL,d1,...,d(L-a)]
# In integer land (big-endian encoding):
#   state = d1*B^(L-1) + d2*B^(L-2) + ... + dL
# After right-shift by a: last a digits become first a digits
@inline function int_translate(state::Int64, pow_shift::Int64, base_shift::Int64)
    # pow_shift = base^(L-a), base_shift = base^a
    tail = state % base_shift
    tail * pow_shift + state ÷ base_shift
end

# Parity (reverse): reverse digit order
@inline function int_reverse(state::Int64, L::Int, base::Int64)
    r = Int64(0)
    v = state
    for _ in 1:L
        r = r * base + v % base
        v ÷= base
    end
    r
end

# Spin-flip: complement each digit (d -> base-1-d)
# For the full integer: base^L - 1 - state
@inline function int_spinflip(state::Int64, maxstate::Int64)
    maxstate - state
end

# ============================================================================
# Prototype: integer shift_canonical! replacement
# ============================================================================

# Identify which action functions are in the AbelianOperator
function identify_actions(g::EDKit.AbelianOperator, L::Int, base::Int64)
    # Probe each generator by applying it to a known state and checking the result
    actions = Symbol[]
    for (i, f) in enumerate(g.f)
        dgt_test = collect(Int64(0):Int64(L-1))  # [0,1,2,...,L-1]
        dgt_copy = copy(dgt_test)
        f(dgt_copy)

        if dgt_copy == circshift(dgt_test, 1)
            push!(actions, :translate_1)
        elseif any(a -> dgt_copy == circshift(dgt_test, a), 2:L-1)
            a = findfirst(a -> dgt_copy == circshift(dgt_test, a), 1:L-1)
            push!(actions, Symbol("translate_$a"))
        elseif dgt_copy == reverse(dgt_test)
            push!(actions, :parity)
        else
            # Check spin-flip: each digit d -> (base-1) - d
            dgt_flip = (base - 1) .- dgt_test
            if dgt_copy == dgt_flip
                push!(actions, :spinflip)
            else
                push!(actions, :unknown)
            end
        end
    end
    actions
end

# Build integer action functions from identified actions
function build_int_actions(actions::Vector{Symbol}, L::Int, base::Int64)
    maxstate = base^L - Int64(1)
    int_fns = Function[]
    for act in actions
        if act == :translate_1
            pow = Int64(base)^(L - 1)
            bshift = Int64(base)
            push!(int_fns, state -> int_translate(state, pow, bshift))
        elseif startswith(string(act), "translate_")
            a = parse(Int, string(act)[11:end])
            pow = Int64(base)^(L - a)
            bshift = Int64(base)^a
            push!(int_fns, state -> int_translate(state, pow, bshift))
        elseif act == :parity
            push!(int_fns, state -> int_reverse(state, L, base))
        elseif act == :spinflip
            push!(int_fns, state -> int_spinflip(state, maxstate))
        else
            error("Unknown action: $act")
        end
    end
    int_fns
end

# Integer version of shift_canonical!
# Returns (Im, phase) -- the minimum 0-based index and the phase
function int_shift_canonical(state0::Int64, g::EDKit.AbelianOperator,
                              int_fns::Vector{Function})
    Im = state0
    best_s = ones(Int, length(g.s))  # identity counters (1-indexed)
    current_s = ones(Int, length(g.s))

    # Iterate through all group elements (same order as AbelianOperator call)
    curr = state0
    for _ in 1:EDKit.order(g)
        # Advance: apply next group element (same logic as ag(dgt))
        for i in eachindex(current_s)
            curr = int_fns[i](curr)
            current_s[i] += 1
            if current_s[i] > g.g[i]
                current_s[i] = 1
            else
                break
            end
        end
        if curr < Im
            Im = curr
            best_s .= current_s
        end
    end

    # Compute phase from best_s
    ph = prod(g.c[i][best_s[i]] for i in eachindex(g.c))
    Im, ph
end

# ============================================================================
# Correctness test: compare integer vs current shift_canonical! for all states
# ============================================================================

function test_correctness(label, b::EDKit.AbelianBasis)
    L = length(b)
    base = Int64(b.B)
    total = Int64(base)^L

    actions = identify_actions(b.G, L, base)
    println("  $label: actions = $actions")

    if :unknown in actions
        println("    SKIPPED: unknown action detected")
        return false
    end

    int_fns = build_int_actions(actions, L, base)

    mismatches = 0
    for s in Int64(1):total
        # Load state into basis
        EDKit.change!(b.dgt, s, base=base)

        # Current implementation
        dgt_save = copy(b.dgt)
        EDKit.init!(b.G)
        Im_orig, g_orig = EDKit.shift_canonical!(b.dgt, b.G; base=base)

        # Get phase from current
        ph_orig = EDKit.phase(g_orig)

        # Integer implementation
        state0 = s - Int64(1)  # 0-based
        Im_int, ph_int = int_shift_canonical(state0, b.G, int_fns)

        # Compare (Im_orig is 1-based, Im_int is 0-based)
        if Im_orig != Im_int + 1
            mismatches += 1
            if mismatches <= 3
                println("    INDEX MISMATCH at s=$s: orig=$Im_orig, int=$(Im_int+1)")
            end
        elseif !(ph_orig ≈ ph_int)
            mismatches += 1
            if mismatches <= 3
                println("    PHASE MISMATCH at s=$s: orig=$ph_orig, int=$ph_int")
            end
        end
    end

    if mismatches == 0
        println("    PASS: all $total states match")
    else
        println("    FAIL: $mismatches / $total mismatches")
    end
    mismatches == 0
end

println("="^70)
println("  CORRECTNESS TESTS")
println("="^70)

# Base=2 tests
for L in [6, 8, 10]
    println("\nL=$L, base=2:")
    test_correctness("translation k=0",
        basis(L=L, base=2, k=0, f=x->sum(x)==L÷2))
    test_correctness("translation k=1",
        basis(L=L, base=2, k=1, f=x->sum(x)==L÷2))
    test_correctness("trans+parity k=0,p=1",
        basis(L=L, base=2, k=0, p=1, f=x->sum(x)==L÷2))
    test_correctness("trans+flip k=0,z=1",
        basis(L=L, base=2, k=0, z=1, f=x->sum(x)==L÷2))
    if L <= 8
        test_correctness("trans+parity+flip k=0,p=1,z=1",
            basis(L=L, base=2, k=0, p=1, z=1, f=x->sum(x)==L÷2))
    end
end

# Base=3 tests
for L in [6, 8]
    println("\nL=$L, base=3:")
    test_correctness("translation k=0",
        basis(L=L, base=3, k=0, f=x->sum(x)==L))
    test_correctness("translation k=1",
        basis(L=L, base=3, k=1, f=x->sum(x)==L))
    test_correctness("trans+parity k=0,p=1",
        basis(L=L, base=3, k=0, p=1, f=x->sum(x)==L))
end

# ============================================================================
# Performance: end-to-end Operator * vector with integer orbit search
# ============================================================================
println("\n\n")
println("="^70)
println("  PERFORMANCE: integer orbit vs current (Operator * vector)")
println("="^70)

# Prototype: patch index(B::AbelianBasis) to use integer orbit
function patched_mul_vec(opt, v, int_fns)
    b = opt.B
    L = length(b)
    base = Int64(b.B)
    ctype = promote_type(eltype(opt), eltype(v))
    target = zeros(ctype, size(opt, 1))

    for j = 1:length(v)
        r = EDKit.change!(b, j)
        C = isone(r) ? v[j] : v[j] / r

        for ti = 1:length(opt.M)
            M_local, I_sites = opt.M[ti], opt.I[ti]
            rows, vals = rowvals(M_local), nonzeros(M_local)
            jj = EDKit.index(b.dgt, I_sites, base=b.B)
            for ii in nzrange(M_local, jj)
                row, val = rows[ii], vals[ii]
                EDKit.change!(b.dgt, I_sites, row, base=b.B)

                # Integer orbit search instead of shift_canonical!
                state0 = EDKit.index(b.dgt, base=b.B) - Int64(1)
                Im_int, ph_int = int_shift_canonical(state0, b.G, int_fns)

                ind = EDKit.binary_search(b.I, Im_int + Int64(1))
                if !iszero(ind)
                    coeff = ph_int * b.R[ind]
                    target[ind] += C * coeff * val
                end
            end
            # Restore subsystem digits
            EDKit.change!(b.dgt, I_sites, jj, base=b.B)
        end
    end
    target
end

function bench_comparison(label, L, base, bfn)
    b = bfn()
    dim = size(b, 1)
    dim == 0 && return

    if base == 2
        Sx = [0 0.5; 0.5 0]; Sy = [0 -0.5im; 0.5im 0]; Sz = [0.5 0; 0 -0.5]
    else
        Sx = [0 1 0; 1 0 1; 0 1 0] / sqrt(2)
        Sy = [0 -1im 0; 1im 0 -1im; 0 1im 0] / sqrt(2)
        Sz = [1.0 0 0; 0 0 0; 0 0 -1.0]
    end
    h = kron(Sx, Sx) + kron(Sy, Sy) + kron(Sz, Sz)
    mats = [sparse(h) for _ in 1:L]
    inds = [[i, mod1(i+1, L)] for i in 1:L]
    opt = operator(mats, [collect(ind) for ind in inds], b)

    v = randn(ComplexF64, dim)

    # Current
    ref = opt * v

    # Integer orbit
    actions = identify_actions(b.G, L, Int64(base))
    if :unknown in actions
        println("  $label: SKIPPED (unknown action)")
        return
    end
    int_fns = build_int_actions(actions, L, Int64(base))
    res = patched_mul_vec(opt, v, int_fns)

    # Verify correctness
    err = norm(res - ref) / norm(ref)
    if err > 1e-10
        println("  $label (dim=$dim): INCORRECT! err=$err")
        return
    end

    b_cur = @benchmark $opt * $v
    # Fresh basis for patched version each time
    b_new = @benchmark patched_mul_vec($opt, $v, $int_fns)

    t_cur = median(b_cur).time
    t_new = median(b_new).time
    println("  $label (dim=$dim): current=$(BenchmarkTools.prettytime(t_cur))  " *
            "integer=$(BenchmarkTools.prettytime(t_new))  " *
            "speedup=$(round(t_cur/t_new, digits=2))x")
end

println("\nBase=2:")
for L in [12, 14, 16]
    bench_comparison("L=$L trans k=0", L, 2,
        () -> basis(L=L, base=2, k=0, f=x->sum(x)==L÷2))
end
for L in [12, 14]
    bench_comparison("L=$L trans+par k=0,p=1", L, 2,
        () -> basis(L=L, base=2, k=0, p=1, f=x->sum(x)==L÷2))
    bench_comparison("L=$L trans+flip k=0,z=1", L, 2,
        () -> basis(L=L, base=2, k=0, z=1, f=x->sum(x)==L÷2))
end

println("\nBase=3:")
for L in [6, 8, 10]
    bench_comparison("L=$L trans k=0", L, 3,
        () -> basis(L=L, base=3, k=0, f=x->sum(x)==L))
end
for L in [6, 8]
    bench_comparison("L=$L trans+par k=0,p=1", L, 3,
        () -> basis(L=L, base=3, k=0, p=1, f=x->sum(x)==L))
end

println("\nDone!")
