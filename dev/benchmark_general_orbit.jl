# Prototype: general-base integer orbit search for specialized translational bases
# Currently only base=2 uses _index_bits; base>=3 uses circshift!+index.
# Test extending integer orbit to all bases.

using EDKit
using LinearAlgebra, SparseArrays
using BenchmarkTools

BenchmarkTools.DEFAULT_PARAMETERS.samples = 10
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 3

# ============================================================================
# General-base integer orbit search prototypes
# ============================================================================

# Translation: right-shift by A positions in base-B digit representation
# circshift!(dgt, A) moves last A elements to front
# Integer: extract last A "digits", prepend them
@inline function int_rotate(state::Int64, pow_shift::Int64, base_shift::Int64)
    tail = state % base_shift
    tail * pow_shift + state ÷ base_shift
end

# Reverse digit order (parity)
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
@inline function int_spinflip(state::Int64, maxstate::Int64)
    maxstate - state  # maxstate = base^L - 1
end

# ============================================================================
# Prototype: TranslationalBasis index for general base
# ============================================================================
function proto_index_trans(b_dgt, b_I, b_R, b_C, b_A, b_B, ncyc)
    I0 = EDKit.index(b_dgt, base=b_B)
    state = I0 - Int64(1)
    Im, M = I0, 0
    L = length(b_dgt)
    A = b_A
    pow_shift = Int64(b_B)^(L - A)
    base_shift = Int64(b_B)^A

    curr = state
    for i in 1:ncyc-1
        curr = int_rotate(curr, pow_shift, base_shift)
        In = curr + Int64(1)
        if isequal(In, I0)
            break
        elseif In < Im
            Im, M = In, i
        end
    end

    i = EDKit.binary_search(b_I, Im)
    iszero(i) && return (zero(eltype(b_C)), one(b_B))
    @inbounds N = b_C[M+1] * b_R[i]
    N, i
end

# ============================================================================
# Prototype: TranslationParityBasis index for general base
# ============================================================================
function proto_index_transparity(b_dgt, b_I, b_R, b_C, b_P, b_A, b_B, ncyc)
    Ia0 = EDKit.index(b_dgt, base=b_B)
    state = Ia0 - Int64(1)
    L = length(b_dgt)
    A = b_A
    pow_shift = Int64(b_B)^(L - A)
    base_shift = Int64(b_B)^A
    base_i = Int64(b_B)

    Ib0 = int_reverse(state, L, base_i) + Int64(1)
    Iam, Ibm = Ia0, Ib0
    R, M = 0, 0

    curr = state
    for i in 1:ncyc-1
        curr = int_rotate(curr, pow_shift, base_shift)
        Ia1 = curr + Int64(1)
        if isequal(Ia1, Ia0)
            break
        elseif Ia1 < Iam
            Iam, R = Ia1, i
        end
        Ib1 = int_reverse(curr, L, base_i) + Int64(1)
        if Ib1 < Ibm
            Ibm, M = Ib1, i
        end
    end

    i, n = if Ibm < Iam
        EDKit.binary_search(b_I, Ibm), b_P * b_C[M+1]
    else
        EDKit.binary_search(b_I, Iam), b_C[R+1]
    end
    iszero(i) && return (zero(eltype(b_C)), one(eltype(b_I)))
    @inbounds n * b_R[i], i
end

# ============================================================================
# Prototype: TranslationFlipBasis index for general base
# ============================================================================
function proto_index_transflip(b_dgt, b_I, b_R, b_C, b_P, b_M, b_A, b_B, ncyc)
    Ia0 = EDKit.index(b_dgt, base=b_B)
    state = Ia0 - Int64(1)
    L = length(b_dgt)
    A = b_A
    pow_shift = Int64(b_B)^(L - A)
    base_shift = Int64(b_B)^A
    maxstate = Int64(b_B)^L - Int64(1)

    Ib0 = int_spinflip(state, maxstate) + Int64(1)
    Iam, Ibm = Ia0, Ib0
    R, M = 0, 0

    curr = state
    for i in 1:ncyc-1
        curr = int_rotate(curr, pow_shift, base_shift)
        Ia1 = curr + Int64(1)
        if isequal(Ia1, Ia0)
            break
        elseif Ia1 < Iam
            Iam, R = Ia1, i
        end
        Ib1 = int_spinflip(curr, maxstate) + Int64(1)
        if Ib1 < Ibm
            Ibm, M = Ib1, i
        end
    end

    i, n = if Ibm < Iam
        EDKit.binary_search(b_I, Ibm), b_P * b_C[M+1]
    else
        EDKit.binary_search(b_I, Iam), b_C[R+1]
    end
    iszero(i) && return (zero(eltype(b_C)), one(eltype(b_I)))
    @inbounds n * b_R[i], i
end

# ============================================================================
# Correctness: verify against existing index(b) for ALL states
# ============================================================================

function test_trans(L, base, k)
    b = TranslationalBasis(L=L, base=base, k=k, f=x->sum(x)==(base-1)*L÷2)
    dim = size(b, 1)
    dim == 0 && return true
    total = base^L
    ncyc = EDKit.ncycle(b)
    mismatches = 0
    for s in Int64(1):Int64(total)
        EDKit.change!(b.dgt, Int64(s), base=Int64(base))
        C_ref, pos_ref = EDKit.index(b)
        EDKit.change!(b.dgt, Int64(s), base=Int64(base))  # restore
        C_new, pos_new = proto_index_trans(b.dgt, b.I, b.R, b.C, b.A, b.B, ncyc)
        if pos_ref != pos_new || !(C_ref ≈ C_new)
            mismatches += 1
            mismatches <= 2 && println("    MISMATCH s=$s: ref=($C_ref,$pos_ref) new=($C_new,$pos_new)")
        end
    end
    pass = mismatches == 0
    println("  TransBasis L=$L base=$base k=$k (dim=$dim): $(pass ? "PASS" : "FAIL ($mismatches/$total)")")
    pass
end

function test_transparity(L, base)
    b = TranslationParityBasis(L=L, base=base, k=0, p=1, f=x->sum(x)==(base-1)*L÷2)
    dim = size(b, 1)
    dim == 0 && return true
    total = base^L
    ncyc = EDKit.ncycle(b)
    mismatches = 0
    for s in Int64(1):Int64(total)
        EDKit.change!(b.dgt, Int64(s), base=Int64(base))
        C_ref, pos_ref = EDKit.index(b)
        EDKit.change!(b.dgt, Int64(s), base=Int64(base))
        C_new, pos_new = proto_index_transparity(b.dgt, b.I, b.R, b.C, b.P, b.A, b.B, ncyc)
        if pos_ref != pos_new || !(C_ref ≈ C_new)
            mismatches += 1
            mismatches <= 2 && println("    MISMATCH s=$s: ref=($C_ref,$pos_ref) new=($C_new,$pos_new)")
        end
    end
    pass = mismatches == 0
    println("  TransParBasis L=$L base=$base k=0,p=1 (dim=$dim): $(pass ? "PASS" : "FAIL ($mismatches/$total)")")
    pass
end

function test_transflip(L, base)
    b = TranslationFlipBasis(L=L, base=base, k=0, p=1, f=x->sum(x)==(base-1)*L÷2)
    dim = size(b, 1)
    dim == 0 && return true
    total = base^L
    ncyc = EDKit.ncycle(b)
    mismatches = 0
    for s in Int64(1):Int64(total)
        EDKit.change!(b.dgt, Int64(s), base=Int64(base))
        C_ref, pos_ref = EDKit.index(b)
        EDKit.change!(b.dgt, Int64(s), base=Int64(base))
        C_new, pos_new = proto_index_transflip(b.dgt, b.I, b.R, b.C, b.P, b.M, b.A, b.B, ncyc)
        if pos_ref != pos_new || !(C_ref ≈ C_new)
            mismatches += 1
            mismatches <= 2 && println("    MISMATCH s=$s: ref=($C_ref,$pos_ref) new=($C_new,$pos_new)")
        end
    end
    pass = mismatches == 0
    println("  TransFlipBasis L=$L base=$base k=0,p=1 (dim=$dim): $(pass ? "PASS" : "FAIL ($mismatches/$total)")")
    pass
end

println("="^70)
println("  CORRECTNESS: general-base integer orbit for specialized bases")
println("="^70)

println("\nBase=2:")
for L in [8, 10, 12]; test_trans(L, 2, 0); end
for L in [8, 10, 12]; test_trans(L, 2, 1); end
for L in [8, 10]; test_transparity(L, 2); end
for L in [8, 10]; test_transflip(L, 2); end

println("\nBase=3:")
for L in [6, 8]; test_trans(L, 3, 0); end
for L in [6, 8]; test_trans(L, 3, 1); end
for L in [6, 8]; test_transparity(L, 3); end
for L in [6, 8]; test_transflip(L, 3); end

println("\nBase=4:")
for L in [4, 6]; test_trans(L, 4, 0); end
for L in [4, 6]; test_transparity(L, 4); end

# ============================================================================
# Performance: end-to-end Operator * vector
# ============================================================================
println("\n\n")
println("="^70)
println("  PERFORMANCE: Operator * vector (current vs integer orbit)")
println("="^70)

function heisenberg_op(L, basis; base=2)
    if base == 2
        Sx = [0 0.5; 0.5 0]; Sy = [0 -0.5im; 0.5im 0]; Sz = [0.5 0; 0 -0.5]
    elseif base == 3
        Sx = [0 1 0; 1 0 1; 0 1 0] / sqrt(2)
        Sy = [0 -1im 0; 1im 0 -1im; 0 1im 0] / sqrt(2)
        Sz = [1.0 0 0; 0 0 0; 0 0 -1.0]
    elseif base == 4
        s = 3/2
        m = [s, s-1, s-2, s-3]
        Sp = diagm(1 => [sqrt(s*(s+1)-m[i]*(m[i]+1)) for i in 1:3])
        Sm = Sp'
        Sx = (Sp + Sm) / 2; Sy = (Sp - Sm) / 2im; Sz = diagm(m)
    end
    h = kron(Sx, Sx) + kron(Sy, Sy) + kron(Sz, Sz)
    mats = [sparse(h) for _ in 1:L]
    inds = [[i, mod1(i+1, L)] for i in 1:L]
    operator(mats, [collect(ind) for ind in inds], basis)
end

# Patch: replace index(b) call in colmn! by calling the prototype
# Can't easily do this without modifying source, so instead compare
# just the orbit search overhead isolated in a tight loop.

println("\nOrbit search cost per call (averaged over all basis states):")

function bench_trans(L, base)
    Nfill = (base - 1) * L ÷ 2
    b = TranslationalBasis(L=L, base=base, k=0, f=x->sum(x)==Nfill)
    dim = size(b, 1)
    dim == 0 && return
    ncyc = EDKit.ncycle(b)
    bI, bR, bC, bA, bB = b.I, b.R, b.C, b.A, b.B

    b_cur = @benchmark begin
        s = 0
        for j in 1:$dim
            EDKit.change!($b, j)
            _, pos = EDKit.index($b)
            s += pos
        end
        s
    end
    b_new = @benchmark begin
        s = 0
        for j in 1:$dim
            EDKit.change!($b, j)
            _, pos = proto_index_trans($b.dgt, $bI, $bR, $bC, $bA, $bB, $ncyc)
            s += pos
        end
        s
    end
    t_cur = median(b_cur).time; t_new = median(b_new).time
    println("  L=$L base=$base (dim=$dim): current=$(BenchmarkTools.prettytime(t_cur/dim))/call  " *
            "integer=$(BenchmarkTools.prettytime(t_new/dim))/call  " *
            "speedup=$(round(t_cur/t_new, digits=2))x")
end

function bench_transparity(L, base)
    Nfill = (base - 1) * L ÷ 2
    b = TranslationParityBasis(L=L, base=base, k=0, p=1, f=x->sum(x)==Nfill)
    dim = size(b, 1)
    dim == 0 && return
    ncyc = EDKit.ncycle(b)
    bI, bR, bC, bP, bA, bB = b.I, b.R, b.C, b.P, b.A, b.B

    b_cur = @benchmark begin
        s = 0
        for j in 1:$dim
            EDKit.change!($b, j)
            _, pos = EDKit.index($b)
            s += pos
        end
        s
    end
    b_new = @benchmark begin
        s = 0
        for j in 1:$dim
            EDKit.change!($b, j)
            _, pos = proto_index_transparity($b.dgt, $bI, $bR, $bC, $bP, $bA, $bB, $ncyc)
            s += pos
        end
        s
    end
    t_cur = median(b_cur).time; t_new = median(b_new).time
    println("  L=$L base=$base (dim=$dim): current=$(BenchmarkTools.prettytime(t_cur/dim))/call  " *
            "integer=$(BenchmarkTools.prettytime(t_new/dim))/call  " *
            "speedup=$(round(t_cur/t_new, digits=2))x")
end

function bench_transflip(L, base)
    Nfill = (base - 1) * L ÷ 2
    b = TranslationFlipBasis(L=L, base=base, k=0, p=1, f=x->sum(x)==Nfill)
    dim = size(b, 1)
    dim == 0 && return
    ncyc = EDKit.ncycle(b)
    bI, bR, bC, bP, bM, bA, bB = b.I, b.R, b.C, b.P, b.M, b.A, b.B

    b_cur = @benchmark begin
        s = 0
        for j in 1:$dim
            EDKit.change!($b, j)
            _, pos = EDKit.index($b)
            s += pos
        end
        s
    end
    b_new = @benchmark begin
        s = 0
        for j in 1:$dim
            EDKit.change!($b, j)
            _, pos = proto_index_transflip($b.dgt, $bI, $bR, $bC, $bP, $bM, $bA, $bB, $ncyc)
            s += pos
        end
        s
    end
    t_cur = median(b_cur).time; t_new = median(b_new).time
    println("  L=$L base=$base (dim=$dim): current=$(BenchmarkTools.prettytime(t_cur/dim))/call  " *
            "integer=$(BenchmarkTools.prettytime(t_new/dim))/call  " *
            "speedup=$(round(t_cur/t_new, digits=2))x")
end

println("\nTranslationalBasis:")
for (L, base) in [(12, 2), (14, 2), (16, 2), (8, 3), (10, 3)]
    bench_trans(L, base)
end

println("\nTranslationParityBasis:")
for (L, base) in [(10, 2), (12, 2), (8, 3)]
    bench_transparity(L, base)
end

println("\nTranslationFlipBasis:")
for (L, base) in [(10, 2), (12, 2), (8, 3)]
    bench_transflip(L, base)
end

println("\nDone!")
