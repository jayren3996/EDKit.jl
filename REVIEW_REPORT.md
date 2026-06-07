# EDKit.jl — Consolidated Review Report

## 1. Executive summary

EDKit.jl is physically correct on its core job: every reviewed basis (onsite, translational, parity/flip, Abelian 2D/3D, spinless-fermion) was numerically verified to reconstruct full spectra to machine precision, the operator kernels agree across paths, and the explicit-buffer thread-safety contract is honored on every internal hot path. The defects cluster into a few recurring themes rather than scattered one-offs: **silent integer overflow** when `dtype`/base/L exceed Int width, **`ComplexF64` eltype carried by bases whose phases are purely real**, and **missing guards on documented-but-unenforced preconditions** (half-filling + inversion, symmetric predicates, sorted spectra). The single most consequential correctness bug is `abelian-1` (fixed-N plus a spin-flip generator off half-filling silently drops basis states and returns a wrong spectrum). Two API-contract bugs also matter broadly: `productstate` is silently wrong for every reduced basis, and the cached-sparse fast path advertised in the docs is unreachable from the standard `mul!` iterative-solver entry point. Most optimizations and all extensions are lower-stakes; the highest-leverage systemic win is decoupling integer width and eltype from hardcoded assumptions across the basis layer.

## 2. Top priorities

| # | Category | Severity | Location | Issue | Fix |
|---|----------|----------|----------|-------|-----|
| 1 | bug | critical | `src/Basis/AbelianBasis.jl:614-625` | Fixed N + inversion generator (z/custom) off half-filling silently drops orbits → wrong spectrum, no error | Gate `use_gosper` on `(!has_inv \|\| 2*Ndigits==L*(base-1))`; else full-scan path, or error per documented restriction |
| 2 | bug | critical | `src/ToolKit.jl:107-113` | `productstate` silently wrong for every reduced basis: drops phase coeff, Float64 output, writes to vector #1 on zero-norm | Capture `(c,I)=index(B,dgt)`; error on `iszero(c)`; allocate `zeros(eltype(B),...)`; `s[I]=c` |
| 3 | bug | high | `src/Basis/AbstractBasis.jl:297-310`, stored `ProjectedBasis.jl:202-203,242` | `index` accumulator inherits *base* type, so `dtype=Int32` or large-L base>2 silently stores corrupt/negative representatives | Decouple accumulator width from base; add constructor check `base^L-1 ≤ typemax(dtype)` |
| 4 | bug | high | `src/Operator.jl:426-451,480,505` | `mul!` (3-/5-arg) never consults the sparse cache → `sparse!` benefit unreachable from block/iterative solvers; matrix `mul!` also single-threaded | In matrix `mul!` delegate to `mul!(target,S,m,α,β)` when cached |
| 5 | bug | high | `src/Basis/TranslationalBasis.jl:352`, `TranslationalParityBasis.jl:227`, `TranslationalFlipBasis.jl:210` | `order()` returns L not orbit size L/a → mis-normalizes symmetrizer/embedding/MPS by factor a for a>1 unit cells | Return `ncycle(b)` (and `2*ncycle(b)` for parity/flip) |
| 6 | bug | high | `src/Basis/ParityBasis.jl:46`, `FlipBasis.jl:44`, `ParityFlipBasis.jl:75` | Predicate `f` checked only on representative → non-symmetric `f` silently builds over-counted, invalid basis | Re-check `f` on flipped/reversed orbit members; reject non-`f`-closed orbits |
| 7 | bug | high | `src/Basis/ParityBasis.jl:105`, `FlipBasis.jl:100`, `ParityFlipBasis.jl:121,131` | Bases inherit `eltype==ComplexF64` despite ±1 real phases → type-unstable `index`, complex matrices for real Hamiltonians | Add `eltype(...)=Float64` (mirror `TranslationalParityBasis.jl:25`) |
| 8 | bug | high | `src/Basis/AbelianBasis.jl:648,664,681,698,705` | Non-default `dtype` crashes int/Gosper construction (`Is=Int[]`, `base::Int64`) → MethodError | Thread `dtype` into `Is`/`dgt`/`change!`/`base` arg and fallback `_abelian_select` |
| 9 | bug | high | `src/Schmidt.jl:105-113,142` | `entropy()` silently ignores `cutoff` for Rényi α∉{0,1}; noise singular values corrupt result for 0<α<1 | Give `renyi_entropy` a `cutoff` and forward from line 111 |
| 10 | bug | high | `src/Basis/AbelianBasis.jl:759-761` | 2-arg `index(B,dgt)` mutates shared `B.G` odometer → not thread-safe, violates contract | `index(B,dgt)=index(B,dgt,_shallow_workspace(B.G))` |

High-leverage optimizations and extensions follow in Sections 4–5; the next tier (`expm-1`, `operator-1`, `abelian-4/5`, `linearmap-1`) is covered there.

## 3. Bugs & correctness

### Critical

**`abelian-1` — Fixed N + inversion generator off half-filling drops states** (`src/Basis/AbelianBasis.jl:614-625`, root `:615`)
`use_gosper` only checks `base==2 && 0≤Ndigits≤L`, ignoring whether any generator carries an inversion mask. The Gosper candidate list contains only fixed-weight integers, but `check_min_int`/`shift_canonical_int` apply the inversion (weight-w → weight-(L-w)), so canonical representatives of the other weight are never enumerated. Verified: `basis(L=6,N=2,z=1)` gives dim 5 + 5 = 10 vs the true 15, and the z-resolved spectrum even contains spurious eigenvalues. The docstring says "N and z are only compatible at half-filling" but nothing enforces it. **Fix:** compute `has_inv=any(any,G.inv)`; require `(!has_inv || 2*Ndigits==L*(base-1))` for `use_gosper`, else full-scan; or error when N fixed + inv present + `2N≠L`.

**`productstate-1` — silently wrong for every reduced/symmetry-adapted basis** (`src/ToolKit.jl:107-113`)
Advertised to work for "symmetry-reduced bases alike" but only correct for `AbstractOnsiteBasis`. Three failures: (1) does `I=index(B)[2]; s[I]=1`, discarding the orbit phase/normalization coeff; (2) `zeros(size(B,1))` is Float64, cannot store a complex phase; (3) an out-of-sector config returns `(0,1)` so it writes unit amplitude on basis vector #1 with no error. **Fix:** `c,I=index(B,collect(v)); iszero(c)&&error(...); s=zeros(eltype(B),size(B,1)); s[I]=c`. (Also resolves `productstate-2` by using a local buffer.)

### High

**Integer-width overflow family** (`bug-1`, `bug-2`; related low-severity `bug-3`, `abelian-7`, uncertain `parityflip-3`)
- `bug-1` (`AbstractBasis.jl:297-310`; stored `ProjectedBasis.jl:202-203,242`): the `index` accumulator `N=zero(T)` binds to the *base* type, so `dtype=Int32` with base>2, or base=2 at L≥31, silently stores wrapped/negative representatives used by all later binary searches. Verified `ProjectedBasis(Int32;L=20,N=0,base=3)` stores `Int32[-808182895]`. In-sector states still round-trip (same overflowing function on both ends), so there is no crash — only objectively wrong indices and possible sort-order/collision corruption.
- `bug-2` (`ProjectedBasis.jl:36-54`, lines 41,44,47): `_binary_fixed_weight_indices` uses `one(T)<<L`; at L=63 Int64 this is `typemin`, so the n==L branch returns a negative index and (more commonly) the loop guard `state<limit` fails immediately, yielding an **empty** basis for any 1≤n≤L-1. Verified `ProjectedBasis(L=63,N=1,base=2)` → size (0,0).
- **Shared fix:** decouple accumulator width from base; guard `L < 8*sizeof(T)-1` and add a constructor-time check `base^L-1 ≤ typemax(dtype)` that errors loudly. This single guard covers bug-1, bug-2, and the low-severity `bug-3` (`AbstractBasis.jl:255-256`, `TensorBasis` size overflows to 0/negative for L≥63) and `abelian-7` (`AbelianBasis.jl:554-569`, Gosper hard-capped at L≤62 while Benes compiles for L≤64).

**`trans-1` — `order()` returns L instead of orbit size L/a** (`TranslationalBasis.jl:352`, `TranslationalParityBasis.jl:227`, `TranslationalFlipBasis.jl:210`)
The orbit has `len=L/a` copies, but `order()` returns L. Consumed by `basis_embedding` (`LinearMap.jl:62,74`) and `mps2vec` (`ITensorsKit.jl:90`). Invisible at a=1; for a>1 every symmetrizer/embedding/MPS amplitude is off by 1/a. Verified: symmetrizer squared-norm ratio is 0.25 at a=2, 0.0625 at a=4 (=1/a²); the fix restores isometry. **Fix:** `order(b::TranslationalBasis)=ncycle(b)`, `2*ncycle(b)` for the parity/flip variants. Add an a>1 isometry regression test.

**`parityflip-2` — projective predicate checked only on representative** (`ParityBasis.jl:46`, `FlipBasis.jl:44`, `ParityFlipBasis.jl:75`)
`judge.F(dgt)` is evaluated once on the representative, then the whole orbit is accepted. Valid only if `f` is symmetry-invariant. The half-filling asserts (`FlipBasis.jl:77`, `ParityFlipBasis.jl:92`) cover only N-sectors, not arbitrary `f`. Verified: `FlipBasis(L=4,f=dgt->dgt[1]==0)` gives p=±1 dims 8+8=16 vs the true f-constrained dim 8; `index` maps f-false flipped states into the basis with no error. **Fix:** after picking the representative, also require `f` on the change!-decoded flipped/reversed members, else `(false,0.0)`.

**`parityflip-1` — ComplexF64 eltype despite ±1 real phases** (`ParityBasis.jl:105`, `FlipBasis.jl:100`, `ParityFlipBasis.jl:121,131`)
None override `eltype`, so they inherit `ComplexF64` (`AbstractBasis.jl:121`); their sibling `TranslationParityBasis.jl:25` already overrides to Float64. Consequences: `index` is type-unstable (`Union{Tuple{ComplexF64,Int}, Tuple{Float64,Int}}`) for Parity/Flip, and real symmetric Hamiltonians assemble as ComplexF64 (2× memory, complex BLAS/eig). **Fix:** add `eltype(...)=Float64`; change the `ParityFlipBasis` literal `0.0` to `zero(eltype(b))`. (Same theme as `abelian-5` below.)

**`operator-2` — `mul!` never uses the cached sparse matrix** (`Operator.jl:426-451,480,505`)
Only `*(opt,m)` and `mul(opt,m)` consult `_cached_sparse`; none of the `mul!` methods do, and the matrix `mul!` is single-threaded. The `sparse!` docstring promises 10–1000× for "applying H to a block of eigenvectors / inside an iterative solver," but block/LOBPCG solvers calling the standard `mul!` get neither cached SpMM nor threading. *(Scope note: KrylovKit single-vector `eigsolve` applies on vectors via `*(opt,v)`, which is matrix-free by design even when cached; the real gap is the matrix `mul!`.)* **Fix:** in matrix `mul!` (and the 5-arg form), `S=_cached_sparse(opt); S!==nothing && return mul!(target,S,m,α,β)` before the matrix-free fallback. This also resolves `operator-5`.

**`abelian-2` — non-default dtype crashes int/Gosper paths** (`AbelianBasis.jl:648,664,681,698,705`)
`_abelian_select_int`/`_abelian_select_gosper` hardcode `Is=Int[]` and `change!(...;base=Int(2))`, so `AbelianBasis{Ti}` cannot unify `Vector{Int32}` dgt with `Vector{Int64}` I → MethodError. Verified `basis(Int32;L=6,N=3,k=0)` throws. Verifiers note the fix must be broader than the suggestion: also convert the `base` argument to `dtype` at the `:639` constructor call (the `B::Ti` field), and fix the fallback `_abelian_select` (`:724`) which infers `Int` from the `1:base^L` range. **Fix:** thread `dtype` through `Is`/`dgt`/`change!`, the `base` arg, and the fallback range.

**`abelian-3` — 2-arg `index(B,dgt)` mutates shared `B.G`** (`AbelianBasis.jl:759-761`)
`index(B,dgt)=index(B,dgt,B.G)` passes the shared odometer; `shift_canonical_int`/`shift_canonical!` `init!` and mutate `g.s` in place. Concurrent calls (even with per-thread dgt) race and return wrong pairs — verified 10574 mismatches under 4 threads. The internal mul!/schmidt paths are safe (they use `_shallow_workspace`); only the public 2-arg API and `index_nocheck` (`LinearMap.jl:15`) are exposed. **Fix:** `index(B,dgt)=index(B,dgt,_shallow_workspace(B.G))`.

**`schmidt-1` — `entropy()` ignores `cutoff` for Rényi α∉{0,1}** (`Schmidt.jl:105-113,142`)
`renyi_entropy(s,α)` neither accepts nor applies `cutoff`, unlike the Shannon and Rényi-0 branches. For 0<α<1, a noise singular value of 1e-8 perceptibly corrupts the result. Verified: `entropy(s,α=0.5,cutoff=1e-10)` returns 1.3862944111 vs `log(4)=1.3862943611`, unchanged when `cutoff` is tightened. **Fix:** `renyi_entropy(s,α;cutoff=1e-20)=log(sum(si^α for si in s if si>cutoff))/(1-α)`; forward from line 111.

**`expm-1` — `expm`/`expv` fixed low-order Taylor, no scaling-and-squaring** (`ToolKit.jl:58-91`)
Degree-10 Taylor truncation with no argument-norm reduction; silently returns >100% error for `‖A‖≳5` (ordinary for an ED Hamiltonian × time step). Verified relerr 1.14 at `‖A‖=5`, ~1944 at 10. **Severity note:** verifiers split high/medium — the functions are not exported, have no internal callers, and the docs route real dynamics to `timeevolve`/Lindblad, so the blast radius is a user reaching for a labeled quick utility. **Fix:** add scaling-and-squaring (`s=ceil(log2(opnorm(A)))`, Taylor of `A/2^s`, square s times), or at minimum a loud `opnorm`-based warning/error above the convergence threshold.

### Medium

**`schmidt-2` — `renyi_entropy` assumes normalized input** (`Schmidt.jl:142,105-113`)
`log(sum(s.^α))/(1-α)` is the Rényi entropy only when `sum(s)==1`; unnormalized input (reachable via the documented `EDKit.entropy` and the Lindblad `entropy(::DensityMatrix)`, since `densitymatrix` does not normalize) silently gives wrong values (even negative). The α≈1 instability claim is overstated — the numerator →0 linearly with (1-α), so error is negligible until `|1-α|≲1e-13`. **Fix:** normalize over above-cutoff entries inside `renyi_entropy`; optionally widen the von-Neumann dispatch to `abs(α-1)<1e-8`.

**`operator-1` — sparse cache keyed on `objectid(opt)`, goes stale on in-place mutation** (`Operator.jl:306,354-358,370`)
`objectid` of an immutable struct is identity-derived, not content-derived; mutating `opt.M[i].nzval` in place leaves `*(opt,m)`/`mul(opt,m)` returning the old cached matrix with no error (verified). The public API never mutates `opt.M` in place (scalar ops build fresh Operators), so triggering requires reaching into internals. **Fix:** store `(content_hash, S)` and miss on hash mismatch, or loudly document that `sparse!` snapshots contents.

**`bug-2` (ITensors) — `tebd4` output incompatible with `tebd_n!`; no −i in time step** (`TEBD.jl:163-191,185`, docstring `153,158`)
Line 185 is `exp(h[l]*sweep.τ*τ)` with no −im, while the docstring calls τ the "physical time step." A user passing a real `dt` silently gets non-unitary imaginary-time evolution (the dev notebook passes `-1im*dt`). **Correction to original framing:** the "no working stepper" claim is overstated — `tebd4`'s flat `Vector{ITensor}` is consumable by ITensors' own `apply(G,psi)`; it just bypasses `tebd_n!`'s SVD truncation. **Fix:** clarify the docstring (operator is `exp(τ·h)`, real-time needs `τ=-i·dt`) and/or add an EDKit-level driver with truncation control, or an `evolve=:real` keyword.

### Low (real but narrow)

- **`bug-4`** (`AbstractBasis.jl:360`): `change!(dgt::AbstractVector{T}, ind::T; base::T=2)` ties buffer/index/base to one type → MethodError on a mismatched scratch buffer. Loosen to `ind::Integer`/`base::Integer` with internal conversion (matches the docstring at `:346`). *(The second cited location `:392` is mischaracterized — `ind` there is already `::Integer`; only `base` is coupled.)*
- **`parityflip-4`** (`ParityBasis.jl`, `FlipBasis.jl`, `ParityFlipBasis.jl`): only reduced bases lacking `copy`; `copy(b::ParityBasis)` is a MethodError. This is an API-consistency gap, **not** a thread-safety bug (Schmidt/mul already use local `similar` buffers). Add the three `copy` methods mirroring the others.
- **`operator-6`** (`Operator.jl:63-77`): `iszero(mats[i])&&continue` can produce a silent zero-term Operator; the real cost is the O(num²) `findfirst` dedup. Use a `Dict{Vector{Int},Int}` for O(1) dedup.
- **`gapratio-1`** (`ToolKit.jl:15-29`): no sorted check; exactly-degenerate levels give a silently-wrong 0/0→1.0. Add a `sorted` keyword + guard and an optional degeneracy `tol`; document single-sector use. (`meangapratio` inherits it.)
- **`qim-1`** (`QIM.jl:22-31`): docstring claims `½⟨{hᵢ,hⱼ}⟩−⟨hᵢ⟩⟨hⱼ⟩` but code computes `Re⟨hᵢψ|hⱼψ⟩−…`, equal only for Hermitian `hᵢ`. Result is still real-symmetric (the `Symmetric` wrap is safe), just not the documented form for non-Hermitian input. Document the Hermitian-only assumption + `ishermitian` guard, or compute the explicit anticommutator form.

### Uncertain (verifiers split)

- **`te-1`** (`TimeEvolution.jl:392-393`, mirrored `LindbladEvolution.jl:384`): defect-monitor sample count hard-capped at 513, so for `ω·τ` large the monitor can be sampled coarser than its own Nyquist requirement and `_max_valid_interval` may accept a tol-violating interval. One verifier confirmed (worked example `ω≈20,τ≈100` exceeds even bare Nyquist at `ω·τ≳1608`); another found empirically the monitor is a smooth band-limited trig sum with a 4× safety factor, so no tol-violating peak is missed in practice (worst-case 1.22× peak underestimate at `ω·τ≈8000`). **Assessment:** real error-control hole at extreme bandwidth×time, but likely low-severity in realistic regimes. If addressed: shrink the candidate interval so `ω·τ` stays resolvable, or derive the cap from `m_max` and warn.
- **`parityflip-3`** (`FlipBasis.jl:20`, `ParityFlipBasis.jl:21`): `M::Int`/`MAX` ignore `dtype`. One verifier confirmed silent corruption on the default Int64 path at L≥63; another showed `dtype=Int128` fails *loudly* (InexactError at construction), and the base=2 overflow is a generic `base^L` issue unrelated to the field type. **Assessment:** at most a cosmetic consistency tweak (make `M::Ti`), unreachable in practice (2^63 states).
- **`o1`** (`FermionOperator.jl:116-121`): the "dense 2^span matrix" mechanism is **refuted** — the local matrix is already built sparse (the comment documents this was done to avoid the dense OOM). The real residual is the `2^span+1` Int64 `colptr` (~2 GiB at span 28), so long-range single-bond terms still scale memory with 2^span. The endpoint-only + parity-coefficient fix is directionally valid but the severity/OOM framing is wrong.
- **`optimization-2`** (`AbstractBasis.jl:483-488`), **`operator-7`** (`Operator.jl:183-189`), **`o2`** (`FermionOperator.jl:305-310`), **`te-3`** (`TimeEvolution.jl:309-323`): all technically accurate but immaterial — cold construction code or branches dominated by far larger costs (binary-search cache misses, O(m²·N) Lanczos, JW kron chains). Treat as parked micro-notes, not actionable without profiling. *(Note: `operator-7`'s reasoning "size(opt,2) alone exceeds 1e6 for L≥20" is wrong — the L=20 half-filled dim is 184,756; it's the product that overshoots.)*

## 4. Optimizations

**`abelian-5` — AbelianBasis always ComplexF64 even for real sectors** (`AbelianBasis.jl:584-590`; falls back to `AbstractBasis.jl:121`)
For k=0, p=±1, z=±1 sectors all characters are real, yet operator assembly/diagonalization is forced complex (2× memory, ~2× FLOPs). **Benefit:** halves cost on the most common (k=0, parity/flip) workflows. **Change:** the existing second type parameter already records this — `eltype(::AbelianBasis{Ti,Tg}) where {Ti,Tg} = Tg<:Real ? Float64 : ComplexF64` (mirrors `TranslationalBasis.jl:31`). Pairs with `parityflip-1` as one eltype-correctness sweep.

**`optimization-1` — base>2 fixed-charge enumeration wastes 70–95% of iterations** (`ProjectedBasis.jl:195-205`; same pattern `TranslationalBasis.jl:180-182`)
`multiexponents(L,N)+filter` discards out-of-range compositions: measured 69.5% waste at L6/N6/base3, 82.8% at L8/N8, growing to ~95%+. **Benefit:** ~10–20× faster construction of large qudit bases (one-time cost, not per-multiply). **Change:** use a bounded mixed-radix composition generator that never emits an out-of-range digit. *(The secondary "predicate fast-path" point is minor — the branch already short-circuits `f` via `isnothing(f)||`; the decode+index is inherent to base>2.)*

**`linearmap-1` — symmetrizer / DoubleBasis-action scale as dense base^L** (`LinearMap.jl:59-77,167-193,210-212`)
`basis_embedding` allocates a dense ComplexF64 `base^L × dim` matrix though the embedding has ≤1 nonzero per row; unusable around L~18-20 (16 MB/column at L=20). This is a user-facing utility path (the standard `trans_inv_operator(M,inds,DoubleBasis(...))` workflow never calls it), but the docstring recommends the `symmetrizer(DoubleBasis(B,Bfull))` pattern. **Benefit:** follows reduced-basis size instead of base^L. **Change:** build via a sparse triplet accumulator and parameterize element type on `promote_type(eltype(B1),eltype(B2))`. *(Note: Parity/Flip have ComplexF64 eltype, so the real-savings cases are Tensor/Projected and the Float64 reduced bases.)*

**`abelian-4` — per-call heap allocation in the index hot path** (`AbelianBasis.jl:514,439-440`)
`shift_canonical_int`/`shift_canonical!` allocate `ms=g.s[:]` (and `tmp=similar(dgt)`) per call, once per nonzero per column. **Benefit:** removes GC pressure from large sparse builds (the headline "~570 B/call" is overstated — single-generator is 64 B/call; it grows with the generator count). **Change:** carry preallocated `ms`/`tmp` in the workspace already threaded through `colmn!`→`index`.

**`operator-4` — 5-arg `mul!` single-threaded while `mul()` threads** (`Operator.jl:429-433,447-451,463-495`)
The standard `mul!` entry point leaves all but one core idle. **Benefit:** near-linear speedup for matrix/block `mul!` (single-vector Lanczos benefits less and needs a min-work guard like `sparse`'s `min_cols_per_thread`). **Change:** factor the threaded accumulate-and-reduce body out of `mul` into a helper using per-thread buffers of `eltype(target)`. Subsumed by the `operator-2` fix.

Lower-value (`trans-3`, `parityflip-5`, `expm-2`, `schmidt-3`, `opt-1` ITensors): real but minor — see the JSON. `parityflip-5` (`ParityBasis.jl:96-107`) specializing reflection to a base=2 bit-reversal is a constant-factor shave, not the claimed O(1). `opt-1` (`PauliBasis.jl`) is cold dev tooling with no `src/` callers.

## 5. Extensions to increase versatility

**Spinful / multi-species fermions** (`e1`, `SpinlessFermionBasis.jl:13-17`; `fermion_operator` requires base==2 at `:155`) — *enables Hubbard/t-J/multi-orbital.* Add a 2L-mode basis with a fixed global JW ordering, `(site,spin)→mode` constructors, `(N↑,N↓)` sectors, and a Hubbard helper. Documented as deferred (`docs/.../fermions.md:240`).

**Symmetry-resolved fermion basis** (`e2`, `FermionOperator.jl:170-179`) — *enables per-momentum/parity fermion ED instead of full fixed-N.* The guard correctly rejects naive JW-on-permute embedding; add a momentum-resolved basis baking the boundary sign `(-1)^(N-1)` into `R` and the allowed-k selection, then relax the guard.

**Steady-state Lindblad solver** (`lindblad-6`, `LindbladEvolution.jl:138-151`) — *computes ρ_ss directly instead of hand-built fixed points or long propagation.* Add `steadystate(A)`: Krylov eigensolver for the eigenvalue nearest 0 (consuming the matrix-free `LiouvillianMap`), reshape, Hermitize, trace-normalize. Note: needs a KrylovKit/Arpack dependency (not currently in `Project.toml`) or a dense fallback.

**Backward / two-sided time evolution** (`te-4`, `TimeEvolution.jl:556-577`) — *enables OTOC / Heisenberg-picture correlators.* The machinery is time-reversible for Hermitian H; only the `cis(-τλ)` sign (`:333,347`) and the `t<0`/monotonicity guards (`:565,588,590,621,623`) block it. Allow signed τ with `|τ|` in the monitor and a per-cache direction flag.

**Mixed local dimensions** (`extension-1`, `AbstractBasis.jl:229-258,297-404`) — *enables cavity-QED/impurity/lattice-gauge (qubit+truncated-boson, spin-½+spin-1).* All digit↔index math bakes in one scalar base across every basis. Add a mixed-radix basis type (`base::AbstractVector`) with Horner/divrem per-site, leaving scalar fast paths untouched.

**RDM / mutual-information helper** (`schmidt-4`, `Schmidt.jl:144-168,189-201`) — *enables observables on A, negativity, I(A:B).* Every `schmidt` already returns S; add `rdm(v,Ainds,b;...)=S*S'` and `mutual_information(...)`, reusing all per-basis dispatch.

Lower-priority: qudit Pauli/Lindblad generalization (`ext-3`, `PauliBasis.jl:11-29` — `GELLMANN` already ships at `Operator.jl:830-841`); n-body fermion strings + `:right` convention (`e3`, `FermionOperator.jl:81`); a hermitian-hopping helper avoiding the L=2 double-count (`e4`, `FermionOperator.jl:287-292`); per-jump rates / ħ in Lindblad (`lindblad-7` — mostly a documentation note); SFF/number-variance/unfolding diagnostics (`diagnostics-1`, `ToolKit.jl:1-41`). Missing test coverage: Pauli MPS/MPO round-trips and `density_expect` (`ext-1` — note a round-trip test exists in `test/TensorTest.jl:198-221` but is **not** run by `Pkg.test()`), and TEBD Trotter accuracy vs exact (`ext-2`).

## 6. Cross-cutting themes

1. **Integer width is decoupled from the basis everywhere it matters.** `dtype` is a documented constructor argument but the index accumulator (`bug-1`), Gosper enumerator (`bug-2`), `TensorBasis` size (`bug-3`), Abelian construction (`abelian-2`), Abelian Gosper ceiling (`abelian-7`), and Flip `M` fields (`parityflip-3`) all hardcode `Int`/base-type arithmetic and overflow silently. **Systemic win:** one constructor-time check `base^L-1 ≤ typemax(dtype)` (errors loudly) plus threading `dtype` into the enumerators kills most of this family.

2. **`ComplexF64` eltype carried by real-phase bases.** Parity/Flip/ParityFlip (`parityflip-1`) and AbelianBasis (`abelian-5`) inherit `ComplexF64` despite real ±1 phases, causing type instability and 2× memory/FLOPs; `linearmap-1`'s embedding hardcodes ComplexF64 too. **Systemic win:** one eltype-correctness sweep mirroring `TranslationalBasis.jl:31`/`TranslationalParityBasis.jl:25`.

3. **Documented preconditions are unenforced.** Half-filling + inversion (`abelian-1`, critical), symmetric predicate `f` (`parityflip-2`), sorted/non-degenerate spectra (`gapratio-1`), Hermitian QIM operators (`qim-1`), and the `tebd4` time-step sign (ITensors `bug-2`) are all stated in docstrings but never checked — turning user mistakes into silent wrong answers. **Systemic win:** add guards/asserts at the documented boundaries.

4. **The cached-sparse path is under-wired.** `sparse!` is advertised as a transparent accelerator but is unreachable from `mul!` (`operator-2`), absent from in-place paths (`operator-5`), single-threaded where `mul` threads (`operator-4`), and stale on mutation (`operator-1`). **Systemic win:** route `mul!` through `_cached_sparse` once and reuse the threaded helper.

5. **`order()`/normalization for a>1 unit cells** (`trans-1`) is a self-contained correctness gap exercised only by an existing documented feature.

## 7. Suggested next steps

1. **Fix the two criticals first:** `abelian-1` (guard `use_gosper` / error off half-filling) and `productstate-1` (capture coeff, use `eltype(B)`, error on zero-norm). Both are silent-wrong-answer bugs on documented APIs.
2. **Do the integer-width sweep** (`bug-1/2/3`, `abelian-2/7`): add a constructor-time overflow check that errors loudly, and thread `dtype` through the Abelian enumerators and `base` argument. Add regression tests for `dtype=Int32` base>2 and L≥63 base=2.
3. **Do the eltype-correctness sweep** (`parityflip-1`, `abelian-5`): add `eltype=Float64`/`Tg`-based overrides; verify real Hamiltonians assemble real and `index` is type-stable.
4. **Fix `trans-1`** (`order()=ncycle`) with an a>1 isometry regression test, and **`parityflip-2`** (re-check `f` on orbit members).
5. **Wire the cached-sparse path into `mul!`** (`operator-2`, subsuming `operator-4`/`operator-5`); add an in-place 5-arg `mul!` delegating to `mul!(target,S,m,α,β)`.
6. **Fix `abelian-3`** (workspace in 2-arg `index`) and the entropy `cutoff`/normalization bugs (`schmidt-1`, `schmidt-2`).
7. **Lower-risk cleanups:** `expm-1` scaling-and-squaring (or a loud warning), `operator-1` content-hash key, the `tebd4` docstring/sign, `bug-4` signature, the three `copy` methods, `gapratio-1`/`qim-1` guards.
8. **Investigate `te-1`** (defect-monitor cap) — reproduce at large `ω·τ` to decide between a real fix and a documentation/efficiency note.
9. **Then schedule extensions** by demand: spinful fermions and the steady-state solver are the highest-value additions, followed by symmetry-resolved fermions, backward evolution, and the RDM helper.
10. **Add the missing correctness tests** (Pauli round-trips into `Pkg.test()`, TEBD vs exact, a>1 isometry, `dtype` overflow) so these regressions are caught in CI — which currently runs no Julia tests.