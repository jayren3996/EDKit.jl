#----------------------------------------------------------------------------------------------------
# TEBD1
#----------------------------------------------------------------------------------------------------
export tebd_n!
"""
    tebd_n!(ψ, G1, G2; cutoff=1e-14, maxdim=30)

n-site TEBD.
"""
function tebd_n!(ψ::MPS, G1::Vector, G2::Vector; cutoff::Real=1e-14, maxdim::Integer=30)
    L = length(ψ)
    s = siteinds(ψ)
    n = length(inds(G1[1])) ÷ 2
    orthogonalize!(ψ, 1)
    #------------------------------------------------------------
    # First block
    block = ITensor(1.0)
    for i in 1:n
        block *= ψ[i]
    end
    wf = G1[1] * block |> noprime! 
    U, S, V = svd(wf, s[1]; cutoff, maxdim)
    ψ[1], T = U, S * V
    link = commonind(ψ[1], T) 
    #------------------------------------------------------------
    # moving right
    for i in 2:L-n
        wf = G1[i] * (T * ψ[i+n-1]) |> noprime! 
        U, S, V = svd(wf, link, s[i]; cutoff, maxdim)
        ψ[i], T = U, S * V
        link = commonind(ψ[i], T)
    end
    #------------------------------------------------------------
    # dealing with the last block
    wf = apply(G2[1], G1[L-n+1]) * (T * ψ[L]) |> noprime!
    U, S, V = svd(wf, link, s[L-n+1:L-1]...; cutoff, maxdim)
    T, ψ[L] = U * S, V
    #------------------------------------------------------------
    # moving left
    for j in 2:L-n
        i = L-j-n+2
        wf = G2[j] * (ψ[i] * T) |> noprime! 
        link = commonind(ψ[i-1], ψ[i])
        U, S, V = svd(wf, commonind(ψ[i-1], ψ[i]), s[i:i+n-2]...; cutoff, maxdim)
        T, ψ[i+n-1] = U * S, V
    end
    #------------------------------------------------------------
    # dealing with the first block
    wf = G2[L-n+1] * (ψ[1] * T) |> noprime!
    for i in n:-1:2
        U, S, V = svd(wf, s[1:i-1]...; cutoff, maxdim)
        wf, ψ[i] = U * S, V
    end
    ψ[1] = wf
    #------------------------------------------------------------
    # restore gauge
    ψ.llim = 0
    ψ.rlim = 2
    normalize!(ψ)
end

#----------------------------------------------------------------------------------------------------
# TEBD4
#----------------------------------------------------------------------------------------------------
struct TrotterSweep
    τ::Float64
    step::Int64
    offset::Int64
end
#----------------------------------------------------------------------------------------------------
const Order4 = begin # Bartel & Zhang, Annals of Physics 418, 168165 (2020) Eq. 32
    τ2 = 1/(4-4^(1/3))
    τ2b = 1/(4-4^(1/3))/2
    τ3 = (1-4/(4-4^(1/3)))
    τ4 = (1-3/(4-4^(1/3)))/2
    TS = TrotterSweep
    [
        TS(τ2b, 2, 1),
        TS(τ2, -2, 0),
        TS(τ2,  2, 1),
        TS(τ2, -2, 0),
        TS(τ4,  2, 1),
        TS(τ3, -2, 0),
        TS(τ4,  2, 1),
        TS(τ2, -2, 0),
        TS(τ2,  2, 1),
        TS(τ2, -2, 0),
        TS(τ2b, 2, 1)
    ]
end
#----------------------------------------------------------------------------------------------------
const Order4n3 = begin # Bartel & Zhang, Annals of Physics 418, 168165 (2020) Eq. 53
    u = BigFloat("0.095968145884398107402")
    q1 = BigFloat("0.43046123580897338276")
    r1 = BigFloat("-0.075403897922216340661")
    q2 = BigFloat("-0.12443549678124729963")
    r2 = BigFloat("0.5") - (u + r1)
    u1 = u
    v1 = q1 - u1
    u2 = r1 - v1
    v2 = q2 - u2
    u3 = r2 - v2 # = 1/2 -q2 - q1

    TS = TrotterSweep
    [
        TS(u,    3, 0),
        TS(u,   -3, 1),
        TS(q1  , 3, 2),
        TS(v1  ,-3, 1),
        TS(r1  , 3, 0),
        TS(u2  ,-3, 1),
        TS(q2  , 3, 2),
        TS(v2  ,-3, 1),
        TS(r2  , 3, 0),
        TS(u3  ,-3, 1),
        TS(2*u3, 3, 2),
        TS(u3  ,-3, 1),
        TS(r2  , 3, 0),
        TS(v2  ,-3, 1),
        TS(q2  , 3, 2),
        TS(u2  ,-3, 1),
        TS(r1  , 3, 0),
        TS(v1  ,-3, 1),
        TS(q1  , 3, 2),
        TS(u,   -3, 1),
        TS(u,    3, 0)
    ]
end
#----------------------------------------------------------------------------------------------------
export tebd4
function tebd4(
    h::Vector{<:AbstractMatrix},
    s::Vector{<:Index},
    τ::Number
)
    n = round(Int64, log(space(s[1]), size(h,1)) )
    SW = if n == 2
        Order4
    elseif n == 3
        Order4n3
    else
        "n should be 2 or 3, instead got $n."
    end

    N = length(s)
    gates = ITensor[]
    for sweep in SW
        step = abs(sweep.step)
        ran = (1+sweep.offset):step:N-step+1
        ran = sweep.step < 0 ? reverse(ran) : ran
        
        for l in ran
            U = exp(h[l] * sweep.τ * τ)
            G = mat2op(U, s[l:l+step-1]...)
            push!(gates, G)
        end
    end
    gates
end



