#----------------------------------------------------------------------------------------------------
# TEBD
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

