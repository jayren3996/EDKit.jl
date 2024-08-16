function fdaoe(s::AbstractVector, l::Integer, γ::Real)
    L = length(s)
    κ = exp(-2γ)
    #=----------------------------------------------------------
    Wᴵᴵ = I(l+1)⁺ ⊕ (∑|n⟩⟨n+1| + κ|l⟩⟨l|)⁻
    Wᶻᶻ = (∑|n⟩⟨n+1| + κ|l+1⟩⟨l+1|)⁺ ⊕ I(l)⁻ 
    Wˣˣ = Wʸʸ = [0 A; B 0], where
        A = I(l)⁺⁻ + κ|l+1⟩⁺⟨l|⁻
        B = ∑|n⟩⁺⟨n+1|⁻
    ----------------------------------------------------------=#
    Wi = zeros(2l+1, 2l+1, 4, 4)
    Wi[1:l+1, 1:l+1, 1, 1] = I(l+1)
    Wi[l+2:end, l+2:end, 1, 1] = diagm(1 => ones(l-1))
    Wi[end, end, 1, 1] = κ

    Wi[1:l+1, 1:l+1, 4, 4] = diagm(1 => ones(l))
    Wi[l+1, l+1, 4, 4] = κ 
    Wi[l+2:end, l+2:end, 4, 4] = I(l)

    A = diagm(l+1, l, ones(l))
    A[end, end] = κ 
    B = diagm(l, l+1, 1 => ones(l))
    Wi[1:l+1, l+2:end, 2, 2] = A
    Wi[1:l+1, l+2:end, 3, 3] = A
    Wi[l+2:end, 1:l+1, 2, 2] = B 
    Wi[l+2:end, 1:l+1, 3, 3] = B 
    #=----------------------------------------------------------
    L = ⟨1|
    R = ∑ₙ|n⟩
    ----------------------------------------------------------=#
    W1 = Wi[1,:,:,:]
    WL = sum(Wi[:,i,:,:] for i in axes(Wi, 2))

    links = [Index(2l+1, "Link,l=$i") for i in 1:L-1]
    D = MPO(L)
    D[1] = ITensor(W1, links[1], s[1]', s[1])
    for i in 2:L-1 
        D[i] = ITensor(Wi, links[i-1], links[i], s[i]', s[i])
    end
    D[L] = ITensor(WL, links[L-1], s[L]', s[L])
    D
end
