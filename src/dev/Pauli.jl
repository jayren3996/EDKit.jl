include("../EDKit.jl")
using Main.EDKit, LinearAlgebra, Test, Profile, ITensors, ITensorMPS

function get_el(psi, el)
    s = siteinds(psi)
    V = ITensor(1.0)
    for j in eachindex(psi)
        V *= (psi[j]*state(s[j], el[j]))
    end
    scalar(V)
end

function randH(sites)
    L = length(sites)
    os = OpSum()
    for i in 1:L 
        os += randn(), "Sx", i, "Sx", mod(i,L)+1
        os += randn(), "Sx", i, "Sy", mod(i,L)+1
        os += randn(), "Sx", i, "Sz", mod(i,L)+1
        os += randn(), "Sy", i, "Sx", mod(i,L)+1
        os += randn(), "Sy", i, "Sy", mod(i,L)+1
        os += randn(), "Sy", i, "Sz", mod(i,L)+1
        os += randn(), "Sz", i, "Sx", mod(i,L)+1
        os += randn(), "Sz", i, "Sy", mod(i,L)+1
        os += randn(), "Sz", i, "Sz", mod(i,L)+1
    end
    MPO(os, sites)
end

L = 5
s = siteinds("S=1/2", L)
S = siteinds("Pauli", L)

H = randH(s)

begin
    ψ = randomMPS(s)
    Hψ = apply(H, ψ)
    A = mps2pmps(Hψ, S)

    O = mps2pmps(ψ, S)
    PH = mpo2pmpo(H, S)
    B = apply(PH, O)
end;
mps2vec(A)-mps2vec(B) |> norm


