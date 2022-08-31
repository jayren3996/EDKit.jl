#---------------------------------------------------------------------------------------------------
# Quantum Circuits
#---------------------------------------------------------------------------------------------------
"""
    MonitoredCircuit{O1<:Operator, O2<:Operator, R<:Real}

Structure representing the monitored dynamics.

Properties:
-----------
- `H` : Effective non-Hermitian Hamiltonian.
- `L` : List of Lindblad jump operators.
- `dt`: Decimated time interval.
"""
struct MonitoredCircuit{O1<:Operator, O2<:Operator, R<:Real}
    H::O1
    L::Vector{O2}
    dt::R
end

export monitoredcircuit
"""
    monitoredcircuit(H, L, I, dt)

Constructor for `MonitoredCircuit`.

Inputs:
-------
- `H` : Hamiltonian. 
- `L` : List of matrices representing jump operators.
- `I` : List of vectors storing the sites on with jump operators act.
- `dt`: Decimated time interval.

Outputs:
--------
- `mc`: MonitoredCircuit.
"""
function monitoredcircuit(
    H::Operator, 
    L::Vector{<:AbstractMatrix}, 
    I::AbstractVector{<:AbstractVector{<:Integer}}, 
    dt::Real
)
    L = map(1:length(L)) do i
        H += operator(-0.5im * L[i]' * L[i], I[i], H.B) 
        operator(L[i], I[i], H.B)
    end
    MonitoredCircuit(H, L, dt)
end

"""
Action of `MonitoredCircuit` on vector.
"""
function *(mc::MonitoredCircuit, v::Vector)
    v = expv(mc.H, v, order=5, λ=-1im * mc.dt) |> normalize
    for l in mc.L 
        lv = l * v
        nlv = norm(lv)
        rand() < nlv^2 * mc.dt || continue
        v = lv / nlv
    end
    v
end

