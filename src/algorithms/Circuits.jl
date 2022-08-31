#---------------------------------------------------------------------------------------------------
# Quantum Circuits
#---------------------------------------------------------------------------------------------------
struct MonitoredCircuit{O1<:Operator, O2<:Operator, R<:Real}
    H::O1
    L::Vector{O2}
    dt::R
end

function monitoredcircuits(
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

