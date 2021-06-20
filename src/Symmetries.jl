"""
    cyclebits!(dgt::AbstractVector{<:Integer})

Right-shift digits (in place).
"""
function cyclebits!(dgt::AbstractVector{<:Integer})
    p = dgt[end]
    for i = length(dgt):-1:2
        dgt[i] = dgt[i-1]
    end
    dgt[1] = p
    dgt
end

"""
    translation_check(dgt::AbstractVector{<:Integer}, I0::Integer, base::Integer)

Check if the digits is of minimum index under cycling.

Inputs:
-------
- `dgt  : Digits.
- `I0`  : Index of the `dgt`.
- `base`: Base.

Outputs:
--------
- `Q`: Whether `dgt` is minimum.
- `R`: Periodicity of the `dgt`, return 0 if !Q.
"""
function translation_check(dgt::AbstractVector{<:Integer}, I0::Integer, base::Integer)
    Q, R = true, length(dgt)
    cyclebits!(dgt)
    for i=1:length(dgt)-1
        In = index(dgt, base=base)
        if In < I0
            Q = false
            R = 0
            change!(dgt, I0, base=base)
            break
        elseif In == I0
            R = i
            break
        end
        cyclebits!(dgt)
    end
    Q, R
end

"""
    translation_check_min!(dgt::AbstractVector{<:Integer}, Im::Integer, base::Integer)

Given an index `Im`, cycle the digits to see:

1. Whether `Im` is the smallest index.
2. Whether the shifted digits have the same index as `Im`.
3. If the shifted digits have the same index as `Im`, find the cycling number.

Check if the digits is of minimum index under cycling.

Inputs:
-------
- `dgt  : Digits.
- `Im`  : Minimum index.
- `base`: Base.
- `R`   : (Optional) Periodicity of `dgt`, `R = length(dgt)` by default.

Outputs:
--------
- `Qm`: Whether `Im` is minimum index.
- `Qs`: Whether `dgt` can be shifted to have index `Im`.
- `M` : Cycling number by which `dgt` is shifted to have index `Im`.
"""
function translation_check_min!(dgt::AbstractVector{<:Integer}, Im::Integer, base::Integer, R::Integer=length(dgt))
    Qm, Qs, M = true, false, 0
    for i=0:R-1
        In = index(dgt, base=base)
        if In < Im
            Qm, Qs = false, false
            break
        elseif In == Im && !Qs
            Qs = true
            M = i
            break
        end
        cyclebits!(dgt)
    end
    Qm, Qs, M
end

"""
translation_index(dgt::AbstractVector{<:Integer}, base::Integer)

Given a digits, return the minimum index under shifting.

Inputs:
-------
- `dgt  : Digits.
- `base`: Base.

Outputs:
--------
- `Im`: Minimum index.
- `M` : Cycling number by which `dgt` is shifted to have index `Im`.
"""
function translation_index(dgt::AbstractVector{<:Integer}, base::Integer)
    I0 = index(dgt, base=base)
    Im, M = I0, 0
    cyclebits!(dgt)
    for i=1:length(dgt)-1
        In = index(dgt, base=base)
        if In == I0
            break
        elseif In < Im
            Im, M = In, i
        end
        cyclebits!(dgt)
    end
    Im, M
end

"""
    spinflip(v::AbstractVector{<:Integer}, base::Integer)

Flip spins Sz on each site.
"""
function spinflip(v::AbstractVector{<:Integer}, base::Integer)
    vf = Vector{eltype(v)}(undef, length(v))
    for i = 1:length(vf)
        vf[i] = base - v[i] - 1
    end
    vf
end

"""
    translation_parity_check(parity, dgt::AbstractVector{<:Integer}, I0::Integer, base::Integer)

Check whether a state is a valid representing state, and whether it is reflect-translation-invariant.

Inputs:
-------
- `parity`: Parity function, can be spatio-reflection or spin-reflection.
- `dgt`   : Digits.
- `I0`    : Index of the `dgt`.
- `base`  : Base.

Outputs:
--------
- `Qm`: Whether `dgt` is a representing vector.
- `Qp`: Whether `dgt` is reflect-translation-invariant.
- `R` : Minimum R that Tᴿ⋅|dgt⟩ = |dgt⟩.
- `T` : Minimum M that Tᴹ⋅P⋅|dgt⟩ = |dgt⟩, 0 if it is not reflect-translation-invariant.
"""
function translation_parity_check(parity, dgt::AbstractVector{<:Integer}, I0::Integer, base::Integer)
    Q, R = translation_check(dgt, I0, base)
    Qm, Qp, M = if Q
        translation_check_min!(parity(dgt), I0, base, R)
    else
        false, false, 0
    end
    Qm, Qp, R, M
end

"""
    translation_parity_index(parity, dgt::AbstractVector{<:Integer}, base::Integer)

Return the index of a digits, and also tell whether the minimum index is in the other parity sector, and return the number of shifts.

Inputs:
-------
- `parity`: Parity function, can be spatio-reflection or spin-reflection.
- `dgt`   : Digits.
- `base`  : Base.

Outputs:
--------
- `Qp`: Whether the minimum index is in the other parity sector.
- `Im`: Minimum index.
- `M` : Minimum M that Tᴹ⋅|dgt⟩ = |Im⟩ or Tᴹ⋅P⋅|dgt⟩ = |Im⟩.
"""
function translation_parity_index(parity, dgt::AbstractVector{<:Integer}, base::Integer)
    Im1, T1 = translation_index(dgt, base)
    Im2, T2 = translation_index(parity(dgt), base)
    Qp, Im, M = if Im2 < Im1 
        true, Im2, T2
    else
        false, Im1, T1
    end
    Qp, Im, M
end
