# Focused binary_search benchmark - preventing dead code elimination.
using BenchmarkTools
BenchmarkTools.DEFAULT_PARAMETERS.samples = 30
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 2

function binary_search_current(list::AbstractVector{<:Integer}, i::Integer)
    l::Int = 1
    r::Int = length(list)
    c::Int = (l + r) ÷ 2
    while true
        t = @inbounds list[c]
        (i < t) ? (r = c - 1) : (i > t) ? (l = c + 1) : break
        (l > r) ? (c = 0; break) : (c = (l + r) ÷ 2)
    end
    c
end

@inline function binary_search_stdlib(list::AbstractVector{<:Integer}, i::Integer)
    idx = searchsortedfirst(list, i)
    @inbounds (idx <= length(list) && list[idx] == i) ? idx : 0
end

for dim in [246, 810, 3432, 12870]
    sorted = sort(unique(rand(1:dim*10, dim*2)))[1:dim]
    targets = [rand(sorted) for _ in 1:1000]

    # Force result to be used to prevent dead code elimination
    b1 = @benchmark begin
        s = 0
        for t in $targets
            s += binary_search_current($sorted, t)
        end
        s
    end

    b2 = @benchmark begin
        s = 0
        for t in $targets
            s += binary_search_stdlib($sorted, t)
        end
        s
    end

    # Verify correctness
    for t in targets
        @assert binary_search_current(sorted, t) == binary_search_stdlib(sorted, t) "Mismatch at $t"
    end

    t1 = median(b1).time
    t2 = median(b2).time
    println("dim=$dim: current=$(BenchmarkTools.prettytime(t1))  " *
            "stdlib=$(BenchmarkTools.prettytime(t2))  " *
            "speedup=$(round(t1/t2, digits=2))x")
end
