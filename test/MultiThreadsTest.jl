#-------------------------------------------------------------------------------------------------------------------------
# Test Translational PXP
#-------------------------------------------------------------------------------------------------------------------------
@testset "Multi-threads" begin
    L, k, p = 28, 0, 1
    mat = begin
        P = Diagonal([1, 1, 1, 0, 1, 1, 0, 0])
        X = [0 1; 1 0]
        P * kron(I(2), X, I(2)) * P
    end
    pxpf(v::Vector{Int}) = all(v[i]==0 || v[mod(i, length(v))+1]==0 for i=1:length(v))
    println("--------------------------------------")
    print("Single-threads:")
    @time bs = translationalbasis(pxpf, k, L, threaded=false)
    print("Multi-threads :")
    @time bm = translationalbasis(pxpf, k, L, threaded=true)
    @test bs.I == bm.I
    @test norm(bs.R-bm.R) ≈ 0.0
end

#-------------------------------------------------------------------------------------------------------------------------
# Test Translational Parity PXP
#-------------------------------------------------------------------------------------------------------------------------
@testset "Multi-threads" begin
    L, k, p = 28, 0, 1
    mat = begin
        P = Diagonal([1, 1, 1, 0, 1, 1, 0, 0])
        X = [0 1; 1 0]
        P * kron(I(2), X, I(2)) * P
    end
    pxpf(v::Vector{Int}) = all(v[i]==0 || v[mod(i, length(v))+1]==0 for i=1:length(v))
    println("--------------------------------------")
    print("Single-threads:")
    @time bs = translationparitybasis(pxpf, k, p, L, threaded=false)
    print("Multi-threads :")
    @time bm = translationparitybasis(pxpf, k, p, L, threaded=true)
    @test bs.I == bm.I
    @test norm(bs.R-bm.R) ≈ 0.0
end
