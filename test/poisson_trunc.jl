@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

using MarkovChains
@testset "poisson_trunc" begin
    λ = 20.0
    tol = 1e-6
    ltp, rtp, probs = poisson_trunc_point(λ, tol)
    @test ltp == 3
    @test rtp == 45
    @test abs(1.0 - sum(probs)) < tol
    λ = 400.0
    ltp, rtp, probs = poisson_trunc_point(λ, tol)
    @test abs(1.0 - sum(probs)) < tol
end
