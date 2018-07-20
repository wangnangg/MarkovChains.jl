@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

importall MarkovChains
@testset "fintime_solve" begin
    @testset "irreducible1" begin
        chain = ContMarkovChain()
        n0 = add_state!(chain)
        n1 = add_state!(chain)
        n2 = add_state!(chain)
        n3 = add_state!(chain)
        add_transition!(chain, n0, n1, 1.0)
        add_transition!(chain, n1, n2, 1.0)
        add_transition!(chain, n2, n3, 1.0)
        add_transition!(chain, n3, n2, 3.0)
        add_transition!(chain, n2, n1, 2.0)
        add_transition!(chain, n1, n0, 1.0)
        init_prob = sparsevec([1, 2], [0.1, 0.9])
        res = fintime_solve_prob(chain, init_prob, 0.0)
        @test res ≈ [0.1, 0.9, 0.0, 0.0]

        res = fintime_solve_prob(chain, init_prob, 1.0)
        @test res[1] > 0.1 && res[1] < 0.375
        @test res[2] < 0.9 && res[2] > 0.375
        @test res[3] < 0.1875
        @test res[4] < 0.0625

        res = fintime_solve_prob(chain, init_prob, 20.0)
        @test res ≈ [0.375, 0.375, 0.1875, 0.0625]
    end
end