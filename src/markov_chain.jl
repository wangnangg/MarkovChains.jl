using .Graphs

export ContMarkovChain, add_state!, add_transition!, state_count

struct ContMarkovChain
    state_graph::Digraph
    trans_rates::Vector{Float64}
    function ContMarkovChain()
        new(Digraph(), Vector{Float64}())
    end
end

add_state!(chain::ContMarkovChain) = add_node!(chain.state_graph)

function add_transition!(chain::ContMarkovChain, src::Integer, dst::Integer, rate::Real)
    idx = add_edge!(chain.state_graph, src, dst)
    push!(chain.trans_rates, rate)
    return idx
end

state_count(chain) = node_count(chain.state_graph)
