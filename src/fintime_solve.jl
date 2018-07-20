
function max_out_rate(chain)
    max = 0
    for src in nodes(chain.state_graph)
        total_rate = 0.0
        for edge in out_edges(chain.state_graph, src)
            rate = chain.trans_rates[edge]
            total_rate += rate
        end
        if total_rate > max
            max = total_rate
        end
    end
    return max
end

function trans_prob_matrix(chain, unif_rate)
    rows = Vector{Int}()
    cols = Vector{Int}()
    vals = Vector{Float64}()
    for src in nodes(chain.state_graph)
        total_rate = 0.0
        for edge in out_edges(chain.state_graph, src)
            dst = dst_node(chain.state_graph, edge)
            rate = chain.trans_rates[edge]
            total_rate += rate
            push!(rows, dst)
            push!(cols, src)
            push!(vals, rate / unif_rate)
        end
        push!(rows, src)
        push!(cols, src)
        self_prob = 1.0 - total_rate / unif_rate
        @assert self_prob >= 0
        push!(vals, self_prob)
    end
    dim = node_count(chain.state_graph)
    sparse(rows, cols, vals, dim, dim)
end

function trans_rate_matrix(chain)
    rows = Vector{Int}()
    cols = Vector{Int}()
    vals = Vector{Float64}()
    for src in nodes(chain.state_graph)
        total_rate = 0.0
        for edge in out_edges(chain.state_graph, src)
            dst = dst_node(chain.state_graph, edge)
            rate = chain.trans_rates[edge]
            total_rate += rate
            push!(rows, dst)
            push!(cols, src)
            push!(vals, rate)
        end
        push!(rows, src)
        push!(cols, src)
        push!(vals, -total_rate)
    end
    dim = node_count(chain.state_graph)
    sparse(rows, cols, vals, dim, dim)
end

export fintime_solve_prob, trans_rate_matrix, trans_prob_matrix, state_prob, state_cumtime

include("poisson_trunc.jl")

function state_prob(sol::Vector{Float64}, state::Integer)
    return sol[state]
end

function state_cumtime(sol::Vector{Float64}, state::Integer)
    return sol[state]
end

function fintime_solve_prob(chain::ContMarkovChain, init_prob::SparseVector, time::Real; unif_rate_factor=1.05, tol=1e-6, ss_check_interval=10)
    @assert unif_rate_factor >= 1.0
    unif_rate = max_out_rate(chain) * unif_rate_factor
    P = trans_prob_matrix(chain, unif_rate)

    prob = fill(0.0, state_count(chain))
    for i in 1:length(init_prob.nzind)
        prob[init_prob.nzind[i]] = init_prob.nzval[i]
    end

    ltp, rtp = poisson_trunc_point(unif_rate * time, tol)
    weight, wsum = poisson_term_weight(unif_rate * time, ltp, rtp, tol)

    checkpoint = copy(prob)
    prob_old = copy(prob)

    for k in 1:ltp
        prob, prob_old = prob_old, prob
        A_mul_B!(prob, P, prob_old)
        if k % ss_check_interval == 0
            checkpoint -= prob
            diff = maximum(abs, checkpoint)
            checkpiont = prob
            if  diff < tol
                return prob
            end
        end
    end
    summed_term = 0.0
    prob_t = fill(0.0, length(prob))
    for k in ltp:rtp
        p_term = weight[k - ltp + 1] / wsum
        summed_term += p_term
        prob_t += p_term * prob

        prob, prob_old = prob_old, prob
        A_mul_B!(prob, P, prob_old)
        if (k - ltp + 1) % ss_check_interval == 0
            checkpoint -= prob
            diff = maximum(abs, checkpoint)
            checkpiont = prob
            if diff < tol
                prob_t += prob * (1.0 - summed_term)
                return prob_t
            end
        end
    end

    return prob_t
end

