function find_rtp(λ, tol, m)
    if λ < 400.0
        λ = 400.0
    end
    a_λ = (1.0 + 1.0 / λ) * exp(1.0 / 16.0) * √(2.0)

    r = 0
    k = 3.0
    while true
        d_k_λ = 1.0 / (1.0 - exp((-2.0 / 9.0) * (k * √(2.0 * λ) + 1.5)))
        if a_λ * d_k_λ * exp(-k * k / 2.0) / (k * √(2.0 * π)) < tol / 2.0
            r = Int(ceil(m + k * √(2.0 * λ) + 1.5))
            break
        end
        k += 1.0
    end
    return r, k
end

function poisson_trunc_point(λ, tol)
    underflow_th = 1e-30
    overflow_th = 1e30

    m = floor(λ)
    rtp, k_r = find_rtp(λ, tol, m) 
    ltp = 0
    if λ < 25.0
        return ltp, rtp
    end

    b_λ = (1.0 + 1.0 / λ) * exp(1.0 / (8.0 * λ))
    k = 1.0
    while true
        if b_λ * exp(-k * k / 2.0) / (k * √(2.0 * π)) < tol / 2.0
            break
        end
        k += 1.0
    end
    ltp = Int(floor(m - k * √(λ) - 1.5))
    if ltp < 0
        ltp = 0
    end

    #overflow check omitted for now
    bound = 0.0
    c_m = (1.0 / √(2.0 * π * m)) * exp((m - λ - 1) / (12.0 * m))
    if λ >= 400.0
        k_hat = k_r * √(2.0) + 3.0 / (2.0 * √(λ))
        bound = c_m * exp(-(k_hat + 1) * (k_hat + 1) / 2.0)
        if overflow_th * bound / (1.0e10 * (rtp - ltp)) < underflow_th
            throw(OverflowError())
        end
    end

    k_tilde = k + 3.0 / (2.0 * √(λ))

    if k_tilde > 0 && k_tilde <= √(λ) / 2.0
        bound = c_m * exp(-k_tilde * k_tilde / 2.0 - -k_tilde * k_tilde * k_tilde / (3.0 * √(λ)))
    elseif k_tilde <= √(Float64(m + 1))
        bound = max(c_m * (1.0 - k_tilde / √(Float64(m + 1)))^(k_tilde * √(Float64(m + 1))), exp(-λ))
        λ
        exp(-λ)
    end

    bound
    overflow_th * bound / (1.0e10 * (rtp - ltp))
    underflow_th
    if bound != 0.0 && overflow_th * bound / (1.0e10 * (rtp - ltp)) < underflow_th
        throw(OverflowError())
    end
    return ltp, rtp
end


function poisson_term_weight(λ, ltp::Integer, rtp::Integer, tol)
    underflow_th = 1e-30
    overflow_th = 1e30

    weight = Vector{Float64}(rtp - ltp + 1)
    m = Int(floor(λ))
    weight[m - ltp + 1] = overflow_th / (1e10 * (rtp - ltp))
    for j in m:-1:ltp + 1
        weight[j - ltp] = (j / λ) * weight[j - ltp + 1]
    end

    if λ < 400.0
        if rtp > 600
            throw(OverflowError())
        end
        for j in m:rtp - 1
            qty = λ / (j + 1)
            if weight[j - ltp + 1] > underflow_th / qty
                weight[j - ltp + 2] = qty * weight[j - ltp + 1]
            else
                rtp = j
            end
        end
    else
        for j in m:rtp - 1
            weight[j - ltp + 2] = weight[j - ltp + 1] * (λ / (j + 1))
        end
    end

    s = 1
    t = length(weight)
    wsum =  0
    while s < t
        if weight[s] <= weight[t]
            wsum += weight[s]
            s += 1
        else
            wsum += weight[t]
            t -= 1
        end
    end
    wsum += weight[s]
    return weight, wsum
end