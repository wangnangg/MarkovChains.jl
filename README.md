# MarkovChains.jl

[![Build Status](https://travis-ci.org/wangnangg/MarkovChains.jl.svg?branch=master)](https://travis-ci.org/wangnangg/MarkovChains.jl)
[![Coverage Status](https://coveralls.io/repos/github/wangnangg/MarkovChains.jl/badge.svg?branch=master)](https://coveralls.io/github/wangnangg/MarkovChains.jl?branch=master)

This pacakge provides functions to solve continuous/discrete time Markov chains for state probablities or accumulated sojourn times at a certain time point, including time infinity.

# Example
## continuous time Markov chains

The following example is about solving a 4 states birth-death chain at time infinity.

```julia
importall MarkovChains
chain = ContMarkovChain()
n0 = add_state!(chain)
n1 = add_state!(chain)
n2 = add_state!(chain)
n3 = add_state!(chain)
add_transition!(chain, n0, n1, 1.0) #transition from n0 to n1 with rate = 1.0
add_transition!(chain, n1, n2, 1.0)
add_transition!(chain, n2, n3, 1.0)
add_transition!(chain, n3, n2, 3.0)
add_transition!(chain, n2, n1, 2.0)
add_transition!(chain, n1, n0, 1.0)
init_prob = sparsevec([1], [1.0])
res = inftime_state(chain, init_prob)
@show cumtime = collect(map(state -> state_cumtime(res, state), n0:n3))
# cumtime = [inf, inf, inf, inf]
@show prob = collect(map(state -> state_prob(res, state), n0:n3))
# prob = [0.375, 0.375, 0.1875, 0.0625]
```

## discrete time Markov chains
