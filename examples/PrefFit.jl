# Import libraries.
using Turing, Plots, StatsPlots, Random
pyplot()

# Iterate from having seen 0 observations to 100 observations.
Ns = 0:1000;

# Draw data from a Bernoulli distribution, i.e. draw heads or tails.
Random.seed!(12)
data = rand(Normal(25.0, 10.0), last(Ns))
histogram(data)

# Declare our Turing model.
@model temppref(y) = begin
    # Our prior belief about the probability of heads in a coin.
    m ~ Uniform(20, 30)
    v ~ Uniform(5, 15)
    # The number of observations.
    N = length(y)
    for n in 1:N
        # Heads or tails of a coin are drawn from a Bernoulli distribution.
        y[n] ~ Normal(m, v)
    end
end;

# Settings of the Hamiltonian Monte Carlo (HMC) sampler.
iterations = 1000
ϵ = 0.05
τ = 10

# Start sampling.
#chain = sample(temppref(data), HMC(iterations, ϵ, τ));

Nsamples = 2000
Nadapt = 1000
δ = .85
num_chains = 2

chains = mapreduce(c -> sample(temppref(data), NUTS(Nsamples,Nadapt, δ)), chainscat, 1:num_chains)
subset = chains[Nadapt:Nsamples, :, :]
plot(subset[[:m, :v]])
histogram(subset[[:m, :v]])
describe(subset)
gelmandiag(subset)
gewekediag(subset)
