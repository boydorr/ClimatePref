using PhyloNetworks
using Distributions
using JLD
using JuliaDB
using JuliaDBMeta
using DataFrames
using StatsBase
using Statistics
using GLM

tree = readTopology("Qian2016.tree")
tipnames = tipLabels(tree)
tip_names = join.(split.(tipnames, "_"), " ")
top_common_names = JLD.load("Common_species_names.jld", "spp_names")
missing_species = setdiff(tip_names, top_common_names)
for i in eachindex(missing_species)
    deleteleaf!(tree, join(split(missing_species[i], " "), "_"))
end

temp_trait = ParamsBM(25, 10.0)
sim1 = simulate(tree, temp_trait)
trait1 = DataFrame(temp = sim1[:Tips], tipNames = tipLabels(sim1))

lambda = phyloNetworklm(@formula(temp ~ 1), trait1, tree, model="lambda")
lambda_estim(lambda)

samp = sample(1:5_000, 250, replace = false)
trait1[:temp_shift] = trait1[:temp]
trait1[samp, :temp_shift] .+= 10.0

lambda = phyloNetworklm(@formula(temp ~ 1), trait1, tree, model="lambda")
lambda_estim(lambda)

# Add individual records
trait_individual = by(trait1, :tipNames, temp = :temp => x -> rand(Normal(x[1], 0.1), 100))
trait_individual = DataFrames.flatten(trait_individual, :temp)
samp = sample(1:nrow(trait_individual), round(Int64, nrow(trait_individual) * 0.05), replace = false)
trait_individual[!, :temp_shift] = trait_individual[!, :temp]
trait_individual[samp, :temp_shift] .+= 10.0
trait_mean = combine(groupby(trait_individual, :tipNames), :temp => mean)

lambda = phyloNetworklm(@formula(temp_mean ~ 1), trait_mean, tree, model="lambda")
lambda_estim(lambda)

# Overall loop
reps = 10; shift = 10.0; ind_var = 1.0; prop_to_shift = [0.01, 0.05, 0.1, 0.2]
lambda_mat = zeros(Float64, reps, length(prop_to_shift))
lambda_shift_mat = zeros(Float64, reps, length(prop_to_shift))

for r in 1:reps
    for p in eachindex(prop_to_shift)
        temp_trait = ParamsBM(25, 10.0)
        sim1 = simulate(tree, temp_trait)
        trait1 = DataFrame(temp = sim1[:Tips], tipNames = tipLabels(sim1))
        trait_individual = by(trait1, :tipNames, temp = :temp => x -> rand(Normal(x[1], ind_var), 100))
        trait_individual = DataFrames.flatten(trait_individual, :temp)
        trait_mean = by(trait_individual, :tipNames, :temp => mean)
        lambda = phyloNetworklm(@formula(temp_mean ~ 1), trait_mean, tree, model="lambda")
        lambda_mat[r, p] = lambda_estim(lambda)

        samp = sample(1:nrow(trait_individual), round(Int64, nrow(trait_individual) * prop_to_shift[p]), replace = false)
        trait_individual[!, :temp_shift] = trait_individual[!, :temp]
        trait_individual[samp, :temp_shift] .+= shift
        trait_mean = by(trait_individual, :tipNames, :temp_shift => mean)
        lambda = phyloNetworklm(@formula(temp_shift_mean ~ 1), trait_mean, tree, model="lambda")
        lambda_shift_mat[r, p] = lambda_estim(lambda)
    end
    println("Repeat $r")
end