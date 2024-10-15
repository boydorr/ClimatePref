# SPDX-License-Identifier: BSD-2-Clause

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

lambda = phyloNetworklm(@formula(temp~1), trait1, tree, model = "lambda")
lambda_estim(lambda)

samp = sample(1:5_000, 250, replace = false)
trait1[:temp_shift] = trait1[:temp]
trait1[samp, :temp_shift] .+= 10.0

lambda = phyloNetworklm(@formula(temp~1), trait1, tree, model = "lambda")
lambda_estim(lambda)

# Add individual records
trait_individual = by(trait1, :tipNames,
                      temp = :temp => x -> rand(Normal(x[1], 0.1), 100))
trait_individual = DataFrames.flatten(trait_individual, :temp)
samp = sample(1:nrow(trait_individual),
              round(Int64, nrow(trait_individual) * 0.2), replace = false)
trait_individual[!, :temp_shift] = trait_individual[!, :temp]
trait_individual[samp, :temp_shift] .+= 10.0
trait_mean = by(trait_individual, :tipNames, :temp => mean)

lambda = phyloNetworklm(@formula(temp_mean~1), trait_mean, tree,
                        model = "lambda")
lambda_estim(lambda)

# Overall loop
reps = 10;
shift = 10.0;
ind_var = 1.0;
prop_to_shift = [0.05, 0.1, 0.25, 0.5];
lambda_mat = zeros(Float64, reps, length(prop_to_shift))
lambda_shift_mat = zeros(Float64, reps, length(prop_to_shift))

for r in 1:reps
    for p in eachindex(prop_to_shift)
        # Read tree and clip
        tree = readTopology("Qian2016.tree")
        for i in eachindex(missing_species)
            deleteleaf!(tree, join(split(missing_species[i], " "), "_"))
        end
        # Simulate trait
        temp_trait = ParamsBM(25, 10.0)
        sim1 = simulate(tree, temp_trait)
        trait1 = DataFrame(temp = sim1[:Tips], tipNames = tipLabels(sim1))
        # Simulate individual records per trait
        trait_individual = by(trait1, :tipNames,
                              temp = :temp => x -> rand(Normal(x[1], ind_var),
                                                        1_000))
        trait_individual = DataFrames.flatten(trait_individual, :temp)
        trait_mean = by(trait_individual, :tipNames, :temp => mean)
        lambda = phyloNetworklm(@formula(temp_mean~1), trait_mean, tree,
                                model = "lambda")
        lambda_mat[r, p] = lambda_estim(lambda)
        println(lambda_estim(lambda))
        # Shift proportion of records
        for i in unique(trait_individual[!, :tipNames])
            pos = findall(trait_individual[!, :tipNames] .== i)
            samp = sample(pos, round(Int64, length(pos) * prop_to_shift[p]),
                          replace = false)
            trait_individual[samp, :temp] .+= rand(Normal(shift, 0.5 * shift))
        end
        trait_mean = by(trait_individual, :tipNames, :temp => mean)
        lambda_shift = phyloNetworklm(@formula(temp_mean~1), trait_mean, tree,
                                      model = "lambda")
        lambda_shift_mat[r, p] = lambda_estim(lambda_shift)
        println(lambda_estim(lambda_shift))
    end
    println("Repeat $r")
end

JLD.save("lambda_shift.jld", "lambda", lambda_shift_mat)

lambda_shift_mat = JLD.load("examples/lambda_shift.jld", "lambda")

mean_shift = mean(lambda_shift_mat, dims = 1)[1, :]
using Plots
plot(prop_to_shift, mean_shift, grid = false, label = "Shifted",
     ylab = "Pagel's lambda", xlab = "Proportion records shifted",
     bottom_margin = 20.0 * Plots.px)
hline!([mean_shift[1]], label = "Raw")
Plots.pdf("Lambda_shift.pdf")
