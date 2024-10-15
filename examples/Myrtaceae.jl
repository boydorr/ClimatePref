# SPDX-License-Identifier: BSD-2-Clause

# 20. Phylogenetic analysis of plant climate preferences at different taxonomic levels
@everywhere begin
    using PhyloNetworks
    using GLM
    using JuliaDB
    using JuliaDBMeta
    using Unitful
    using ClimatePref
    using ClimatePref.Unitful
    using StatsBase
    using JLD
    using OnlineStats
    using Unitful.DefaultSymbols
    using Statistics
    using DataFrames
    using Unitful
    using Unitful.DefaultSymbols
    using OnlineStats
    using StatsBase
    using ClimatePref
    using CSV
end

# Load all records with taxonomic information
cera_taxo = JuliaDB.load("CERA_taxo")

# Filter for species in the Myrtaceae family
family = filter(c -> c.family == "Myrtaceae", cera_taxo)
family = collect(family)
JLD.save("Myrtaceae.jld", "family", family)
family = @transform family {tmin = ustrip(:tmin),
                            tmax = ustrip(:tmax),
                            tmean = ustrip(:tmean),
                            trng = ustrip(:trng),
                            stl1 = ustrip(:stl1mean),
                            stl2 = ustrip(:stl2mean),
                            stl3 = ustrip(:stl3mean),
                            stl4 = ustrip(:stl4mean),
                            swvl1 = ustrip(:swvl1mean),
                            swvl2 = ustrip(:swvl2mean),
                            swvl3 = ustrip(:swvl3mean),
                            swvl4 = ustrip(:swvl4mean),
                            ssr = ustrip(:ssrmean),
                            tp = ustrip(:tpmean)}
CSV.write("Myrtaceae.csv", family)

# Load tree
tree = readTopology("Qian2016.tree")
tipnames = tipLabels(tree)
tip_names = join.(split.(tipnames, "_"), " ")
spp_names = collect(JuliaDB.select(family, :species))
cross_species = spp_names ∩ tip_names

# Filter for species actually in the tree
family_fil = filter(g -> g.species in cross_species, family)
missing_species = setdiff(tip_names, cross_species)
for i in eachindex(missing_species)
    deleteleaf!(tree, join(split(missing_species[i], " "), "_"))
end

# Group gbif data by Species ID and get mean and percentiles of data
phylo_traits = @groupby family_fil :species {tmin = mean(ustrip(:tmin)),
                                             tmax = mean(ustrip(:tmax)),
                                             tmean = mean(ustrip(:tmean)),
                                             trng = mean(ustrip(:trng)),
                                             stl1 = mean(ustrip(:stl1mean)),
                                             stl2 = mean(ustrip(:stl2mean)),
                                             stl3 = mean(ustrip(:stl3mean)),
                                             stl4 = mean(ustrip(:stl4mean)),
                                             swvl1 = mean(ustrip(:swvl1mean)),
                                             swvl2 = mean(ustrip(:swvl2mean)),
                                             swvl3 = mean(ustrip(:swvl3mean)),
                                             swvl4 = mean(ustrip(:swvl4mean)),
                                             ssr = mean(ustrip(:ssrmean)),
                                             tp = mean(ustrip(:tpmean))}

# Add in tip names to data and save
phylo_traits = @transform phylo_traits {tipNames = join.(split.(:species, " "),
                                                         "_")}

# Convert to dataframe
dat = DataFrame(collect(phylo_traits))

# Fit lambda models and save
lambdas = fitLambdas(tree, dat, names(dat)[2:(end - 1)], 0.1)
JLD.save("Lambdas_family_Myrtaceae.jld", "lambdas", lambdas)
