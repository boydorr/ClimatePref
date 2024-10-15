# SPDX-License-Identifier: BSD-2-Clause

# 12. Reconstruct climate preferences using ten-fold cross validation

using JuliaDB
using Unitful
using ClimatePref
using ClimatePref.Unitful
using PhyloNetworks
using JLD
using DataFrames
using Statistics
using Distances

### RAW DATA ###

mins = [
    197.0K,
    197.0K,
    197.0K,
    0K,
    197.0K,
    197.0K,
    197.0K,
    197.0K,
    0.0m^3,
    0.0m^3,
    0.0m^3,
    0.0m^3,
    0.0J / m^2,
    0.0m
]
maxs = [
    320.0K,
    320.0K,
    320.0K,
    80K,
    320.0K,
    320.0K,
    320.0K,
    320.0K,
    1.0m^3,
    1.0m^3,
    1.0m^3,
    1.0m^3,
    3.0e7J / m^2,
    0.1m
]

# Read in tip data and tree
phylo_traits = JuliaDB.load("Phylo_traits")
tree = readTopology("Qian2016.tree")

# Cut tree down to common species
top_common_names = JLD.load("Common_species_names.jld", "spp_names")
tipnames = tipLabels(tree)
tip_names = join.(split.(tipnames, "_"), " ")
missing_species = setdiff(tip_names, top_common_names)
for i in eachindex(missing_species)
    deleteleaf!(tree, join(split(missing_species[i], " "), "_"))
end

# Filter traits to common species
phylo_traits_filter = filter(p -> p.tipNames in join.(split.(top_common_names,
                                                             " "), "_"),
                             phylo_traits)
dat = DataFrame(collect(phylo_traits_filter))

# Reconstruct values for different climate vars
climate_vars = [
    :tmin,
    :tmax,
    :tmean,
    :trng,
    :stl1,
    :stl2,
    :stl3,
    :stl4,
    :swvl1,
    :swvl2,
    :swvl3,
    :swvl4,
    :ssr,
    :tp
]
recon = map(climate_var -> fitMissings(tree, dat, climate_var, false),
            climate_vars)

# Convert result to dataframe and save
recon_mat = Array{Union{Float64, Missing}}(undef, 5000, length(recon))
for i in 1:length(recon)
    recon_mat[:, i] .= recon[i]
end
recon_dat = DataFrame(recon_mat)
names!(recon_dat, climate_vars)

tipName = [n.name for n in tree.leaf]
recon_dat[:tipName] = tipName

JLD.save("AncRecon_raw.jld", "traits", recon_dat)

recon_dat = JLD.load("AncRecon_raw.jld", "traits")
anc_corr_raw = [cor(recon_dat[c], dat[c]) for c in climate_vars]
JLD.save("Corr_raw.jld", "corr", anc_corr_raw)

# Randomise tips and fit again

recon = map(climate_var -> fitMissings(tree, dat, climate_var, true),
            climate_vars)
recon_mat = Array{Union{Float64, Missing}}(undef, 5000, length(recon))
for i in 1:length(recon)
    recon_mat[:, i] .= recon[i]
end
recon_dat = DataFrame(recon_mat)
names!(recon_dat, climate_vars)

tipName = [n.name for n in tree.leaf]
recon_dat[:tipName] = tipName

JLD.save("RandAncRecon_raw.jld", "traits", recon_dat)

anc_corr_raw = [cor(recon_dat[c], dat[c]) for c in climate_vars]
JLD.save("Corr_rand_raw.jld", "corr", anc_corr_raw)

### EFFORT ADJUSTED DATA ###

# Read in tip data and tree
phylo_traits = JuliaDB.load("Phylo_traits_adjust")
tree = readTopology("Qian2016.tree")

# Cut tree down to common species
top_common_names = JLD.load("Common_species_names.jld", "spp_names")
tipnames = tipLabels(tree)
tip_names = join.(split.(tipnames, "_"), " ")
missing_species = setdiff(tip_names, top_common_names)
for i in eachindex(missing_species)
    deleteleaf!(tree, join(split(missing_species[i], " "), "_"))
end

# Filter traits to common species
phylo_traits_filter = filter(p -> p.tipNames in join.(split.(top_common_names,
                                                             " "), "_"),
                             phylo_traits)
dat = DataFrame(collect(phylo_traits_filter))

# Reconstruct values for different climate vars
climate_vars = [
    :tmin,
    :tmax,
    :tmean,
    :trng,
    :stl1,
    :stl2,
    :stl3,
    :stl4,
    :swvl1,
    :swvl2,
    :swvl3,
    :swvl4,
    :ssr,
    :tp
]
recon = map(climate_var -> fitMissings(tree, dat, climate_var, false),
            climate_vars)

# Convert result to dataframe and save
recon_mat = Array{Union{Float64, Missing}}(undef, 5000, length(recon))
for i in 1:length(recon)
    recon_mat[:, i] .= recon[i]
end
recon_dat = DataFrame(recon_mat)
names!(recon_dat, climate_vars)

tipName = [n.name for n in tree.leaf]
recon_dat[:tipName] = tipName

JLD.save("AncRecon_adjust.jld", "traits", recon_dat)

anc_corr_adjust = [cor(recon_dat[c], dat[c]) for c in climate_vars]
JLD.save("Corr_adjust.jld", "corr", anc_corr_adjust)

# Randomise tips and fit again

recon = map(climate_var -> fitMissings(tree, dat, climate_var, true),
            climate_vars)
recon_mat = Array{Union{Float64, Missing}}(undef, 5000, length(recon))
for i in 1:length(recon)
    recon_mat[:, i] .= recon[i]
end
recon_dat = DataFrame(recon_mat)
names!(recon_dat, climate_vars)

tipName = [n.name for n in tree.leaf]
recon_dat[:tipName] = tipName

JLD.save("RandAncRecon_adjust.jld", "traits", recon_dat)

anc_corr_adjust = [cor(recon_dat[c], dat[c]) for c in climate_vars]
JLD.save("Corr_rand_adjust.jld", "corr", anc_corr_adjust)

### CLIMATE- EFFORT ADJUSTED DATA ###

# Read in tip data and tree
phylo_traits = JuliaDB.load("Phylo_traits_adjust2")
tree = readTopology("Qian2016.tree")

# Cut tree down to common species
top_common_names = JLD.load("Common_species_names.jld", "spp_names")
tipnames = tipLabels(tree)
tip_names = join.(split.(tipnames, "_"), " ")
missing_species = setdiff(tip_names, top_common_names)
for i in eachindex(missing_species)
    deleteleaf!(tree, join(split(missing_species[i], " "), "_"))
end

# Filter traits to common species
phylo_traits_filter = filter(p -> p.tipNames in join.(split.(top_common_names,
                                                             " "), "_"),
                             phylo_traits)
dat = DataFrame(collect(phylo_traits_filter))

# Reconstruct values for different climate vars
climate_vars = [
    :tmin,
    :tmax,
    :tmean,
    :trng,
    :stl1,
    :stl2,
    :stl3,
    :stl4,
    :swvl1,
    :swvl2,
    :swvl3,
    :swvl4,
    :ssr,
    :tp
]
recon = map(climate_var -> fitMissings(tree, dat, climate_var, false),
            climate_vars)

# Convert result to dataframe and save
recon_mat = Array{Union{Float64, Missing}}(undef, 5000, length(recon))
for i in 1:length(recon)
    recon_mat[:, i] .= recon[i]
end
recon_dat = DataFrame(recon_mat)
names!(recon_dat, climate_vars)

tipName = [n.name for n in tree.leaf]
recon_dat[:tipName] = tipName

JLD.save("AncRecon_adjust2.jld", "traits", recon_dat)

anc_corr_adjust2 = [cor(recon_dat[c], dat[c]) for c in climate_vars]
JLD.save("Corr_adjust2.jld", "corr", anc_corr_adjust2)

# Randomise tips and fit again

recon = map(climate_var -> fitMissings(tree, dat, climate_var, true),
            climate_vars)
recon_mat = Array{Union{Float64, Missing}}(undef, 5000, length(recon))
for i in 1:length(recon)
    recon_mat[:, i] .= recon[i]
end
recon_dat = DataFrame(recon_mat)
names!(recon_dat, climate_vars)

tipName = [n.name for n in tree.leaf]
recon_dat[:tipName] = tipName

JLD.save("RandAncRecon_adjust2.jld", "traits", recon_dat)

anc_corr_adjust2 = [cor(recon_dat[c], dat[c]) for c in climate_vars]
JLD.save("Corr_rand_adjust2.jld", "corr", anc_corr_adjust2)

### PLOT CORR RESULTS AS HEATMAP ###
using JLD
using Plots
import Plots.px
pyplot()

files = ["Corr_raw.jld", "Corr_adjust.jld", "Corr_adjust2.jld"]
corrs = map(files) do f
    return JLD.load(f, "corr")
end
x = ["Raw", "Effort", "Effort + Climate"]
y = [
    "tmin",
    "tmax",
    "tmean",
    "stl1",
    "stl2",
    "stl3",
    "swvl1",
    "swvl2",
    "swvl3",
    "swvl4",
    "ssr",
    "tp"
]
corrs = hcat(corrs...)
h = heatmap(x, y, corrs, seriescolor = :RdBu, clim = (-1, 1), layout = (1, 2),
            subplot = 1, title = "Imputed data", size = (1000, 700),
            bottom_margin = 20px, left_margin = 20px, right_margin = 20px,
            top_margin = 30px, tickfontsize = 12, colorbar = :none)

files = ["Corr_rand_raw.jld", "Corr_rand_adjust.jld", "Corr_rand_adjust2.jld"]
corrs = map(files) do f
    return JLD.load(f, "corr")
end
x = ["Raw", "Effort", "Effort + Climate"]
y = [
    "tmin",
    "tmax",
    "tmean",
    "stl1",
    "stl2",
    "stl3",
    "swvl1",
    "swvl2",
    "swvl3",
    "swvl4",
    "ssr",
    "tp"
]
corrs = hcat(corrs...)
h = heatmap!(x, y, corrs, seriescolor = :RdBu, clim = (-1, 1), subplot = 2,
             title = "Randomised \n imputed data", tickfontsize = 12)
png(h, "Correlation_heatmap.png")
