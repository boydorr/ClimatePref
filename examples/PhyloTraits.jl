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
@everywhere using Unitful
@everywhere using Unitful.DefaultSymbols
@everywhere using OnlineStats
@everywhere using StatsBase

# Load tree
tree = readTopology("Qian2016.tree")
tipnames = tipLabels(tree)
tip_names = join.(split.(tipnames, "_"), " ")

# Load GBIF data
gbif = JuliaDB.load("CERA_JOIN_SIMPLE")
spp = collect(JuliaDB.select(gbif, :SppID))
numspp = unique(spp)

# Load Species names
spp_names = JLD.load("Species_names.jld", "spp_names")
spp_ids = JLD.load("Species_names.jld", "spp_ids")
sppdict = Dict(zip(spp_ids, spp_names))
iddict = Dict(zip(spp_names, spp_ids))
spp_names = [sppdict[i] for i in numspp]

# Get top 5000 most common species
cross_species = spp_names ∩ tip_names
cross_ids = [iddict[i] for i in cross_species]
sorted_counts = countmap(spp)
sorted_counts = filter(kv -> kv.first in cross_ids, sorted_counts)
sorted_counts = sort(sorted_counts, byvalue=true, rev=true)
top_common_ids = collect(keys(sorted_counts))[1:5000]
top_common_names = [sppdict[i] for i in top_common_ids]
#JLD.save("Common_species_names.jld", "spp_names", top_common_names)
#p = bar(collect(values(sorted_counts))[1:5000], grid = false)
#png(p, "GBIF_counts.png")

# Filter GBIF data for common species and delete from tree
gbif_fil = filter(g->g.SppID in cross_ids, gbif)
missing_species = setdiff(tip_names, top_common_names)
for i in eachindex(missing_species)
    deleteleaf!(tree, join(split(missing_species[i], " "), "_"))
end

### LAMBDAS FOR RAW DATA ###

# Group gbif data by Species ID and get mean and percentiles of data
phylo_traits = @groupby gbif_fil :SppID {tmin = mean(ustrip(:tmin)), tmin10 = percentile(ustrip(:tmin), 10), tmax = mean(ustrip(:tmax)), tmax90 = percentile(ustrip(:tmax), 90), tmean = mean(ustrip(:tmean)), trng = mean(ustrip(:trng)), stl1 = mean(ustrip(:stl1mean)), stl2 = mean(ustrip(:stl2mean)), stl3 = mean(ustrip(:stl3mean)), stl4 = mean(ustrip(:stl4mean)), swvl1 = mean(ustrip(:swvl1mean)), swvl2 = mean(ustrip(:swvl2mean)), swvl3 = mean(ustrip(:swvl3mean)), swvl4 = mean(ustrip(:swvl4mean)), ssr = mean(ustrip(:ssrmean)), tp = mean(ustrip(:tpmean))}

# Add in tip names to data and save
trait_ids = collect(JuliaDB.select(phylo_traits, :SppID))
new_cross_species = [sppdict[i] for i in trait_ids]
phylo_traits = pushcol(phylo_traits, :tipNames, join.(split.(new_cross_species, " "), "_"))
JuliaDB.save(phylo_traits, "Phylo_traits_new")

# Filter for common species
phylo_traits_filter = filter(p -> p.tipNames in join.(split.(top_common_names, " "), "_"), phylo_traits)

# Convert to dataframe
dat = DataFrame(collect(phylo_traits_filter))

# Fit lambda models and save
lambdas = fitLambdas(tree, dat, names(dat)[2:end-1])
JLD.save("Lambdas_common.jld", "lambdas", lambdas)


### LAMBDAS FOR EFFORT-ADJUSTED DATA ###

mins = [197.0K, 197.0K, 197.0K, 0K, 197.0K, 197.0K, 197.0K, 197.0K, 0.0m^3, 0.0m^3, 0.0m^3, 0.0m^3, 0.0J/m^2, 0.0m]
maxs = [320.0K, 320.0K, 320.0K, 80K, 320.0K, 320.0K, 320.0K, 320.0K, 1.0m^3, 1.0m^3, 1.0m^3, 1.0m^3, 3.0e7J/m^2, 0.1m]

# Load EVi and gbif counts
total_evi_counts = JLD.load("Total_evi_counts.jld", "total")
total_gbif_counts = JLD.load("Total_gbif_counts.jld", "total")

# Adjustment - remove NaNs and Infs
adjustment = total_gbif_counts ./ total_evi_counts
adjustment[isnan.(adjustment)] .= 1
adjustment[isinf.(adjustment)] .= 1

# Apply adjustment to data grouped by Species ID
phylo_traits_adj = @groupby gbif_fil :SppID {tmin = adjust(uconvert.(K, :tmin), adjustment[:, 1], mins[1], maxs[1]),tmax = adjust(uconvert.(K, :tmax), adjustment[:, 2], mins[2], maxs[2]), tmean = adjust(uconvert.(K, :tmean), adjustment[:, 3], mins[3], maxs[3]), trng = adjust(uconvert.(K, :trng), adjustment[:, 4], mins[4], maxs[4]), stl1 = adjust(uconvert.(K, :stl1mean), adjustment[:, 5], mins[5], maxs[5]), stl2 = adjust(uconvert.(K, :stl2mean), adjustment[:, 6], mins[6], maxs[6]), stl3 = adjust(uconvert.(K, :stl3mean), adjustment[:, 7], mins[7], maxs[7]), stl4 = adjust(uconvert.(K, :stl4mean), adjustment[:, 8], mins[8], maxs[8]), swvl1 = adjust(:swvl1mean, adjustment[:, 9], mins[9], maxs[9]), swvl2 =  adjust(:swvl2mean, adjustment[:, 10], mins[10], maxs[10]), swvl3 =  adjust(:swvl3mean, adjustment[:, 11], mins[11], maxs[11]), swvl4 =  adjust(:swvl4mean, adjustment[:, 12], mins[12], maxs[12]), ssr =  adjust(:ssrmean, adjustment[:, 13], mins[13], maxs[13]), tp =  adjust(:tpmean, adjustment[:, 14], mins[14], maxs[14])}

# Add tip names and save data
trait_ids = collect(JuliaDB.select(phylo_traits_adj, :SppID))
new_cross_species = [sppdict[i] for i in trait_ids]
phylo_traits_adj = pushcol(phylo_traits_adj, :tipNames, join.(split.(new_cross_species, " "), "_"))
JuliaDB.save(phylo_traits_adj, "Phylo_traits_adjust")

# Filter for common species
phylo_traits_filter = filter(p -> p.tipNames in join.(split.(top_common_names, " "), "_"), phylo_traits_adj)

# Convert to Dataframe
dat = DataFrame(collect(phylo_traits_filter))

# Fit lambda models and save
lambdas = fitLambdas(tree, dat, names(dat)[2:end-1])
JLD.save("Lambdas_EVI_adjust.jld", "lambdas", lambdas)


### LAMBDAS FOR CLIMATE-EFFORT ADJUSTED DATA ###

# Load all adjustment counts: EVI, GBIF, CERA
total_evi_counts = JLD.load("Total_evi_counts.jld", "total")
total_gbif_counts = JLD.load("Total_gbif_counts.jld", "total")
total_cera_counts = JLD.load("Total_cera_counts.jld", "total")

# Make adjustment for effort and divide by CERA counts
adjustment = total_gbif_counts ./ (total_evi_counts .* total_cera_counts)
adjustment[isnan.(adjustment)] .= 1
adjustment[isinf.(adjustment)] .= 1

# Apply adjustment to data grouped by Species ID
phylo_traits_adj2 = @groupby gbif_fil :SppID {tmin = adjust(uconvert.(K, :tmin), adjustment[:, 1], mins[1], maxs[1]),tmax = adjust(uconvert.(K, :tmax), adjustment[:, 2], mins[2], maxs[2]), tmean = adjust(uconvert.(K, :tmean), adjustment[:, 3], mins[3], maxs[3]), trng = adjust(uconvert.(K, :trng), adjustment[:, 4], mins[4], maxs[4]), stl1 = adjust(uconvert.(K, :stl1mean), adjustment[:, 5], mins[5], maxs[5]), stl2 = adjust(uconvert.(K, :stl2mean), adjustment[:, 6], mins[6], maxs[6]), stl3 = adjust(uconvert.(K, :stl3mean), adjustment[:, 7], mins[7], maxs[7]), stl4 = adjust(uconvert.(K, :stl4mean), adjustment[:, 8], mins[8], maxs[8]), swvl1 = adjust(:swvl1mean, adjustment[:, 9], mins[9], maxs[9]), swvl2 =  adjust(:swvl2mean, adjustment[:, 10], mins[10], maxs[10]), swvl3 =  adjust(:swvl3mean, adjustment[:, 11], mins[11], maxs[11]), swvl4 =  adjust(:swvl4mean, adjustment[:, 12], mins[12], maxs[12]), ssr =  adjust(:ssrmean, adjustment[:, 13], mins[13], maxs[13]), tp =  adjust(:tpmean, adjustment[:, 14], mins[14], maxs[14])}

# Add tip names and save data
phylo_traits_adj2 = pushcol(phylo_traits_adj2, :tipNames, join.(split.(new_cross_species, " "), "_"))
JuliaDB.save(phylo_traits_adj2, "Phylo_traits_adjust2")

# Filter for common species
phylo_traits_filter = filter(p -> p.tipNames in join.(split.(top_common_names, " "), "_"), phylo_traits_adj2)

# Convert to Dataframe
dat = DataFrame(collect(phylo_traits_filter))

# Fit lambda models and save
lambdas = fitLambdas(tree, dat, names(dat)[2:end-1])
JLD.save("Lambdas_temp_adjust.jld", "lambdas", lambdas)


### PLOT LAMBDA RESULTS AS HEATMAP ###
using JLD
using Plots
import Plots.px
pyplot()
lambdas_1 = JLD.load("data/Lambdas_common.jld", "lambdas")
lambdas_2 = JLD.load("data/Lambdas_EVI_adjust.jld", "lambdas")
lambdas_3 = JLD.load("data/Lambdas_temp_adjust.jld", "lambdas")
subset = [[1,3, 5]; collect(7:16)]
subset2 = [collect(1:3); collect(5:14)]
x = ["Raw", "Effort", "Effort + \n Climate"]
y = ["tmin", "tmax", "tmean", "stl1", "stl2", "stl3", "stl4", "swvl1", "swvl2", "swvl3", "swvl4", "ssr", "tp"]
lambdas = hcat(lambdas_1[subset], lambdas_2[subset2], lambdas_3[subset2])
heatmap(y, x, transpose(lambdas), seriescolor = :Blues, colorbar = :legend, legend = :top, size = (900, 200), guidefontsize = 12, tickfontsize = 12, xrotation = 90, clim = (0, 1), colorbar_title = "λ")
Plots.pdf("plots/Lambda_heatmap.pdf")


lambdas_1 = JLD.load("data/Lambdas_raw_continent.jld", "lambdas")
lambdas_2 = JLD.load("data/Lambdas_effort_continent.jld", "lambdas")
lambdas_3 = JLD.load("data/Lambdas_climate_continent.jld", "lambdas")
subset = [collect(1:3); collect(5:14)]
x = ["Raw", "Effort", "Effort + Climate"]
y = ["tmin", "tmax", "tmean", "stl1", "stl2", "stl3", "stl4", "swvl1", "swvl2", "swvl3", "swvl4", "ssr", "tp"]
lambdas = hcat(lambdas_1[subset], lambdas_2[subset], lambdas_3[subset])
heatmap(x, y, lambdas, seriescolor = :Blues, clim = (0, 1), guidefontsize =10, tickfontsize=10, size = (500, 500), top_margin = 20px, bottom_margin = 20px, colorbar_title = "λ")
Plots.pdf("plots/Lambda_continent_heatmap.pdf")

files = ["data/Corr_raw.jld", "data/Corr_adjust.jld", "data/Corr_adjust2.jld"]
corrs = map(files) do f
    JLD.load(f, "corr")
end
x = ["Raw", "Effort", "Effort \n + \n Climate"]
y = ["tmin", "tmax", "tmean", "stl1", "stl2", "stl3", "swvl1", "swvl2", "swvl3", "swvl4", "ssr", "tp"]
corrs = hcat(corrs ...)
h = heatmap(x, y, corrs, seriescolor = :RdBu, clim = (-1, 1), layout = (1,2), subplot =1, title = "Imputed data", size = (700, 500),bottom_margin=20px, left_margin=20px, right_margin=20px, top_margin=30px, tickfontsize = 12, colorbar = :none)

files = ["data/Corr_rand_raw.jld", "data/Corr_rand_adjust.jld", "data/Corr_rand_adjust2.jld"]
corrs = map(files) do f
    JLD.load(f, "corr")
end
x = ["Raw", "Effort", "Effort \n + \n Climate"]
y = ["tmin", "tmax", "tmean", "stl1", "stl2", "stl3", "swvl1", "swvl2", "swvl3", "swvl4", "ssr", "tp"]
corrs = hcat(corrs ...)
h = heatmap!(x, y, corrs, seriescolor = :RdBu, clim = (-1, 1), subplot =2, title = "Randomised \n imputed data", tickfontsize = 12)
Plots.pdf(h, "plots/Correlation_heatmap.pdf")


### PLOT EXAMPLE OF SPECIES ADJUSTMENT ###
using JuliaDB
using JLD
using Plots
using OnlineStats
using Unitful
using Unitful.DefaultSymbols
using MyUnitful
pyplot()

gbif = JuliaDB.load("GBIF_JOIN/CERA_JOIN_SIMPLE")
solanum_dulcamara = filter(g -> g.SppID == 336226, gbif)
tmin = collect(select(solanum_dulcamara, :tmin))
edges = range(197K, stop = 320K, length = 1000)
h = Hist(edges)
fit!(h, uconvert.(K, tmin))
JLD.save("solanum_dulcamara_counts.jld", "total", h.counts)

total_evi_counts = JLD.load("Total_evi_counts.jld", "total")
total_gbif_counts = JLD.load("Total_gbif_counts.jld", "total")
total_cera_counts = JLD.load("Total_cera_counts.jld", "total")
solanum_dulcamara = JLD.load("solanum_dulcamara_counts.jld", "total")

# Make adjustment for effort and divide by CERA counts
adjustment1 = total_gbif_counts ./ total_evi_counts
adjustment2 = total_gbif_counts ./ (total_evi_counts .* total_cera_counts)

edges = range(197K, stop = 320K, length = 1000)
tmin = ustrip.(collect(uconvert.(°C, edges)))
bar(tmin[418:815], solanum_dulcamara[418:815], layout = (3, 1), subplot = 1, size = (600, 800), grid = false, fillcolor=1, linecolor =:match, legend = false)
bar!(tmin[418:815], solanum_dulcamara[418:815] .* adjustment1[418:815, 1], layout = (3, 1), subplot = 2, size = (600, 800), grid = false, fillcolor = 2, linecolor =:match, legend = false)
bar!(tmin[418:815], solanum_dulcamara[418:815] .* adjustment2[418:815, 1], layout = (3, 1), subplot = 3, size = (600, 800), grid = false, fillcolor = 3, linecolor =:match, legend = false)



### LAMBDAS FOR PERCENTILES ###

mins = [197.0K, 0.0m]
maxs = [320.0K, 0.1m]

# Load EVi and gbif counts
total_evi_counts = JLD.load("Total_evi_counts.jld", "total")
total_gbif_counts = JLD.load("Total_gbif_counts.jld", "total")

# Adjustment - remove NaNs and Infs
adjustment = total_gbif_counts ./ total_evi_counts
adjustment[isnan.(adjustment)] .= 1
adjustment[isinf.(adjustment)] .= 1

# Apply adjustment to data grouped by Species ID
phylo_traits_tmin = @groupby gbif_fil :SppID {tmin10 = adjust_percentile(uconvert.(K, :tmin), adjustment[:, 1], mins[1], maxs[1], 10), tmin20 = adjust_percentile(uconvert.(K, :tmin), adjustment[:, 1], mins[1], maxs[1], 20), tmin30 = adjust_percentile(uconvert.(K, :tmin), adjustment[:, 1], mins[1], maxs[1], 30), tmin40 = adjust_percentile(uconvert.(K, :tmin), adjustment[:, 1], mins[1], maxs[1], 40), tmin50 = adjust_percentile(uconvert.(K, :tmin), adjustment[:, 1], mins[1], maxs[1], 50), tmin60 = adjust_percentile(uconvert.(K, :tmin), adjustment[:, 1], mins[1], maxs[1], 60), tmin70 = adjust_percentile(uconvert.(K, :tmin), adjustment[:, 1], mins[1], maxs[1], 70), tmin80 = adjust_percentile(uconvert.(K, :tmin), adjustment[:, 1], mins[1], maxs[1], 80), tmin90 = adjust_percentile(uconvert.(K, :tmin), adjustment[:, 1], mins[1], maxs[1], 90)}

phylo_traits_tmax = @groupby gbif_fil :SppID {tmax10 = adjust_percentile(uconvert.(K, :tmax), adjustment[:, 1], mins[1], maxs[1], 10), tmax20 = adjust_percentile(uconvert.(K, :tmax), adjustment[:, 1], mins[1], maxs[1], 20), tmax30 = adjust_percentile(uconvert.(K, :tmax), adjustment[:, 1], mins[1], maxs[1], 30), tmax40 = adjust_percentile(uconvert.(K, :tmax), adjustment[:, 1], mins[1], maxs[1], 40), tmax50 = adjust_percentile(uconvert.(K, :tmax), adjustment[:, 1], mins[1], maxs[1], 50), tmax60 = adjust_percentile(uconvert.(K, :tmax), adjustment[:, 1], mins[1], maxs[1], 60), tmax70 = adjust_percentile(uconvert.(K, :tmax), adjustment[:, 1], mins[1], maxs[1], 70), tmax80 = adjust_percentile(uconvert.(K, :tmax), adjustment[:, 1], mins[1], maxs[1], 80), tmax90 = adjust_percentile(uconvert.(K, :tmax), adjustment[:, 1], mins[1], maxs[1], 90)}

phylo_traits_tp = @groupby gbif_fil :SppID {tp10 = adjust_percentile(uconvert.(mm, :tpmean), adjustment[:, 14], mins[2], maxs[2], 10), tp20 = adjust_percentile(uconvert.(mm, :tpmean), adjustment[:, 14], mins[2], maxs[2], 20), tp30 = adjust_percentile(uconvert.(mm, :tpmean), adjustment[:, 14], mins[2], maxs[2], 30), tp40 = adjust_percentile(uconvert.(mm, :tpmean), adjustment[:, 14], mins[2], maxs[2], 40), tp50 = adjust_percentile(uconvert.(mm, :tpmean), adjustment[:, 14], mins[2], maxs[2], 50), tp60 = adjust_percentile(uconvert.(mm, :tpmean), adjustment[:, 14], mins[2], maxs[2], 60), tp70 = adjust_percentile(uconvert.(mm, :tpmean), adjustment[:, 14], mins[2], maxs[2], 70), tp80 = adjust_percentile(uconvert.(mm, :tpmean), adjustment[:, 14], mins[2], maxs[2], 80), tp90 = adjust_percentile(uconvert.(mm, :tpmean), adjustment[:, 14], mins[2], maxs[2], 90)

# Add tip names and save data
trait_ids = collect(JuliaDB.select(phylo_traits_tmin, :SppID))
new_cross_species = cross_species[indexin(trait_ids, cross_ids)]
phylo_traits_tmin = pushcol(phylo_traits_tmin, :tipNames, join.(split.(new_cross_species, " "), "_"))
phylo_traits_tmax = pushcol(phylo_traits_tmax, :tipNames, join.(split.(new_cross_species, " "), "_"))
phylo_traits_tp = pushcol(phylo_traits_tp, :tipNames, join.(split.(new_cross_species, " "), "_"))

# Filter for common species
phylo_traits_tmin = filter(p -> p.tipNames in join.(split.(top_common_names, " "), "_"), phylo_traits_tmin)
phylo_traits_tmax = filter(p -> p.tipNames in join.(split.(top_common_names, " "), "_"), phylo_traits_tmax)
phylo_traits_tp = filter(p -> p.tipNames in join.(split.(top_common_names, " "), "_"), phylo_traits_tp)

# Convert to Dataframe
dat1 = DataFrame(collect(phylo_traits_tmin))
dat2 = DataFrame(collect(phylo_traits_tmax))
dat3 = DataFrame(collect(phylo_traits_tp))

# Fit lambda models and save
lambdas = fitLambdas(tree, dat1, names(dat1)[2:end-1])
JLD.save("Lambdas_percentiles_tmin.jld", "lambdas", lambdas)

lambdas = fitLambdas(tree, dat2, names(dat2)[2:end-1])
JLD.save("Lambdas_percentiles_tmax.jld", "lambdas", lambdas)

lambdas = fitLambdas(tree, dat3, names(dat3)[2:end-1])
JLD.save("Lambdas_percentiles_tp.jld", "lambdas", lambdas)
