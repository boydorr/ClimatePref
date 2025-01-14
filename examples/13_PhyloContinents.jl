# SPDX-License-Identifier: BSD-2-Clause

# 13. Phylogenetic analysis adjusted on a continent level
using RCall
using JuliaDB
using JuliaDBMeta
using Unitful
using Unitful.DefaultSymbols
using PhyloNetworks
using ClimatePref
using ClimatePref.Unitful
using StatsBase
using JLD
using JLD2
using DataFrames

# Load tree
tree = readTopology("Qian2016.tree")
tipnames = tipLabels(tree)
tip_names = join.(split.(tipnames, "_"), " ")
# Load GBIF data
gbif = JuliaDB.load("CERA_JOIN_SIMPLE")
spp = collect(JuliaDB.select(gbif, :SppID))
numspp = unique(spp)

# Extracting continent values
R"library(raster); library(rgdal)
shapefolder <- 'Continent'; shapefile <- 'Continent_crop'
continent <- readOGR(shapefolder, shapefile)
ext <- extent(-180, 180, -90, 90)
r <- raster(ext, resolution= c(360/481, 180/241))
cont_raster <- rasterize(continent, r)
cont_mat = getValues(cont_raster)
x = coordinates(cont_raster)[, 1]
y = coordinates(cont_raster)[, 2]
"
@rget cont_mat
@rget x
@rget y
cont_tab = table(x, y, cont_mat, names = (:x, :y, :continent))

ref = create_reference(0.75)
refval = extractvalues(x .* °, y .* °, ref)
cont_tab = pushcol(cont_tab, :refval, refval)
JuliaDB.save(cont_tab, "Continents")

# Join GBIF with Continents

continents = JuliaDB.load("Continents")
continents = distribute(continents, 1)
gbif_cont = join(gbif, continents, how = :inner, lkey = :refval, rkey = :refval)
gbif_cont = filter(c -> !ismissing(c.continent), gbif_cont)

# Load Species names
spp_names = JLD.load("Species_names.jld", "spp_names")
spp_ids = JLD.load("Species_names.jld", "spp_ids")
sppdict = Dict(zip(spp_ids, spp_names))
iddict = Dict(zip(spp_names, spp_ids))
spp_names = [sppdict[i] for i in numspp]

# Get top 5000 most common species
cross_species = spp_names ∩ tip_names
cross_ids = [iddict[i] for i in cross_species]
gbif_fil = filter(g -> g.SppID in cross_ids, gbif_cont)

top_common_names = JLD.load("Common_species_names.jld", "spp_names")

missing_species = setdiff(tip_names, top_common_names)
for i in eachindex(missing_species)
    deleteleaf!(tree, join(split(missing_species[i], " "), "_"))
end

spp_dict = Dict(zip(cross_ids, cross_species))

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

@everywhere function adjust_cont(x, adj, min, max)
    edges = range(min, stop = max, length = 1000)
    h = Hist(edges)
    fit!(h, x)
    counts = h.counts .* adj
    return counts
end

function adjustData(gbif::JuliaDB.DIndexedTable, spp_dict::Dict,
                    continent::Int64, filter_names::Vector{String},
                    adjustment::Array{Float64, 2}, mins, maxs)
    # Filter for continent
    gbif1 = filter(g -> g.continent == continent, gbif)

    # Apply adjustment to data grouped by Species ID
    phylo_traits = @groupby gbif1 :SppID {tmin = adjust_cont(uconvert.(K,
                                                                       :tmin),
                                                             adjustment[:, 1],
                                                             mins[1], maxs[1]),
                                          tmax = adjust_cont(uconvert.(K,
                                                                       :tmax),
                                                             adjustment[:, 2],
                                                             mins[2], maxs[2]),
                                          tmean = adjust_cont(uconvert.(K,
                                                                        :tmean),
                                                              adjustment[:, 3],
                                                              mins[3], maxs[3]),
                                          trng = adjust_cont(uconvert.(K,
                                                                       :trng),
                                                             adjustment[:, 4],
                                                             mins[4], maxs[4]),
                                          stl1 = adjust_cont(uconvert.(K,
                                                                       :stl1mean),
                                                             adjustment[:, 5],
                                                             mins[5], maxs[5]),
                                          stl2 = adjust_cont(uconvert.(K,
                                                                       :stl2mean),
                                                             adjustment[:, 6],
                                                             mins[6], maxs[6]),
                                          stl3 = adjust_cont(uconvert.(K,
                                                                       :stl3mean),
                                                             adjustment[:, 7],
                                                             mins[7], maxs[7]),
                                          stl4 = adjust_cont(uconvert.(K,
                                                                       :stl4mean),
                                                             adjustment[:, 8],
                                                             mins[8], maxs[8]),
                                          swvl1 = adjust_cont(:swvl1mean,
                                                              adjustment[:, 9],
                                                              mins[9], maxs[9]),
                                          swvl2 = adjust_cont(:swvl2mean,
                                                              adjustment[:, 10],
                                                              mins[10],
                                                              maxs[10]),
                                          swvl3 = adjust_cont(:swvl3mean,
                                                              adjustment[:, 11],
                                                              mins[11],
                                                              maxs[11]),
                                          swvl4 = adjust_cont(:swvl4mean,
                                                              adjustment[:, 12],
                                                              mins[12],
                                                              maxs[12]),
                                          ssr = adjust_cont(:ssrmean,
                                                            adjustment[:, 13],
                                                            mins[13], maxs[13]),
                                          tp = adjust_cont(:tpmean,
                                                           adjustment[:, 14],
                                                           mins[14], maxs[14]),
                                          samp_size = length(:refval)}

    # Add in tip names to data and save
    phylo_traits1 = @transform phylo_traits {tipNames = join.(split.(spp_dict[:SppID],
                                                                     " "), "_")}

    # Filter for common species
    phylo_traits_filter = filter(p -> p.tipNames in join.(split.(filter_names,
                                                                 " "), "_"),
                                 phylo_traits1)

    # Convert to dataframe
    return phylo_traits_filter
end

@everywhere using StatsBase
@everywhere using Unitful
@everywhere using Unitful.DefaultSymbols
@everywhere using OnlineStats

@everywhere function countmean(x, min, max)
    edges = range(min, stop = max, length = 1000)
    step = edges[2] - edges[1]
    edges = collect(edges) .+ step / 2
    return mean(edges[1:(end - 1)], weights(x))
end

@everywhere function weightmean(x, y)
    return mapslices(x -> mean(x, weights(y)), hcat(x...), dims = 2)
end
function fitAdjust(gbif::JuliaDB.DIndexedTable, continents::Vector{Int64},
                   spp_dict::Dict, top_common_names::Vector{String}, mins, maxs;
                   climate::Bool = false, raw::Bool = false)
    total_dat = []
    for c in continents[1:end]
        evi_counts = JLD.load("Total_evi_counts_$c.jld", "total")
        gbif_counts = JLD.load("Total_gbif_counts_$c.jld", "total")
        cera_counts = JLD.load("Total_cera_counts_$c.jld", "total")
        adjustment = ifelse(climate, gbif_counts ./ evi_counts,
                            gbif_counts ./ (evi_counts .* cera_counts))
        adjustment[isnan.(adjustment)] .= 1
        adjustment[isinf.(adjustment)] .= 1
        if raw
            fill!(adjustment, 1)
        end
        dat = adjustData(gbif, spp_dict, c, top_common_names, adjustment, mins,
                         maxs)
        if length(total_dat) > 0
            total_dat = merge(total_dat, dat)
        else
            total_dat = dat
        end
    end
    total_dat1 = @groupby total_dat :SppID {tmin = countmean(weightmean(:tmin,
                                                                        :samp_size),
                                                             mins[1], maxs[1]),
                                            tmax = countmean(weightmean(:tmax,
                                                                        :samp_size),
                                                             mins[2], maxs[2]),
                                            tmean = countmean(weightmean(:tmean,
                                                                         :samp_size),
                                                              mins[3], maxs[3]),
                                            trng = countmean(weightmean(:trng,
                                                                        :samp_size),
                                                             mins[4], maxs[4]),
                                            stl1 = countmean(weightmean(:stl1,
                                                                        :samp_size),
                                                             mins[5], maxs[5]),
                                            stl2 = countmean(weightmean(:stl2,
                                                                        :samp_size),
                                                             mins[6], maxs[6]),
                                            stl3 = countmean(weightmean(:stl3,
                                                                        :samp_size),
                                                             mins[7], maxs[7]),
                                            stl4 = countmean(weightmean(:stl4,
                                                                        :samp_size),
                                                             mins[8], maxs[8]),
                                            swvl1 = countmean(weightmean(:swvl1,
                                                                         :samp_size),
                                                              mins[9], maxs[9]),
                                            swvl2 = countmean(weightmean(:swvl2,
                                                                         :samp_size),
                                                              mins[10],
                                                              maxs[10]),
                                            swvl3 = countmean(weightmean(:swvl3,
                                                                         :samp_size),
                                                              mins[11],
                                                              maxs[11]),
                                            swvl4 = countmean(weightmean(:swvl4,
                                                                         :samp_size),
                                                              mins[12],
                                                              maxs[12]),
                                            tp = countmean(weightmean(:tp,
                                                                      :samp_size),
                                                           mins[13], maxs[13]),
                                            ssr = countmean(weightmean(:ssr,
                                                                       :samp_size),
                                                            mins[14], maxs[14]),
                                            tipNames = first(:tipNames)}
    total_dat1 = DataFrame(collect(total_dat1))
    total_dat1[:tipNames] = collect(total_dat1[:tipNames])
    return total_dat1
end

dat1 = fitAdjust(gbif_fil, collect(1:6), spp_dict, top_common_names, mins, maxs,
                 raw = true)
JLD2.@save "Phylo_dat_continent.jld" dat1

dat2 = fitAdjust(gbif_fil, collect(1:6), spp_dict, top_common_names, mins, maxs)
JLD2.@save "Phylo_dat_adjust_continent.jld2" dat2

dat3 = fitAdjust(gbif_fil, collect(1:6), spp_dict, top_common_names, mins, maxs,
                 climate = true, raw = false)
JLD2.@save "Phylo_dat_adjust2_continent.jld2" dat3

lambdas1 = fitLambdas(tree, dat1, names(dat1)[2:(end - 1)])
JLD.save("Lambdas_raw_continent.jld", "lambdas", lambdas1)

lambdas2 = fitLambdas(tree, dat2, names(dat2)[2:(end - 1)])
JLD.save("Lambdas_effort_continent.jld", "lambdas", lambdas2)

lambdas3 = fitLambdas(tree, dat3, names(dat3)[2:(end - 1)])
JLD.save("Lambdas_climate_continent.jld", "lambdas", lambdas3)

using JLD
using Plots
pyplot()
lambdas_1 = JLD.load("data/Lambdas_raw_continent.jld", "lambdas")
lambdas_2 = JLD.load("data/Lambdas_effort_continent.jld", "lambdas")
x = ["Raw", "Effort"]
y = [
    "tmin",
    "tmax",
    "tmean",
    "stl1",
    "stl2",
    "stl3",
    "stl4",
    "swvl1",
    "swvl2",
    "swvl3",
    "swvl4",
    "ssr",
    "tp"
]
subset = [collect(1:3); collect(5:12); 14; 13]
lambdas = hcat(lambdas_1, lambdas_2)
heatmap(y, x, transpose(lambdas[subset, :]),
        seriescolor = cgrad(:bluesreds_r, [0, mean(lambdas), 1]),
        colorbar = :legend, legend = :top, size = (900, 200),
        guidefontsize = 12, tickfontsize = 12, xrotation = 90, clim = (0, 1),
        colorbar_title = "λ")
Plots.pdf("Lambda_continent_heatmap.pdf")
