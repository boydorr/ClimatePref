# SPDX-License-Identifier: BSD-2-Clause

# 16. Perform phylogenetic trait analysis on Worldclim data

using Unitful
using Unitful.DefaultSymbols
using ClimatePref
using ClimatePref.Units
using JuliaDB

prec = load("Worldclim/prec")
srad = load("Worldclim/srad")
tavg = load("Worldclim/tavg")
tmax = load("Worldclim/tmax")
tmin = load("Worldclim/tmin")

function addtabs(wc)
    vars = [srad, tavg, tmax, tmin]
    varnames = [:srad, :tavg, :tmax, :tmin]
    for i in eachindex(varnames)
        wc = pushcol(wc, varnames[i], select(vars[i], :val))
    end
    return wc
end
wc = addtabs(renamecol(prec, :val, :prec))
wc = reindex(wc, :refval)
wc = distribute(wc, 1)

gbif = JuliaDB.load("GBIF_TPL")
ref = create_reference(1 / 12)
x = collect(JuliaDB.select(gbif, :decimallatitude))
y = collect(JuliaDB.select(gbif, :decimallongitude))
refval = extractvalues(y .* °, x .* °, ref)
gbif = pushcol(gbif, :refval, refval)
gbif = reindex(gbif, :refval)

gbif_join = join(gbif, wc, how = :inner, lkey = :refval, rkey = :refval,
                 rselect = (:prec, :srad, :tavg, :tmax, :tmin, :month),
                 lselect = (:SppID, :refval, :date))
save(gbif_join, "GBIF_WC")

gbif_join = filter(g -> !isnan(g.prec), gbif_join)

using JLD
using JuliaDB
using JuliaDBMeta
@everywhere using OnlineStats

# Variable names and corresponding min to max values
vars = [:prec, :srad, :tavg, :tmax, :tmin]
mins = [0.0mm, 0.0kJ / (u"d" * m^2), 200K, 200K, 200K]
maxs = [2600.0mm, 50000.0kJ / (u"d" * m^2), 320K, 320K, 320K]

# Bin gbif data into 1,000 bins per variable
total_gbif_counts = [bin(gbif_join, vars[i], mins[i], maxs[i])
                     for i in eachindex(vars)]
total_gbif_counts = hcat(total_gbif_counts...)
JLD.save("Total_gbif_counts_wc.jld", "total", total_gbif_counts)

wc = filter(w -> !isnan(w.prec), wc)
@everywhere function bin(tab::JuliaDB.DIndexedTable, var::Symbol, min, max)
    x = collect(select(tab, var))
    if unit(x[1]) == °C
        x = uconvert.(K, x)
    end
    edges = range(min, stop = max, length = 1000)
    h = Hist(edges)
    fit!(h, x)
    return h.counts
end

# Variable names and corresponding min to max values
vars = [:prec, :srad, :tavg, :tmax, :tmin]
mins = [0.0mm, 0.0kJ / (u"d" * m^2), 200K, 200K, 200K]
maxs = [2600.0mm, 50000.0kJ / (u"d" * m^2), 320K, 320K, 320K]

total_wc_counts = [bin(wc, vars[i], mins[i], maxs[i]) for i in eachindex(vars)]
total_wc_counts = hcat(total_wc_counts...)
JLD.save("Total_wc_counts.jld", "total", total_wc_counts)

using RCall
R"library(raster)
evi = raster('Evi.tif')
print(evi)
res = res(evi)
evi_mat = evi[]
x = coordinates(evi)[, 1]
y = coordinates(evi)[, 2]
"
@rget evi_mat
@rget x
@rget y
#lon = -180:res[1]:180
#lat = -90:res[2]:90
#x = repeat(lon, length(lat))
#y = vcat(map(x-> fill(x, length(lon)), lat)...)
evi_tab = table(x, y, evi_mat, names = (:x, :y, :evi))

ref = create_reference(1 / 12)
refval = extractvalues(x .* °, y .* °, ref)
evi_tab = pushcol(evi_tab, :refval, refval)
save(evi_tab, "EVI_WC")

@everywhere function bin(x, min, max)
    edges = range(min, stop = max, length = 1000)
    h = Hist(edges)
    fit!(h, x)
    return h.counts
end
wc_counts = @groupby wc :refval {prec = bin(:prec, mins[1], maxs[1]),
                                 srad = bin(:srad, mins[2], maxs[2]),
                                 tavg = bin(:tavg, mins[3], maxs[3]),
                                 tmax = bin(:tmax, mins[4], maxs[4]),
                                 tmin = bin(:tmin, mins[5], maxs[5])}
save(wc_counts, "WC_counts")

wc_counts = JuliaDB.load("WC_counts")
evi_tab = JuliaDB.load("EVI_WC")
evi_refs = @groupby evi_tab :refval {evi_means = mean(:evi)}

evi_refs = distribute(evi_refs, 1)

wc_counts = join(wc_counts, evi_refs, how = :inner, lkey = :refval,
                 rkey = :refval)
wc_counts = @transform wc_counts {prec = :prec .* :evi_means,
                                  srad = :srad .* :evi_means,
                                  tavg = :tavg .* :evi_means,
                                  tmax = :tmax .* :evi_means,
                                  tmin = :tmin .* :evi_means}

function count_convert(tab::JuliaDB.DIndexedTable, var::Symbol)
    counts = collect(select(tab, var))
    mat = Array{Union{Float64, Missing}}(undef, length(tab), length(counts[1]))
    for i in 1:length(counts)
        mat[i, :] .= counts[i]
    end
    return mapslices(x -> sum(skipmissing(x)), mat, dims = 1)[1, :]
end
total_evi_counts = [count_convert(wc_counts, i)
                    for i in [:prec, :srad, :tavg, :tmax, :tmin]]
total_evi_counts = hcat(total_evi_counts...)

# Save EVI-adjusted counts
JLD.save("Total_evi_counts_wc.jld", "total", total_evi_counts)

using Unitful
using Unitful.DefaultSymbols
using ClimatePref
using ClimatePref.Units
using JuliaDB
using PhyloNetworks
using GLM
using JLD
using Statistics
using JuliaDBMeta
using StatsBase
using DataFrames

tree = readTopology("Qian2016.tree")
tipnames = tipLabels(tree)
tip_names = join.(split.(tipnames, "_"), " ")
function fitLambdas(tree::HybridNetwork, dat::DataFrame, vars = Vector{Symbol})
    lambdas = map(vars) do v
        dat[v] = ustrip.(dat[v])
        f = @eval @formula($v~1)
        fitPagel = phyloNetworklm(f, dat, tree, model = "lambda")
        return lambda_estim(fitPagel)
    end
    return lambdas
end

function filterGBIF(file::String)
    gbif = JuliaDB.load(file)
    tree = readTopology("Qian2016.tree")
    tipnames = tipLabels(tree)
    tip_names = join.(split.(tipnames, "_"), " ")
    spp_names = JLD.load("Species_names.jld", "spp_names")
    spp_ids = JLD.load("Species_names.jld", "spp_ids")
    sppdict = Dict(zip(spp_names, spp_ids))
    cross_species = spp_names ∩ tip_names
    cross_ids = [sppdict[x] for x in cross_species]
    gbif_fil = filter(g -> g.SppID in cross_ids, gbif)
    gbif_fil = filter(g -> !isnan(g.prec), gbif_fil)
    return gbif_fil
end

#spp = collect(JuliaDB.select(gbif, :SppID))
#numspp = unique(spp)

spp_names = JLD.load("Species_names.jld", "spp_names")
spp_ids = JLD.load("Species_names.jld", "spp_ids")
sppdict = Dict(zip(spp_names, spp_ids))

# Get top 5000 most common species
cross_species = spp_names ∩ tip_names
cross_ids = [sppdict[x] for x in cross_species]
#sorted_counts = countmap(spp)
#sorted_counts = filter((k, v) -> k in cross_ids, sorted_counts)
#sorted_counts = sort(sorted_counts, byvalue=true, rev=true)
#top_common_ids = collect(keys(sorted_counts))[1:5000]
#top_common_names = spp_names[indexin(top_common_ids, numspp)]

top_common_names = JLD.load("Common_species_names.jld", "spp_names")
missing_species = setdiff(tip_names, top_common_names)
for i in eachindex(missing_species)
    deleteleaf!(tree, join(split(missing_species[i], " "), "_"))
end

gbif_fil = filterGBIF("GBIF_WC")
JuliaDB.save(gbif_fil, "GBIF_FIL_WC")

### LAMBDAS FOR RAW DATA ###

# Group gbif data by Species ID and get mean and percentiles of data
@everywhere using Statistics
@everywhere using Unitful

function meanGBIF!(gbif_fil, cross_species, cross_ids)
    gbif_fil = @groupby gbif_fil :SppID {prec = mean(skipmissing(ustrip(:prec))),
                                         srad = mean(skipmissing(ustrip(:srad))),
                                         tavg = mean(skipmissing(ustrip(:tavg))),
                                         tmax = mean(skipmissing(ustrip(:tmax))),
                                         tmin = mean(skipmissing(ustrip(:tmin)))}

    # Add in tip names to data and save
    trait_ids = collect(JuliaDB.select(gbif_fil, :SppID))
    new_cross_species = cross_species[indexin(trait_ids, cross_ids)]
    gbif_fil = pushcol(gbif_fil, :tipNames,
                       join.(split.(new_cross_species, " "), "_"))
    return gbif_fil
end

gbif_fil = JuliaDB.load("GBIF_FIL_WC")
gbif_fil = meanGBIF!(gbif_fil, cross_species, cross_ids)

JuliaDB.save(gbif_fil, "Phylo_traits_WC")

# Filter for common species
top_common_names = JLD.load("Common_species_names.jld", "spp_names")
phylo_traits = JuliaDB.load("Phylo_traits_WC")
phylo_traits_filter = filter(p -> p.tipNames in join.(split.(top_common_names,
                                                             " "), "_"),
                             phylo_traits)

# Convert to dataframe
dat = DataFrame(collect(phylo_traits_filter))

# Fit lambda models and save
lambdas = fitLambdas(tree, dat, names(dat)[2:(end - 1)])
JLD.save("Lambdas_common_wc.jld", "lambdas", lambdas)

### LAMBDAS FOR EFFORT-ADJUSTED DATA ###

mins = [0.0mm, 0.0kJ / (u"d" * m^2), 200K, 200K, 200K]
maxs = [2600.0mm, 50000.0kJ / (u"d" * m^2), 320K, 320K, 320K]

# Load EVi and gbif counts
total_evi_counts = JLD.load("Total_evi_counts_wc.jld", "total")
total_gbif_counts = JLD.load("Total_gbif_counts_wc.jld", "total")

# Adjustment - remove NaNs and Infs
adjustment = total_gbif_counts ./ total_evi_counts
adjustment[isnan.(adjustment)] .= 1
adjustment[isinf.(adjustment)] .= 1
@everywhere using OnlineStats
@everywhere using Statistics
@everywhere using StatsBase
@everywhere using Unitful
@everywhere using Unitful.DefaultSymbols
@everywhere function adjust(x, adj, min, max)
    edges = range(min, stop = max, length = 1000)
    h = Hist(edges)
    fit!(h, x)
    counts = h.counts .* adj
    step = edges[2] - edges[1]
    edges = collect(edges) .+ step / 2
    return mean(edges[1:(end - 1)], weights(counts))
end

# Apply adjustment to data grouped by Species ID

function meanGBIF(gbif_fil, cross_species, cross_ids, mins, maxs)
    gbif_fil = @groupby gbif_fil :SppID {prec = adjust(:prec, adjustment[:, 1],
                                                       mins[1], maxs[1]),
                                         srad = adjust(:srad, adjustment[:, 2],
                                                       mins[2], maxs[2]),
                                         tavg = adjust(uconvert.(K, :tavg),
                                                       adjustment[:, 3],
                                                       mins[3], maxs[3]),
                                         tmax = adjust(uconvert.(K, :tmax),
                                                       adjustment[:, 4],
                                                       mins[4], maxs[4]),
                                         tmin = adjust(uconvert.(K, :tmin),
                                                       adjustment[:, 5],
                                                       mins[5], maxs[5])}

    # Add in tip names to data and save
    trait_ids = collect(JuliaDB.select(gbif_fil, :SppID))
    new_cross_species = cross_species[indexin(trait_ids, cross_ids)]
    return gbif_fil = pushcol(gbif_fil, :tipNames,
                              join.(split.(new_cross_species, " "), "_"))
end
gbif_fil = JuliaDB.load("GBIF_FIL_WC")
phylo_traits_adj = meanGBIF(gbif_fil, cross_species, cross_ids, mins, maxs)

JuliaDB.save(phylo_traits_adj, "Phylo_traits_adjust_WC")

#phylo_traits_adj = JuliaDB.load("Phylo_traits_adjust_WC")
# Filter for common species
phylo_traits_filter = filter(p -> p.tipNames in join.(split.(top_common_names,
                                                             " "), "_"),
                             phylo_traits_adj)

# Convert to Dataframe
dat = DataFrame(collect(phylo_traits_filter))

# Fit lambda models and save
lambdas = fitLambdas(tree, dat, names(dat)[2:(end - 1)])
JLD.save("Lambdas_EVI_adjust_WC.jld", "lambdas", lambdas)

### LAMBDAS FOR CLIMATE-ADJUSTED DATA ###

mins = [0.0mm, 0.0kJ / (u"d" * m^2), 200K, 200K, 200K]
maxs = [2600.0mm, 50000.0kJ / (u"d" * m^2), 320K, 320K, 320K]

# Load EVi and gbif counts
total_evi_counts = JLD.load("../../sdc/Total_evi_counts_wc.jld", "total")
total_gbif_counts = JLD.load("../../sdc/Total_gbif_counts_wc.jld", "total")
total_wc_counts = JLD.load("../../sdc/Total_wc_counts.jld", "total")

# Make adjustment for effort and divide by wc counts
adjustment = total_gbif_counts ./ (total_evi_counts .* total_wc_counts)
adjustment[isnan.(adjustment)] .= 1
adjustment[isinf.(adjustment)] .= 1

gbif_fil = JuliaDB.load("GBIF_FIL_WC")
phylo_traits_adj2 = meanGBIF(gbif_fil, cross_species, cross_ids, mins, maxs)

JuliaDB.save(phylo_traits_adj2, "Phylo_traits_adjust2_WC")

phylo_traits_adj2 = JuliaDB.load("Phylo_traits_adjust2_WC")
# Filter for common species
phylo_traits_filter = filter(p -> p.tipNames in join.(split.(top_common_names,
                                                             " "), "_"),
                             phylo_traits_adj2)

# Convert to Dataframe
dat = DataFrame(collect(phylo_traits_filter))

# Fit lambda models and save
lambdas = fitLambdas(tree, dat, names(dat)[2:(end - 1)])
JLD.save("Lambdas_temp_adjust_WC.jld", "lambdas", lambdas)

using JLD
using Plots
import Plots.px
pyplot()
lambdas_1 = JLD.load("Lambdas_common_wc.jld", "lambdas")
lambdas_2 = JLD.load("Lambdas_EVI_adjust_WC.jld", "lambdas")
lambdas_3 = JLD.load("Lambdas_temp_adjust_WC.jld", "lambdas")
x = ["Raw", "Effort", "Effort + \n Climate"]
y = ["tp", "ssr", "tmean", "tmax", "tmin"]
lambdas = hcat(lambdas_1, lambdas_2, lambdas_3)
heatmap(y[[5, 3, 4, 2, 1]], x, transpose(lambdas[[5, 3, 4, 2, 1], :]),
        seriescolor = :Blues, colorbar = :legend, legend = :top,
        size = (900, 200), guidefontsize = 12, tickfontsize = 12,
        xrotation = 90, clim = (0, 1), colorbar_title = "λ")
Plots.pdf("plots/Lambda_heatmap_wc.pdf")
