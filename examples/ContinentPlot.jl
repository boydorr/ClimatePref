using RCall
using JuliaDB
using JuliaDBMeta
using Unitful
using Unitful.DefaultSymbols
using PhyloNetworks
using MyUnitful
using ClimatePref
using Simulation
using StatsBase
using JLD
using JLD2
using DataFrames

tree = readTopology("Qian2016.tree")
tipnames = tipLabels(tree)
tip_names = join.(split.(tipnames, "_"), " ")
# Load GBIF data
gbif =  JuliaDB.load("CERA_JOIN_SIMPLE")
spp = collect(JuliaDB.select(gbif, :SppID))
numspp = unique(spp)

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
gbif_fil = filter(g->g.SppID in cross_ids, gbif_cont)

top_common_names = JLD.load("Common_species_names.jld", "spp_names")

missing_species = setdiff(tip_names, top_common_names)
for i in eachindex(missing_species)
    deleteleaf!(tree, join(split(missing_species[i], " "), "_"))
end

spp_dict = Dict(zip(cross_ids, cross_species))

mins = [197.0K, 197.0K, 197.0K, 0K, 197.0K, 197.0K, 197.0K, 197.0K, 0.0m^3, 0.0m^3, 0.0m^3, 0.0m^3, 0.0J/m^2, 0.0m]
maxs = [320.0K, 320.0K, 320.0K, 80K, 320.0K, 320.0K, 320.0K, 320.0K, 1.0m^3, 1.0m^3, 1.0m^3, 1.0m^3, 3.0e7J/m^2, 0.1m]

@everywhere function adjust_cont(x, adj, min, max)
    edges = range(min, stop = max, length = 1000)
    h = Hist(edges)
    fit!(h, x)
    counts = h.counts .* adj
    step = edges[2]-edges[1]
    edges = collect(edges) .+ step/2
    return mean(edges[1:(end-1)], weights(counts))
end

function adjustData(gbif::JuliaDB.DIndexedTable,  spp_dict::Dict, continent::Int64, filter_names::Vector{String}, adjustment::Array{Float64, 2}, mins, maxs)
    # Filter for continent
    gbif1 = filter(g->g.continent == continent, gbif_fil)

    # Apply adjustment to data grouped by Species ID
    phylo_traits = @groupby gbif1 :SppID {tmin = adjust_cont(uconvert.(K, :tmin), adjustment[:, 1], mins[1], maxs[1]),tmax = adjust_cont(uconvert.(K, :tmax), adjustment[:, 2], mins[2], maxs[2]), tmean = adjust_cont(uconvert.(K, :tmean), adjustment[:, 3], mins[3], maxs[3]), trng = adjust_cont(uconvert.(K, :trng), adjustment[:, 4], mins[4], maxs[4]), stl1 = adjust_cont(uconvert.(K, :stl1mean), adjustment[:, 5], mins[5], maxs[5]), stl2 = adjust_cont(uconvert.(K, :stl2mean), adjustment[:, 6], mins[6], maxs[6]), stl3 = adjust_cont(uconvert.(K, :stl3mean), adjustment[:, 7], mins[7], maxs[7]), stl4 = adjust_cont(uconvert.(K, :stl4mean), adjustment[:, 8], mins[8], maxs[8]), swvl1 = adjust_cont(:swvl1mean, adjustment[:, 9], mins[9], maxs[9]), swvl2 =  adjust_cont(:swvl2mean, adjustment[:, 10], mins[10], maxs[10]), swvl3 =  adjust_cont(:swvl3mean, adjustment[:, 11], mins[11], maxs[11]), swvl4 =  adjust_cont(:swvl4mean, adjustment[:, 12], mins[12], maxs[12]), ssr =  adjust_cont(:ssrmean, adjustment[:, 13], mins[13], maxs[13]), tp =  adjust_cont(:tpmean, adjustment[:, 14], mins[14], maxs[14]), samp_size = length(:refval)}

    # Add in tip names to data and save
    phylo_traits1 = @transform phylo_traits {tipNames = join.(split.(spp_dict[:SppID], " "), "_")}

    # Filter for common species
    phylo_traits_filter = filter(p -> p.tipNames in join.(split.(filter_names, " "), "_"), phylo_traits1)

    # Convert to dataframe
    dat = DataFrame(collect(phylo_traits_filter))
    return dat
end
@everywhere using StatsBase
@everywhere using Unitful
@everywhere using Unitful.DefaultSymbols
@everywhere using OnlineStats

function fitAdjust!(gbif::JuliaDB.DIndexedTable, tree::HybridNetwork, continents::Vector{Int64}, spp_dict::Dict, top_common_names::Vector{String}, mins, maxs, filename::String; climate::Bool = false, raw::Bool = false)
    for c in continents[1:end]
        evi_counts = JLD.load("Total_evi_counts_$c.jld", "total")
        gbif_counts = JLD.load("Total_gbif_counts_$c.jld", "total")
        cera_counts = JLD.load("Total_cera_counts_$c.jld", "total")
        adjustment = ifelse(climate, gbif_counts ./ evi_counts, gbif_counts ./ (evi_counts .* cera_counts))
        adjustment[isnan.(adjustment)] .= 1
        adjustment[isinf.(adjustment)] .= 1
        if raw
            fill!(adjustment, 1)
        end
        dat = adjustData(gbif_fil, spp_dict, c, top_common_names, adjustment, mins, maxs)
        dat[:tipNames] = collect(dat[:tipNames])
        JLD2.@save filename * "$c.jld" dat
    end
end

fitAdjust!(gbif_fil, tree, collect(1:6), spp_dict, top_common_names, mins, maxs, "Continent_dat_raw", raw = true)
fitAdjust!(gbif_fil, tree, collect(1:6), spp_dict, top_common_names, mins, maxs, "Continent_dat_adjust")
fitAdjust!(gbif_fil, tree, collect(1:6), spp_dict, top_common_names,  mins, maxs, "Continent_dat_adjust2", climate = true, raw = false)

for i in 1:6
    dat1 = JLD.load("Continent_dat_raw$i.jld", "dat")
    lambdas1 = fitLambdas(tree, dat1, names(dat1)[2:end-1])
    JLD.save("Lambdas_raw_continent$i.jld", "lambdas", lambdas1)

    dat2 = JLD.load("Continent_dat_adjust$i.jld", "dat")
    lambdas2 = fitLambdas(tree, dat2, names(dat2)[2:end-1])
    JLD.save("Lambdas_effort_continent$i.jld", "lambdas", lambdas2)

    dat3 = JLD.load("Continent_dat_adjust2$i.jld", "dat")
    lambdas3 = fitLambdas(tree, dat3, names(dat3)[2:end-1])
    JLD.save("Lambdas_climate_continent$i.jld", "lambdas", lambdas3)
end

using JLD
using Plots
pyplot()

heatmap(seriescolor = :Blues, size = (1000, 500), guidefontsize = 12, tickfontsize = 12, xrotation = 90, clim = (0, 1), colorbar_title = "λ", layout = (@layout [a b c; d e f]), link = :both, margin = 5*Plots.mm, titlefontsize = 14)
continents = ["Africa", "South America", "North America", "Europe", "Asia", "Australasia"]
for i in 1:6
    lambdas_1 = JLD.load("data/Lambdas_raw_continent$i.jld", "lambdas")
    lambdas_2 = JLD.load("data/Lambdas_effort_continent$i.jld", "lambdas")
    lambdas_3 = JLD.load("data/Lambdas_climate_continent$i.jld", "lambdas")

    x = ["Raw", "Effort", "Effort + \n Climate"]
    y = ["tmin", "tmax", "tmean", "stl1", "stl2", "stl3", "stl4", "swvl1", "swvl2", "swvl3", "swvl4", "ssr", "tp"]
    legend = ifelse((i == 3) | (i == 6), true, false)
    yax = ifelse((i == 1) | (i == 4), true, false)
    xax = ifelse(i > 3, true, false)
    subset = [collect(1:3); collect(5:14)]
    lambdas = hcat(lambdas_1, lambdas_2, lambdas_3)
    display(heatmap!(y, x, transpose(lambdas[subset, :]), seriescolor = :Blues, colorbar = :legend, guidefontsize = 12, tickfontsize = 12, xrotation = 90, clim = (0, 1), colorbar_title = "λ", subplot = i, title = continents[i], legend = legend, xaxis = xax, yaxis = yax))
end
Plots.pdf("plots/Lambda_continent_heatmap.pdf")

using JuliaDB
using JuliaDBMeta
using Unitful
using DataFrames
using CSV
dat = load("Phylo_traits_adjust")
dat = collect(dat)
dat = @transform dat {tmin = ustrip.(uconvert.(°C, :tmin)), tmax = ustrip.(uconvert.(°C, :tmax)), tmean = ustrip.(uconvert.(°C, :tmean)), trng = ustrip.(uconvert.(°C, :trng)), stl1 = ustrip.(uconvert.(°C, :stl1)), stl2 = ustrip.(uconvert.(°C, :stl2)), stl3 = ustrip.(uconvert.(°C, :stl3)), stl4 = ustrip.(uconvert.(°C, :stl4)), swvl1 = ustrip.(:swvl1), swvl2 = ustrip.(:swvl2), swvl3 = ustrip.(:swvl3), swvl4 = ustrip.(:swvl4), ssr = ustrip.(:ssr), tp = ustrip.(uconvert.(mm, :tp))}
dat2 = DataFrame(dat)
CSV.write("Phylo_traits.csv", dat2)

continent = JLD.load("Phylo_dat_continent.jld", "dat")
names = JLD.load("Common_species_names.jld", "spp_names")
names = join.(split.(names, " "), "_")

names = DataFrame(name = names[indexin(continent[:tipNames], names)], cont = continent[:c])
CSV.write("Common_species_names.csv", names)
