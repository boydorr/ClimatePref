using PhyloNetworks
using GLM
using DataFrames
using StatsBase
using Statistics
using Unitful
using Random
using JuliaDB
using JuliaDBMeta
using RCall
using OnlineStats
using Unitful.DefaultSymbols

"""
    adjust(x, adj, min, max)

    Function to bin data into even sized bins and make adjustment, returning the mean value of count bins.

"""
function adjust(x, adj, min, max)
    edges = range(min, stop = max, length = 1000)
    h = Hist(edges)
    fit!(h, x)
    counts = h.counts .* adj
    step = edges[2]-edges[1]
    edges = collect(edges) .+ step/2
    return mean(edges[1:(end-1)], weights(counts))
end

"""
    adjust_percentile(x, adj, min, max, perc)

    Function to bin data into even sized bins and make adjustment, returning the percentile value of count bins.

"""
function adjust_percentile(x, adj, min, max, perc)
    edges = range(min, stop = max, length = 1000)
    h = Hist(edges)
    fit!(h, x)
    counts = h.counts .* adj
    step = edges[2]-edges[1]
    edges = collect(edges) .+ step/2
    ind = findlast(cumsum(counts) .<= (perc*sum(counts)/100))
    if isnothing(ind)
        ind = 1
    end
    return edges[1:(end-1)][ind]
end


"""
    fitLambdas(tree::HybridNetwork, dat::DataFrame, vars = Vector{Symbol})

    Function to fit Pagel's lambda models to a tree of type `HybridNetwork` for a range of variables, `vars`, in a dataframe of traits, `dat`. Returns a vector of lambda values, one per trait.
"""
function fitLambdas(tree::HybridNetwork, dat::DataFrame, vars = Vector{Symbol})
    lambdas = map(vars) do v
        dat[v] = ustrip.(dat[v])
        f = @eval @formula($v ~ 1)
        fitPagel = phyloNetworklm(f, dat, tree, model="lambda")
        lambda_estim(fitPagel)
    end
    return lambdas
end

"""
    fitMissings(tree::HybridNetwork, dat::DataFrame, climate_var::Symbol)

    Function to section the data into tenths and reconstruct missing values to create a fully reconstructed dataframe.
"""
function fitMissings(tree::HybridNetwork, dat::DataFrame, climate_var::Symbol, shuffle::Bool = false)
    if shuffle
        @warn "Tips have been shuffled randomly!"
        shuffle!(dat[climate_var])
    end
    dat[climate_var] = ustrip.(dat[climate_var])
    allowmissing!(dat)
    vals = collect(1:nrow(dat))
    traits = fill(0.0, length(vals))
    numSamples = round(Int64, nrow(dat) * 0.1)
    for i in 1:10
        smp = sort(sample(vals, numSamples, replace = false))
        vals = setdiff(vals, smp)
        dat2 = dat[[climate_var, :tipNames]]
        tipName = [n.name for n in tree.leaf]
        tipNumber = [n.number for n in tree.leaf]
        tipOrder = indexin(dat2[:tipNames], tipName)
        dat2[smp, climate_var] = missing
        anc = ancestralStateReconstruction(dat2, tree)
        dat2[:tipNums] = tipNumber[tipOrder]
        traits[smp] = anc.traits_nodes[indexin(dat2[:tipNums][smp], anc.NodeNumbers)]
    end
    return traits
end

"""
    rawData(gbif::JuliaDB.DIndexedTable, spp_dict::Dict, continent::Int64, filter_names::Vector{String})

    Ready raw data (no adjustments) for Lambda models.
"""
function rawData(gbif::JuliaDB.DIndexedTable, spp_dict::Dict, continent::Int64, filter_names::Vector{String})
    # Filter for continent
    gbif = filter(g->g.continent == continent, gbif)

    # Group by species
    phylo_traits = @groupby africa :SppID {tmin = mean(ustrip(:tmin)), tmin10 = percentile(ustrip(:tmin), 10), tmax = mean(ustrip(:tmax)), tmax90 = percentile(ustrip(:tmax), 90), tmean = mean(ustrip(:tmean)), trng = mean(ustrip(:trng)), stl1 = mean(ustrip(:stl1mean)), stl2 = mean(ustrip(:stl2mean)), stl3 = mean(ustrip(:stl3mean)), stl4 = mean(ustrip(:stl4mean)), swvl1 = mean(ustrip(:swvl1mean)), swvl2 = mean(ustrip(:swvl2mean)), swvl3 = mean(ustrip(:swvl3mean)), swvl4 = mean(ustrip(:swvl4mean)), ssr = mean(ustrip(:ssrmean)), tp = mean(ustrip(:tpmean)), samp_size = length(:refval)}

    # Add in tip names to data and save
    trait_ids = collect(JuliaDB.select(phylo_traits, :SppID))
    new_cross_species = [spp_dict[i] for i in trait_ids]
    phylo_traits = pushcol(phylo_traits, :tipNames, join.(split.(new_cross_species, " "), "_"))

    # Filter for common species
    phylo_traits_filter = filter(p -> p.tipNames in join.(split.(filter_names, " "), "_"), phylo_traits)

    # Convert to dataframe
    #dat = DataFrame(collect(phylo_traits_filter))
    return phylo_traits_filter
end

"""
    adjustData(gbif::JuliaDB.DIndexedTable,  spp_dict::Dict, continent::Int64, filter_names::Vector{String}, adjustment::Array{Float64, 2}, mins, maxs)

    Ready adjusted data (EVI or EVI/climate) for Lambda models.
"""
function adjustData(gbif::JuliaDB.DIndexedTable,  spp_dict::Dict, continent::Int64, filter_names::Vector{String}, adjustment::Array{Float64, 2}, mins, maxs)
    # Filter for continent
    gbif = filter(g->g.continent == continent, gbif)

    # Apply adjustment to data grouped by Species ID
    phylo_traits = @groupby gbif :SppID {tmin = adjust(uconvert.(K, :tmin), adjustment[:, 1], mins[1], maxs[1]),tmax = adjust(uconvert.(K, :tmax), adjustment[:, 2], mins[2], maxs[2]), tmean = adjust(uconvert.(K, :tmean), adjustment[:, 3], mins[3], maxs[3]), trng = adjust(uconvert.(K, :trng), adjustment[:, 4], mins[4], maxs[4]), stl1 = adjust(uconvert.(K, :stl1mean), adjustment[:, 5], mins[5], maxs[5]), stl2 = adjust(uconvert.(K, :stl2mean), adjustment[:, 6], mins[6], maxs[6]), stl3 = adjust(uconvert.(K, :stl3mean), adjustment[:, 7], mins[7], maxs[7]), stl4 = adjust(uconvert.(K, :stl4mean), adjustment[:, 8], mins[8], maxs[8]), swvl1 = adjust(:swvl1mean, adjustment[:, 9], mins[9], maxs[9]), swvl2 =  adjust(:swvl2mean, adjustment[:, 10], mins[10], maxs[10]), swvl3 =  adjust(:swvl3mean, adjustment[:, 11], mins[11], maxs[11]), swvl4 =  adjust(:swvl4mean, adjustment[:, 12], mins[12], maxs[12]), ssr =  adjust(:ssrmean, adjustment[:, 13], mins[13], maxs[13]), tp =  adjust(:tpmean, adjustment[:, 14], mins[14], maxs[14]), samp_size = length(:refval)}

    # Add in tip names to data and save
    trait_ids = collect(JuliaDB.select(phylo_traits, :SppID))
    new_cross_species = [spp_dict[i] for i in trait_ids]
    phylo_traits = pushcol(phylo_traits, :tipNames, join.(split.(new_cross_species, " "), "_"))

    # Filter for common species
    phylo_traits_filter = filter(p -> p.tipNames in join.(split.(filter_names, " "), "_"), phylo_traits)

    # Convert to dataframe
    #dat = DataFrame(collect(phylo_traits_filter))
    return phylo_traits_filter
end

"""
    extractContinents(filename::String)

    Function to extract continent level data from JuliaDB files (note that the files must have a reference value column). A shapefile must be provided of the continent values.
"""
function extractContinents(filename::String, shapefolder::String, shapefile::String)
    # Load data
    data = JuliaDB.load(filename)
    # Extract reference values
    refvals = collect(JuliaDB.select(data, :refval))
    x, y = convert_coords(refvals, 481)
    # Put into R and extract details from continent raster
    @rput shapefile; @rput shapefolder
    R"library(raster); library(rgdal)
    continent <- readOGR(shapefolder, shapefile)
    ext <- extent(-180, 180, -90, 90)
    r <- raster(ext, resolution= c(360/481, 180/241))
    cont_raster <- rasterize(continent, r)
    vals = getValues(cont_raster, format = 'matrix')
    "
    @rget vals
    [vals[x[i], y[i]] for i in length(x)]
    # Add column of continent
    data_cont = pushcol(data, :continent, vals)
    # Filter out missing values (sea, etc.)
    data_filter = JuliaDB.filter(c-> !ismissing(c.continent), data_cont)
    return data_filter
end
