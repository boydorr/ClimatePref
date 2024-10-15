# SPDX-License-Identifier: BSD-2-Clause

# 10. Bin GBIF data
using ClimatePref
using JuliaDB
using JuliaDBMeta
using JLD
@everywhere using OnlineStats

gbif = JuliaDB.load("CERA_JOIN_SIMPLE")

# Function to bin GBIF data into even sized bins for variable and min/max values
function bin(tab::JuliaDB.DIndexedTable, var::Symbol, min, max)
    x = collect(JuliaDB.select(tab, var))
    if unit(x[1]) == °C
        x = uconvert.(K, x)
    end
    edges = range(min, stop = max, length = 1000)
    h = Hist(edges)
    fit!(h, x)
    return h.counts
end

# Variable names and corresponding min to max values
vars = [
    :tmin,
    :tmax,
    :tmean,
    :trng,
    :stl1mean,
    :stl2mean,
    :stl3mean,
    :stl4mean,
    :swvl1mean,
    :swvl2mean,
    :swvl3mean,
    :swvl4mean,
    :ssrmean,
    :tpmean
]
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

# Bin gbif data into 1,000 bins per variable
total_gbif_counts = [bin(gbif, vars[i], mins[i], maxs[i])
                     for i in eachindex(vars)]
total_gbif_counts = hcat(total_gbif_counts...)
JLD.save("Total_gbif_counts.jld", "total", total_gbif_counts)

# Bin data per continent
gbif = JuliaDB.load("CERA_JOIN_SIMPLE")
continents = JuliaDB.load("Continents")
continents = distribute(continents, 1)
gbif_cont = join(gbif, continents, how = :inner, lkey = :refval, rkey = :refval)
gbif_cont = filter(c -> !ismissing(c.continent), gbif_cont)

function binContinents(gbif::JuliaDB.DIndexedTable, continent::Int64, mins,
                       maxs, vars)
    gbif = filter(c -> c.continent == continent, gbif)
    total_gbif_counts = [bin(gbif, vars[i], mins[i], maxs[i])
                         for i in eachindex(vars)]
    return total_gbif_counts = hcat(total_gbif_counts...)
end

function mapContinents!(gbif::JuliaDB.DIndexedTable, continents::Vector{Int64},
                        mins, maxs, vars)
    map(continents) do cont
        total_gbif_counts = binContinents(gbif, cont, mins, maxs, vars)
        JLD.save("Total_gbif_counts_$cont.jld", "total", total_gbif_counts)
        return print(cont)
    end
end
mapContinents!(gbif_cont, collect(1:6), mins, maxs, vars)

# Extract species names from TPL and save
tpl = ReadTPL("TPL")
spp = collect(select(gbif, :SppID))
numspp = unique(spp)

tab = filter(t -> t.SppID in numspp, tpl)
tpl_species = collect(select(tab, :species))
tpl_id = collect(select(tab, :SppID))
spp_names = tpl_species[indexin(numspp, tpl_id)]
JLD.save("GBIF_JOIN/Species_names.jld", "spp_names", spp_names)
