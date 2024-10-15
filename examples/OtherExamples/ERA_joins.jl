# SPDX-License-Identifier: BSD-2-Clause

using JuliaDB
using ClimatePref
using Unitful
using ClimatePref.Units

# Load GBIF and add in unique ID column
gbif = load("GBIF_TPL")
gbif = pushcol(gbif, :UID, 1:length(gbif))

# Create date column
mth = collect(select(gbif, :month))
yr = collect(select(gbif, :year))
dt = (mth .* month) .+ (yr .* year)
dt = uconvert.(year, dt)
gbif = pushcol(gbif, :date, dt)

# Extract coordinates and compare to reference grid
lat = collect(select(gbif, :decimallatitude))
lon = collect(select(gbif, :decimallongitude))
ref = create_reference(0.75)
refval = extractvalues(lon .* °, lat .* °, ref)
gbif = pushcol(gbif, :refval, refval)
small_gbif = select(gbif, (:UID, :SppID, :date, :refval))
save(small_gbif, "Small_GBIF")

# Import era data and refine
era = load("ECMWF/era_int_all")
era = @transform era {newssr = :ssr / m^4}
era = popcol(era, :ssr)
era = renamecol(era, :newssr, :ssr)
yr = select(era, :year) .* year
mth = select(era, :month) .* month
function era_extend(era::IndexedTable)
    for i in 0:11
        dt = yr .+ mth
        dt .+= (i * month)
        dt = uconvert.(year, dt)
        era = pushcol(era, Symbol("date$i"), dt)
    end
    return era
end
era_date = era_extend(era)
era_date = renamecol(era_date, :refid => :refval)
save(era_date, "ECMWF/era_int")

# Join ERA with gbif (once for every twelve months)
gbif = load("GBIF/Small_GBIF")
era = load("ECMWF/era_int")
era = distribute(era, 1)

function era_gbif_join(era::JuliaDB.DIndexedTable, gbif::JuliaDB.DIndexedTable,
                       i::Int64)
    date_var = Symbol("date$i")
    bif = renamecol(gbif, colnames(gbif)[3] => date_var)
    gbif_era = join(bif, era, how = :inner, lkey = (:refval, date_var),
                    rkey = (:refval, date_var),
                    rselect = (:year, :refval, :x, :y, :month, :ssr, :stl1,
                               :t2m, :tp, date_var))
    gbif_era = renamecol(gbif_era, date_var, :date)
    gbif_era = reindex(gbif_era, (:SppID, :refval, :date))
    return gbif_era
end

function era_gbif_join(era::JuliaDB.DIndexedTable, gbif::JuliaDB.DIndexedTable,
                       filename::String)
    for i in 0:11
        if isfile(filename)
            gbif_era = load(filename)
            new_join = era_gbif_join(era, gbif, i)
            gbif_era = merge(gbif_era, new_join)
            save(gbif_era, filename)
        else
            gbif_era = era_gbif_join(era, gbif, i)
            save(gbif_era, filename)
        end
        print(i, "\n")
    end
end

function era_gbif_join(era::JuliaDB.DIndexedTable, gbif::JuliaDB.DIndexedTable,
                       filename::String)
    for i in 0:11
        gbif_era = era_gbif_join(era, gbif, i)
        save(gbif_era, filename * "_$i")
        print(i, "\n")
    end
end
era_gbif_join(era, gbif, "ECMWF/GBIF_ERA")

function era_load_merge(filename::String)
    for i in 0:11
        if isfile(filename)
            gbif_era = load(filename)
            new_merge = load(filename * "_$i")
            gbif_era = merge(gbif_era, new_join)
            save(gbif_era, filename)
        else
            gbif_era = load(filename * "_$i")
            save(gbif_era, filename)
        end
        print(i, "\n")
    end
end
era_load_merge("ECMWF/GBIF_ERA")

era = load("ECMWF/era_int_all")

# New idea
@everywhere function vec2array(v)
    r = length(v)
    c = length(v[1])
    a = Array{Int64}(r, c)
    for i in 1:r, j in 1:c
        a[i, j] = v[i][j]
    end
    return a
end
variables = [:stl1, :t2m, :tp, :ssr]

splits = Dict(zip(variables,
                  [collect(-60:1:50), collect(-60:1:50), collect(-80:1:40),
                      collect(0:0.1:4), collect(0:0.5:30)]))
for j in 2:7
    @everywhere func(y) = fit(Histogram, ustrip.(y),
                              splits[worldclim[j]]).weights
    bins = @sync @parallel (merge) for i in eachindex(genera)
        genus = JuliaDB.load(string("Genera/Worldclim/", genera[i]))
        spp = unique(select(genus, :species))
        spp = spp[spp .!= ""]
        if length(spp) == 0
            bin = AxisArray(Array{Int64, 2}(0, 30),
                            Axis{:Species}(spp),
                            Axis{:Bins}(splits[worldclim[j]][2:end]))
        else
            res = groupby(:func => func, genus, :species, select = worldclim[j])
            res = filter(p -> p.species != "", res)
            bin = AxisArray(vec2array(select(res, :func)),
                            Axis{:Species}(collect(select(res, :species))),
                            Axis{:Bins}(splits[worldclim[j]][2:end]))
        end
        bin
    end
    JLD.save(string("Binned_data_", worldclim[j], ".jld"),
             string("Binned_data_", worldclim[j]), bins)
end
