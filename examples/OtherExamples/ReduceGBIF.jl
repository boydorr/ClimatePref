# SPDX-License-Identifier: BSD-2-Clause

addprocs(20)
@everywhere using JuliaDB
@everywhere using ClimatePref
@everywhere using Unitful
@everywhere using Unitful.DefaultSymbols
@everywhere function na_mean(x::Float32, y::Float32)
    if typeof(x) <: DataValues.DataValue
        x = x.value
    end
    if typeof(y) <: DataValues.DataValue
        y = y.value
    end
    if isnan(x) & isnan(y)
        return NaN32 * °C
    else
        vec = [x, y]
        return mean(vec[.!isnan.(vec)])
    end
end
@everywhere function na_mean(x::Float64, y::Float64)
    if typeof(x) <: DataValues.DataValue
        x = x.value
    end
    if typeof(y) <: DataValues.DataValue
        y = y.value
    end
    if isnan(x) & isnan(y)
        return NaN * °C
    else
        vec = [x, y]
        return mean(vec[.!isnan.(vec)])
    end
end
dir = "Documents/gbif/Genera/Worldclim/"
for i in searchdir(dir, "")
    genus = load(string(dir, i))
    if i == "Aa"
        tab = groupreduce(na_mean, genus, :species, select = :tavg)
    elseif length(genus) > 0
        red_genus = groupreduce(na_mean, genus, :species, select = :tavg)
        if typeof(select(red_genus, 2)) <: Array{NamedTuples.NamedTuple, 1}
            red_genus = table(select(red_genus, 1),
                              map(x -> collect(select(red_genus, 2)[x])[1],
                                  1:length(red_genus)),
                              names = [:species, :na_mean])
        end
        tab = merge(red_genus, tab)
    end
end
save(tab, "Documents/gbif/Genera/Grouped_genera")

@everywhere bioclim_names = ["annual_tavg", "mean_diurnal_range",
    "isothermality",
    "temp_seasonality", "tmax_warmestm", "tmin_coldestm", "annual_trange",
    "tavg_wettestq", "tavg_driestq", "tavg_warmestq", "tavg_coldestq",
    "annual_prec",
    "prec_wettestm", "prec_driestm", "prec_seasonality", "prec_wettestq",
    "prec_driestq",
    "prec_warmestq", "prec_coldestq"]
dir = "Genera/Bioclim/"
tab = @sync @parallel (merge) for i in searchdir(dir, "")
    genus = load(string(dir, i))
    if length(genus) > 0
        red_genus = groupby(:annual_tavg => na_mean, genus, :species,
                            select = :annual_tavg)
        for j in 2:19
            tab = groupby(Symbol.(bioclim_names)[j] => na_mean, genus,
                          :species, select = Symbol.(bioclim_names)[j])
            red_genus = pushcol(red_genus, Symbol.(bioclim_names)[j],
                                select(tab, Symbol.(bioclim_names)[j]))
        end
    end
    red_genus
end
save(tab, "Documents/gbif/Genera/Grouped_genera")
@everywhere function na_mean(vec)
    if typeof(vec) <: DataValues.DataValueArray
        vec = vec.values
    elseif typeof(vec) <: DataValues.DataValue
        vec = vec.value
    end
    return mean(vec[.!isnan.(vec)])
end
