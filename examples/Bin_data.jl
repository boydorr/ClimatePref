addprocs(8)
@everywhere using StatsBase
@everywhere using Unitful
@everywhere using JuliaDB
@everywhere using AxisArrays
using ClimatePref
using JLD
@everywhere import Unitful.ustrip
@everywhere ustrip(x::DataValues.DataValue) = ustrip(x.value)
@everywhere function vec2array(v)
    r = length(v)
    c = length(v[1])
    a = Array{Int64}(r,c)
    for i in 1:r, j in 1:c
        a[i,j] = v[i][j]
    end
    a
end
dir = "../Worldclim/Monthly"
worldclim = map(searchdir(dir, "")) do str
    Symbol(str)
end
dir = "Genera/Worldclim/Monthly"
genera = searchdir(dir, "")


splits = Dict(zip(worldclim, [collect(0:50:3000), collect(0:500:70000),
    collect(-70:1:50), collect(-70:1:50), collect(-80:1:40),
    collect(0:0.1:4), collect(0:0.5:30)]))
for j in 2:7
@everywhere func(y) = fit(Histogram, ustrip.(y), splits[worldclim[j]]).weights
    bins = @sync @parallel (merge) for i in eachindex(genera)
        genus = JuliaDB.load(string("Genera/Worldclim/", genera[i]))
        spp = unique(select(genus, :species))
        spp = spp[spp .!= ""]
        if length(spp) == 0
            bin = AxisArray(Array{Int64, 2}(0, 30),
            Axis{:Species}(spp), Axis{:Bins}(splits[worldclim[j]][2:end]))
        else
        res = groupby(:func => func, genus, :species, select=worldclim[j])
        res = filter(p->p.species != "", res)
        bin = AxisArray(vec2array(select(res, :func)),
        Axis{:Species}(collect(select(res, :species))), Axis{:Bins}(splits[worldclim[j]][2:end]))
    end
    bin
    end
    JLD.save(string("Binned_data_",worldclim[j],".jld"),
     string("Binned_data_", worldclim[j]), bins)
end
j=1
for j in 2:7
@everywhere func(y) = fit(Histogram, ustrip.(y), splits[worldclim[j]], closed=:right).weights
    bins = @sync @parallel (merge) for i in eachindex(genera)
        genus = JuliaDB.load(string("Genera/Worldclim/Monthly/", genera[i]))
        spp = unique(select(genus, :species))
        spp = spp[spp .!= ""]
        if length(spp) == 0
            bin = AxisArray(Array{Int64, 2}(0, length(splits[worldclim[j]][2:end])),
            Axis{:Species}(spp), Axis{:Bins}(splits[worldclim[j]][2:end]))
        else
        res = groupby(:func => func, genus, :species, select=worldclim[j])
        res = filter(p->p.species != "", res)
        bin = AxisArray(vec2array(select(res, :func)),
        Axis{:Species}(collect(select(res, :species))), Axis{:Bins}(splits[worldclim[j]][2:end]))
    end
    bin
    end
    JLD.save(string("Binned_data_monthly_",worldclim[j],".jld"),
     string("Binned_data_monthly_", worldclim[j]), bins)
end

function detect_anomaly(bins::AxisArray{Int64, 2}, gap::Int64)
    anomalies = Array{Int64, 1}(0)
    for i in 1:size(bins, 1)
        row = bins[i, :]
        if sum(row) > 0
        start = minimum(find(row .> 0))[1]
        ending = length(row) - minimum(find(reverse(row) .> 0))[1] + 1
        row = row[start:ending]
        count = 0
        for j in eachindex(row)
            if count >= gap
                push!(anomalies, i)
                break
            end
            if row[j] == 0
                count = count + 1
            else
                count = 0
            end
        end
    end
end
    return bins.axes[1].val[anomalies]
end

end
JLD.save("Binned_data.jld", "Binned_data", bins)
