#### CODE FOR EVI-ADJUSTED BINS ####
using ClimatePref
using JuliaDB
using JuliaDBMeta
using JLD
using RCall
@everywhere using OnlineStats
@everywhere using Unitful
@everywhere using Unitful.DefaultSymbols

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
evi_tab = table(x, y, evi_mat, names = (:x, :y, :evi))

ref = create_reference(0.75)
refval = extractvalues(x .* °, y .* °, ref)
evi_tab = pushcol(evi_tab, :refval, refval)
save(evi_tab, "EVI")

# Min and max of bins
@everywhere mins = [197.0K, 197.0K, 197.0K, 0K, 197.0K, 197.0K, 197.0K, 197.0K, 0.0m^3, 0.0m^3, 0.0m^3, 0.0m^3, 0.0J/m^2, 0.0m]
@everywhere maxs = [320.0K, 320.0K, 320.0K, 80K, 320.0K, 320.0K, 320.0K, 320.0K, 1.0m^3, 1.0m^3, 1.0m^3, 1.0m^3, 3.0e7J/m^2, 0.1m]

# Function to bin data into equal sized bins of length 1,000
@everywhere function bin(x, min, max)
    edges = range(min, stop = max, length = 1000)
    h = Hist(edges)
    fit!(h, x)
    return h.counts
end

cera_simple = JuliaDB.load("CERA_simple")

# Bin each variable of cera_simple
cera_counts = @groupby cera_simple :refval {tmin = bin(uconvert.(K, :tmin), mins[1], maxs[1]),tmax = bin(uconvert.(K, :tmax), mins[2], maxs[2]), tmean = bin(uconvert.(K, :tmean), mins[3], maxs[3]), trng = bin(uconvert.(K, :trng), mins[4], maxs[4]), stl1 = bin(uconvert.(K, :stl1mean), mins[5], maxs[5]), stl2 = bin(uconvert.(K, :stl2mean), mins[6], maxs[6]), stl3 = bin(uconvert.(K, :stl3mean), mins[7], maxs[7]), stl4 = bin(uconvert.(K, :stl4mean), mins[8], maxs[8]), swvl1 = bin(:swvl1mean, mins[9], maxs[9]), swvl2 =  bin(:swvl2mean, mins[10], maxs[10]), swvl3 =  bin(:swvl3mean, mins[11], maxs[11]), swvl4 =  bin(:swvl4mean, mins[12], maxs[12]), ssr =  bin(:ssrmean, mins[13], maxs[13]), tp =  bin(:tpmean, mins[14], maxs[14])}
JuliaDB.save(cera_counts, "CERA_counts")

cera_counts = JuliaDB.load("CERA_counts")

# Load EVI data and group into means by location
evi_tab = JuliaDB.load("EVI")
evi_refs = @groupby evi_tab :refval {evi_means = mean(:evi)}

# Filter EVI data to those locations that are in CERA
cr = collect(select(cera_counts, :refval))
evi_ref = filter(e -> e.refval in cr, evi_refs)
evi_means = select(evi_ref, :evi_means)

# Add EVI data to CERA and multiply counts
cera_counts = pushcol(cera_counts, :evi, evi_means)
cera_counts = @transform cera_counts {tmin = :tmin .* :evi, tmax = :tmax .* :evi, tmean = :tmean .* :evi, trng = :trng .* :evi, stl1 = :stl1 .* :evi, stl2 = :stl2 .* :evi, stl3 = :stl3 .* :evi, stl4 = :stl4 .* :evi, swvl1 = :swvl1 .* :evi, swvl2 =  :swvl2 .* :evi, swvl3 =  :swvl3 .* :evi, swvl4 =  :swvl4 .* :evi, ssr =  :ssr .* :evi, tp =  :tp .* :evi}

# Function to sum up all counts across the different rows of the dataframe
function count_convert(tab::JuliaDB.DIndexedTable, var::Symbol)
    counts = collect(JuliaDB.select(tab, var))
    mat = Array{Union{Float64, Missing}}(undef, length(tab), length(counts[1]))
    for i in 1:length(counts)
        mat[i, :] .= counts[i]
    end
    return mapslices(x -> sum(skipmissing(x)), mat, dims = 1)[1, :]
end
total_evi_counts = [count_convert(cera_counts, i) for i in [:tmin, :tmax, :tmean, :trng, :stl1, :stl2, :stl3, :stl4, :swvl1, :swvl2, :swvl3, :swvl4, :ssr, :tp]]
total_evi_counts = hcat(total_evi_counts...)

# Save EVI-adjusted counts
JLD.save("Total_evi_counts.jld", "total", total_evi_counts)


## Join cera and evi data with continents and filter missing values
cera_simple = JuliaDB.load("CERA_simple")
continents = JuliaDB.load("Continents")
evi_tab = JuliaDB.load("GBIF_JOIN/EVI")
evi_cont = join(evi_tab, continents, how = :inner, lkey = :refval, rkey = :refval)
evi_cont = filter(e -> !ismissing(e.continent), evi_cont)
continents = distribute(continents, 1)
cera_cont = join(cera_simple, continents, how = :inner, lkey = :refval, rkey = :refval)
cera_cont = filter(e -> !ismissing(e.continent), cera_cont)

function binContinents(cera::JuliaDB.DIndexedTable, evi::JuliaDB.IndexedTable, continent::Int64, mins, maxs)
    cera = filter(c -> c.continent == continent, cera)
    cera_counts = @groupby cera :refval {tmin = bin(uconvert.(K, :tmin), mins[1], maxs[1]),tmax = bin(uconvert.(K, :tmax), mins[2], maxs[2]), tmean = bin(uconvert.(K, :tmean), mins[3], maxs[3]), trng = bin(uconvert.(K, :trng), mins[4], maxs[4]), stl1 = bin(uconvert.(K, :stl1mean), mins[5], maxs[5]), stl2 = bin(uconvert.(K, :stl2mean), mins[6], maxs[6]), stl3 = bin(uconvert.(K, :stl3mean), mins[7], maxs[7]), stl4 = bin(uconvert.(K, :stl4mean), mins[8], maxs[8]), swvl1 = bin(:swvl1mean, mins[9], maxs[9]), swvl2 =  bin(:swvl2mean, mins[10], maxs[10]), swvl3 =  bin(:swvl3mean, mins[11], maxs[11]), swvl4 =  bin(:swvl4mean, mins[12], maxs[12]), ssr =  bin(:ssrmean, mins[13], maxs[13]), tp =  bin(:tpmean, mins[14], maxs[14])}

    evi_tab = filter(e -> e.continent == continent, evi)
    evi_refs = @groupby evi_tab :refval {evi_means = mean(:evi)}

    cr = collect(select(cera_counts, :refval))
    evi_ref = filter(e -> e.refval in cr, evi_refs)
    evi_means = select(evi_ref, :evi_means)

    cera_counts = pushcol(cera_counts, :evi, evi_means)
    cera_counts = @transform cera_counts {tmin = :tmin .* :evi, tmax = :tmax .* :evi, tmean = :tmean .* :evi, trng = :trng .* :evi, stl1 = :stl1 .* :evi, stl2 = :stl2 .* :evi, stl3 = :stl3 .* :evi, stl4 = :stl4 .* :evi, swvl1 = :swvl1 .* :evi, swvl2 =  :swvl2 .* :evi, swvl3 =  :swvl3 .* :evi, swvl4 =  :swvl4 .* :evi, ssr =  :ssr .* :evi, tp =  :tp .* :evi}

    total_evi_counts = [count_convert(cera_counts, i) for i in [:tmin, :tmax, :tmean, :trng, :stl1, :stl2, :stl3, :stl4, :swvl1, :swvl2, :swvl3, :swvl4, :ssr, :tp]]
    return total_evi_counts = hcat(total_evi_counts...)
end

function mapContinents!(cera::JuliaDB.DIndexedTable, evi::JuliaDB.IndexedTable, continents::Vector{Int64}, mins, maxs)
    map(continents) do cont
        total_evi_counts = binContinents(cera, evi, cont, mins, maxs)
        JLD.save("Total_evi_counts_$cont.jld", "total", total_evi_counts)
        print(cont)
    end
end
mapContinents!(cera_cont, evi_cont, collect(1:6))
