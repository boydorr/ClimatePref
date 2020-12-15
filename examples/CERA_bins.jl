#### CODE FOR CERA BINS ####
include("GIT/Chapter4/src/Chapter4.jl")
using .Chapter4
using JuliaDB
using JuliaDBMeta
@everywhere using OnlineStats

cera_simple = JuliaDB.load("CERA_simple")

# Function to bin CERA data into even sized bins for variable and min/max values
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
vars = [:tmin, :tmax, :tmean, :trng, :stl1mean, :stl2mean, :stl3mean, :stl4mean, :swvl1mean, :swvl2mean, :swvl3mean, :swvl4mean, :ssrmean, :tpmean]
mins = [197.0K, 197.0K, 197.0K, 0K, 197.0K, 197.0K, 197.0K, 197.0K, 0.0m^3, 0.0m^3, 0.0m^3, 0.0m^3, 0.0J/m^2, 0.0m]
maxs = [320.0K, 320.0K, 320.0K, 80K, 320.0K, 320.0K, 320.0K, 320.0K, 1.0m^3, 1.0m^3, 1.0m^3, 1.0m^3, 3.0e7J/m^2, 0.1m]

total_cera_counts = [bin(cera_simple, vars[i], mins[i], maxs[i]) for i in eachindex(vars)]
total_cera_counts = hcat(total_cera_counts...)
JLD.save("Total_cera_counts.jld", "total", total_cera_counts)


cera_simple = JuliaDB.load("CERA_simple")
continents = JuliaDB.load("Continents")
continents = distribute(continents, 1)
cera_cont = join(cera_simple, continents, how = :inner, lkey = :refval, rkey = :refval)
cera_cont = filter(c -> !ismissing(c.continent), cera_cont)

function binContinents(cera::JuliaDB.DIndexedTable, continent::Int64, mins, maxs, vars)
    cera = filter(c -> c.continent == continent, cera)
    total_cera_counts = [bin(cera, vars[i], mins[i], maxs[i]) for i in eachindex(vars)]
    return total_cera_counts = hcat(total_cera_counts...)
end

function mapContinents!(cera::JuliaDB.DIndexedTable, continents::Vector{Int64}, mins, maxs, vars)
    map(continents) do cont
        total_cera_counts = binContinents(cera, cont, mins, maxs, vars)
        JLD.save("Total_cera_counts_$cont.jld", "total", total_cera_counts)
        print(cont)
    end
end
mapContinents!(cera_cont, collect(1:6), mins, maxs, vars)
