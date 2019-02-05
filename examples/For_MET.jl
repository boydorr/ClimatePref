## TEST DATA FROM ERA INTERIM - TEMPERATURE AT 2M ##

using NetCDF
using Unitful
using Unitful.DefaultSymbols
using AxisArrays
using MyUnitful
using ClimatePref
using JuliaDB
using Plots
using JLD
pyplot()

include("Met_funcs.jl")

## Load world temperature data downloaded from ERA interim ##
folder = "data/"
file = "era_interim_moda_"

times = [collect(1980year:1month:(1990year - 1.0month)),
    collect(1990year:1month:(2000year - 1.0month)),
    collect(2000year:1month:(2010year - 1.0month)),
    collect(2010year:1month:(2018year- 1.0month))]

temp = extractERA(folder, file, "t2m", times)
plot(temp, 2012year + December)

# Load GBIF data in order to extract t2m values from GPS locations
#sol = JuliaDB.load("/Users/claireh/Documents/PhD/Data/GBIF/Solanum/output")
#sol = pushcol(sol, :id, collect(1:215440))
#coords = select(sol, (:decimallatitude, :decimallongitude))
#x = select(coords, :decimallongitude); y = select(coords, :decimallatitude)
#years = select(sol, :year)
#vals = extractvalues(x * °, y * °, years, temp, 1980, 2000)
#sol = pushcol(sol, :val, ustrip.(vals))
#JuliaDB.save(sol, "/Users/claireh/Documents/PhD/Data/GBIF/Solanum/era_output")

sol = JuliaDB.load("/Users/claireh/Documents/PhD/Data/GBIF/Solanum/era_output")

## Plot histograms of temperature values for subset of species ##
spp_names = ["Solanum dulcamara", "Solanum nigrum", "Solanum americanum","Solanum parvifolium"]
getprofile(spp_names, sol, "Temperature °C", (2,2))


numSpecies = 8; grd = (1,1); totalK = 250000.0kJ/km^2; area = 4.0km^2; req= 20.0kJ; individuals=10000
eco = create_eco(numSpecies, grd, area, totalK, req, individuals)
abun = runsim(eco, 2000months, 100)
plot_abun(abun, 10)

for i in [1, 2, 4, 8, 16]
    numSpecies = i; grd = (1,1); totalK = 250000.0kJ/km^2; area = 4.0km^2; req= 20.0kJ; individuals=10000
    eco = create_eco(numSpecies, grd, area, totalK, req, individuals)
    abun = runsim(eco, 2000months, 10)
    if i == 1
        display(plot([1/i], [mean(abun[:,:,1500:end,:])], seriestype = :scatter, label = "", ylim = (0, 5*10^4), xlim = (0,1.1), grid = false, ylab = "Average abundance at eqm", xlab= "1 / number of species"))
    else
        display(plot!([1/i], [mean(abun[:,:,1500:end,:])], seriestype = :scatter, label = ""))
    end
end

dir = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/wc"
srad = extractworldclim(joinpath(dir, "wc2.0_5m_srad"))
srad.array = srad.array[-10° .. 60°, 35° .. 80°]
meansrad = mean(srad.array[.!isnan.(srad.array)])
totalK = uconvert(kJ/m^2, meansrad * month)

numSpecies = 10000; grd = (10,10); req= 10.0kJ; individuals=1000000; area = 100.0*km^2; totalK = 100000.0kJ/km^2
eco = create_eco(numSpecies, grd, area, totalK, req, individuals)
abun = runsim(eco, 200months, 1)
JLD.save("10000species.jld", "abun", abun)
plot_abun(abun, 1:20, 1)

numSpecies = 5000; grd = (10,10); req= 10.0kJ; individuals=125000; area = 50.0*km^2; totalK = 25000.0kJ/km^2
eco = create_eco(numSpecies, grd, area, totalK, req, individuals)
abun = runsim(eco, 50years, 1)
JLD.save("5000species.jld", "abun", abun)
plot_abun(abun, 1:20, 1)

numSpecies = 200000; grd = (10,10); req= 10.0kJ; individuals=5000000; area = (250.0*km)^2; totalK = 1000.0kJ/km^2
eco = create_eco(numSpecies, grd, area, totalK, req, individuals)
abun = runsim(eco, 2years, 1)
JLD.save("200000species.jld", "abun", abun)
plot_abun(abun, 1:20, 1)
