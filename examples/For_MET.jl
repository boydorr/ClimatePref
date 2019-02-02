############### TEST DATA FROM ERA INTERIM - TEMPERATURE AT 2M #################
############# LOADED AND VALUES EXTRACTED FOR SUBSET OF GBIF DATA ##############

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

### Load world temperature data downloaded from ERA interim ###
folder = "data/"
file = "era_interim_moda_"

times = [collect(1980year:1month:(1990year - 1.0month)),
    collect(1990year:1month:(2000year - 1.0month)),
    collect(2000year:1month:(2010year - 1.0month)),
    collect(2010year:1month:(2018year- 1.0month))]

temp = extractERA(folder, file, "t2m", times)
plot(temp, 2000year+11month)

sol = JuliaDB.load("/Users/claireh/Documents/PhD/Data/GBIF/Solanum/era_output")

# Try plotting histograms of values for subset of species
spp_names = ["Solanum dulcamara", "Solanum nigrum", "Solanum americanum","Solanum parvifolium"]
getprofile(spp_names, sol, "Temperature °C", (2,2))


# Load GBIF data in order to extract t2m values from GPS locations
#sol = JuliaDB.load("/Users/claireh/Documents/PhD/Data/GBIF/Solanum/output")
#sol = pushcol(sol, :id, collect(1:215440))
#coords = select(sol, (:decimallatitude, :decimallongitude))
#x = select(coords, :decimallongitude); y = select(coords, :decimallatitude)
#years = select(sol, :year)
#vals = extractvalues(x * °, y * °, years, temp, 1980, 2000)
#sol = pushcol(sol, :val, ustrip.(vals))
#JuliaDB.save(sol, "/Users/claireh/Documents/PhD/Data/GBIF/Solanum/era_output")
