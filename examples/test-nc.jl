using NetCDF
using Unitful
using Unitful.DefaultSymbols
using AxisArrays
using myunitful

# Load example dataset downloaded from EMCWF
ncinfo("/Users/claireh/Downloads/test.nc")

test = extractERA("/Users/claireh/Downloads/test.nc", "t2m",
    collect(1.0month:1month:2month))

### Load world temperature data downloaded from ERA interim ###

# Basic info about data
ncinfo("data/era_interim_moda_1990")
# Extract data for t2m parameter - temperature at 2m
dir = "data/era_interim_moda_1990"
tempax = extractERA(dir, "t2m", collect(1.0month:1month:10year))

# Collect lats and lons from array
lon = ustrip(axes(tempax.array, 1).val)
lat = ustrip(axes(tempax.array, 2).val)
# Strip array of units to plot
tempsC = ustrip(Array(tempax.array))
# Plot array as pdf with world shapefile on top
using RCall
@rput tempsC
@rput lat; @rput lon
R"pdf(file='plots/testERAint.pdf', paper = 'a4r', height= 8.27, width=11.69 )
library(rgdal)
library(fields)
world = readOGR('data/ne_10m_land/ne_10m_land.shp', layer='ne_10m_land')
image.plot(lon, lat, tempsC[,,1], xlab='', ylab='')
plot(world, add = T)
dev.off()"

# Load data for land cover
file = "data/World.tif"
world = extractfile(file)
wrld = ustrip.(world)
@rput wrld
R"image(wrld)"

# Use to mask out oceans from t2m array and plot as PDF
mask = isnan.(world)
tempsC[find(mask)] = NaN
@rput tempsC
@rput lat; @rput lon
R"pdf(file='plots/testERAint_mask.pdf', paper = 'a4r', height= 8.27, width=11.69 )
library(rgdal)
library(fields)
world = readOGR('data/ne_10m_land/ne_10m_land.shp', layer='ne_10m_land')
image.plot(lon, lat, tempsC[,,1], xlab='', ylab='')
plot(world, add = T)
dev.off()"


# Load GBIF data in order to extract t2m values from GPS locations
using JuliaDB
sol = load("/Users/claireh/Documents/PhD/Data/GBIF/Solanum/output")
sol = pushcol(sol, :id, collect(1:215440))
coords = select(sol, (:decimallatitude, :decimallongitude))
x = select(coords, :decimallongitude); y = select(coords, :decimallatitude)
years = select(sol, :year)
vals = extractvalues(x * °, y * °, years, tempax, 1990)

# Try plotting values for subset of species
spp_names = ["Solanum dulcamara", "Solanum nigrum", "Solanum americanum",
"Solanum parvifolium"]
for i in spp_names
spp = filter(p-> p[:species] == i, sol)
ids = select(spp, :id)
subvals = ustrip.(vals[ids])
res = vcat(subvals...)
res = res[!isnan.(res)]
using RCall
@rput res
@rput i
R"library(ggplot2);library(cowplot);library(gridExtra)
q = qplot(as.vector(res), geom='histogram',
            xlab ='Average temperature (deg C)')
pdf(file=paste('plots/era_',i,'.pdf', sep=''), paper = 'a4r', height= 8.27, width=11.69 )
print(q)
dev.off()
"
end
# Compare directly with worldclim
dir = "/Users/claireh/Documents/PhD/Data/Worldclim/wc2.0_5m/wc2.0_5m_tavg"
tavg = extractworldclim(dir)
vals = extractvalues(x * °, y * °, tavg, 1month:1month:12month)
for i in spp_names
spp = filter(p-> p[:species] == i, sol)
ids = select(spp, :id)
subvals = ustrip.(vals[ids, :])
res = vcat(subvals...)
res = res[!isnan.(res)]
using RCall
@rput res
@rput i
R"library(ggplot2);library(cowplot);library(gridExtra)
q = qplot(as.vector(res), geom='histogram',
            xlab ='Average temperature (deg C)')
pdf(file=paste('plots/worldclim_',i,'.pdf', sep=''), paper = 'a4r', height= 8.27, width=11.69 )
print(q)
dev.off()
"
end
