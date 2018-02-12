using NetCDF
using Unitful
using Unitful.DefaultSymbols
using AxisArrays
using myunitful

ncinfo("/Users/claireh/Downloads/test.nc")

x = ncread("/Users/claireh/Downloads/test.nc", "t2m")

y = x * 1.0

y[y .≈ -32767] = NaN
temps = y .* 0.0009362609303530759K .+ 270.9726153656286K

ncinfo("data/era_interim_moda_1990")
lat = reverse(ncread("data/era_interim_moda_1990", "latitude"))
lon = ncread("data/era_interim_moda_1990", "longitude")
twomtemp = ncread("data/era_interim_moda_1990", "t2m")
twomtemp = twomtemp * 1.0

splitval = find(lon .== 180)[1]
firstseg = collect((splitval+1):size(twomtemp,1))
secondseg = collect(1:splitval)
twomtemp = twomtemp[vcat(firstseg ,secondseg), :, :]
lon[firstseg] = 180 - lon[firstseg]
lon = lon[vcat(reverse(firstseg) ,secondseg)]

twomtemp[twomtemp .≈ -32767] = NaN
temps = twomtemp .* 0.0017312391138308897K .+ 270.9726153656286K
tempsC = uconvert.(°C, temps)
tempsC = ustrip.(tempsC)[:,end:-1:1, :]
step1 = 0.75°
step2 = 180°/241
tempax = AxisArray(uconvert.(°C, temps),
                       Axis{:longitude}(lon * °),
                       Axis{:latitude}(lat * °),
                       Axis{:time}(1month:1month:10year))
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

file = "data/World.tif"
world = extractfile(file)
wrld = ustrip.(world)
@rput wrld
R"image(wrld)"

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

using JuliaDB
sol = load("/Users/claireh/Documents/PhD/Data/GBIF/Solanum/output")
coords = select(sol, (:decimallatitude, :decimallongitude))
x = select(coords, :decimallongitude); y = select(coords, :decimallatitude)
yr = select(sol, :year); mn = select(sol, :month)
vals = map(x * °, y * °, year) do x, y, yr
    extractvalues(x, y, tempax, yr)
end
