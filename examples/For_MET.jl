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
using Diversity
using Diversity.ShortNames
using Diversity.ShortNames: γ
pyplot()

include("Met_funcs.jl")

## Load world temperature data downloaded from ERA interim ##
folder = "data/"
file = "era_interim_moda_"

times = [collect(1980year:1month:(1990year - 1.0month)),
    collect(1990year:1month:(2000year - 1.0month)),
    collect(2000year:1month:(2010year - 1.0month)),
    collect(2010year:1month:(2018year- 1.0month))]

temp = readERA(folder, file, "t2m", times)
plot(temp, 2002year + November)

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

## 10,000 species for 200 months

numSpecies = 10000; grd = (10,10); req= 10.0kJ; individuals=1000000; area = 100.0*km^2; totalK = 100000.0kJ/km^2
#eco = create_eco(numSpecies, grd, area, totalK, req, individuals)
#abun = runsim(eco, 200months, 1)
#JLD.save("10000species.jld", "abun", abun)
abun = JLD.load("data/10000species.jld", "abun")

# Plot abundance over time for first twenty species
plot_abun(abun, grd, 1:20, 1)

# Plot species richness over time
heatmap(mapslices(x->sum(x.>0), abun, dims = 1)[1,:,:,1], xlab = "Time (months)", ylab = "Space")

# Animation abundance over time for first species
for i in 1:5:200
    display(heatmap(reshape(abun[1,:,i,1], 10, 10), clim = (0, 10)))
end
# Animation beta bar over time
for i in 1:10:200
    meta = Metacommunity(abun[:,:,i,1])
    b = reshape(subdiv(β̄(meta), 1)[:diversity], 10, 10)
    B = metadiv(β̄(meta), 1)[:diversity]
    if i ==1
        display(heatmap(b, layout = (2, 1), clim = (1, 4), title = "Subcommunity β̄"))
    else
        display(heatmap!(b, subplot = 1,  clim = (1, 4)))
    end
    display(plot!([i], B, seriestype = :scatter, subplot = 2, xlim = (0, 200), ylim = (1, 4), grid = false, label = "", title = "Metacommunity β̄"))
end
# Animation rho bar over time
for i in 1:10:200
    meta = Metacommunity(abun[:,:,i,1])
    b = reshape(subdiv(ρ̄(meta), 1)[:diversity], 10, 10)
    B = metadiv(ρ̄(meta), 1)[:diversity]
    if i ==1
        display(heatmap(b, layout = (2, 1), clim = (0.2, 0.6), title = "Subcommunity ρ̄"))
    else
        display(heatmap!(b, subplot = 1,  clim = (0.2, 0.6)))
    end
    display(plot!([i], B, seriestype = :scatter, subplot = 2, xlim = (0, 200), ylim = (0, 1), grid = false, label = "", title = "Metacommunity ρ̄"))
end

## 5,000 species for 50 years
numSpecies = 5000; grd = (10,10); req= 10.0kJ; individuals=125000; area = 50.0*km^2; totalK = 25000.0kJ/km^2
#eco = create_eco(numSpecies, grd, area, totalK, req, individuals)
#abun = runsim(eco, 50years, 1)
#JLD.save("5000species.jld", "abun", abun)
abun = JLD.load("data/5000species.jld", "abun")

# Plot abundance over time for twenty species
plot_abun(abun, grd, 740:760, 1)

# Plot species richness over time
heatmap(mapslices(x->sum(x.>0), abun, dims = 1)[1,:,:,1], xlab = "Time (months)", ylab = "Space")

# Animation abundance over time for first species
for i in 1:10:600
    display(heatmap(reshape(abun[1,:,i,1], 10, 10), clim = (0, 10)))
end

# Animation alpha bar and gamma
for i in 1:20:600
    meta = Metacommunity(abun[:,:,i,1])
    b = reshape(subdiv(ᾱ(meta), 1)[:diversity], 10, 10)
    B = metadiv(γ(meta), 1)[:diversity]
    if i ==1
        display(heatmap(b, layout = (2, 1), clim = (1, 1000), title = "Subcommunity ᾱ"))
    else
        display(heatmap!(b, subplot = 1,  clim = (1, 1000)))
    end
    display(plot!([i], B, seriestype = :scatter, subplot = 2, xlim = (0, 600), ylim = (1, 5000), grid = false, label = "", color = :1, title = "Metacommunity γ"))
end

## 200,000 species over 625,000km² area
numSpecies = 200000; grd = (10,10); req= 10.0kJ; individuals=5000000; area = (250.0*km)^2; totalK = 1000.0kJ/km^2
#eco = create_eco(numSpecies, grd, area, totalK, req, individuals)
#abun = runsim(eco, 2years, 1)
#JLD.save("200000species.jld", "abun", abun)
abun = JLD.load("data/200000species.jld", "abun")

# Plot abundance over time for twenty species
plot_abun(abun, grd, 1:20, 1)

# Plot species richness over time
heatmap(mapslices(x->sum(x.>0), abun, dims = 1)[1,:,:,1], xlab = "Time (months)", ylab = "Space")

## 200,000km² with 3000 species for 100 years (30 by 10)
