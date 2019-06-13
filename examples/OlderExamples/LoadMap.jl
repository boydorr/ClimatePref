addprocs(8)
using JuliaDB
using ClimatePref
using RCall

gbif = load("GBIF_worldclim")
rf = collect(columns(gbif)[:refval])

ref = create_reference(1/12)
fill!(ref.array, 0)

for i in eachindex(rf)
    ref.array[rf[i]] += 1
end
ra = Array(ref.array)

@rput ra
R" library(raster); library(fields); library(viridis);
   ra = raster(ra); print(ra); png('gbif_plot.png', height = 4321, width = 2161)
   par(mar=c(0,0,0,0));image(log(1+ra), axes = F, xlab= '', ylab = '', col = c('black',viridis(100)[10:100]))
   dev.off()"
