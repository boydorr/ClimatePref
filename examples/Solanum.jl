using Unitful
using AxisArrays
using ClimatePref
using myunitful

import Unitful: °, °C, mm
import ArchGDAL
import Base.read
const AG = ArchGDAL

using JuliaDB
#sol = loadtable("/Users/claireh/Documents/PhD/Data/GBIF/Solanum",
#       indexcols = [:phylum, :class, :order, :family, :genus, :species, :scientificname],
#       type_detect_rows =5000,
#       colparsers=Dict(:datasetkey=>String,
#                       :occurrenceid=>String,
#                       :gbifid=>String,
#                       :locality=>String,
#                       :publishingorgkey=>String,
#                       :taxonkey=>String,
#                       :institutioncode=>String,
#                       :catalognumber=>String,
#                       :recordnumber=>String))
#save(gbif, "/Users/claireh/Documents/PhD/Data/GBIF/Solanum/output")
sol = load("/Users/claireh/Documents/PhD/Data/GBIF/Solanum/output")

coords = select(sol, (:decimallatitude, :decimallongitude))

dir = "/Users/claireh/Documents/PhD/Data/Worldclim/wc2.0_5m"
tavg = extractworldclim(joinpath(dir, "wc2.0_5m_tavg"))
tmax = extractworldclim(joinpath(dir, "wc2.0_5m_tmax"))
prec = extractworldclim(joinpath(dir, "wc2.0_5m_prec"))
srad = extractworldclim(joinpath(dir, "wc2.0_5m_srad"))
vapr = extractworldclim(joinpath(dir, "wc2.0_5m_vapr"))
wind = extractworldclim(joinpath(dir, "wc2.0_5m_wind"))
bio = extractbioclim(joinpath(dir, "wc2.0_5m_bio"))
x = select(coords, :decimallongitude); y = select(coords, :decimallatitude)
vals_tavg = extractvalues(x * °, y * °, tavg, 1month:1month:12month)
vals_tmax = extractvalues(x * °, y * °, tmax, 1month:1month:12month)
vals_prec = extractvalues(x * °, y * °, prec, 1month:1month:12month)
vals_srad = extractvalues(x * °, y * °, srad, 1month:1month:12month)
vals_vapr = extractvalues(x * °, y * °, vapr, 1month:1month:12month)
vals_wind = extractvalues(x * °, y * °, wind, 1month:1month:12month)
vals_bio = extractvalues(x * °, y * °, bio, 1:1:19)

na_mean(x) = mean(x[.!isnan.(x)])
na_max(x) = maximum(x[.!isnan.(x)])
na_min(x) = minimum(x[.!isnan.(x)])
na_var(x) = var(x[.!isnan.(x)])
sol = pushcol(sol, :id, collect(1:215440))
names = map(searchdir(dir, "wc2.0_5m_")) do str
    Symbol(split(str, "wc2.0_5m_")[2])
end

## Run through 4 tester species and plot histograms of 6 worldclim
## variables for each
spp_names = ["Solanum dulcamara", "Solanum nigrum", "Solanum americanum",
"Solanum parvifolium"]
for i in spp_names
    # Select species and the id of the rows
    spp = filter(p-> p[:species] == i, sol)
    ids = select(spp, :id)
    # Select values for row ids from each variable
    res_prec = ustrip(vals_prec[ids])
    res_srad = ustrip(vals_srad[ids])
    res_tavg = ustrip(vals_tavg[ids])
    res_tmax = ustrip(vals_tmax[ids])
    res_vapr = ustrip(vals_vapr[ids])
    res_wind = ustrip(vals_wind[ids])
    # Combine into one array
    res = [res_prec res_srad res_tavg res_tmax res_vapr res_wind]
    using RCall
    @rput res
    @rput i
    # Plot in R
    R"library(ggplot2);library(cowplot);library(gridExtra)
    res[, 1][res[ , 1]==-32768] = NA
    q1 = qplot(as.vector(res[, 1]), geom='histogram',
                xlab ='Average precipitation (mm)')
    q2 = qplot(as.vector(res[, 2]), geom='histogram',
                xlab ='Solar radiation (kJ m-2 day-1)')
    q3 = qplot(as.vector(res[, 3]), geom='histogram',
                xlab ='Average temperature (deg C)')
    q4 = qplot(as.vector(res[, 4]), geom='histogram',
                xlab ='Maximum temperature (deg C)')
    q5 = qplot(as.vector(res[, 5]), geom='histogram',
                xlab ='Water vapour pressure (kPa)')
    q6 = qplot(as.vector(res[, 6]), geom='histogram',
                xlab ='Wind speed (m s-1)')
    pdf(file=paste('plots/new',i,'.pdf', sep=''), paper = 'a4r', height= 8.27, width=11.69 )
    grid.arrange(q1, q2, q3, q4, q5, q6, nrow=2)
    dev.off()
    "
end

using RCall
tav = ustrip(tavg.array)
@rput tav
R"library(fields);
for (i in 1:12){
    #jpeg(paste('plots/tavg', i, '.jpeg',sep=''), height= 595, width=842)
    image.plot(tav[,,i])
    #dev.off()
    }"

pts = [x[1:4],y[1:4]]
thisstep = ustrip(step(axes(tavg.array, 1).val))
@rput pts; @rput thisstep
R"library(fields);
#jpeg('plots/points.jpeg', height= 595, width=842)
image.plot(seq(-180, 180, by=thisstep),
seq(-90, 90, by=thisstep), tavg[,,1], xlab='', ylab='')
points(pts[[1]], pts[[2]], pch=20)
#dev.off()
"
