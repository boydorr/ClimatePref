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

dir = "/Users/claireh/Documents/PhD/Data/Worldclim/wc2.0_5m/"
folders = searchdir(dir, "wc2.0_5m")
datasets = extractfolders(dir, folders)
x = select(coords, :decimallongitude); y = select(coords, :decimallatitude)
vals = extractvalues(x * °, y * °, datasets)

na_mean(x) = mean(x[.!isnan.(x)])
na_max(x) = maximum(x[.!isnan.(x)])
na_min(x) = minimum(x[.!isnan.(x)])
na_var(x) = var(x[.!isnan.(x)])
sol = pushcol(sol, :id, collect(1:215440))
names = map(folders) do str
    Symbol(split(str, "wc2.0_5m_")[2])
end
for i in 1:length(names)
    sol = pushcol(sol, names[i] , mapslices(na_mean, vals[i], 2)[:,1])
end
spp_names = ["Solanum dulcamara", "Solanum nigrum", "Solanum americanum",
"Solanum parvifolium"]
for i in spp_names
spp = filter(p-> p[:species] == i, sol)
ids = select(spp, :id)
res = map( x-> ustrip(vals[x][ids]), 1:7)
using RCall
@rput res
@rput i
R"library(ggplot2);library(cowplot);library(gridExtra)
res[[2]][res[[2]]==-32768] = NA
q2 = qplot(as.vector(res[[2]]), geom='histogram',
            xlab ='Average precipitation (mm)')
q3 = qplot(as.vector(res[[3]]), geom='histogram',
            xlab ='Solar radiation (kJ m-2 day-1)')
q4 = qplot(as.vector(res[[4]]), geom='histogram',
            xlab ='Average temperature (°C)')
q5 = qplot(as.vector(res[[5]]), geom='histogram',
            xlab ='Maximum temperature (°C)')
q6 = qplot(as.vector(res[[6]]), geom='histogram',
            xlab ='Water vapour pressure (kPa)')
q7 = qplot(as.vector(res[[7]]), geom='histogram',
            xlab ='Wind speed (m s-1)')
pdf(file=paste('plots/',i,'.pdf', sep=''), paper = 'a4r', height= 8.27, width=11.69 )
grid.arrange(q2, q3, q4, q5, q6, q7, nrow=2)
dev.off()
"
end
using RCall
tavg = datasets[4]
tavg = ustrip(tavg)
@rput tavg
R"library(fields);
for (i in 1:12){
    jpeg(paste('plots/tavg', i, '.jpeg',sep=''), height= 595, width=842)
    image.plot(tavg[,,i])
    dev.off()}"
pts = [x[1:4],y[1:4]]
thisstep = ustrip(step(axes(datasets[1], 1).val))
@rput pts; @rput thisstep
R"library(fields);
jpeg('plots/points.jpeg', height= 595, width=842)
image.plot(seq(-180, 180, by=thisstep),
seq(-90, 90, by=thisstep), tavg[,,1], xlab='', ylab='')
points(pts[[1]], pts[[2]], pch=20)
dev.off()"
