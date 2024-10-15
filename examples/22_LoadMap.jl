using JuliaDB
using ClimatePref
using Unitful
using Unitful.DefaultSymbols
using RCall

gbif = load("GBIF_TPL_new")
lat = collect(select(gbif, :decimallatitude))
lon = collect(select(gbif, :decimallongitude))
ref = create_reference(0.08983153)
exclude = (lat .* ° .<= ref.array.axes[2].val[end]) .& (lon .* ° .< ref.array.axes[1].val[end])
refval = extractvalues(lon[exclude] .* °, lat[exclude] .* °, ref)
fill!(ref.array, 0)

for i in eachindex(refval)
   ref.array[refval[i]] += 1
end
ra = transpose(Array(ref.array))
ra = ra[end:-1:1, :]
@rput ra

R" library(raster); library(fields); library(viridis);library(rgdal);library(rasterVis)
    world = readOGR('../WorldShp/', 'world')
    evi = raster('../Evi.tif')
    ra = raster(ra, xmn = extent(evi)[1], xmx = extent(evi)[2], ymn = extent(evi)[3], ymx = extent(evi)[4], crs = crs(evi))
    print(ra)
    evi[evi[] == 0] = NA
    #evi[(evi[]) == 0] = NA
    effort = ra/evi
    effort = aggregate(effort, fact = 4, fun = mean)
    print(effort)
    effort = projectRaster(effort, crs = ('+proj=moll +lon_0=0 +x_0=0 +y_0=0'))
    world = spTransform(world, CRS = ('+proj=moll +lon_0=0 +x_0=0 +y_0=0'))
    l = levelplot(log(1 + effort),
          margin = list(FUN = 'sum'), colorkey=list(space='bottom'),
          par.settings=viridisTheme) +
  layer(sp.polygons(world, lwd=1))
  png('Figure_1.png', width = 11, height = 8, units = 'in', res = 300)
  #pdf('Figure_1.pdf', paper = 'a4r', width = 11, height = 8, dpi = 300)
    print(l)
  dev.off()
  "