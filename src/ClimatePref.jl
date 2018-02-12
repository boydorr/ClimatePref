module ClimatePref

include("DataTypes.jl")
export Worldclim, Bioclim, ERA

include("GDAL.jl")
export read, searchdir, extractworldclim, extractbioclim, extractERA, extractvalues


end
