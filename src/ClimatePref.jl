module ClimatePref

include("DataTypes.jl")
export Worldclim, Bioclim, ERA, Reference

include("GDAL.jl")
export read, searchdir, extractworldclim, extractbioclim, extractERA, extractvalues

include("Phylo_models.jl")
export Brownian, Lambda, fitBrownian, fitLambda, varcovar

end
