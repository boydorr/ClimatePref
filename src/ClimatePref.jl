module ClimatePref

include("DataTypes.jl")
export Worldclim, Bioclim, ERA, Reference

include("GDAL.jl")
export read, searchdir, extractworldclim, extractbioclim, extractERA,
 extractvalues, extractfile

include("Phylo_models.jl")
export Brownian, Lambda, fitBrownian, fitLambda, varcovar

include("drop_tip.jl")
export getinternalnodes, drop_tip!

end
