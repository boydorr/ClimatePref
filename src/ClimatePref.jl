module ClimatePref

include("DataTypes.jl")
export Worldclim, Bioclim, ERA, CERA, Reference

include("GDAL.jl")
export read, searchdir, extractworldclim, extractbioclim, extractERA, extractCERA,
 extractvalues, extractfile

include("Tools.jl")
export create_reference, gardenmask, genus_worldclim_average,
    genus_worldclim_monthly, upresolution

include("Phylo_models.jl")
export Brownian, Lambda, fitBrownian, fitLambda, varcovar

include("drop_tip.jl")
export getinternalnodes, drop_tip!

end
