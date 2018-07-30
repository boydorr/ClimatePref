using ClimatePref
using AxisArrays
using MyUnitful
using Unitful

ClimatePref.extractCERA("data", "cera_20C", "t2m")



test = ustrip.(cera1900.array[:,:,1])
using RCall
@rput test
R"library(viridis);library(fields)
    image.plot(test, col = magma(20))"
