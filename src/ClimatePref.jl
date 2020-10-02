module ClimatePref

using Requires
function __init__()
    @require PyCall="438e738f-606a-5dbb-bf0a-cddfbfd45ab0" begin
        println("Creating ECMWF interface ...")
        include("ERA_interim_tools.jl")
        export retrieve_era_interim
        include("ECMWF_tools.jl")
        export retrieve_ECMWF
    end
end

module Units

import Unitful
using Unitful: @unit
day = Unitful.d
week = Unitful.wk
@unit month "month" Month 2.628e6 * Unitful.s false
@unit year "year" Year 31536000 * Unitful.s false

const days = day
const weeks = week
const months = month
const years = year
const January = 0month
const February = 1month
const March = 2months
const April = 3months
const May = 4months
const June = 5months
const July = 6months
const August = 7months
const September = 8months
const October = 9months
const November = 10months
const December = 11months

const localunits = Unitful.basefactors
function __init__()
    merge!(Unitful.basefactors, localunits)
    Unitful.register(Units)
end

export day, days, week, weeks, month, months, year, years, Rates,
January, February, March, April, May, June, July, August,
September, October, November, December

end

include("ClimateTypes.jl")
export Worldclim, Bioclim, ERA, CERA, Reference

include("ReadData.jl")
export read, searchdir, readworldclim, readbioclim, readERA, readCERA, readfile, readCHELSA

include("ReadGBIF.jl")
export ReadGBIF

include("ReadTPL.jl")
export ReadTPL

include("ExtractClimate.jl")
export extractvalues

include("DataCleaning.jl")
export create_reference, gardenmask, genus_worldclim_average,
    genus_worldclim_monthly, upresolution, downresolution, mask

include("Conversion.jl")
export worldclim_to_DB, era_to_DB, CHELSA_to_DB

include("Plotting.jl")
export getprofile

include("PhyloModels.jl")
export Brownian, Lambda, fitBrownian, fitLambda, varcovar

end
