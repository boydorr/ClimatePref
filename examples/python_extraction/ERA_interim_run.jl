using ClimatePref
using PyCall

yearStart = 1979
yearEnd = 2018

## Reanalysis products
# E.g. temperature at 2m above sea level and soil temperature at level 1.
tempat2m = "167.128" # use param names from ECMWF
soiltemp_lvl1 = "139.128"
retrieve_era_interim(tempat2m, yearStart, yearEnd,
    "era_int_temp2m")
retrieve_era_interim(soiltemp_lvl1, yearStart, yearEnd,
    "era_int_soiltemp1")
retrieve_era_interim(soiltemp_lvl2, yearStart, yearEnd,
    "era_int_soiltemp2")
retrieve_era_interim(soiltemp_lvl3, yearStart, yearEnd,
    "era_int_soiltemp3")
retrieve_era_interim(soiltemp_lvl4, yearStart, yearEnd,
    "era_int_soiltemp4")
## Forecast products
# E.g. total precipitation and net solar radiation at the surface.
# For these I need to access a separate stream and model type.

totalprec = "228.128"
surfacenetsolar = "176.128"
stream = "mdfa" # monthly means of daily forecast accumulations
type ="fc"
retrieve_era_interim(totalprec, yearStart, yearEnd, "era_int_totalprec",
    stream = stream, modeltype = type)
retrieve_era_interim(surfacenetsolar, yearStart, yearEnd, "era_int_netsolar", stream = stream, modeltype = type)
