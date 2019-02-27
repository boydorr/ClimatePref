yearStart = 1979
yearEnd = 1979

## Reanalysis products
# E.g. temperature at 2m above sea level and soil temperature at level 1.
tempat2m = "167.128" # use param names from ECMWF
soiltemp_lvl1 = "139.128"
retrieve_era_interim(tempat2m, yearStart, yearEnd,
    "era_int_temp2m")
retrieve_era_interim(soiltemp_lvl1, yearStart, yearEnd,
    'era_int_soiltemp1')

## Forecast products
# E.g. total precipitation and net solar radiation at the surface.
# For these I need to access a separate stream and model type.

totalprec = "228.128"
surfacenetsolar = "176.128"
stream = "mdfa" # monthly means of daily forecast accumulations
type ="fc"
retrieve_era_interim(totalprec, yearStart, yearEnd,
    stream, type, filename = "era_int_totalprec")
retrieve_era_interim(surfacenetsolar, yearStart, yearEnd,
    stream, type, filename = "era_int_netsolar")
