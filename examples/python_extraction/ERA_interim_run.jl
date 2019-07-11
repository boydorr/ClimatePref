using ClimatePref
using PyCall

yearStart = 1979
yearEnd = 2018

## Reanalysis products
# E.g. temperature at 2m above sea level and soil temperature at level 1.
tempat2m = "167.128" # use param names from ECMWF
soiltemp_lvl1 = "139.128"
soiltemp_lvl2 = "170.128"
soiltemp_lvl3 = "183.128"
soiltemp_lvl4 = "236.128"
soilwater_lv1 = "39.128"
soilwater_lv2 = "40.128"
soilwater_lv3 = "41.128"
soilwater_lv4 = "42.128"
retrieve_ECMWF(tempat2m, yearStart, yearEnd,
    "era_int_temp2m")
retrieve_ECMWF(soiltemp_lvl1, yearStart, yearEnd,
    "era_int_soiltemp1")
retrieve_ECMWF(soiltemp_lvl2, yearStart, yearEnd,
    "era_int_soiltemp2")
retrieve_ECMWF(soiltemp_lvl3, yearStart, yearEnd,
    "era_int_soiltemp3")
retrieve_ECMWF(soiltemp_lvl4, yearStart, yearEnd,
    "era_int_soiltemp4")
retrieve_ECMWF(soilwater_lv1, yearStart, yearEnd,
    "era_int_soilwater1")
retrieve_ECMWF(soilwater_lv2, yearStart, yearEnd,
    "era_int_soilwater2")
retrieve_ECMWF(soilwater_lv3, yearStart, yearEnd,
    "era_int_soilwater3")
retrieve_ECMWF(soilwater_lv4, yearStart, yearEnd,
    "era_int_soilwater4")
## Forecast products
# E.g. total precipitation and net solar radiation at the surface.
# For these I need to access a separate stream and model type.

totalprec = "228.128"
surfacenetsolar = "176.128"
stream = "mdfa" # monthly means of daily forecast accumulations
type ="fc"
retrieve_ECMWF(totalprec, yearStart, yearEnd, "era_int_totalprec",
    stream = stream, modeltype = type)
retrieve_ECMWF(surfacenetsolar, yearStart, yearEnd, "era_int_netsolar", stream = stream, modeltype = type)

# CERA 20C

yearStart = 1901
yearEnd = 2010

## Reanalysis products
# E.g. temperature at 2m above sea level and soil temperature at level 1.
tempat2m = "167.128" # use param names from ECMWF
soiltemp_lvl1 = "139.128"
soiltemp_lvl2 = "170.128"
soiltemp_lvl3 = "183.128"
soiltemp_lvl4 = "236.128"
soilwater_lv1 = "39.128"
soilwater_lv2 = "40.128"
soilwater_lv3 = "41.128"
soilwater_lv4 = "42.128"
retrieve_ECMWF(tempat2m, yearStart, yearEnd,
    "cera_20c_temp2m", eclass = "ep", dataset = "cera20c", stream = "edmo")
retrieve_ECMWF(soiltemp_lvl1, yearStart, yearEnd,
    "cera_20c_soiltemp1", eclass = "ep", dataset = "cera20c", stream = "edmo")
retrieve_ECMWF(soiltemp_lvl2, yearStart, yearEnd,
    "cera_20c_soiltemp2", eclass = "ep", dataset = "cera20c", stream = "edmo")
retrieve_ECMWF(soiltemp_lvl3, yearStart, yearEnd,
    "cera_20c_soiltemp3", eclass = "ep", dataset = "cera20c", stream = "edmo")
retrieve_ECMWF(soiltemp_lvl4, yearStart, yearEnd,
    "cera_20c_soiltemp4", eclass = "ep", dataset = "cera20c", stream = "edmo")
retrieve_ECMWF(soilwater_lv1, yearStart, yearEnd,
    "cera_20c_soilwater1", eclass = "ep", dataset = "cera20c", stream = "edmo")
retrieve_ECMWF(soilwater_lv2, yearStart, yearEnd,
    "cera_20c_soilwater2", eclass = "ep", dataset = "cera20c", stream = "edmo")
retrieve_ECMWF(soilwater_lv3, yearStart, yearEnd,
    "cera_20c_soilwater3", eclass = "ep", dataset = "cera20c", stream = "edmo")
retrieve_ECMWF(soilwater_lv4, yearStart, yearEnd,
    "cera_20c_soilwater4", eclass = "ep", dataset = "cera20c", stream = "edmo")
## Forecast products
# E.g. total precipitation and net solar radiation at the surface.
# For these I need to access a separate stream and model type.

totalprec = "228.128"
surfacenetsolar = "176.128"
#stream = "mdfa" # monthly means of daily forecast accumulations
type ="fc"
retrieve_ECMWF(totalprec, yearStart, yearEnd, "cera_20c_totalprec",
    modeltype = type, eclass = "ep", dataset = "cera20c", stream = "edmo")
retrieve_ECMWF(surfacenetsolar, yearStart, yearEnd, "cera_20c_surfacenetsolar", modeltype = type, eclass = "ep", dataset = "cera20c", stream = "edmo")
