using ClimatePref
using MyUnitful
using JuliaDB
folder = "."
files = ["era_int_netsolar", "era_int_soiltemp1", "era_int_temp2m", "era_int_totalprec"]
params = ["ssr", "stl1", "t2m", "tp"]
times = [collect(1980year:1month:(1990year - 1.0month)),
    collect(1990year:1month:(2000year - 1.0month)),
    collect(2000year:1month:(2010year - 1.0month)),
    collect(2010year:1month:(2018year- 1.0month))]
for i in eachindex(files)
    era = readERA(folder, files[i], params[i], times)
    eraDB = era_to_DB(era)
    save(eraDB, files[i])
end
