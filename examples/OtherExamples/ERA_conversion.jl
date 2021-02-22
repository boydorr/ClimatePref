using ClimatePref
using ClimatePref.Units
using JuliaDB
folder = "."
files = ["era_int_netsolar", "era_int_soiltemp1", "era_int_temp2m", "era_int_totalprec"]
params = ["ssr", "stl1", "t2m", "tp"]
times = [collect(1979year:1month:(1980year - 1.0month)),
    collect(1980year:1month:(1990year - 1.0month)),
    collect(1990year:1month:(2000year - 1.0month)),
    collect(2000year:1month:(2010year - 1.0month)),
    collect(2010year:1month:(2019year- 1.0month))]
for i in eachindex(files)
    era = readERA(folder, files[i], params[i], times)
    eraDB = era_to_DB(era)
    save(eraDB, files[i])
end

folder = "../../ECMWF"
files = ["era_int_netsolar_", "era_int_soiltemp1_", "era_int_temp2m_", "era_int_totalprec_"]
params = ["ssr", "stl1", "t2m", "tp"]
times = [collect(1979year:1month:(1980year - 1.0month)),
    collect(1980year:1month:(1990year - 1.0month)),
    collect(1990year:1month:(2000year - 1.0month)),
    collect(2000year:1month:(2010year - 1.0month)),
    collect(2010year:1month:(2019year- 1.0month))]
for i in eachindex(files)
    era = readERA(folder, files[i], params[i], times)
    gbif = extractvalues(gbif, era, Symbol(params[i]))
    save(eraDB, files[i])
end

using ClimatePref
using ClimatePref.Units
using JuliaDB
folder = "../../ECMWF"
files = ["era_int_netsolar", "era_int_soiltemp1", "era_int_temp2m", "era_int_totalprec"]
params = [:ssr, :stl1, :t2m, :tp]
era = load(files[1])
