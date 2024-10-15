# SPDX-License-Identifier: BSD-2-Clause

# 5. Convert ERA data to JuliaDB tables
using ClimatePref
using ClimatePref.Units
using JuliaDB

folder = "."
files = [
    "era_int_netsolar_",
    "era_int_temp2m_",
    "era_int_totalprec_",
    "era_int_soiltemp1_",
    "era_int_soiltemp2_",
    "era_int_soiltemp3_",
    "era_int_soiltemp4_",
    "era_int_soilwater1_",
    "era_int_soilwater2_",
    "era_int_soilwater3_",
    "era_int_soilwater4_"
]
params = [
    "ssr",
    "t2m",
    "tp",
    "stl1",
    "stl2",
    "stl3",
    "stl4",
    "swvl1",
    "swvl2",
    "swvl3",
    "swvl4"
]
times = [collect((1979year):(1month):(1980year - 1.0month)),
    collect((1980year):(1month):(1990year - 1.0month)),
    collect((1990year):(1month):(2000year - 1.0month)),
    collect((2000year):(1month):(2010year - 1.0month)),
    collect((2010year):(1month):(2019year - 1.0month))]
for i in eachindex(files)
    era = readERA(folder, files[i], params[i], times)
    eraDB = era_to_DB(era)
    save(eraDB, files[i] * "new")
end

# Combine all data together
function combine_era(era::IndexedTable)
    for i in eachindex(files)[2:end]
        new_era = load(files[i] * "new")
        vals = select(new_era, :val)
        era = pushcol(era, Symbol(params[i]), vals)
    end
    return era
end

era = load(files[1] * "new")
era = renamecol(era, :val, Symbol(params[1]))
era = combine_era(era)
era = select(era,
             (:year, :refval, :x, :y, :month, :stl1, :stl2, :stl3, :stl4,
              :swvl1, :swvl2, :swvl3, :swvl4, :ssr, :t2m, :tp))
era = reindex(era, (:year, :refval))
save(era, "era_int_all")
