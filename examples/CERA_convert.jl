# 6. Convert CERA data into JuliaDB tables
# and combine with era
using ClimatePref
using ClimatePref.Units
using Unitful
using JuliaDB
using JuliaDBMeta
using AxisArrays

# Save cera files as JuliaDB tables
folder = "."
files = ["cera_20c_soiltemp1_", "cera_20c_soiltemp2_", "cera_20c_soiltemp3_", "cera_20c_soiltemp4_", "cera_20c_soilwater1_", "cera_20c_soilwater2_", "cera_20c_soilwater3_", "cera_20c_soilwater4_", "cera_20c_surfacenetsolar_", "cera_20c_temp2m_", "cera_20c_totalprec_"]
params = ["stl1","stl2","stl3","stl4","swvl1","swvl2","swvl3","swvl4","ssr", "t2m", "tp"]
for i in eachindex(files)[1:end]
    cera = readCERA(folder, files[i], params[i])
    AxisArrays.axes(cera.array)[3].val .-= 1month
    ceraDB = era_to_DB(cera)
    save(ceraDB, files[i] * "new")
end

function combine_cera(cera::IndexedTable)
    for i in eachindex(files)[2:end]
        new_cera = load(files[i] * "new")
        vals = select(new_cera, :val)
        cera = pushcol(cera, Symbol(params[i]), vals)
    end
    return cera
end

cera = load(files[1] * "new")
cera = renamecol(cera, :val, Symbol(params[1]))
cera = combine_cera(cera)
cera = select(cera, (:x, :y, :month, :year, :refval, :stl1, :stl2, :stl3, :stl4, :swvl1, :swvl2, :swvl3, :swvl4, :ssr, :t2m, :tp))
cera = reindex(cera, (:year, :refval))
save(cera, "cera_all_new")

# Merge with ERA interim for combined data
era = load("../ECMWF/era_int_all")
cera = load("cera_all_new")
cera_filter = filter(c-> c.year < 1979, cera)
cera_era = merge(cera_filter, era)
save(cera_era, "CERA_ERA_new")

# Create date columns (one for each month of data being extracted)
yr = select(cera_era, :year) .* year
mth = select(cera_era, :month) .* month
function era_extend(era::IndexedTable)
    for i in 0:11
        dt = yr .+ mth
        dt .+= (i * month)
        dt = uconvert.(year, dt)
        era = pushcol(era, Symbol("date$i"), dt)
    end
    return era
end
cera_era_date = era_extend(cera_era)
save(cera_era_date, "CERA_ERA_new")
