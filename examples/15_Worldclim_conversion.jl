# 15. Convert worldclim data to JuliaDB format
using Unitful
using AxisArrays
using ClimatePref
using ClimatePref.Units
using JuliaDB

dir = "/Users/claireh/Documents/PhD/Data/Worldclim/wc2.0_5m"
tavg = readworldclim(joinpath(dir, "wc2.0_5m_tavg"))
tmax = readworldclim(joinpath(dir, "wc2.0_5m_tmax"))
tmin = readworldclim(joinpath(dir, "wc2.0_5m_tmin"))
prec = readworldclim(joinpath(dir, "wc2.0_5m_prec"))
srad = readworldclim(joinpath(dir, "wc2.0_5m_srad"))
vapr = readworldclim(joinpath(dir, "wc2.0_5m_vapr"))
wind = readworldclim(joinpath(dir, "wc2.0_5m_wind"))
bio = readbioclim(joinpath(dir, "wc2.0_5m_bio"))


worldclim_to_DB(tavg)


using ClimatePref
using ClimatePref.Units
using JuliaDB
folder = "wc2.0_5m"
files = ["wc2.0_5m_tavg", "wc2.0_5m_tmax", "wc2.0_5m_tmin", "wc2.0_5m_prec", "wc2.0_5m_srad", "wc2.0_5m_vapr", "wc2.0_5m_wind"]
params = ["tavg", "tmax", "tmin", "prec", "srad", "vapr", "wind"]
for i in eachindex(files)
    wc = readworldclim(joinpath(folder, files[i]))
    wcDB = worldclim_to_DB(wc)
    save(wcDB, joinpath("Worldclim", params[i]))
end
