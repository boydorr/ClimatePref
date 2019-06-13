using Unitful
using AxisArrays
using ClimatePref
using MyUnitful
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
