# SPDX-License-Identifier: BSD-2-Clause

using ClimatePref
using Statistics
using JuliaDB
# Download CHELSA data from online server
# for var in ["prec", "tmean", "tmax", "tmin"]
#     for i in 1979:2013
#         for j in 1:12
#             j = @sprintf "%02d" j
#             download("https://www.wsl.ch/lud/chelsa/data/timeseries/$var/CHELSA_$var" *"_$i" *"_$j" *"_V1.2.1.tif", "/home/claireh/Documents/CHELSA/ts/$var/CHELSA_$var" *"_$i" *"_$j.tif")
#         end
#     end
# end

# Read in data and save to JuliaDB
dir = "CHELSA/ts"
varlist = ["prec", "tmax", "tmean", "tmin"]
aggfun = [sum, mean, mean, mean]
for i in eachindex(varlist)
    chelsa = readCHELSA(joinpath(dir, varlist[i]), varlist[i], res = 10,
                        fn = aggfun[i])
    chDB = ClimatePref.CHELSA_to_DB(chelsa)
    JuliaDB.save(chDB, joinpath(dir, varlist[i] * "DB"))
end
