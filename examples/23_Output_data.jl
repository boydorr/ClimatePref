# SPDX-License-Identifier: BSD-2-Clause

using PhyloNetworks
using GLM
using JuliaDB
using JuliaDBMeta
using Unitful
using ClimatePref
using ClimatePref.Unitful
using StatsBase
using JLD
using OnlineStats
using Unitful.DefaultSymbols
using Statistics
using DataFrames
using CSV
@everywhere using Unitful
@everywhere using Unitful.DefaultSymbols
@everywhere using OnlineStats
@everywhere using StatsBase
@everywhere using ClimatePref

cera_taxo = JuliaDB.load("CERA_taxo")

mins = [
    197.0K,
    197.0K,
    197.0K,
    0K,
    197.0K,
    197.0K,
    197.0K,
    197.0K,
    0.0m^3,
    0.0m^3,
    0.0m^3,
    0.0m^3,
    0.0J / m^2,
    0.0m
]
maxs = [
    320.0K,
    320.0K,
    320.0K,
    80K,
    320.0K,
    320.0K,
    320.0K,
    320.0K,
    1.0m^3,
    1.0m^3,
    1.0m^3,
    1.0m^3,
    3.0e7J / m^2,
    0.1m
]

# Load EVi and gbif counts
total_evi_counts = JLD.load("Total_evi_counts.jld", "total")
total_gbif_counts = JLD.load("Total_gbif_counts.jld", "total")

# Adjustment - remove NaNs and Infs
adjustment = total_gbif_counts ./ total_evi_counts
adjustment[isnan.(adjustment)] .= 1
adjustment[isinf.(adjustment)] .= 1

cera_taxo = @transform cera_taxo {swvl1mean = (x -> if (x < 0.0m^3)
                                                   return 0.0m^3
                                               else
                                                   return x
                                               end)(:swvl1mean),
                                  swvl2mean = (x -> if (x < 0.0m^3)
                                                   return 0.0m^3
                                               else
                                                   return x
                                               end)(:swvl2mean),
                                  swvl3mean = (x -> if (x < 0.0m^3)
                                                   return 0.0m^3
                                               else
                                                   return x
                                               end)(:swvl3mean),
                                  swvl4mean = (x -> if (x < 0.0m^3)
                                                   return 0.0m^3
                                               else
                                                   return x
                                               end)(:swvl4mean)}

# Apply adjustment to data grouped by Species
phylo_traits = @groupby cera_taxo :SppID {species = first(:species),
                                          genus = first(:genus),
                                          family = first(:family),
                                          order = first(:order),
                                          class = first(:class),
                                          phylum = first(:phylum),
                                          numrecords = length(:UID),
                                          tmin = adjust(uconvert.(K, :tmin),
                                                        adjustment[:, 1],
                                                        mins[1], maxs[1]),
                                          tmax = adjust(uconvert.(K, :tmax),
                                                        adjustment[:, 2],
                                                        mins[2], maxs[2]),
                                          tmean = adjust(uconvert.(K, :tmean),
                                                         adjustment[:, 3],
                                                         mins[3], maxs[3]),
                                          trng = adjust(uconvert.(K, :trng),
                                                        adjustment[:, 4],
                                                        mins[4], maxs[4]),
                                          stl1 = adjust(uconvert.(K, :stl1mean),
                                                        adjustment[:, 5],
                                                        mins[5], maxs[5]),
                                          stl2 = adjust(uconvert.(K, :stl2mean),
                                                        adjustment[:, 6],
                                                        mins[6], maxs[6]),
                                          stl3 = adjust(uconvert.(K, :stl3mean),
                                                        adjustment[:, 7],
                                                        mins[7], maxs[7]),
                                          stl4 = adjust(uconvert.(K, :stl4mean),
                                                        adjustment[:, 8],
                                                        mins[8], maxs[8]),
                                          swvl1 = adjust(:swvl1mean,
                                                         adjustment[:, 9],
                                                         mins[9], maxs[9]),
                                          swvl2 = adjust(:swvl2mean,
                                                         adjustment[:, 10],
                                                         mins[10], maxs[10]),
                                          swvl3 = adjust(:swvl3mean,
                                                         adjustment[:, 11],
                                                         mins[11], maxs[11]),
                                          swvl4 = adjust(:swvl4mean,
                                                         adjustment[:, 12],
                                                         mins[12], maxs[12]),
                                          ssr = adjust(:ssrmean,
                                                       adjustment[:, 13],
                                                       mins[13], maxs[13]),
                                          tp = adjust(:tpmean,
                                                      adjustment[:, 14],
                                                      mins[14], maxs[14])}

tab = collect(phylo_traits)
tab = @transform tab {tmin = ustrip(:tmin),
                      tmax = ustrip(:tmax),
                      tmean = ustrip(:tmean),
                      trng = ustrip(:trng),
                      stl1 = ustrip(:stl1),
                      stl2 = ustrip(:stl2),
                      stl3 = ustrip(:stl3),
                      stl4 = ustrip(:stl4),
                      swvl1 = ustrip(:swvl1),
                      swvl2 = ustrip(:swvl2),
                      swvl3 = ustrip(:swvl3),
                      swvl4 = ustrip(:swvl4),
                      ssr = ustrip(:ssr),
                      tp = ustrip(:tp)}
CSV.write("Global_plant_climate_envelopes.csv", tab)
