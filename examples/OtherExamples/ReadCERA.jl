# SPDX-License-Identifier: BSD-2-Clause

using ClimatePref
using AxisArrays
using ClimatePref.Units
using Unitful

ClimatePref.readCERA("data", "cera_20C", "t2m")

test = ustrip.(cera1900.array[:, :, 1])
using RCall
@rput test
R"library(viridis);library(fields)
    image.plot(test, col = magma(20))"
