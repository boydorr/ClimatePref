# SPDX-License-Identifier: BSD-2-Clause

using ClimatePref
using ClimatePref.Units
using AxisArrays
using Test

wc = AxisArray(fill(1.0, 100, 100, 12))
@test_nowarn Worldclim(wc)

era = AxisArray(fill(1.0, 100, 100, 12), Axis{:lat}(1:100), Axis{:lon}(1:100),
                Axis{:time}((1month):(1month):(12months)))
@test_nowarn ERA(era)
@test_nowarn CERA(era)
