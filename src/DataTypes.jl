using AxisArrays
using Unitful
using MyUnitful
using RecipesBase
import AxisArrays.axes
abstract type AbstractClimate end

mutable struct Worldclim <: AbstractClimate
    array::AxisArray
    function Worldclim(array::AxisArray)
        size(array, 3) == 12 ||
            error("There should be 12 months of data for worldclim")
        new(array)
    end
end

mutable struct Bioclim <: AbstractClimate
    array::AxisArray
    function Bioclim(array::AxisArray)
        size(array, 3) == 19 ||
            error("There should 19 climate variables for bioclim")
        new(array)
    end
end

mutable struct ERA <: AbstractClimate
    array::AxisArray
    function ERA(array::AxisArray)
        typeof(collect(axes(array, 3).val)[1])<: Unitful.Time ||
            error("Third dimension of array must be time")
        new(array)
    end
end

@recipe function f(era::ERA, time::Unitful.Time)
    tm = ustrip.(uconvert(year, time))
    yr = floor(Int64, tm)
    ind = round(Int64, (tm - yr)/(1/12))
    typeof(ind) <: Int64 || error("NO")
    mnth = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"][ind+1]
    A = transpose(ustrip.(era.array[:, :, time]))
    x = -179.25:0.75:179.25
    y = -89.25:0.75:90
    seriestype  :=  :heatmap
    grid --> false
    title --> "$yr $mnth"
    x, y, A
end

mutable struct CERA <: AbstractClimate
    array::AxisArray
    function CERA(array::AxisArray)
        typeof(collect(axes(array, 3).val)[1])<: Unitful.Time ||
            error("Third dimension of array must be time")
        new(array)
    end
end

mutable struct Reference <: AbstractClimate
    array::AxisArray
end



function TestERA()
    dir = dirname(pathof(ClimatePref)) * "/../test/Testdata/TestERA"
    data = extractERA(dir, "t2m", collect(1.0month:1month:10year))
    data.array = data.array[-10° .. 60°, 35° .. 80°, :]
    return data
end

function TestWorldclim()
    dir = dirname(pathof(ClimatePref)) * "/../test/Testdata/TestWorldclim/"
    data = extractworldclim(joinpath(dir, "wc2.0_5m_srad"))
    return data
end
