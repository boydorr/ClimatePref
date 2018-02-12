#using BinDeps

#@BinDeps.setup

#function version_check(name, handle)
#    fptr = Libdl.dlsym_e(handle, :GDALVersionInfo)
#    if fptr == C_NULL   #lookup failure
#        return false
#    end
#    versionptr = ccall(fptr,Cstring,(Cstring,),"RELEASE_NAME")
#    versionstring = unsafe_string(versionptr)
#    gdalversion = convert(VersionNumber, versionstring)
#    gdalversion >= v"2.1.0"
#end
#Libdl.dlopen("/home/claireh/.julia/v0.6/Conda/deps/usr/lib/libgdal.so", Libdl.RTLD_GLOBAL)
#libgdal = library_dependency("libgdal",
#                             aliases=["gdal","gdal201", "gdal202",
#                             "gdal_w32","gdal_w64","libgdal.20",
#                             "libgdal.a", "libgdal.la",
#                             "libgdal.so", "libgdal.so.20",
#                             "libgdal.so.20.3.2"],
#                             validate=version_check)
#using Conda
#Conda.add_channel("mychannel")
#provides(Conda.Manager, "gdal", libgdal, installed_libpath="/home/claireh/.julia/v0.6/Conda/deps/usr/lib/")
#@BinDeps.install Dict(:libgdal => :libgdal)
 #Only needs to be done once

#Pkg.clone("https://github.com/visr/GDAL.jl.git")
#Pkg.build("GDAL")
#Pkg.clone("https://github.com/yeesian/ArchGDAL.jl.git")
#Pkg.build("ArchGDAL")

using Unitful
using Unitful.DefaultSymbols
using myunitful
using AxisArrays

import Unitful: °, °C, mm
import ArchGDAL
import Base.read
const AG = ArchGDAL

vardict = Dict("bio" => NaN, "prec" => mm, "srad" => u"kJ"* u"m"^-2 * day^-1,
"tavg" => °C, "tmax" => °C, "tmin" => °C, "vapr" => u"kPa", "wind" => u"m" * u"s"^-1)
"""
    read(f, filename)

Function to read raster file into julia.
"""
function read(f, filename)
    return AG.registerdrivers() do
        AG.read(filename) do dataset
            f(dataset)
        end
    end
end
"""
    searchdir(path,key)

Function to search a directory `path` using a given `key` string.
"""
searchdir(path,key) = filter(x->contains(x,key), readdir(path))
"""
    extractfile(dir::String)

Function to import a selected file from a path string.
"""
function extractfile(dir::String)
    txy = [Float64, Int64(1), Int64(1)]

    read(file) do dataset
        #txy[1] = AG.getdatatype(AG.getband(dataset, 1))
        txy[2] = AG.width(AG.getband(dataset, 1))
        txy[3] = AG.height(AG.getband(dataset, 1))
        print(dataset)
    end

    a = Array{txy[1], 2}(txy[2], txy[3])
    read(file) do dataset
        bd = AG.getband(dataset, 1);
        AG.read!(bd, a);
    end;
    lat, long = size(a, 1), size(a, 2);
    step1 = 360.0° / lat;
    step2 = 180.0° / long;

    world = AxisArray(a[:, long:-1:1],
                           Axis{:latitude}(-180.0°:step1:(180.0° - step1 / 2)),
                           Axis{:longitude}(-90.0°:step2:(90.0°-step2/2)));

    if txy[1] <: AbstractFloat
        world[world .== world[Axis{:latitude}(0°),
                                             Axis{:longitude}(step2/2)]] *= NaN;
    end;
    world
end
"""
    extractfolder(dir::String)

Function to extract all raster files from a specified folder directory,
and convert into an axis array.
"""
function extractfolder(dir::String)
    files = map(searchdir(dir, ".tif")) do files
        joinpath(dir, files)
    end
    txy = [Float32, Int32(1), Int32(1)];

    read(files[1]) do dataset
        txy[1] = AG.getdatatype(AG.getband(dataset, 1))
        txy[2] = AG.width(AG.getband(dataset, 1))
        txy[3] = AG.height(AG.getband(dataset, 1))
        print(dataset)
    end

    numfiles = length(files)
    b = Array{txy[1], 3}(Int64(txy[2]), Int64(txy[3]), numfiles);
    map(eachindex(files)) do count
    a = Array{txy[1], 2}(txy[2], txy[3]);
    read(files[count]) do dataset
        bd = AG.getband(dataset, 1);
        AG.read!(bd, a);
    end;
    b[:, :, count] = a
    end
    lat, long = size(b, 1), size(b, 2);
    variable = split(dir, "wc2.0_5m_")[2]
    unit = vardict[variable]
    step = 180.0° / long;
    if numfiles == 12
        thirdaxis = Axis{:time}(1month:1month:12month)
    else
        thirdaxis = Axis{:var}(1:1:numfiles)
        unit = 1.0
    end
    world = AxisArray(b[:, long:-1:1, :] * unit,
                           Axis{:latitude}(-180.0°:step:(180.0° - step / 2)),
                           Axis{:longitude}(-90.0°:step:(90.0°-step/2)),
                           thirdaxis);

    if txy[1] <: AbstractFloat
        world[world .== world[Axis{:latitude}(0°),
                                             Axis{:longitude}(0°),
                                             thirdaxis(thirdaxis[1])]] *= NaN;
    end;
    world
end
"""
    extractfolders(dir::String, folders::Array{String, 1})

Function to extract all raster files from several specified folder directories,
`folders`, and convert into an axis array.
"""
function extractfolders(dir::String, folders::Array{String, 1})
    extract = map(folders) do folder
        extractfolder(joinpath(dir, folder))
    end
end


"""
    extractvalues(x::Vector{typeof(1.0°)},y::Vector{typeof(1.0°)},
     array::AxisArray, dim::Unitful.Time)

Function to extract values from a climate array, at specified x, y locations at
a specific time, `dim`.
"""

function extractvalues(x::Vector{typeof(1.0°)},y::Vector{typeof(1.0°)},
   array::AxisArray, dim::Unitful.Time)
   all(x .<= 180.0) && all(x .>= -180.0) ||
   error("X coordinate is out of bounds")
   all(y .< 90.0) && all(y .> -90.0) ||
   error("Y coordinate is out of bounds")
   return map((i, j) -> array[(i-step/2)..(i+step/2),
                              (j-step/2)..(j+step/2), dim][1], x, y)
end
"""
    extractvalues(x::Vector{typeof(1.0°)},y::Vector{typeof(1.0°)},
       array::AxisArray, dim::StepRange{typeof(1month)})

Function to extract values from a climate array, at specified x, y locations and
time period, `dim`.
"""

function extractvalues(x::Vector{typeof(1.0°)},y::Vector{typeof(1.0°)},
   array::AxisArray, dim::StepRange{typeof(1month)})
   all(x .<= 180.0) && all(x .>= -180.0) ||
   error("X coordinate is out of bounds")
   all(y .< 90.0) && all(y .> -90.0) ||
   error("Y coordinate is out of bounds")
   res = map((i, j) -> array[(i-step/2)..(i+step/2),
                              (j-step/2)..(j+step/2),
                              start(dim)..last(dim)][1,1,:], x, y)
   return transpose(hcat(res...))
end
"""
    extractvalues(x::Vector{typeof(1.0°)},y::Vector{typeof(1.0°)},
       array::typeof(files))

Function to extract values from several climate arrays, at specified x, y
locations, over all axes.
"""
function extractvalues(x::Vector{typeof(1.0°)},y::Vector{typeof(1.0°)},
   array::AbstractArray)
   thisstep = step(axes(array[1], 1).val)
   all(x .<= 180.0) && all(x .>= -180.0) ||
   error("X coordinate is out of bounds")
   all(y .< 90.0) && all(y .> -90.0) ||
   error("Y coordinate is out of bounds")
   map(array) do array
    dim = axes(array, 3).val
    res = map((i, j) -> array[(i-thisstep/2)..(i+thisstep/2),
                              (j-thisstep/2)..(j+thisstep/2),
                              start(dim)..last(dim)][1,1,:], x, y)
    transpose(hcat(res...))
    end
end
