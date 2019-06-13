using Unitful
using AxisArrays
using ClimatePref
using myunitful
using Iterators

import Unitful: °, °C, mm
import ArchGDAL
import Base.read
import IndexedTables.rows
import IndexedTables.Columns
rows(x::AbstractArray) = x
function Columns(cols::AbstractArray...; names::Union{Vector,Tuple{Vararg{Any}},Void}=nothing)
    if isa(names, Void) || any(x->!(x isa Symbol), names)
        Columns{IndexedTables.eltypes(typeof(cols)),typeof(cols)}(cols)
    else
        dt = eval(:(@NT($(names...)))){map(eltype, cols)...}
        ct = eval(:(@NT($(names...)))){map(typeof, cols)...}
        Columns{dt,ct}(ct(cols...))
    end
end
const AG = ArchGDAL

# Load in worldclim and bioclim datasets
dir = "/Users/claireh/Documents/PhD/Data/Worldclim/wc2.0_5m"
tavg = readworldclim(joinpath(dir, "wc2.0_5m_tavg"))
tmax = readworldclim(joinpath(dir, "wc2.0_5m_tmax"))
tmin = readworldclim(joinpath(dir, "wc2.0_5m_tmin"))
prec = readworldclim(joinpath(dir, "wc2.0_5m_prec"))
srad = readworldclim(joinpath(dir, "wc2.0_5m_srad"))
vapr = readworldclim(joinpath(dir, "wc2.0_5m_vapr"))
wind = readworldclim(joinpath(dir, "wc2.0_5m_wind"))
bio = readbioclim(joinpath(dir, "wc2.0_5m_bio"))
using JuliaDB

# Create a reference array the same size as the climate data
gridsize = step(axes(tavg.array, 1).val)
refarray = AxisArray(Array{Int64, 2}(4321, 2161),
                     Axis{:longitude}(-180.0°:gridsize:180.0°),
                     Axis{:latitude}(-90.0°:gridsize:90.0°))
refarray[1:length(refarray)]= collect(1:length(refarray))
ref = ClimatePref.Reference(refarray)

# Collect together all of the worldclim data
worldclim_collect = [prec, srad, tavg, tmax, tmin, vapr, wind]
# Average over the 12 months of data
worldclim_mean = map(x -> mapslices(mean, x.array, 3)[:,:,1], worldclim_collect)
# Read in the names (note that this starts with bioclim)
worldclim_names = map(searchdir(dir, "wc2.0_5m_")) do str
    split(str, "wc2.0_5m_")[2]
end

# Generate all x and y coordinate pairs for the entire grid
x = collect(axes(worldclim_collect[1].array, 1).val)
y = collect(axes(worldclim_collect[1].array, 2).val)
expandedxy = collect(product(x, y))
newx = map(x-> x[1], expandedxy)
newy = map(x-> x[2], expandedxy)

# Loop through the averaged datasets, creating a juliaDB table of all
# possible grid locations and the climate at these points.
# Save them separately to be loaded in again for faster use.
for i in 2:length(worldclim_mean)
    values = worldclim_mean[i][1:end]
    worldclim_tab = table(newx, newy, values,
        names = [:x, :y, Symbol(worldclim_names[i+1])])
    if (i == 1)
        coords = hcat(select(worldclim_tab, :x), select(worldclim_tab, :y))
        vals = extractvalues(coords[:, 1], coords[:, 2], ref)
    end
    worldclim_tab = pushcol(worldclim_tab, :refval, vals)
    save(worldclim_tab, string("data/Worldclim/", worldclim_names[i+1]))
end
# Can store values for each month as an array (with some tweaks to the JuliaDB code),
# but this is very slow - get in contact with JuliaDB people
len = length(worldclim_collect[1].array[:,:,1])
for i in 1:length(worldclim_collect)
    reshaped_array = reshape(worldclim_collect[i].array, len, 12)
    values = mapslices(x-> [x], reshaped_array, 2)
    worldclim_tab = table(newx, newy, values,
        names = [:x, :y, Symbol(worldclim_names[i+1])])
    if (i == 1)
        coords = hcat(select(worldclim_tab, :x), select(worldclim_tab, :y))
        vals = extractvalues(coords[:, 1], coords[:, 2], ref)
    end
    worldclim_tab = pushcol(worldclim_tab, :refval, vals)
    save(worldclim_tab, string("data/Worldclim/", worldclim_names[i+1]))
end

# Load in cleaned Solanum data
sol = load("data/Clean_data")
# Select coordinates from solanum data
coords_sol = hcat(select(sol, :decimallatitude), select(sol, :decimallongitude))
sol = pushcol(sol, :refval, extractvalues(coords_sol[:, 2] * °, coords_sol[:, 1] * °, ref))
# loop through worldclim variables and join based on reference grid square
for i in 1:7
    worldclim_dat = load(string("data/Worldclim/", worldclim_names[i+1]))
    sol = join(sol, worldclim_dat,  how=:left, lkey=:refval, rkey =:refval, rselect=Symbol(worldclim_names[i+1]))
end
# Save dataset as solanum worldclim data
save(sol, "data/Worldclim/sol_worldclim")

# Do the same for 19 bioclim variables!
# Create name and unit vectors
bioclim_names = ["annual_tavg", "mean_diurnal_range", "isothermality",
"temp_seasonality", "tmax_warmestm", "tmin_coldestm", "annual_trange",
"tavg_wettestq", "tavg_driestq", "tavg_warmestq", "tavg_coldestq", "annual_prec",
"prec_wettestm", "prec_driestm", "prec_seasonality", "prec_wettestq", "prec_driestq",
"prec_warmestq", "prec_coldestq"]
units = [°C, °C, 1, °C, °C, °C, °C, °C, °C, °C, °C,
u"mm", u"mm",u"mm",u"mm",u"mm",u"mm",u"mm",u"mm",]

# Create x, y, grid locations
x = collect(axes(bio.array[:,:,1], 1).val)
y = collect(axes(bio.array[:,:,1], 2).val)
expandedxy = collect(product(x, y))
newx = map(x-> x[1], expandedxy)
newy = map(x-> x[2], expandedxy)

# Loop through each variable and append to juliaDB table
for i in 1:length(bioclim_names)
    if (i == 1)
        bioclim_tab = table(newx, newy, bio.array[:,:,i][1:end] * units[i], names = [:x, :y, Symbol(bioclim_names[i])])
    else
        bioclim_tab = pushcol(bioclim_tab, Symbol(bioclim_names[i]), bio.array[:,:,i][1:end] * units[i])
    end
end
# Extract the reference values at the different grid locations
coords = hcat(select(bioclim_tab, :x), select(bioclim_tab, :y))
bioclim_tab = pushcol(bioclim_tab, :refval,
extractvalues(coords[:, 1], coords[:, 2], ref))
# Save as bioclim file
save(bioclim_tab, string("data/Bioclim/bio"))

sol = load("data/Clean_data")
# Select coordinates from solanum data
coords_sol = hcat(select(sol, :decimallatitude), select(sol, :decimallongitude))
sol = pushcol(sol, :refval, extractvalues(coords_sol[:, 2] * °, coords_sol[:, 1] * °, ref))
bioclim_dat = load("data/Bioclim/bio")
sol = join(sol, bioclim_dat,  how=:left, lkey=:refval, rkey =:refval,
    rselect=Tuple(map(x -> Symbol(x), bioclim_names)))
save(sol, "data/Bioclim/sol_bioclim")
