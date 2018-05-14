addprocs(6)
using Unitful
using AxisArrays
using ClimatePref
using myunitful
using Iterators
@everywhere using JuliaDB
using Unitful.DefaultSymbols


function create_reference(gridsize::Float64)
    x = 360 * (1/gridsize) + 1
    y = 180 * (1/gridsize) + 1
    gridsize = gridsize * °
    refarray = AxisArray(Array{Int64, 2}(Int(x), Int(y)),
                         Axis{:longitude}(-180.0°:gridsize:180.0°),
                         Axis{:latitude}(-90.0°:gridsize:90.0°))
    refarray[1:length(refarray)]= collect(1:length(refarray))
    ref = Reference(refarray)
end

function gardenmask(occ::IndexedTables.NextTable, gard::IndexedTables.NextTable,
     masksize::Float64)
    coords = hcat(select(gard, :Latitude), select(gard, :Longitude))
    ref = create_reference(masksize)

    # Add in refval column to garden info of which reference grid square it falls in
    gard = pushcol(gard, :refval, extractvalues(coords[:, 2] * °, coords[:, 1] * °, ref))

    coords_occ = hcat(select(occ, :decimallatitude), select(occ, :decimallongitude))
    occ = pushcol(occ, :refval, extractvalues(coords_occ[:, 2] * °, coords_occ[:, 1] * °, ref))
    # Use anti-join to filter out those that have the same reference value
    occ = join(occ, gard, how=:anti, lkey=:refval, rkey =:refval)
    occ = popcol(occ, :refval)
    return occ
end

function genus_clean(genus::IndexedTables.NextTable)
    gardens = load("data/garden_table")
    genus = gardenmask(genus, gardens, 0.02)
    return genus
end


@everywhere function genus_worldclim(genus::IndexedTables.NextTable)
    dir = "/Users/claireh/Documents/PhD/Data/Worldclim/wc2.0_5m"
    genus = genus_clean(genus)
    ref = create_reference(1/12)
    coords = hcat(select(genus, :decimallatitude), select(genus, :decimallongitude))
    genus = pushcol(genus, :refval, extractvalues(coords[:, 2] * °,
     coords[:, 1] * °, ref))
    # loop through worldclim variables and join based on reference grid square
    worldclim_names = map(searchdir(dir, "wc2.0_5m_")) do str
        split(str, "wc2.0_5m_")[2]
    end
    for i in 1:7
        worldclim_dat = load(string("data/Worldclim/", worldclim_names[i+1]))
        genus = join(genus, worldclim_dat,  how=:left, lkey=:refval, rkey =:refval, rselect=Symbol(worldclim_names[i+1]))
    end
    genus = popcol(genus, :refval)
    return genus
end

sol = load("/Users/claireh/Documents/PhD/Data/GBIF/Solanum/output")
types = map(x-> fieldtype(eltype(sol), x), fieldnames(eltype(sol)))
names = fieldnames(eltype(sol))
dir = "/Users/claireh/Documents/PhD/Data/GBIF/final/"
genera = searchdir(dir, ".csv")

@sync @parallel for i in genera[1:end]
    genus = loadtable(string(dir, i),
        indexcols = [:phylum, :class, :order, :family, :genus, :species, :scientificname],
        colparsers=Dict(zip(names,types)), distributed=false)
    genus = genus_worldclim(genus)
    save(genus, string("data/Genera/Worldclim/", split(i, ".csv")[1]))
end





# Load in cleaned Solanum data
#sol = load("data/Clean_data")
# Select coordinates from solanum data
#coords_sol = hcat(select(sol, :decimallatitude), select(sol, :decimallongitude))
#sol = pushcol(sol, :refval, extractvalues(coords_sol[:, 2] * °, coords_sol[:, 1] * °, ref))
# loop through worldclim variables and join based on reference grid square
#for i in 1:7
#    worldclim_dat = load(string("data/Worldclim/", worldclim_names[i+1]))
#    sol = join(sol, worldclim_dat,  how=:left, lkey=:refval, rkey =:refval, rselect=Symbol(worldclim_names[i+1]))
#end
# Save dataset as solanum worldclim data
#save(sol, "data/Worldclim/sol_worldclim")
