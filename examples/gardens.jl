using JuliaDB
using DataFrames
using CSV
using Unitful
using Unitful.DefaultSymbols
using AxisArrays

gardens = loadtable("/Users/claireh/Documents/PhD/Data/GBIF/garden",
    colparsers=Dict(:gardenid=>Int64,
                    :institution=>String,
                    :city=>String,
                    :state=>String,
                    :country_name=>String,
                    :InstitutionType=>String))
coords = hcat(select(gardens, :Latitude), select(gardens, :Longitude))

# Create reference array
gridsize = 0.02°
refarray = AxisArray(Array{Int64, 2}(18001, 9001), Axis{:longitude}(-180.0°:gridsize:180.0°),
                     Axis{:latitude}(-90.0°:gridsize:90.0°))
refarray[1:length(refarray)]= collect(1:length(refarray))
ref = Reference(refarray)

# Add in refval column to garden info of which reference grid square it falls in
gardens = pushcol(gardens, :refval, extractvalues(coords[:, 2] * °, coords[:, 1] * °, ref))

# Do the same for solanum species
sol = load("/Users/claireh/Documents/PhD/Data/GBIF/Solanum/output")
sol = pushcol(sol, :id, collect(1:215440))
coords_sol = hcat(select(sol, :decimallatitude), select(sol, :decimallongitude))
sol = pushcol(sol, :refval, extractvalues(coords_sol[:, 2] * °, coords_sol[:, 1] * °, ref))
# Use anti-join to filter out those that have the same reference value
sol = join(sol, gardens, how=:anti, lkey=:refval, rkey =:refval)

function gardenmask!(occ::JuliaDB.Table, gard::JuliaDB.Table, masksize::Float64)
    coords = hcat(select(gard, :Latitude), select(gard, :Longitude))
    coords2 = hcat(select(sol, :decimallatitude), select(sol, :decimallongitude))
    removals = Array{Array{Int64, 1}, 1}(length(gard))
    for i in 1:length(gard)
        res = filter(p-> (p.decimallatitude > (coords[i, 1] - masksize)) &
         (p.decimallatitude <(coords[i, 1] + masksize)) &
         (p.decimallongitude > (coords[i,2] - masksize)) &
         (p.decimallongitude < (coords[i,2] + masksize)),
        occ)
        removals[i] = select(res, :id)
    end
    rem = vcat(removals[!isempty.(removals)]...)

end
function gardenmask(occ::IndexedTables.NextTable, gard::IndexedTables.NextTable,
     masksize::Float64)
    coords = hcat(select(gard, :Latitude), select(gard, :Longitude))
    coords2 = hcat(select(occ, :decimallatitude), select(occ, :decimallongitude))
    removals = Array{Array{Int64, 1}, 1}(length(gard))
    for i in 1:length(gard)
        removals[i] = find((coords2[:,1] .> (coords[i, 1] - masksize)) .&
         (coords2[:,1] .< (coords[i, 1] + masksize)) .&
         (coords2[:,2] .> (coords[i,2] - masksize)) .&
         (coords2[:,2] .< (coords[i,2] + masksize)))
    end
    rem = vcat(removals[.!isempty.(removals)]...)
    tab = table(rem,collect(1:length(rem)), names = [:id, :test])
    occ = join(occ, tab, how=:anti, lkey=:id, rkey =:id)
    numrows = length(rem)
    info("$numrows rows removed")
    return occ
end
rows(gardens)[find(!isempty.(removals))]
sol = gardenmask(sol, gardens, 0.01)
# Check locality info against gardens
info_check = DataFrame(gard = Vector{String}(3441),
                locality = Vector{Vector{Vector{String}}}(3441))
for i in 1:length(removals)
    info_check[:gard][i] = rows(gardens)[i][:institution]
    locs = map(x -> select(filter(p->p.id == x, sol), :locality), removals[i])
    info_check[:locality][i] = locs
end
info_check = info_check[.!isempty.(info_check[:,2]), :]
CSV.write("data/garden_removal.csv", info_check)
