"""
    create_reference(gridsize::Float64)

Function to create a reference grid array of type `Reference`.
"""
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
"""
    gardenmask(occ::IndexedTables.NextTable, gard::IndexedTables.NextTable,
     masksize::Float64)

Function to mask botanic garden locations found in `gard` from occurrence
records found in `occ`, with an optional size for the mask, `masksize`.
"""
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
"""
    genus_clean(genus::IndexedTables.NextTable)

Function to clean occurrence data of botanic garden information.
"""
function genus_clean(genus::IndexedTables.NextTable)
    gardens = load("data/garden_table")
    genus = gardenmask(genus, gardens, 0.02)
    return genus
end

"""
    genus_worldclim(genus::IndexedTables.NextTable)

Function to clean occurrence data of botanic garden information, and
join with worldclim data.
"""
function genus_worldclim_average(genus::IndexedTables.NextTable,
    worldclim::IndexedTables.NextTable)
    worldclim_names = [:prec, :srad, :tavg, :tmax, :tmin, :vapr, :wind]
    genus = genus_clean(genus)
    ref = create_reference(1/12)
    coords = hcat(select(genus, :decimallatitude), select(genus, :decimallongitude))
    genus = pushcol(genus, :refval, extractvalues(coords[:, 2] * °,
     coords[:, 1] * °, ref))
    genus = join(genus, worldclim_dat,  how=:left, lkey=:refval, rkey =:refval,
        rselect=Tuple(Symbol.(worldclim_names)))
    genus = popcol(genus, :refval)
    return genus
end
function genus_worldclim_monthly(genus::IndexedTables.NextTable,
    worldclim::IndexedTables.NextTable)
    worldclim_names = [:prec, :srad, :tavg, :tmax, :tmin, :vapr, :wind]
    genus = genus_clean(genus)
    ref = create_reference(1/12)
    coords = hcat(select(genus, :decimallatitude), select(genus, :decimallongitude))
    genus = pushcol(genus, :refval, extractvalues(coords[:, 2] * °,
     coords[:, 1] * °, ref))
    genus = join(genus, worldclim_dat,  how=:left, lkey=:refval, rkey =:refval,
        rselect=Tuple(Symbol.(worldclim_names)))
    genus = popcol(genus, :refval)
    genus = filter(p-> !isnull(p.tavg), genus)
    if length(genus)== 0
     len = 0
    else
     len = Int(length(genus)/12)
    end
    months = repmat(collect(1:12), len)
    genus = pushcol(genus, :worldclim_month, months)
    return genus
end
