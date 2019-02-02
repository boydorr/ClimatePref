using IndexedTables
using AxisArrays
using Unitful.DefaultSymbols
using Plots
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
function gardenmask(occ::IndexedTable, gard::IndexedTable,
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
function genus_clean(genus::IndexedTable)
    gardens = load("data/garden_table")
    genus = gardenmask(genus, gardens, 0.02)
    return genus
end

"""
    genus_worldclim_average(genus::IndexedTables.NextTable)

Function to clean occurrence data of botanic garden information, and
join with worldclim data.
"""
function genus_worldclim_average(genus::IndexedTable,
    worldclim::IndexedTable)
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
"""
    genus_worldclim_monthly(genus::IndexedTables.NextTable)

Function to clean occurrence data of botanic garden information, and
join with monthly worldclim data.
"""
function genus_worldclim_monthly(genus::IndexedTable,
    worldclim::IndexedTable)
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

"""
    upresolution(genus::IndexedTables.NextTable)

Function to clean occurrence data of botanic garden information.
"""
function upresolution(era::ERA, rescale::Int64)
    array = upresolution(era.array, rescale)
    return ERA(array)
end
function upresolution(wc::Worldclim, rescale::Int64)
    array = upresolution(wc.array, rescale)
    return Worldclim(array)
end
function upresolution(bc::Bioclim, rescale::Int64)
    array = upresolution(bc.array, rescale)
    return Bioclim(array)
end
function upresolution(aa::AxisArray{T, 3} where T, rescale::Int64)
    grid = size(aa)
    grid = (grid[1] .* rescale, grid[2] .* rescale, grid[3])
    array = Array{typeof(aa[1]), 3}(grid)
    map(1:grid[3]) do time
        for x in 1:size(aa, 1)
            for y in 1:size(aa, 2)
        array[(rescale*x-(rescale-1)):(rescale*x),
            (rescale*y-(rescale - 1)):(rescale*y), time] = aa[x, y, time]
            end
        end
    end
    lon = aa.axes[1].val
    smallstep = (lon[2] - lon[1]) / rescale
    if lon[1] == -180°
        newlon = collect(lon[1]:smallstep:(lon[end]+smallstep))
    else
        newlon = collect((lon[1] -smallstep):smallstep:lon[end])
    end
    lat = aa.axes[2].val
    smallstep = (lat[2] - lat[1]) / rescale
    if lat[1] == -90°
        newlat = collect(lat[1]:smallstep:(lat[end]+smallstep))
    else
        newlat = collect((lat[1]-smallstep):smallstep:lat[end])
    end
    return AxisArray(array,
        Axis{:longitude}(newlon),
        Axis{:latitude}(newlat),
        Axis{:time}(aa.axes[3].val))
end
function upresolution(aa::AxisArray{T, 2} where T, rescale::Int64)
    grid = size(aa) .* rescale
    array = Array{typeof(aa[1]), 2}(grid)
    for x in 1:size(aa, 1)
        for y in 1:size(aa, 2)
            array[(rescale*x-(rescale-1)):(rescale*x),
            (rescale*y-(rescale - 1)):(rescale*y)] = aa[x, y]
        end
    end
    lon = aa.axes[1].val
    smallstep = (lon[2] - lon[1]) / rescale
    if lon[1] == -180°
        newlon = collect(lon[1]:smallstep:(lon[end]+smallstep))
    else
        newlon = collect((lon[1] -smallstep):smallstep:lon[end])
    end
    lat = aa.axes[2].val
    smallstep = (lat[2] - lat[1]) / rescale
    if lat[1] == -90°
        newlat = collect(lat[1]:smallstep:(lat[end]+smallstep))
    else
        newlat = collect((lat[1]-smallstep):smallstep:lat[end])
    end
    return AxisArray(array,
        Axis{:longitude}(newlon),
        Axis{:latitude}(newlat))
end
import Plots: px
function getprofile(spp_names::Array{String, 1}, data::IndexedTable, variable_name::String, dims::Tuple{Int64, Int64} = (1,1))
    for i in eachindex(spp_names)
        (dims[1] * dims[2] == length(spp_names) || dims == (1,1)) || error("Dimensions not big enough for number of species")
        spp = filter(p-> p[:species] == spp_names[i], data)
        vals = select(spp, :val)
        res = vcat(vals...)
        res = res[.!isnan.(res)]
        sp = ifelse(dims == (1,1), 1, i)
        lg = ifelse(dims == (1,1), :left, :none)
        title = ifelse(dims == (1,1), "", spp_names[i])
        if i == 1
            display(histogram(res, bins = -20:2:30, grid = false, xlabel = variable_name, label=spp_names[i], layout = dims, legend= lg, top_margin = 20px,
            bottom_margin = 20px, title = title))
        else
            display(histogram!(res, bins = -20:2:30, grid = false, xlabel = variable_name, label=spp_names[i], subplot = sp, legend= lg, top_margin = 20px,
            bottom_margin = 20px, title = title))
        end
    end
end
