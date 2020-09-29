using JuliaDB
using ClimatePref
using Unitful
using ClimatePref.Units

# Load GBIF and add in unique ID column
gbif = load("GBIF_TPL")
gbif = pushcol(gbif, :UID, 1:length(gbif))

# Create date column
mth = collect(select(gbif, :month))
yr = collect(select(gbif, :year))
dt = (mth .* month) .+ (yr .* year)
dt = uconvert.(year, dt)
gbif = pushcol(gbif, :date, dt)

# Extract coordinates and compare to reference grid
lat = collect(select(gbif, :decimallatitude))
lon = collect(select(gbif, :decimallongitude))
ref = create_reference(0.75)
refval = extractvalues(lon .* °, lat .* °, ref)
gbif = pushcol(gbif, :refval, refval)
small_gbif = select(gbif, (:UID, :SppID, :date, :refval))
save(small_gbif, "Small_GBIF")

# Import era data and refine
era = load("ECMWF/era_int_all")
yr = select(era, :year) .* year
mth = select(era, :month) .* month
function era_extend(era::IndexedTable)
    for i in 0:11
        dt = yr .+ mth
        dt .+= (i * month)
        dt = uconvert.(year, dt)
        era = pushcol(era, Symbol("date$i"), dt)
    end
    return era
end
era_date = era_extend(era)
era_date = renamecol(era_date, :refid => :refval)
save(era_date, "ECMWF/era_int")

# Join ERA with gbif (once for every twelve months)
gbif = load("GBIF/Small_GBIF")
era = load("ECMWF/era_int")
era = distribute(era, 1)

function era_gbif_join(era::JuliaDB.DNextTable, gbif::JuliaDB.DNextTable, date_var::Symbol)
    bif = renamecol(gbif, colnames(gbif)[3] => date_var)
    gbif_era = join(bif, era, how = :inner, lkey = (:refval, date_var), rkey = (:refval, date_var), rselect = (:year, :refval, :x, :y, :month, :ssr, :stl1, :t2m, :tp, date_var))
    gbif_era = renamecol(gbif_era, date_var, :date)
    return gbif_era
end

function era_gbif_join(era::JuliaDB.DNextTable, gbif::JuliaDB.DNextTable)
    gbif_era = era_gbif_join(era, gbif, :date0)
    for i in 1:11
        new_join = era_gbif_join(era, gbif, Symbol("date$i"))
        gbif_era = merge(gbif_era, new_join)
    end
    return gbif_era
end

joined_data = era_gbif_join(era, gbif)
