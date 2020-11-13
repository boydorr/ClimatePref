using JuliaDB


# Dictionary to match each gbif heading to the correct type in Julia - the auto-match does not seem to work very well in JuliaDB with complex strings.
coldict = Dict(:phylum => String, :class => String, :order => String, :family => String, :genus => String, :species => String, :scientificname => String, :gbifid => String, :datasetkey => String, :occurrenceid => String, :kingdom => String, :infraspecificepithet => String, :taxonrank => String, :countrycode => String, :locality => String, :publishingorgkey => String, :decimallatitude => Union{Missing, Float64}, :decimallongitude => Union{Missing, Float64}, :coordinateuncertaintyinmeters => Union{Missing, Float64}, :coordinateprecision => Union{Missing, Float64}, :elevation => Union{Missing, Float64}, :elevationaccuracy => Union{Missing, Float64}, :depth => Union{Missing, Float64}, :depthaccuracy => Union{Missing, Float64}, :eventdate => String, :day => Union{Missing, Int64}, :month => Union{Missing, Int64}, :year => Union{Missing, Int64}, :taxonkey => String, :specieskey => String, :basisofrecord => String, :institutioncode => String, :collectioncode => String, :catalognumber => String, :recordnumber => String, :identifiedby => String, :license => String, :rightsholder => String, :recordedby => String, :typestatus => String, :establishmentmeans => String, :lastinterpreted => String, :mediatype => String, :issue => String)

"""
    ReadGBIF(folder::String)

Function to read in data from GBIF in the form of a folder housing CSVs of each Genera.
"""
function ReadGBIF(folder::String)
    gbif = loadtable(folder,
    indexcols = [:phylum, :class, :order, :family, :genus, :species, :scientificname], colparsers = coldict)
    return gbif
end
