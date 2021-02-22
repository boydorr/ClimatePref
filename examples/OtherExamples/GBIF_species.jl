addprocs(8)
@everywhere using JuliaDB
gbif = load("/Users/claireh/Documents/PhD/Data/GBIF/test/output")
gbif = pushcol(gbif, :id, collect(1:25331))

tb = groupby(length, gbif, :species)
groupreduce(+, gbif, :species, select=:id)

save(tb, "gbif_species")


using JuliaDB
test = loadtable("/Users/claireh/Documents/PhD/Data/GBIF/test",
       indexcols = [:phylum, :class, :order, :family, :genus, :species, :scientificname],
       type_detect_rows =5000,
       colparsers=Dict(:datasetkey=>String,
                       :occurrenceid=>String,
                       :gbifid=>String,
                       :locality=>String,
                       :publishingorgkey=>String,
                       :taxonkey=>String,
                       :institutioncode=>String,
                       :catalognumber=>String,
                       :recordnumber=>String))
save(test, "/Users/claireh/Documents/PhD/Data/GBIF/test/output")

using JuliaDB
test = load("/Users/claireh/Documents/PhD/Data/GBIF/test/output")
test = pushcol(test, :id, collect(1:25331))
res = filter((x->x.phylum ==""), test)
collect(select(res, :id))
result = groupby(length, test, :phylum)
test2 = filter((x->x.species !=""), test)
save(test2, "/Users/claireh/Documents/PhD/Data/GBIF/test/output2")
test2 = load("/Users/claireh/Documents/PhD/Data/GBIF/test/output2")
result = groupby(length, test2, :species)

addprocs(8)
using JuliaDB
test = loadtable("/Users/claireh/Documents/PhD/Data/GBIF/test_copy",
       indexcols = [:phylum, :class, :order, :family, :genus, :species, :scientificname],
       type_detect_rows =5000,
       colparsers=Dict(:datasetkey=>String,
                       :occurrenceid=>String,
                       :gbifid=>String,
                       :locality=>String,
                       :publishingorgkey=>String,
                       :taxonkey=>String,
                       :institutioncode=>String,
                       :catalognumber=>String,
                       :recordnumber=>String))
save(test,  "/Users/claireh/Documents/PhD/Data/GBIF/test_copy/output")

test = load("/Users/claireh/Documents/PhD/Data/GBIF/test_copy/output")
test = pushcol(test, :id, collect(1:25331))
result = groupby(length, test, :genus)


addprocs(8)
using JuliaDB
gbif = load("GBIF_worldclim")
tavg = load("tavg")
dtavg = distribute(tavg, 8)
gbif_tavg = join(gbif, dtavg, how=:left, lkey=:refval, rkey=:refval, rselect=:tavg)

if typeof(select(red_genus, 2))<: Array{NamedTuples.NamedTuple, 1}
    red_genus = table(select(red_genus, 1),
    map(x-> collect(select(red_genus,2)[x])[1], 1:length(red_genus)),
    names = [:species, :na_mean])
end
function na_mean(x, y)
    if typeof(x) <: DataValues.DataValue x = x.value end
    if typeof(y) <: DataValues.DataValue y = y.value end
    if isnan(x) & isnan(y)
        z = NaN32 *°C
    else
        vec = [x, y]
        z = return mean(vec[.!isnan.(vec)])
    end
    return collect(z)[1]
end
