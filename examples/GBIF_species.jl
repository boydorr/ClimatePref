addprocs(8)
@everywhere using JuliaDB
gbif = load("output2")

tb = groupby(length, gbif, :species)

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
