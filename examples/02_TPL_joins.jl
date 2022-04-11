# 2. Script to clean GBIF occurrences by plant list
using JuliaDB
using ClimatePref

# Load plant list dataset (from source text files generated from http://www.theplantlist.org)
tpl = readTPL("TPL")
save(tpl, "ThePlantList")

# Load plant list - needs garbage collector workaround to avoid Julia crashing - issue described here: https://github.com/JuliaComputing/MemPool.jl/issues/26
Base.GC.enable(false)
tpl = load("ThePlantList")
specieslist = collect(select(tpl, :species))
Base.GC.enable(true)
# Load gbif (make sure reindexed for quicker joins)
gbif = load("gbif/full_data/Era_GBIF")
# Join gbif/tpl to include only species in both tables
gbif_tpl = join(gbif, tpl, how = :inner, lkey = :species, rkey =:species)
save(gbif_tpl, "GBIF_TPL")

# If joins take up too much space on device - this is an alternative option using filter
gbif = filter(g -> g.species in specieslist, gbif)
save(gbif_tpl, "GBIF_TPL")
