# 11. Phylogenetic analysis of plant climate preferences
using PhyloNetworks
using GLM
using JuliaDB
using JuliaDBMeta
using Unitful
using ClimatePref
using ClimatePref.Unitful
using StatsBase
using JLD
using OnlineStats
using Unitful.DefaultSymbols
using Statistics
using DataFrames
@everywhere using Unitful
@everywhere using Unitful.DefaultSymbols
@everywhere using OnlineStats
@everywhere using StatsBase

# Load tree
tree = readTopology("Qian2016.tree")
tipnames = tipLabels(tree)
tip_names = join.(split.(tipnames, "_"), " ")

# Get record IDs of records without citizen science
gbif = JuliaDB.load("/home/claireh/Documents/gbif/full_data/GBIF_filtered/")
gbif = pushcol(gbif, :UID, 1:length(gbif))
top = JuliaDB.loadtable("Top_institutions.csv")
citizen_science = collect(JuliaDB.select(filter(t->t.Citizen_science == "TRUE", top), :Institution))
gbif_fil = filter(g -> g.institutioncode ∉ citizen_science, gbif)
recordID = unique(collect(JuliaDB.select(gbif_fil, :UID)))
JLD.save("Record_ID_no_citizen_science.jld", "recordID", recordID)

# Filter cera records for these
cera_records = JuliaDB.load("CERA_JOIN_SIMPLE")
recordID = JLD.load("Record_ID_no_citizen_science.jld", "recordID")
cera_fil = filter(c -> c.UID ∈ recordID, cera_records)
# Load Species names
spp_names = JLD.load("Species_names.jld", "spp_names")
spp_ids = JLD.load("Species_names.jld", "spp_ids")
sppdict = Dict(zip(spp_ids, spp_names))
iddict = Dict(zip(spp_names, spp_ids))
spp_names = [sppdict[i] for i in numspp]

# Get top 5000 most common species
cross_species = spp_names ∩ tip_names
cross_ids = [iddict[i] for i in cross_species]
sorted_counts = countmap(spp)
sorted_counts = filter(kv -> kv.first in cross_ids, sorted_counts)
sorted_counts = sort(sorted_counts, byvalue=true, rev=true)
top_common_ids = collect(keys(sorted_counts))[1:5000]
top_common_names = [sppdict[i] for i in top_common_ids]
#JLD.save("Common_species_names.jld", "spp_names", top_common_names)
#p = bar(collect(values(sorted_counts))[1:5000], grid = false)
#png(p, "GBIF_counts.png")

# Filter GBIF data for common species and delete from tree
gbif_fil = filter(g->g.SppID in cross_ids, gbif)
missing_species = setdiff(tip_names, top_common_names)
for i in eachindex(missing_species)
    deleteleaf!(tree, join(split(missing_species[i], " "), "_"))
end
cera_fil = filter(c-> c.UID ∈ recordID, gbif_fil)

### LAMBDAS FOR RAW DATA ###

# Group gbif data by Species ID and get mean and percentiles of data
phylo_traits = @groupby gbif_fil :SppID {tmin = mean(ustrip(:tmin)), tmin10 = percentile(ustrip(:tmin), 10), tmax = mean(ustrip(:tmax)), tmax90 = percentile(ustrip(:tmax), 90), tmean = mean(ustrip(:tmean)), trng = mean(ustrip(:trng)), stl1 = mean(ustrip(:stl1mean)), stl2 = mean(ustrip(:stl2mean)), stl3 = mean(ustrip(:stl3mean)), stl4 = mean(ustrip(:stl4mean)), swvl1 = mean(ustrip(:swvl1mean)), swvl2 = mean(ustrip(:swvl2mean)), swvl3 = mean(ustrip(:swvl3mean)), swvl4 = mean(ustrip(:swvl4mean)), ssr = mean(ustrip(:ssrmean)), tp = mean(ustrip(:tpmean))}
