using JuliaDB
using Phylo
using Compat

# Load tree data and select names of species
t = open(io -> parsenewick(io, NamedPolytomousTree), "data/Qian2016.tree")
tree_names = filter(name -> contains(name, r"Solanum"), collect(nodenamefilter(isleaf, t)))
# Clean spp names to remove _ between Genus and Species
tree_names = join.(split.(tree_names, "_"), " ")

# Load cleaned Solanum data
sol = load("data/Clean_data")
occ_names = unique(select(sol, :species))

# Find which species cross over
cross_species = occ_names ∩ tree_names

# Only 179 species cross over out of ~1000.

# Loop through each species, collect environmental data and assign to tree tip
map(cross_species) do spp_name
    clean_dat = filter(p-> p.species == spp_name, sol)
    select(clean_dat, :tavg)
