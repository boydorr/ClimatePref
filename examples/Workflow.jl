using JuliaDB
using Phylo
using Compat


IndexedTables.NextTable() = table([])
# Load tree data and select names of species
t = open(io -> parsenewick(io, PolytomousTree{LeafInfo, IndexedTables.NextTable}),
        "data/Qian2016.tree")
tree_names = collect(nodenamefilter(isleaf, t))
sol_names = filter(name -> contains(name, r"Solanum"), collect(nodenamefilter(isleaf, t)))
# Clean spp names to remove _ between Genus and Species
species_names = join.(split.(tree_names, "_"), " ")

using Unitful
# Load cleaned Solanum data
sol = load("data/Worldclim/sol_worldclim")
occ_names = unique(select(sol, :species))

# Find which species cross over
cross_species = occ_names ∩ species_names
cross_species_names = join.(split.(cross_species, " "), "_")
# Only 179 species cross over out of ~1000.

drop_tip!(t, cross_species_names)

using RCall
@rput cross_species_names
R" library(ape)
    tree = read.newick("Qian2016.tree")
    names_tree = tree$tip.label
    keep_tips = match(cross_species_names, names_tree)
    drop_tips = names_tree[-keep_tips]
    newtree = drop.tip(tree, drop_tips, trim.internal = T)
    write.tree(newtree, "Solanum_tree.tree")
"

IndexedTables.NextTable() = table([])
t = open(io -> parsenewick(io, PolytomousTree{LeafInfo, IndexedTables.NextTable}),
        "data/Solanum_tree.tree")
cross_species_names = collect(collect(nodenamefilter(isleaf, t)))
cross_species = join.(split.(cross_species_names, "_"), " ")
using Unitful
# Load cleaned Solanum data
sol = load("data/Worldclim/sol_worldclim")
occ_names = unique(select(sol, :species))
# Loop through each species, collect environmental data and assign to tree tip
map(cross_species, cross_species_names) do spp_name, tree_name
    clean_dat = filter(p-> p.species == spp_name, sol)
    #tavg = ustrip.(select(clean_dat, :tavg))
    setnoderecord!(t, tree_name, clean_dat)
end

