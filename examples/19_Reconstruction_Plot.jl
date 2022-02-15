using JLD
using JuliaDB
using DataFrames
using Plots
using Unitful.DefaultSymbols
using StatsBase
using ClimatePref
using ClimatePref.Unitful
using PhyloNetworks

# Load raw data vs reconstructed
raw_traits = JuliaDB.load("data/Phylo_traits_adjust")
recon_traits = JLD.load("data/AncRecon_adjust.jld", "traits")
top_common_names = JLD.load("data/Common_species_names.jld", "spp_names")
raw_traits_filter = filter(p -> p.tipNames in join.(split.(top_common_names, " "), "_"), raw_traits)

old = collect(JuliaDB.select(raw_traits_filter, :tmin))
new = recon_traits[!, :tmin]

# Plot against one another
scatter(old ./ K, new, grid = false, xlabel = "Raw minimum temperature (K)",
       ylabel = "Reconstructed minimum temperature (K)", ms = 2, label =  "", markercolour = :grey, markerstrokecolour = false)
plot!(250:300, 250:300, label = "", colour = :black)
Plots.pdf("examples/Reconstructedvsrawtmin.pdf")


old = collect(JuliaDB.select(raw_traits_filter, :tmax))
new = recon_traits[!, :tmax]

# Plot against one another
scatter(old ./ K, new, grid = false, xlabel = "Raw maximum temperature (K)",
       ylabel = "Reconstructed maximum temperature (K)", ms = 2, label =  "", markercolour = :grey, markerstrokecolour = false)
plot!(270:310, 270:310, label = "", colour = :black)
Plots.pdf("examples/Reconstructedvsrawtmax.pdf")

# Calculate root mean squared deviation
rmsd(old ./K, new)


## Plot ancestral trait reconstruction for tmin
# Read in tip data and tree
phylo_traits = JuliaDB.load("examples/Phylo_traits_new")
tree = readTopology("examples/Qian2016.tree")

# Cut tree down to common species
top_common_names = JLD.load("examples/Common_species_names.jld", "spp_names")
tipnames = tipLabels(tree)
tip_names = join.(split.(tipnames, "_"), " ")
missing_species = setdiff(tip_names, top_common_names)
for i in eachindex(missing_species)
    deleteleaf!(tree, join(split(missing_species[i], " "), "_"))
end

# Filter traits to common species
phylo_traits_filter = filter(p -> p.tipNames in join.(split.(top_common_names, " "), "_"), phylo_traits)
dat = DataFrame(collect(phylo_traits_filter))

# Reconstruct tmin
anc = ancestralStateReconstruction(dat[:, [:tipNames, :tmin]], tree)
ancInt = expectations(anc)
anc = ancestralStateReconstruction(dat[:, [:tipNames, :ssr]], tree)
ancInt[!, :ssr] = expectations(anc)[!, :condExpectation]
anc = ancestralStateReconstruction(dat[:, [:tipNames, :tp]], tree)
ancInt[!, :tp] = expectations(anc)[!, :condExpectation]
# Save results
JLD.save("examples/AncRecon_dat.jld", "anc", ancInt)
writeTopology(tree, "examples/Qian_filtered.newick")

using Phylo
using JLD
tree2 = open(parsenewick, "examples/Qian_filtered.newick")
ancInt = JLD.load("examples/AncRecon_dat.jld", "anc")
nodeindices = map(a -> PhyloNetworks.getIndexNode(a, tree), ancInt[:, :nodeNumber])
nodes = tree2.nodes[nodeindices]
for i in 1:nrow(ancInt)
    setnodedata!(tree2, nodes[i], "tmin",  ancInt[i, :condExpectation]-273.15)
    setnodedata!(tree2, nodes[i], "ssr",  ancInt[i, :ssr]/1000)
    setnodedata!(tree2, nodes[i], "tp",  ancInt[i, :tp] * 1000)
end
sort!(tree2)
plot(tree2, line_z = "tmin", lw = 1, layout = (@layout [a b c]), size = (1600, 800), 
title = "A", title_loc = :left, colorbar_title = "Minimum temperature (°C)", guidefontsize = 16, titlefontsize = 18,
tickfontsize = 14, linecolor = :matter)
plot!(tree2, line_z = "tp", lw = 1, subplot = 2, linecolor = :haline, size = (1600, 800), 
title = "B", title_loc = :left, colorbar_title = "Total precipitation (mm)", guidefontsize = 16, titlefontsize = 18,
tickfontsize = 14)
plot!(tree2, line_z = "ssr", lw = 1, subplot = 3, linecolor = :cividis, size = (1600, 800), 
title = "C", title_loc = :left, colorbar_title = "Solar radiation (kJ/m²)", guidefontsize = 16, titlefontsize = 18,
tickfontsize = 14)

Plots.pdf("examples/Ancestral_trait_recon_sep.pdf")