using JuliaDB
using Phylo
using Compat
using ClimatePref


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
    tree = read.newick('Qian2016.tree')
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

# Inference step
using Optim
using Distributions

# Create fake tree
n = 10
dist = Ultrametric{PolytomousTree{LeafInfo, Vector{Float64}}}(n)
tree = rand(dist)

function varcovar(tree::AbstractTree)
    tips = collect(nodenamefilter(isleaf, tree))
    root = collect(nodenamefilter(isroot, tree))[1]
    V = zeros(Float64, length(tips), length(tips))
    for i in 1:(length(tips) - 1)
        for j in i+1:length(tips)
            V[i, i] =  distance(tree, root, tips[i])
            V[j, j]= V[i,i]
            inter = getancestors(tree, tips[i]) ∩ getancestors(tree, tips[j])
            common = indmax(map(x-> distance(tree, root, x), inter))
            V[i, j] = distance(tree, root, inter[common])
            V[j, i] = V[i, j]
        end
    end
    return V
end

V = varcovar(tree)
traits = rand(Normal(0, 5), n)
#tips = collect(nodenamefilter(isleaf, tree))
#for i in eachindex(traits)
#    setnoderecord!(tree, tips[i], [traits[i]])
#end
O = ones(n)

LL(x) = 1/2 * (n * log(2π) + log(abs(det(x[1] * V))) +
transpose(traits - x[2] * O) * inv(x[1] * V) * (traits - x[2] * O))

result = optimize(LL, [0.1, 0.1])
opts = Optim.minimizer(result)
H = ForwardDiff.hessian(LL, opts)
se = sqrt.(diag(inv(H)))


opt1 = inv(transpose(O) * inv(V) * O) * (transpose(O) * inv(V) * traits)
opt2 = transpose(traits - opt1 * O) * inv(V) * (traits - opt1 * O)/n


res = fitBrownian(tree, traits)

na_mean(x) = mean(x[.!isnan.(x)])
# Fit to real tree

soltraits = map(x -> na_mean(select(getnoderecord(t, x), :tavg)), cross_species_names)
res = fitBrownian(t, ustrip.(soltraits))

fitLambda(tree, traits)

fitLambda(t, ustrip.(soltraits))

# Test against R code
include(joinpath(Pkg.dir("Phylo"), "src", "rcall.jl"));
R"library(ape)"
rt = rcall(:rtree, 10)
jt = NamedTree(rt)
traits = rand(Normal(0, 5), 10)
fitBrownian(jt, traits)
@rput jt; @rput traits
R"library(geiger)
names(traits) = jt$tip.label
fitContinuous(jt, traits, model='BM', control = list(method ='Nelder-Mead'))"

fitLambda(jt, traits)
    R"library(geiger)
    names(traits) = jt$tip.label
    fitContinuous(jt, traits, model='lambda', , control = list(method ='Nelder-Mead'))"
# Test real tree against R code
sol_traits = ustrip.(soltraits)
@rput cross_species_names; @rput sol_traits
R" library(ape)
    tree = read.tree('./data/Qian2016.tree')
    names_tree = tree$tip.label
    keep_tips = match(cross_species_names, names_tree)
    drop_tips = names_tree[-keep_tips]
    names(sol_traits) = cross_species_names
    newtree = drop.tip(tree, drop_tips, trim.internal = T)
    fitContinuous(newtree, sol_traits, model = 'BM')
"
@rget newtree
fitBrownian(newtree, sol_traits)

    R" library(ape)
        tree = read.tree('./data/Qian2016.tree')
        names_tree = tree$tip.label
        keep_tips = match(cross_species_names, names_tree)
        drop_tips = names_tree[-keep_tips]
        names(sol_traits) = cross_species_names
        newtree = drop.tip(tree, drop_tips, trim.internal = T)
        fitContinuous(newtree, sol_traits, model = 'lambda')
    "
@rget newtree
fitLambda(newtree, sol_traits)


using Phylo
using Compat
using ClimatePref
using JuliaDB

IndexedTables.NextTable() = table([])
# Load tree data and select names of species
t = open(io -> parsenewick(io, PolytomousTree{LeafInfo, IndexedTables.NextTable}),
        "data/Qian2016.tree")
tree_names = collect(nodenamefilter(isleaf, t))
species_names = join.(split.(tree_names, "_"), " ")

gbif_names = load("data/gbif_species")
gbif_species = keys(collect(gbif_names), :species)

cross_species = gbif_species ∩ species_names
cross_species_names = join.(split.(cross_species, " "), "_")
species_names = table(cross_species, names = [:names])
save(species_names, "joint_names")


addprocs(8)
using JuliaDB
# Find which species cross over
gbif = load("data/output_new2")
names = load("data/joint_names")
selection = filter(x -> x.species in select(names, :names), gbif)
save(selection, "gbif_tree")


using Phylo
using Compat
using ClimatePref
using JuliaDB
using FileIO
using Unitful
gbif = load("data/Grouped_genera")
datavals = find(map(x->typeof(x)<: DataValues.DataValue, select(gbif, :na_mean)))
map(x-> select(gbif, :na_mean)[x] = select(gbif, :na_mean)[x].value, datavals)
gbif = groupreduce(na_mean, gbif, :species, select = :na_mean)
save(gbif, "data/Grouped_genera2")

gbif = JuliaDB.load("data/Grouped_genera2")
spp_names = select(gbif, :species)
spp_names = join.(split.(spp_names, " "), "_")

IndexedTables.NextTable() = table([])
# Load tree data and select names of species
t = open(io -> parsenewick(io, PolytomousTree{LeafInfo, IndexedTables.NextTable}),
        "data/Qian2016.tree")
tree_names = getleafnames(t)
cross_species_names_un = tree_names ∩ spp_names
cross_species_names = join.(split.(cross_species_names_un, "_"), " ")
drop_tip!(t, cross_species_names_un)

newtab = filter(p -> p.species in cross_species_names, gbif)
newtab = pushcol(newtab, :tavg, ustrip.(select(newtab, :na_mean)))
newtab = popcol(newtab, :na_mean)
newtab = pushcol(newtab, :species_name, join.(split.(select(newtab, :species), " "), "_"))
newtab = popcol(newtab, :species)
FileIO.save("data/Tree_annotations.csv", newtab)

map(cross_species_names, cross_species_names_un) do spp_name, tree_name
    clean_dat = filter(p -> p.species == spp_name, gbif)
    setnoderecord!(t, tree_name, clean_dat)
end

matches = map(x-> find(x .== select(gbif, :species)), cross_species_names)
traits = ustrip.(select(gbif, :na_mean)[vcat(matches...)])
res = fitBrownian(t, traits)

addprocs(20)
@everywhere using Phylo
@everywhere using IndexedTables
tips = collect(nodenamefilter(isleaf, t))
root = collect(nodenamefilter(isroot, t))[1]
V = SharedArray(zeros(Float64, length(tips), length(tips)))
@sync @parallel for i in 1:(length(tips) - 1)
    for j in (i+1):length(tips)
        V[i, i] =  distance(t, root, tips[i])
        V[j, j]= V[i,i]
        inter = getancestors(t, tips[i]) ∩ getancestors(t, tips[j])
        common = indmax(map(x-> distance(t, root, x), inter))
        V[i, j] = distance(t, root, inter[common])
        V[j, i] = V[i, j]
    end
end

t = open(io -> parsenewick(io, PolytomousTree{LeafInfo, Float64}),
        "../GBIF_tree.tree")
traits = CSV.read("../Tree_annotations.csv")
traits[:,2] = join.(split.(traits[:,2], " "), "_")
ord = map(x-> find(traits[:,2] .== x)[1], tree_names)
traits = traits[ord, :]
@everywhere d = SharedArray{Float64,1}(Int64(nnodes*(nnodes+1)/2))
@everywhere ta = SharedArray{Float64,1}(Int64(nnodes*(nnodes+1)/2))
for i in 1:(nnodes-1)
    d[((i-1)*(nnodes-1)):((i)*(nnodes-1))] = map(x-> distance(t, tree_names[i], x), tree_names[i+1:end])
    ta[((i-1)*(nnodes-1)+i):((i)*(nnodes-1))] = abs(traits[i, 1] .- traits[(i+1):end, 1])
    n = n + 1
end
for i in 1:(nnodes-1)
    if i == 1
      d = map(x-> distance(t, tree_names[i], x), tree_names[i+1:end])
      ta = abs(traits[i, 1] .- traits[(i+1):end, 1])
    else
      d = [d; map(x-> distance(t, tree_names[i], x), tree_names[i+1:end])]
      ta = [ta; abs(traits[i, 1] .- traits[(i+1):end, 1])]
    end
end
map(x-> distance(tree_names[i], x), tree_names[i+1:end])

found_species = unique(select(total, :species))
setdiff(species_names, found_species)
