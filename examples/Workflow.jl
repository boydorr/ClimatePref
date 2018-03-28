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
