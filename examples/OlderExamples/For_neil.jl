addprocs(8)
using JuliaDB

sol = load("/Users/claireh/Documents/PhD/Data/GBIF/Solanum/output")
types = map(x-> fieldtype(eltype(sol), x), fieldnames(eltype(sol)))
names = fieldnames(eltype(sol))
dict = Dict(zip(names,types))
dict[:decimallatitude] = DataValues.DataValue{Float64}
dict[:decimallongitude] = DataValues.DataValue{Float64}
dict[:year] = DataValues.DataValue{Int64}
#dir = "/Users/claireh/Documents/PhD/Data/For_neil/finalL"
#liliales = loadtable(dir,
#        indexcols = [:phylum, :class, :order, :family, :genus, :species, :scientificname],
#        colparsers=dict)
#save(liliales, "/Users/claireh/Documents/PhD/Data/For_neil/Liliales")
liliales = load("/Users/claireh/Documents/PhD/Data/For_neil/Liliales")
species = loadtable("/Users/claireh/Documents/PhD/Data/For_neil/testspecies.csv",
distributed = false)
gen = collect(select(species, :Genus))
spp = collect(select(species, :Species))
species_names = map((x,y) -> string(x," ", y), gen, spp)
without_syn = filter(x-> x.species in species_names, liliales)
tab = collect(without_syn)
using FileIO
FileIO.save("/Users/claireh/Documents/PhD/Data/For_neil/matches.csv", tab)


synonyms = loadtable("/Users/claireh/Documents/PhD/Data/For_neil/testsynonyms.csv",
distributed = false, colparsers = Dict(:Synonym => String))
syns = collect(select(synonyms, :Synonym))
extras = map(x->string(split(x, " ")[1], " ", split(x, " ")[2], " ", split(x, " ")[3], " ", split(x, " ")[4]), syns)
extras2 = map(x->string(split(x, " ")[1], " ", split(x, " ")[2], " ", split(x, " ")[3]), syns)
Synonym = map(x->split(x, " Synonym ")[1], syns)
all_synonyms = vcat(Synonym, extras, extras2)

#dir = "/Users/claireh/Documents/PhD/Data/For_neil/finalA"
#asparagales = loadtable(dir,
#        indexcols = [:phylum, :class, :order, :family, :genus, :species, :scientificname],
#        colparsers=dict)
#save(asparagales, "/Users/claireh/Documents/PhD/Data/For_neil/Asparagales")
asparagales = load("/Users/claireh/Documents/PhD/Data/For_neil/Asparagales")
without_syn = filter(x-> x.species in species_names, asparagales)
total = merge(tab, collect(without_syn))
using FileIO
FileIO.save("/Users/claireh/Documents/PhD/Data/For_neil/matches_without_synonyms.csv", total)


with_syn = filter(x-> x.scientificname in all_synonyms, liliales)
syn_tab = collect(with_syn)
with_syn = filter(x-> x.scientificname in all_synonyms, asparagales)
total_syn = merge(syn_tab, collect(with_syn))
FileIO.save("/Users/claireh/Documents/PhD/Data/For_neil/matches_with_synonyms.csv", total_syn)
