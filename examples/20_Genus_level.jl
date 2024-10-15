# SPDX-License-Identifier: BSD-2-Clause

# 20. Phylogenetic analysis of plant climate preferences at different taxonomic levels
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
@everywhere using ClimatePref

# Get record IDs of records without citizen science
gbif = JuliaDB.load("/home/claireh/Documents/gbif/full_data/GBIF_filtered/")
gbif = pushcol(gbif, :UID, 1:length(gbif))
taxonomy = unique(collect(JuliaDB.select(gbif,
                                         (:species, :genus, :family, :order,
                                          :class, :phylum, :UID))))
taxonomy = table(taxonomy)
JuliaDB.save(taxonomy, "Taxonomy")

# Filter cera records for these
cera_records = JuliaDB.load("CERA_JOIN_SIMPLE")
#recordID = JLD.load("Record_ID_no_citizen_science.jld", "recordID")
taxonomy = JuliaDB.load("Taxonomy")
taxonomy = distribute(taxonomy, 1)
cera_taxo = JuliaDB.join(cera_records, taxonomy, how = :left, rkey = :UID,
                         lkey = :UID)
JuliaDB.save(cera_taxo, "CERA_taxo")

cera_taxo = JuliaDB.load("CERA_taxo")
taxonomy = JuliaDB.load("Taxonomy")
taxonomy = distribute(taxonomy, 1)
taxo_fil = filter(t -> t.species ∈ tip_names, taxonomy)
big_genera = @groupby taxonomy :genus {nspp = length(unique(:species))}
big_genera = filter(b -> b.nspp > 50, big_genera)

big_genera = JuliaDB.load("Big_genera")
## Genus loop ##
genera = unique(collect(JuliaDB.select(big_genera, :genus)))
for i in eachindex(genera)
    if isfile("Lambdas_genus_$(genera[i]).jld")
        continue
    end
    # Load tree
    tree = readTopology("Qian2016.tree")
    tipnames = tipLabels(tree)
    tip_names = join.(split.(tipnames, "_"), " ")

    genus = filter(c -> c.genus == genera[i], cera_taxo)
    spp_names = collect(JuliaDB.select(genus, :species))
    cross_species = spp_names ∩ tip_names

    genus_fil = filter(g -> g.species in cross_species, genus)
    missing_species = setdiff(tip_names, cross_species)
    for i in eachindex(missing_species)
        deleteleaf!(tree, join(split(missing_species[i], " "), "_"))
    end

    # Group gbif data by Species ID and get mean and percentiles of data
    phylo_traits = @groupby genus_fil :species {tmin = mean(ustrip(:tmin)),
                                                tmax = mean(ustrip(:tmax)),
                                                tmean = mean(ustrip(:tmean)),
                                                trng = mean(ustrip(:trng)),
                                                stl1 = mean(ustrip(:stl1mean)),
                                                stl2 = mean(ustrip(:stl2mean)),
                                                stl3 = mean(ustrip(:stl3mean)),
                                                stl4 = mean(ustrip(:stl4mean)),
                                                swvl1 = mean(ustrip(:swvl1mean)),
                                                swvl2 = mean(ustrip(:swvl2mean)),
                                                swvl3 = mean(ustrip(:swvl3mean)),
                                                swvl4 = mean(ustrip(:swvl4mean)),
                                                ssr = mean(ustrip(:ssrmean)),
                                                tp = mean(ustrip(:tpmean))}

    # Add in tip names to data and save
    phylo_traits = @transform phylo_traits {tipNames = join.(split.(:species,
                                                                    " "), "_")}

    # Convert to dataframe
    dat = DataFrame(collect(phylo_traits))

    # Fit lambda models and save
    lambdas = fitLambdas(tree, dat, names(dat)[2:(end - 1)], 0.1)
    JLD.save("Lambdas_genus_$(genera[i]).jld", "lambdas", lambdas)

    mins = [
        197.0K,
        197.0K,
        197.0K,
        0K,
        197.0K,
        197.0K,
        197.0K,
        197.0K,
        0.0m^3,
        0.0m^3,
        0.0m^3,
        0.0m^3,
        0.0J / m^2,
        0.0m
    ]
    maxs = [
        320.0K,
        320.0K,
        320.0K,
        80K,
        320.0K,
        320.0K,
        320.0K,
        320.0K,
        1.0m^3,
        1.0m^3,
        1.0m^3,
        1.0m^3,
        3.0e7J / m^2,
        0.1m
    ]

    # Load EVi and gbif counts
    total_evi_counts = JLD.load("Total_evi_counts.jld", "total")
    total_gbif_counts = JLD.load("Total_gbif_counts.jld", "total")

    # Adjustment - remove NaNs and Infs
    adjustment = total_gbif_counts ./ total_evi_counts
    adjustment[isnan.(adjustment)] .= 1
    adjustment[isinf.(adjustment)] .= 1

    # Apply adjustment to data grouped by Species ID
    phylo_traits_adj = @groupby genus_fil :species {tmin = adjust(uconvert.(K,
                                                                            :tmin),
                                                                  adjustment[:,
                                                                             1],
                                                                  mins[1],
                                                                  maxs[1]),
                                                    tmax = adjust(uconvert.(K,
                                                                            :tmax),
                                                                  adjustment[:,
                                                                             2],
                                                                  mins[2],
                                                                  maxs[2]),
                                                    tmean = adjust(uconvert.(K,
                                                                             :tmean),
                                                                   adjustment[:,
                                                                              3],
                                                                   mins[3],
                                                                   maxs[3]),
                                                    trng = adjust(uconvert.(K,
                                                                            :trng),
                                                                  adjustment[:,
                                                                             4],
                                                                  mins[4],
                                                                  maxs[4]),
                                                    stl1 = adjust(uconvert.(K,
                                                                            :stl1mean),
                                                                  adjustment[:,
                                                                             5],
                                                                  mins[5],
                                                                  maxs[5]),
                                                    stl2 = adjust(uconvert.(K,
                                                                            :stl2mean),
                                                                  adjustment[:,
                                                                             6],
                                                                  mins[6],
                                                                  maxs[6]),
                                                    stl3 = adjust(uconvert.(K,
                                                                            :stl3mean),
                                                                  adjustment[:,
                                                                             7],
                                                                  mins[7],
                                                                  maxs[7]),
                                                    stl4 = adjust(uconvert.(K,
                                                                            :stl4mean),
                                                                  adjustment[:,
                                                                             8],
                                                                  mins[8],
                                                                  maxs[8]),
                                                    swvl1 = adjust(:swvl1mean,
                                                                   adjustment[:,
                                                                              9],
                                                                   mins[9],
                                                                   maxs[9]),
                                                    swvl2 = adjust(:swvl2mean,
                                                                   adjustment[:,
                                                                              10],
                                                                   mins[10],
                                                                   maxs[10]),
                                                    swvl3 = adjust(:swvl3mean,
                                                                   adjustment[:,
                                                                              11],
                                                                   mins[11],
                                                                   maxs[11]),
                                                    swvl4 = adjust(:swvl4mean,
                                                                   adjustment[:,
                                                                              12],
                                                                   mins[12],
                                                                   maxs[12]),
                                                    ssr = adjust(:ssrmean,
                                                                 adjustment[:,
                                                                            13],
                                                                 mins[13],
                                                                 maxs[13]),
                                                    tp = adjust(:tpmean,
                                                                adjustment[:,
                                                                           14],
                                                                mins[14],
                                                                maxs[14])}
    phylo_traits_adj = @transform phylo_traits_adj {tipNames = join.(split.(:species,
                                                                            " "),
                                                                     "_")}

    # Convert to dataframe
    dat = DataFrame(collect(phylo_traits_adj))

    # Fit lambda models and save
    lambdas = fitLambdas(tree, dat, names(dat)[2:(end - 1)], 0.1)
    JLD.save("Lambdas_effort_genus_$(genera[i]).jld", "lambdas", lambdas)
end

using Plots
using JLD
using ClimatePref
using DataFrames
genera = ClimatePref.searchdir("data/genus", "Lambdas_genus")
lambdas = fill(0.0, length(genera), 14)
for i in eachindex(genera)
    lambda = JLD.load(joinpath("data/genus", genera[i]), "lambdas")
    lambdas[i, :] .= lambda
end

families = ClimatePref.searchdir("data/family", "Lambdas_family")
lambdas_f = fill(0.0, length(families), 14)
for i in eachindex(families)
    lambda = JLD.load(joinpath("data/family", families[i]), "lambdas")
    lambdas_f[i, :] .= lambda
end

orders = ClimatePref.searchdir("data/order", "Lambdas_order")
lambdas_o = fill(0.0, length(orders), 14)
for i in eachindex(orders)
    lambda = JLD.load(joinpath("data/order", orders[i]), "lambdas")
    lambdas_o[i, :] .= lambda
end

classes = ClimatePref.searchdir("data/class", "Lambdas_class")
lambdas_c = fill(0.0, length(classes), 14)
for i in eachindex(classes)
    lambda = JLD.load(joinpath("data/class", classes[i]), "lambdas")
    lambdas_c[i, :] .= lambda
end

lambdas_k = JLD.load("data/Lambdas_common_full.jld", "lambdas")

tmin = DataFrame(λ = lambdas[:, 1], type = fill("genus", length(genera)))
tmin = [tmin;
        DataFrame(λ = lambdas_f[:, 1], type = fill("family", length(families)))]
tmin = [tmin;
        DataFrame(λ = lambdas_o[:, 1], type = fill("order", length(orders)))]
tmin = [tmin;
        DataFrame(λ = lambdas_c[:, 1], type = fill("class", length(classes)))]
tmin = [tmin; DataFrame(λ = lambdas_k[1], type = "kingdom")]

scatter(tmin[!, :type], tmin[!, :λ], zcolor = tmin[!, :λ],
        m = cgrad([:red, :white, :blue], [0, 0.5, 1]), label = "")
plot!(["genus", "family", "order", "class"],
      [mean(lambdas), mean(lambdas_f), mean(lambdas_o), mean(lambdas_c)],
      colour = :black, label = "")
plot!(["class", "kingdom"], [mean(lambdas_c), lambdas_k[1]], linestyle = :dash,
      colour = :black, label = "")
Plots.pdf("examples/scatterphylo.pdf")

cera_taxo = JuliaDB.load("CERA_taxo")
taxonomy = JuliaDB.load("Taxonomy")
taxonomy = distribute(taxonomy, 12)
taxo_fil = filter(t -> t.species ∈ tip_names, taxonomy)
big_family = @groupby taxo_fil :family {nspp = length(unique(:species))}
big_family = filter(b -> b.nspp > 50, big_family)
JuliaDB.save(big_family, "Big_family")

big_family = JuliaDB.load("Big_family")
## Family loop ##
families = unique(collect(JuliaDB.select(big_family, :family)))
for i in eachindex(families)
    if isfile("Lambdas_family_$(families[i]).jld")
        continue
    end
    # Load tree
    tree = readTopology("Qian2016.tree")
    tipnames = tipLabels(tree)
    tip_names = join.(split.(tipnames, "_"), " ")

    family = filter(c -> c.family == families[i], cera_taxo)
    spp_names = collect(JuliaDB.select(family, :species))
    cross_species = spp_names ∩ tip_names

    family_fil = filter(g -> g.species in cross_species, family)
    missing_species = setdiff(tip_names, cross_species)
    for i in eachindex(missing_species)
        deleteleaf!(tree, join(split(missing_species[i], " "), "_"))
    end

    # Group gbif data by Species ID and get mean and percentiles of data
    phylo_traits = @groupby family_fil :species {tmin = mean(ustrip(:tmin)),
                                                 tmax = mean(ustrip(:tmax)),
                                                 tmean = mean(ustrip(:tmean)),
                                                 trng = mean(ustrip(:trng)),
                                                 stl1 = mean(ustrip(:stl1mean)),
                                                 stl2 = mean(ustrip(:stl2mean)),
                                                 stl3 = mean(ustrip(:stl3mean)),
                                                 stl4 = mean(ustrip(:stl4mean)),
                                                 swvl1 = mean(ustrip(:swvl1mean)),
                                                 swvl2 = mean(ustrip(:swvl2mean)),
                                                 swvl3 = mean(ustrip(:swvl3mean)),
                                                 swvl4 = mean(ustrip(:swvl4mean)),
                                                 ssr = mean(ustrip(:ssrmean)),
                                                 tp = mean(ustrip(:tpmean))}

    # Add in tip names to data and save
    phylo_traits = @transform phylo_traits {tipNames = join.(split.(:species,
                                                                    " "), "_")}

    # Convert to dataframe
    dat = DataFrame(collect(phylo_traits))

    # Fit lambda models and save
    try
        lambdas = fitLambdas(tree, dat, names(dat)[2:(end - 1)], 0.1)
        JLD.save("Lambdas_family_$(families[i]).jld", "lambdas", lambdas)
    catch
        continue
    end

    mins = [
        197.0K,
        197.0K,
        197.0K,
        0K,
        197.0K,
        197.0K,
        197.0K,
        197.0K,
        0.0m^3,
        0.0m^3,
        0.0m^3,
        0.0m^3,
        0.0J / m^2,
        0.0m
    ]
    maxs = [
        320.0K,
        320.0K,
        320.0K,
        80K,
        320.0K,
        320.0K,
        320.0K,
        320.0K,
        1.0m^3,
        1.0m^3,
        1.0m^3,
        1.0m^3,
        3.0e7J / m^2,
        0.1m
    ]

    # Load EVi and gbif counts
    total_evi_counts = JLD.load("Total_evi_counts.jld", "total")
    total_gbif_counts = JLD.load("Total_gbif_counts.jld", "total")

    # Adjustment - remove NaNs and Infs
    adjustment = total_gbif_counts ./ total_evi_counts
    adjustment[isnan.(adjustment)] .= 1
    adjustment[isinf.(adjustment)] .= 1

    # Apply adjustment to data grouped by Species ID
    phylo_traits_adj = @groupby family_fil :species {tmin = adjust(uconvert.(K,
                                                                             :tmin),
                                                                   adjustment[:,
                                                                              1],
                                                                   mins[1],
                                                                   maxs[1]),
                                                     tmax = adjust(uconvert.(K,
                                                                             :tmax),
                                                                   adjustment[:,
                                                                              2],
                                                                   mins[2],
                                                                   maxs[2]),
                                                     tmean = adjust(uconvert.(K,
                                                                              :tmean),
                                                                    adjustment[:,
                                                                               3],
                                                                    mins[3],
                                                                    maxs[3]),
                                                     trng = adjust(uconvert.(K,
                                                                             :trng),
                                                                   adjustment[:,
                                                                              4],
                                                                   mins[4],
                                                                   maxs[4]),
                                                     stl1 = adjust(uconvert.(K,
                                                                             :stl1mean),
                                                                   adjustment[:,
                                                                              5],
                                                                   mins[5],
                                                                   maxs[5]),
                                                     stl2 = adjust(uconvert.(K,
                                                                             :stl2mean),
                                                                   adjustment[:,
                                                                              6],
                                                                   mins[6],
                                                                   maxs[6]),
                                                     stl3 = adjust(uconvert.(K,
                                                                             :stl3mean),
                                                                   adjustment[:,
                                                                              7],
                                                                   mins[7],
                                                                   maxs[7]),
                                                     stl4 = adjust(uconvert.(K,
                                                                             :stl4mean),
                                                                   adjustment[:,
                                                                              8],
                                                                   mins[8],
                                                                   maxs[8]),
                                                     swvl1 = adjust(:swvl1mean,
                                                                    adjustment[:,
                                                                               9],
                                                                    mins[9],
                                                                    maxs[9]),
                                                     swvl2 = adjust(:swvl2mean,
                                                                    adjustment[:,
                                                                               10],
                                                                    mins[10],
                                                                    maxs[10]),
                                                     swvl3 = adjust(:swvl3mean,
                                                                    adjustment[:,
                                                                               11],
                                                                    mins[11],
                                                                    maxs[11]),
                                                     swvl4 = adjust(:swvl4mean,
                                                                    adjustment[:,
                                                                               12],
                                                                    mins[12],
                                                                    maxs[12]),
                                                     ssr = adjust(:ssrmean,
                                                                  adjustment[:,
                                                                             13],
                                                                  mins[13],
                                                                  maxs[13]),
                                                     tp = adjust(:tpmean,
                                                                 adjustment[:,
                                                                            14],
                                                                 mins[14],
                                                                 maxs[14])}
    phylo_traits_adj = @transform phylo_traits_adj {tipNames = join.(split.(:species,
                                                                            " "),
                                                                     "_")}

    # Convert to dataframe
    dat = DataFrame(collect(phylo_traits_adj))

    # Fit lambda models and save
    try
        lambdas = fitLambdas(tree, dat, names(dat)[2:(end - 1)], 0.1)
        JLD.save("Lambdas_effort_family_$(families[i]).jld", "lambdas", lambdas)
    catch
        continue
    end
end

cera_taxo = JuliaDB.load("CERA_taxo")
taxonomy = JuliaDB.load("Taxonomy")
taxonomy = distribute(taxonomy, 12)
taxo_fil = filter(t -> t.species ∈ tip_names, taxonomy)
big_order = @groupby taxo_fil :order {nspp = length(unique(:species))}
big_order = filter(b -> b.nspp > 50, big_order)
JuliaDB.save(big_order, "Big_order")

big_order = JuliaDB.load("Big_order")
## order loop ##
orders = unique(collect(JuliaDB.select(big_order, :order)))
for i in eachindex(orders)
    if isfile("Lambdas_order_$(orders[i]).jld")
        continue
    end
    # Load tree
    tree = readTopology("Qian2016.tree")
    tipnames = tipLabels(tree)
    tip_names = join.(split.(tipnames, "_"), " ")

    order = filter(c -> c.order == orders[i], cera_taxo)
    spp_names = collect(JuliaDB.select(order, :species))
    cross_species = spp_names ∩ tip_names

    order_fil = filter(g -> g.species in cross_species, order)
    missing_species = setdiff(tip_names, cross_species)
    for i in eachindex(missing_species)
        deleteleaf!(tree, join(split(missing_species[i], " "), "_"))
    end

    # Group gbif data by Species ID and get mean and percentiles of data
    phylo_traits = @groupby order_fil :species {tmin = mean(ustrip(:tmin)),
                                                tmax = mean(ustrip(:tmax)),
                                                tmean = mean(ustrip(:tmean)),
                                                trng = mean(ustrip(:trng)),
                                                stl1 = mean(ustrip(:stl1mean)),
                                                stl2 = mean(ustrip(:stl2mean)),
                                                stl3 = mean(ustrip(:stl3mean)),
                                                stl4 = mean(ustrip(:stl4mean)),
                                                swvl1 = mean(ustrip(:swvl1mean)),
                                                swvl2 = mean(ustrip(:swvl2mean)),
                                                swvl3 = mean(ustrip(:swvl3mean)),
                                                swvl4 = mean(ustrip(:swvl4mean)),
                                                ssr = mean(ustrip(:ssrmean)),
                                                tp = mean(ustrip(:tpmean))}

    # Add in tip names to data and save
    phylo_traits = @transform phylo_traits {tipNames = join.(split.(:species,
                                                                    " "), "_")}

    # Convert to dataframe
    dat = DataFrame(collect(phylo_traits))

    # Fit lambda models and save
    try
        lambdas = fitLambdas(tree, dat, names(dat)[2:(end - 1)], 0.1)
        JLD.save("Lambdas_order_$(orders[i]).jld", "lambdas", lambdas)
    catch
        continue
    end

    mins = [
        197.0K,
        197.0K,
        197.0K,
        0K,
        197.0K,
        197.0K,
        197.0K,
        197.0K,
        0.0m^3,
        0.0m^3,
        0.0m^3,
        0.0m^3,
        0.0J / m^2,
        0.0m
    ]
    maxs = [
        320.0K,
        320.0K,
        320.0K,
        80K,
        320.0K,
        320.0K,
        320.0K,
        320.0K,
        1.0m^3,
        1.0m^3,
        1.0m^3,
        1.0m^3,
        3.0e7J / m^2,
        0.1m
    ]

    # Load EVi and gbif counts
    total_evi_counts = JLD.load("Total_evi_counts.jld", "total")
    total_gbif_counts = JLD.load("Total_gbif_counts.jld", "total")

    # Adjustment - remove NaNs and Infs
    adjustment = total_gbif_counts ./ total_evi_counts
    adjustment[isnan.(adjustment)] .= 1
    adjustment[isinf.(adjustment)] .= 1

    # Apply adjustment to data grouped by Species ID
    phylo_traits_adj = @groupby order_fil :species {tmin = adjust(uconvert.(K,
                                                                            :tmin),
                                                                  adjustment[:,
                                                                             1],
                                                                  mins[1],
                                                                  maxs[1]),
                                                    tmax = adjust(uconvert.(K,
                                                                            :tmax),
                                                                  adjustment[:,
                                                                             2],
                                                                  mins[2],
                                                                  maxs[2]),
                                                    tmean = adjust(uconvert.(K,
                                                                             :tmean),
                                                                   adjustment[:,
                                                                              3],
                                                                   mins[3],
                                                                   maxs[3]),
                                                    trng = adjust(uconvert.(K,
                                                                            :trng),
                                                                  adjustment[:,
                                                                             4],
                                                                  mins[4],
                                                                  maxs[4]),
                                                    stl1 = adjust(uconvert.(K,
                                                                            :stl1mean),
                                                                  adjustment[:,
                                                                             5],
                                                                  mins[5],
                                                                  maxs[5]),
                                                    stl2 = adjust(uconvert.(K,
                                                                            :stl2mean),
                                                                  adjustment[:,
                                                                             6],
                                                                  mins[6],
                                                                  maxs[6]),
                                                    stl3 = adjust(uconvert.(K,
                                                                            :stl3mean),
                                                                  adjustment[:,
                                                                             7],
                                                                  mins[7],
                                                                  maxs[7]),
                                                    stl4 = adjust(uconvert.(K,
                                                                            :stl4mean),
                                                                  adjustment[:,
                                                                             8],
                                                                  mins[8],
                                                                  maxs[8]),
                                                    swvl1 = adjust(:swvl1mean,
                                                                   adjustment[:,
                                                                              9],
                                                                   mins[9],
                                                                   maxs[9]),
                                                    swvl2 = adjust(:swvl2mean,
                                                                   adjustment[:,
                                                                              10],
                                                                   mins[10],
                                                                   maxs[10]),
                                                    swvl3 = adjust(:swvl3mean,
                                                                   adjustment[:,
                                                                              11],
                                                                   mins[11],
                                                                   maxs[11]),
                                                    swvl4 = adjust(:swvl4mean,
                                                                   adjustment[:,
                                                                              12],
                                                                   mins[12],
                                                                   maxs[12]),
                                                    ssr = adjust(:ssrmean,
                                                                 adjustment[:,
                                                                            13],
                                                                 mins[13],
                                                                 maxs[13]),
                                                    tp = adjust(:tpmean,
                                                                adjustment[:,
                                                                           14],
                                                                mins[14],
                                                                maxs[14])}
    phylo_traits_adj = @transform phylo_traits_adj {tipNames = join.(split.(:species,
                                                                            " "),
                                                                     "_")}

    # Convert to dataframe
    dat = DataFrame(collect(phylo_traits_adj))

    # Fit lambda models and save
    try
        lambdas = fitLambdas(tree, dat, names(dat)[2:(end - 1)], 0.1)
        JLD.save("Lambdas_effort_order_$(orders[i]).jld", "lambdas", lambdas)
    catch
        continue
    end
end

cera_taxo = JuliaDB.load("CERA_taxo")
taxonomy = JuliaDB.load("Taxonomy")
taxonomy = distribute(taxonomy, 12)
taxo_fil = filter(t -> t.species ∈ tip_names, taxonomy)
big_class = @groupby taxo_fil :class {nspp = length(unique(:species))}
big_class = filter(b -> b.nspp > 50, big_class)
JuliaDB.save(big_class, "Big_class")

big_class = JuliaDB.load("Big_class")
## class loop ##
classes = unique(collect(JuliaDB.select(big_class, :class)))
for i in eachindex(classes)
    if isfile("Lambdas_class_$(classes[i]).jld")
        continue
    end
    # Load tree
    tree = readTopology("Qian2016.tree")
    tipnames = tipLabels(tree)
    tip_names = join.(split.(tipnames, "_"), " ")

    class = filter(c -> c.class == classes[i], cera_taxo)
    spp_names = collect(JuliaDB.select(class, :species))
    cross_species = spp_names ∩ tip_names

    class_fil = filter(g -> g.species in cross_species, class)
    missing_species = setdiff(tip_names, cross_species)
    for i in eachindex(missing_species)
        deleteleaf!(tree, join(split(missing_species[i], " "), "_"))
    end

    # Group gbif data by Species ID and get mean and percentiles of data
    phylo_traits = @groupby class_fil :species {tmin = mean(ustrip(:tmin)),
                                                tmax = mean(ustrip(:tmax)),
                                                tmean = mean(ustrip(:tmean)),
                                                trng = mean(ustrip(:trng)),
                                                stl1 = mean(ustrip(:stl1mean)),
                                                stl2 = mean(ustrip(:stl2mean)),
                                                stl3 = mean(ustrip(:stl3mean)),
                                                stl4 = mean(ustrip(:stl4mean)),
                                                swvl1 = mean(ustrip(:swvl1mean)),
                                                swvl2 = mean(ustrip(:swvl2mean)),
                                                swvl3 = mean(ustrip(:swvl3mean)),
                                                swvl4 = mean(ustrip(:swvl4mean)),
                                                ssr = mean(ustrip(:ssrmean)),
                                                tp = mean(ustrip(:tpmean))}

    # Add in tip names to data and save
    phylo_traits = @transform phylo_traits {tipNames = join.(split.(:species,
                                                                    " "), "_")}
    JuliaDB.save(phylo_traits, "Phylo_traits_Magnoliopsida")

    # Convert to dataframe
    dat = DataFrame(collect(phylo_traits))

    # Fit lambda models and save
    try
        lambdas = fitLambdas(tree, dat, names(dat)[2:(end - 1)], 0.1)
        JLD.save("Lambdas_class_$(classes[i]).jld", "lambdas", lambdas)
    catch
        continue
    end

    # mins = [197.0K, 197.0K, 197.0K, 0K, 197.0K, 197.0K, 197.0K, 197.0K, 0.0m^3, 0.0m^3, 0.0m^3, 0.0m^3, 0.0J/m^2, 0.0m]
    # maxs = [320.0K, 320.0K, 320.0K, 80K, 320.0K, 320.0K, 320.0K, 320.0K, 1.0m^3, 1.0m^3, 1.0m^3, 1.0m^3, 3.0e7J/m^2, 0.1m]

    # # Load EVi and gbif counts
    # total_evi_counts = JLD.load("Total_evi_counts.jld", "total")
    # total_gbif_counts = JLD.load("Total_gbif_counts.jld", "total")

    # # Adjustment - remove NaNs and Infs
    # adjustment = total_gbif_counts ./ total_evi_counts
    # adjustment[isnan.(adjustment)] .= 1
    # adjustment[isinf.(adjustment)] .= 1

    # # Apply adjustment to data grouped by Species ID
    # phylo_traits_adj = @groupby class_fil :species {tmin = adjust(uconvert.(K, :tmin), adjustment[:, 1], mins[1], maxs[1]),tmax = adjust(uconvert.(K, :tmax), adjustment[:, 2], mins[2], maxs[2]), tmean = adjust(uconvert.(K, :tmean), adjustment[:, 3], mins[3], maxs[3]), trng = adjust(uconvert.(K, :trng), adjustment[:, 4], mins[4], maxs[4]), stl1 = adjust(uconvert.(K, :stl1mean), adjustment[:, 5], mins[5], maxs[5]), stl2 = adjust(uconvert.(K, :stl2mean), adjustment[:, 6], mins[6], maxs[6]), stl3 = adjust(uconvert.(K, :stl3mean), adjustment[:, 7], mins[7], maxs[7]), stl4 = adjust(uconvert.(K, :stl4mean), adjustment[:, 8], mins[8], maxs[8]), swvl1 = adjust(:swvl1mean, adjustment[:, 9], mins[9], maxs[9]), swvl2 =  adjust(:swvl2mean, adjustment[:, 10], mins[10], maxs[10]), swvl3 =  adjust(:swvl3mean, adjustment[:, 11], mins[11], maxs[11]), swvl4 =  adjust(:swvl4mean, adjustment[:, 12], mins[12], maxs[12]), ssr =  adjust(:ssrmean, adjustment[:, 13], mins[13], maxs[13]), tp =  adjust(:tpmean, adjustment[:, 14], mins[14], maxs[14])}
    # phylo_traits_adj = @transform phylo_traits_adj {tipNames = join.(split.(:species, " "), "_")}

    # # Convert to dataframe
    # dat = DataFrame(collect(phylo_traits_adj))

    # # Fit lambda models and save
    # try 
    # lambdas = fitLambdas(tree, dat, names(dat)[2:end-1], 0.1)
    # JLD.save("Lambdas_effort_class_$(classes[i]).jld", "lambdas", lambdas)
    # catch
    #     continue
    # end
end

cera_taxo = JuliaDB.load("CERA_taxo")
taxonomy = JuliaDB.load("Taxonomy")
taxonomy = distribute(taxonomy, 12)
taxo_fil = filter(t -> t.species ∈ tip_names, taxonomy)
big_phyla = @groupby taxo_fil :phylum {nspp = length(unique(:species))}
big_phyla = filter(b -> b.nspp > 25, big_phyla)
JuliaDB.save(big_phyla, "Big_phyla")

big_phyla = JuliaDB.load("Big_phyla")
## class loop ##
phyla = unique(collect(JuliaDB.select(big_phyla, :phylum)))
for i in eachindex(phyla)
    if isfile("Lambdas_phylum_$(phyla[i]).jld")
        continue
    end
    # Load tree
    tree = readTopology("Qian2016.tree")
    tipnames = tipLabels(tree)
    tip_names = join.(split.(tipnames, "_"), " ")

    phylum = filter(c -> c.phylum == phyla[i], cera_taxo)
    spp_names = collect(JuliaDB.select(phylum, :species))
    cross_species = spp_names ∩ tip_names

    phylum_fil = filter(g -> g.species in cross_species, phylum)
    missing_species = setdiff(tip_names, cross_species)
    for i in eachindex(missing_species)
        deleteleaf!(tree, join(split(missing_species[i], " "), "_"))
    end

    # Group gbif data by Species ID and get mean and percentiles of data
    phylo_traits = @groupby phylum_fil :species {tmin = mean(ustrip(:tmin)),
                                                 tmax = mean(ustrip(:tmax)),
                                                 tmean = mean(ustrip(:tmean)),
                                                 trng = mean(ustrip(:trng)),
                                                 stl1 = mean(ustrip(:stl1mean)),
                                                 stl2 = mean(ustrip(:stl2mean)),
                                                 stl3 = mean(ustrip(:stl3mean)),
                                                 stl4 = mean(ustrip(:stl4mean)),
                                                 swvl1 = mean(ustrip(:swvl1mean)),
                                                 swvl2 = mean(ustrip(:swvl2mean)),
                                                 swvl3 = mean(ustrip(:swvl3mean)),
                                                 swvl4 = mean(ustrip(:swvl4mean)),
                                                 ssr = mean(ustrip(:ssrmean)),
                                                 tp = mean(ustrip(:tpmean))}

    # Add in tip names to data and save
    phylo_traits = @transform phylo_traits {tipNames = join.(split.(:species,
                                                                    " "), "_")}

    # Convert to dataframe
    dat = DataFrame(collect(phylo_traits))

    # Fit lambda models and save
    try
        lambdas = fitLambdas(tree, dat, names(dat)[2:(end - 1)], 0.1)
        JLD.save("Lambdas_phylum_$(phyla[i]).jld", "lambdas", lambdas)
    catch
        continue
    end

    # mins = [197.0K, 197.0K, 197.0K, 0K, 197.0K, 197.0K, 197.0K, 197.0K, 0.0m^3, 0.0m^3, 0.0m^3, 0.0m^3, 0.0J/m^2, 0.0m]
    # maxs = [320.0K, 320.0K, 320.0K, 80K, 320.0K, 320.0K, 320.0K, 320.0K, 1.0m^3, 1.0m^3, 1.0m^3, 1.0m^3, 3.0e7J/m^2, 0.1m]

    # # Load EVi and gbif counts
    # total_evi_counts = JLD.load("Total_evi_counts.jld", "total")
    # total_gbif_counts = JLD.load("Total_gbif_counts.jld", "total")

    # # Adjustment - remove NaNs and Infs
    # adjustment = total_gbif_counts ./ total_evi_counts
    # adjustment[isnan.(adjustment)] .= 1
    # adjustment[isinf.(adjustment)] .= 1

    # # Apply adjustment to data grouped by Species ID
    # phylo_traits_adj = @groupby phylum_fil :species {tmin = adjust(uconvert.(K, :tmin), adjustment[:, 1], mins[1], maxs[1]),tmax = adjust(uconvert.(K, :tmax), adjustment[:, 2], mins[2], maxs[2]), tmean = adjust(uconvert.(K, :tmean), adjustment[:, 3], mins[3], maxs[3]), trng = adjust(uconvert.(K, :trng), adjustment[:, 4], mins[4], maxs[4]), stl1 = adjust(uconvert.(K, :stl1mean), adjustment[:, 5], mins[5], maxs[5]), stl2 = adjust(uconvert.(K, :stl2mean), adjustment[:, 6], mins[6], maxs[6]), stl3 = adjust(uconvert.(K, :stl3mean), adjustment[:, 7], mins[7], maxs[7]), stl4 = adjust(uconvert.(K, :stl4mean), adjustment[:, 8], mins[8], maxs[8]), swvl1 = adjust(:swvl1mean, adjustment[:, 9], mins[9], maxs[9]), swvl2 =  adjust(:swvl2mean, adjustment[:, 10], mins[10], maxs[10]), swvl3 =  adjust(:swvl3mean, adjustment[:, 11], mins[11], maxs[11]), swvl4 =  adjust(:swvl4mean, adjustment[:, 12], mins[12], maxs[12]), ssr =  adjust(:ssrmean, adjustment[:, 13], mins[13], maxs[13]), tp =  adjust(:tpmean, adjustment[:, 14], mins[14], maxs[14])}
    # phylo_traits_adj = @transform phylo_traits_adj {tipNames = join.(split.(:species, " "), "_")}

    # # Convert to dataframe
    # dat = DataFrame(collect(phylo_traits_adj))

    # # Fit lambda models and save
    # try 
    # lambdas = fitLambdas(tree, dat, names(dat)[2:end-1], 0.1)
    # JLD.save("Lambdas_effort_phylum_$(phyla[i]).jld", "lambdas", lambdas)
    # catch
    #     continue
    # end
end
