# SPDX-License-Identifier: BSD-2-Clause

using DataFrames, CSV, JSON, DataStructures

data = CSV.read("data/Global_plant_climate_envelopes.csv");
select!(data, Not(:SppID));
plantlist = JSON.parsefile("data/plant_list_2024-06.json");

species_genus = Dict{String, String}()
genus_family = Dict{String, String}()
synonym_to_species = Dict{String, String}()
species_to_synonyms = Dict{String, Set{String}}()
binomial_id = Dict{String, Int}()
name_id = Dict{String, Int}()
dropped = String[]

df_wfo = DataFrame(rank = String[], binomial = String[], name = String[],
                   species = Union{Nothing, String}[], genus = Union{Nothing, String}[],
                   family = Union{Nothing, String}[], order = Union{Nothing, String}[],
                   phylum = Union{Nothing, String}[], kingdom = Union{Nothing, String}[])
                   
for plant in plantlist
    if haskey(plant, "role_s")
        if plant["role_s"] == "accepted"
            row = (rank = plant["rank_s"], binomial = plant["full_name_string_alpha_s"],
                   name = plant["full_name_string_plain_s"],
                   species = get(plant, "placed_in_species_s", nothing),
                   genus = get(plant, "placed_in_genus_s", nothing),
                   family = get(plant, "placed_in_family_s", nothing),
                   order = get(plant, "placed_in_order_s", nothing),
                   phylum = get(plant, "placed_in_phylum_s", nothing),
                   kingdom = get(plant, "placed_in_kingdom_s", nothing))
            if row.rank != "code" && !haskey(binomial_id, row.binomial)
                push!(df_wfo, row)
                binomial_id[row.binomial] = nrow(df_wfo)
                name_id[row.name] = nrow(df_wfo)
                if row.rank == "species"
                    species_genus[row.binomial] = row.genus
                    species_to_synonyms[row.binomial] = Set{String}()
                elseif plant["rank_s"] == "genus"
                    genus_family[row.genus] = row.family
                end
            end
        end
    end
end

non_species = 0
for plant in plantlist
    if haskey(plant, "role_s")
        if plant["role_s"] == "synonym"
            if plant["rank_s"] == "species"
                accepted = plant["accepted_full_name_string_plain_s"]
                name = plant["full_name_string_alpha_s"]
                if haskey(binomial_id, accepted)
                    row = df_wfo[binomial_id[accepted], :]
                    if row.rank == "species"
                        if name != accepted && !haskey(binomial_id, name)
                            synonym_to_species[name] = accepted
                            push!(species_to_synonyms[accepted], name)
                        end
                    elseif row.rank == "subspecies"
                        label = row.genus * " " * row.species
                        if name != label && !haskey(binomial_id, name)
                            synonym_to_species[name] = label
                            push!(species_to_synonyms[label], name)
                        end
                    else
                        non_species += 1
                    end
                elseif haskey(name_id, accepted)
                    row = df_wfo[name_id[accepted], :]
                    if row.rank == "species"
                        if name != row.binomial && !haskey(binomial_id, name)
                            synonym_to_species[name] = row.binomial
                            push!(species_to_synonyms[row.binomial], name)
                        end
                    elseif row.rank == "subspecies" || row.rank == "variety" || row.rank == "form" || row.rank == "subvariety"
                        label = row.genus * " " * row.species
                        if name != label && !haskey(binomial_id, name)
                            synonym_to_species[name] = label
                            push!(species_to_synonyms[label], name)
                        end
                    else
                        non_species += 1
                    end
                else
                    push!(dropped, accepted)
                end
            end
        end
    end
end

df_climate = copy(data)
df_climate.species = [get(synonym_to_species, s, s) for s in df_climate.species];
sort!(df_climate, :species)
select!(df_climate, Not(:genus))
select!(df_climate, Not(:family))
select!(df_climate, Not(:order))
select!(df_climate, Not(:class))
select!(df_climate, Not(:phylum))
spp = df_climate.species;
duplicates = spp[1:end-1] .== spp[2:end];
for i in findall(duplicates)
    n1 = df_climate[i, :numrecords]
    n2 = df_climate[i+1, :numrecords]
    for col in names(df_climate)[3:end]
        df_climate[i, col] = (df_climate[i, col] * n1 +
                              df_climate[i+1, col] * n2) / (n1 + n2)
    end
    df_climate[i+1, :numrecords] = n1 + n2
end
push!(duplicates, false);
deleterows!(df_climate, duplicates);
@info "Merging $(sum(duplicates)) duplicate species"
@info "$(nrow(df_climate)) accepted species"

CSV.write("data/climate.csv", df_climate)
