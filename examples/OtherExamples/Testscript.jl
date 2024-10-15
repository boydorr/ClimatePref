# SPDX-License-Identifier: BSD-2-Clause

# Add more workers if needed
addprocs(8)

using JuliaDB
# Load data
test = load("/Users/claireh/Documents/PhD/Data/GBIF/test/output")
# Group by species, with count of each
result = groupby(length, test, :species)
# Filter out blank species names
filter((x -> x.species != ""), test)
