# SPDX-License-Identifier: BSD-2-Clause

# 1. Clean GBIF ready for joining with The Plant List
using ClimatePref
using JuliaDB

# Load and save GBIF in between runs to save on disk space

# Read in full gbif repo (broken up into Genus csvs first)
gbif = ReadGBIF("final")
save(gbif, "Full_GBIF")

# Filter out those records missing coordinates
gbif = filter(g -> !ismissing(g.decimallatitude) &
                   !ismissing(g.decimallongitude), gbif)
save(gbif, "Geo_GBIF")

# Extract grid info for each record (used to match to ERA data later)
ref = create_reference(0.75)
gbif = extractvalues(gbif, ref, :refid)
# Filter for blank species and reindex for faster joins
gbif = filter(g -> g.species != "", gbif)
gbif = reindex(gbif, :species)
save(gbif, "Era_GBIF")
