using ClimatePref
using JuliaDB
using JuliaDBMeta

# Load and save GBIF in between runs to save on disk space

# Read in full gbif repo (broken up into Genus csvs first)
gbif = ReadGBIF("final")
save(gbif, "Full_GBIF_new")

# Filter out those records missing coordinates
gbif = filter(g -> !ismissing(g.decimallatitude) & !ismissing(g.decimallongitude), gbif)
save(gbif, "Geo_GBIF_new")

# Extract grid info for each record (used to match to ERA data later)
ref = create_reference(0.75)
gbif = extractvalues(gbif, ref, :refid)
save(gbif, "Era_GBIF_new")

# Load info on gardens and country centroids
gardens = loadtable("gardens.csv", chunks = 1)
centroids = loadtable("Centroids.csv", chunks = 1)
centroids = @transform centroids {Latitude = :y, Longitude = :x}

# Filter them out
gbif = mask(gbif, gardens, 0.02)
gbif = mask(gbif, centroids, 0.02)
save(gbif, "GBIF_mask_new")
