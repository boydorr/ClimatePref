using ClimatePref
using JuliaDB
using JuliaDBMeta
gbif = ReadGBIF("final")
save(gbif, "Full_GBIF_new")

gbif = filter(g -> !ismissing(g.decimallatitude) & !ismissing(g.decimallongitude), gbif)
save(gbif, "Geo_GBIF_new")

ref = create_reference(0.75)
gbif = extractvalues(gbif, ref, :refid)
save(gbif, "Era_GBIF_new")

gardens = loadtable("gardens.csv", chunks = 1)
centroids = loadtable("Centroids.csv", chunks = 1)
centroids = @transform centroids {Latitude = :y, Longitude = :x}

gbif = mask(gbif, gardens, 0.02)
gbif = mask(gbif, centroids, 0.02)
