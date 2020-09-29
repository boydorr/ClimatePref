using ClimatePref
using JuliaDB
gbif = ReadGBIF("final")
save(gbif, "Full_GBIF")

gbif = filter(g -> !ismissing(g.decimallatitude) & !ismissing(g.decimallongitude), gbif)
save(gbif, "Geo_GBIF")

ref = create_reference(0.75)
gbif = extractvalues(gbif, ref, :refid)
save(gbif, "Era_GBIF")

gardens = loadtable("gardens.csv")
centroids = loadtable("Centroids.csv")
