# 3. Post-join cleaning for botanic gardens & country centroids
using JuliaDB
using JuliaDBMeta
using Unitful.DefaultSymbols
using ClimatePref
using ClimatePref.Units

# Load info on gardens and country centroids
gardens = loadtable("gardens.csv", chunks = 1)
centroids = loadtable("Centroids.csv", chunks = 1)
centroids = @transform centroids {Latitude = :y, Longitude = :x}

# Filter them out
# gbif = load("GBIF_TPL")
# gbif = mask(gbif, gardens, 0.02)
# gbif = mask(gbif, centroids, 0.02)
# save(gbif, "GBIF_mask")

# Takes too much memory to do in one function, so do mask in steps, saving out at each turn
gardens = loadtable("gardens.csv", chunks = 1)
centroids = loadtable("Centroids.csv", chunks = 1)
centroids = @transform centroids {Latitude = :y, Longitude = :x}

coords = hcat(collect(select(gardens, :Latitude)), collect(select(gardens, :Longitude)))
ref = create_reference(0.02)

# Add in refval column to garden info of which reference grid square it falls in
gardens = pushcol(gardens, :refval, extractvalues(coords[:, 2] * °, coords[:, 1] * °, ref))
save(gardens, "Gardens")

coords = hcat(collect(select(centroids, :Latitude)), collect(select(centroids, :Longitude)))
ref = create_reference(0.02)

# Add in refval column to garden info of which reference grid square it falls in
centroids = pushcol(centroids, :refval, extractvalues(coords[:, 2] * °, coords[:, 1] * °, ref))
save(centroids, "Centroids")

# Load gbif data and add reference values for cleaning datasets
gbif = load("GBIF_TPL")
ref = create_reference(0.02)
lat = collect(select(gbif, :decimallatitude))
lon = collect(select(gbif, :decimallongitude))
refval = extractvalues(lon .* °, lat .* °, ref)
gbif = pushcol(gbif, :refval, refval)
gbif = reindex(gbif, :refval)
save(gbif, "GBIF_antijoin")

# Use anti-join to filter out those that have the same reference value
gbif = load("GBIF_antijoin")

gardens = load("Gardens")
gbif = join(gbif, gardens, how=:anti, lkey=:refval, rkey =:refval)
gbif = popcol(gbif, :refval)

centroids = load("Centroids")
gbif = join(gbif, centroids, how=:anti, lkey=:refval, rkey =:refval)
gbif = popcol(gbif, :refval)

# Alternative to anti-join using filtering
gbif = load("GBIF_antijoin")
gardens = load("Gardens")
refs = collect(select(gardens, :refval))
gbif = filter(g -> g.refval ∉ refs, gbif)

gbif = load("GBIF_Gardens")
centroids = load("Centroids")
refs = collect(select(centroids, :refval))
gbif = filter(g -> g.refval ∉ refs, gbif)
gbif = popcol(gbif, :refval)
save(gbif, "GBIF_filtered")

# Streamline data to only necessary columns

# Load GBIF and add in unique record ID and Species ID column
gbif = load("GBIF_filtered")
gbif = pushcol(gbif, :UID, 1:length(gbif))
spplist = collect(select(gbif, :species))
uniquespp = sort(unique(spplist))
sppdict = Dict(zip(uniquespp, collect(1:length(uniquespp))))
sppids = [sppdict[sp] for sp in spplist]
gbif = pushcol(gbif, :SppID, sppids)

# Save reversal of species dict for later analyses
sppdict = Dict(zip(collect(1:length(uniquespp)), uniquespp))
JLD.save("Species_names.jld", "spp_names", uniquespp, "spp_ids", collect(1:length(uniquespp)))

# Create date column
mth = collect(select(gbif, :month))
yr = collect(select(gbif, :year))
dt = (mth .* month) .+ (yr .* year)
dt = uconvert.(year, dt)
gbif = pushcol(gbif, :date, dt)

# Extract coordinates and compare to reference grid
lat = collect(select(gbif, :decimallatitude))
lon = collect(select(gbif, :decimallongitude))
ref = create_reference(0.75)
refval = extractvalues(lon .* °, lat .* °, ref)
gbif = pushcol(gbif, :refval, refval)

# Only keep necessary columns for easier joins later (smaller dataset)
small_gbif = select(gbif, (:UID, :SppID, :date, :refval))
small_gbif = reindex(small_gbif, (:refval, :date))
save(small_gbif, "Small_GBIF")
