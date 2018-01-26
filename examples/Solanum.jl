using Unitful
using AxisArrays
using ClimatePref

import Unitful: °, °C, mm
import ArchGDAL
import Base.read
const AG = ArchGDAL

using JuliaDB
gbif = loadtable("/Users/claireh/Documents/PhD/Data/GBIF/Solanum",
       indexcols = [:phylum, :class, :order, :family, :genus, :species, :scientificname],
       type_detect_rows =5000,
       colparsers=Dict(:datasetkey=>String,
                       :occurrenceid=>String,
                       :gbifid=>String,
                       :locality=>String,
                       :publishingorgkey=>String,
                       :taxonkey=>String,
                       :institutioncode=>String,
                       :catalognumber=>String,
                       :recordnumber=>String))
save(gbif, "/Users/claireh/Documents/PhD/Data/GBIF/Solanum/")
