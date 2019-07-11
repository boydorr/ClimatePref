
library("MODISTools")
prods = mt_products()
sites = mt_sites()
bands = mt_bands(product = "MOD13Q1")
subset <- mt_subset(product = "MOD13Q1",
                    lat = 40,
                    lon = -110,
                    band = "250m_16_days_EVI",
                    start = "2004-01-01",
                    end = "2004-03-31",
                    progress = FALSE)


library(MODIS)
world <- list(xmax = -179, xmin = 179, ymin = -89, ymax = 89)
worldEVI <- runGdal(product = "MOD13C2",
        extent = raster::extent(-180, 180, -90, 90),
        begin = "2000.02.01", end = "2018.12.31",
        SDSstring = "01", outDirPath = "~/Documents/EVI")


