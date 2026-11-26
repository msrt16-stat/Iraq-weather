## =====================================================================
## get_boundaries.R -- download the Iraq governorate boundaries (GADM)
##
## GADM's licence does not allow redistribution of its data, so the
## shapefile is not stored in this repository. This script downloads
## GADM version 3.6, level 1 (governorates) for Iraq and writes it to
## Iraq_shp/gadm36_IRQ_1.shp, the path expected by maps_plotting.R.
##
## Usage (from the repository root):   Rscript get_boundaries.R
## Requires: geodata, terra, sf
## Alternative: download "Iraq, shapefile, version 3.6" manually from
## https://gadm.org/download_country_v3.html and unzip into Iraq_shp/.
## =====================================================================

for (p in c("geodata", "terra", "sf")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

dir.create("Iraq_shp", showWarnings = FALSE)
irq <- geodata::gadm(country = "IRQ", level = 1, version = "3.6",
                     path = tempdir())
irq_sf <- sf::st_as_sf(irq)
sf::st_write(irq_sf, "Iraq_shp/gadm36_IRQ_1.shp", delete_layer = TRUE)
message("Wrote Iraq_shp/gadm36_IRQ_1.shp (", nrow(irq_sf), " governorates)")
