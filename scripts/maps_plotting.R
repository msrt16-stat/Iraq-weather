## ============================================================
## Iraq humidity map plotting script (Figures 1-10, 2015-2024)
##
## Produces one PNG per year (2015-2024), each showing:
##   - an IDW-interpolated humidity choropleth over Iraq's governorates,
##     using a FIXED 5-class color scale (identical breakpoints in every
##     year) so that shading is directly comparable across the full
##     ten-year sequence -- see "FIXED COLOR SCALE" note below.
##   - a scatter-pie chart per governorate showing the relative (min-max
##     normalized) contribution of Temperature, Top-of-Atmosphere
##     radiation, Wind Speed, Precipitation, and Soil Moisture.
##   - governorate boundaries, labels, neighboring-country context,
##     north arrow, and scale bar.
##
## *** BUG FIXES applied relative to the version this was reconstructed
## from (see inline "FIX:" comments for exact locations) ***:
##   1. `neighbors` was referenced (via st_point_on_surface(neighbors))
##      but its definition was commented out, so the script could not
##      run at all. Restored (uncommented) below.
##   2. read.xlsx() requires the `openxlsx` package, not `readxl`
##      (which was the only Excel-reading package loaded). Added
##      library(openxlsx).
##   3. The predictor column names used in the yearly summarise() step
##      did not match the actual Excel column names ("precipitation
##      rain" / "soil" vs. the real "Precipitation" / "Soil Moisture").
##      Corrected to match the source file.
##   4. The `lookup` table (City ID -> governorate name, used only to
##      attach a `City` column to the shapefile for potential future
##      use) had governorates 10 onward in the wrong order relative to
##      the canonical City-ID numbering used throughout the rest of the
##      project's R scripts and the paper. This does NOT currently
##      affect the map output, since no downstream step in this script
##      actually uses iraq_sf$City (the choropleth, pie-chart positions,
##      and text labels all come directly from each record's own
##      Latitude/Longitude or from the shapefile's own NAME_1 field).
##      Corrected anyway, as a latent bug that would silently corrupt
##      results if anyone later adds a step that does join on City.
##   5. Removed a redundant unconditional install.packages() call at
##      the top (the conditional installer further down is sufficient
##      and avoids reinstalling the package on every run).
##
## FIXED COLOR SCALE: humidity_breaks/humidity_labels are defined ONCE,
## outside the year loop, using breakpoints derived from the pooled
## 2015-2024 governorate-year humidity distribution (see the paper's
## Methods section and fixed_humidity_breaks.R). Do NOT move this
## definition inside the yearly loop or recompute it per year -- doing
## so was the original bug this fixed: an independently-computed color
## scale for each year makes shading incomparable across years.
## ============================================================

## Set the path to your desired working folder
PROJECT_ROOT <- Sys.getenv("IRAQ_WEATHER_ROOT", unset = getwd())  # set by run_all.R; see README
OUTPUT_DIR   <- file.path(PROJECT_ROOT, "outputs", "10_maps_plotting")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
setwd(OUTPUT_DIR)
getwd()

## ---- 1. Load libraries ----
library(dplyr)
library(readxl)
library(openxlsx)      # FIX 2: required for read.xlsx() below (readxl alone is not enough)
library(sf)
library(ggplot2)
library(scatterpie)
library(ggspatial)
library(ggnewscale)
library(gstat)
library(rnaturalearth)
library(rnaturalearthdata)

## Install rnaturalearthdata only if not already present (do not force a
## reinstall on every run -- FIX 5 removed a redundant duplicate call here)
if (!requireNamespace("rnaturalearthdata", quietly = TRUE)) {
  install.packages("rnaturalearthdata")
}

ocean <- ne_download(
  scale = "medium",
  type = "ocean",
  category = "physical",
  returnclass = "sf"
)

## ---- 2. Read the humidity/climate data ----
## Expects the standard project data file (e.g. "Dalya5x.xlsx"), with
## columns: City, Humidity, Temp., "Top of Atmosphere", "Wind Speed",
## Precipitation, "Soil Moisture", Latitude, Longitude -- 18 governorates
## x 120 monthly rows.
path <- file.path(PROJECT_ROOT, "Dalya5x.xlsx")   # no interactive prompt -- place the data file directly in PROJECT_ROOT
stopifnot("Data file not found -- check DATA_FILE name/location" = file.exists(path))
raw <- read.xlsx(path)

## Diagnostic: print the exact column names as R actually sees them.
## If any step below still fails with a "not found" error, the output
## of this line is the first thing to check -- copy it exactly (quotes
## and all) rather than guessing a spelling.
print(colnames(raw))

## Rename columns POSITIONALLY (1st column -> City, 2nd -> Humidity, ...)
## rather than matching by exact name string. This is a defensive
## approach already used in this project's other scripts (e.g.
## weather_R_FINAL.R, rerun_sr_random_split_75_25.R), and avoids any
## dependence on exact header spelling, capitalization, or hidden
## whitespace. If your file's column ORDER differs from City, Humidity,
## Temp, TOA, Wind, Precip, Soil, Lat, Lon, adjust this vector to match
## -- but the exact text of each original header no longer matters.
stopifnot(ncol(raw) == 9)   # fails loudly and clearly if the column count doesn't match expectations
colnames(raw) <- c("City", "Humidity", "Temp", "TOA", "Wind", "Precip", "Soil", "Lat", "Lon")

## Reconstruct Year/Month: the source file has no explicit date column,
## but each governorate's 120 rows are known to run January 2015
## through December 2024 in file order (verified against the seasonal
## temperature cycle; see the paper's Methods section).
years  <- rep(2015:2024, each = 12)
months <- rep(1:12, times = 10)

raw <- raw %>%
  group_by(City) %>%
  mutate(
    Year = years,
    Month = months
  ) %>%
  ungroup()

## ---- 3. Read the Iraq governorate shapefile and prepare boundaries ----
iraq_sf <- st_read(file.path(PROJECT_ROOT, "Iraq_shp", "gadm36_IRQ_1.shp"))

## Correct two governorate names as spelled in the raw GADM shapefile
## ("Sala ad-Din" missing an "h"; "Al-Muthannia" with an extra "i") to
## their standard spellings, "Salah al-Din" and "Muthanna". Doing this
## once, immediately after reading the shapefile, means every downstream
## use of iraq_sf$NAME_1 -- the lookup/join below, and the map's text
## labels (which are drawn directly from iraq_sf$NAME_1) -- automatically
## uses the corrected spelling, with no separate fix needed elsewhere.
iraq_sf$NAME_1[iraq_sf$NAME_1 == "Sala ad-Din"]    <- "Salah al-Din"
iraq_sf$NAME_1[iraq_sf$NAME_1 == "Al-Muthannia"]   <- "Muthanna"

## Diagnostic only: inspect the shapefile's governorate-name spelling
## (NAME_1), after the correction above, before building `lookup` below,
## since `lookup$NAME_1` must match these strings exactly for the join
## to succeed.
print(data.frame(
  Row = 1:nrow(iraq_sf),
  NAME_1 = iraq_sf$NAME_1
))

## City-ID-to-governorate-name lookup table, using the same canonical
## numbering (1-18) and spelling used throughout the project's R scripts
## and the paper. NOTE: as of this script, iraq_sf$City (created by the
## left_join below) is not actually used downstream -- the choropleth,
## pie-chart positions, and text labels all come from each record's own
## Latitude/Longitude or from the shapefile's own NAME_1 field directly.
## This table is kept correct anyway (FIX 4) so that any future addition
## that does join on City will not silently misassign governorates.
lookup <- data.frame(
  City = 1:18,
  NAME_1 = c(
    "Dihok", "Ninawa", "Arbil", "As-Sulaymaniyah", "At-Ta'mim",
    "Salah al-Din", "Diyala", "Baghdad", "Al-Anbar", "Wasit",
    "Karbala'", "Babil", "An-Najaf", "Al-Qadisiyah", "Dhi-Qar",
    "Maysan", "Muthanna", "Al-Basrah"   # 14-16 verified by point-in-polygon on the data coordinates
  )
)

iraq_sf <- left_join(iraq_sf, lookup, by = "NAME_1")
st_crs(iraq_sf) <- 4326

## Outer national boundary (dissolved governorate polygons), used to
## clip the IDW interpolation surface so it doesn't spill outside Iraq
iraq_outer <- st_union(iraq_sf)

## World boundaries, for neighboring-country context and the ocean layer
world <- ne_countries(scale = "medium", returnclass = "sf")

## FIX 1: this was commented out in the version the script was
## reconstructed from, leaving `neighbors` undefined and causing every
## downstream reference to it (st_point_on_surface(neighbors) below) to
## fail with "object 'neighbors' not found". Restored here.
neighbors <- world %>%
  filter(admin %in% c(
    "Turkey",
    "Iran",
    "Syria",
    "Jordan",
    "Saudi Arabia",
    "Kuwait"
  ))

country_labels <- st_point_on_surface(neighbors)

country_labels <- cbind(
  st_drop_geometry(country_labels),
  st_coordinates(country_labels)
)

## Optional manual nudge of neighboring-country label positions (left
## disabled here as in the original; uncomment and tune if labels
## overlap borders or other map elements in your output)
# country_labels$X <- country_labels$X +
#   c(0.8, -15, 3, 4, 3, 0.5)
# country_labels$Y <- country_labels$Y +
#   c(-1, 0.2, 0.6, 0.8, 0.3, -0.2)

## Min-max normalization, used to bring the five pie-chart predictors
## (originally on very different scales/units) onto a common [0,1]
## range so each governorate's pie slices are comparable
min_max_norm <- function(x) {
  if (all(is.na(x)) || max(x, na.rm = TRUE) == min(x, na.rm = TRUE)) return(x)
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

## ---- FIXED humidity color scale (see header note above) ----
## Breakpoints derived from the pooled 2015-2024 governorate-year
## humidity distribution (Jenks natural breaks; see the paper's Methods
## section and fixed_humidity_breaks.R). Defined ONCE here, outside the
## year loop, and reused identically for every year's map.
humidity_breaks <- c(4.2, 5.0, 5.4, 5.9, 6.3, 6.9)
humidity_labels <- c("4.2 - 5.0", "5.0 - 5.4", "5.4 - 5.9", "5.9 - 6.3", "6.3 - 6.9")

humidity_colors <- c(
  "4.2 - 5.0" = "#D8F0FA",
  "5.0 - 5.4" = "#A7D8EF",
  "5.4 - 5.9" = "#66C2E5",
  "5.9 - 6.3" = "#2596BE",
  "6.3 - 6.9" = "#00688B"
)

output_dir <- getwd()   # already the dedicated outputs/10_maps_plotting subfolder set up above
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

## ---- 4. Build a standard interpolation grid over Iraq ----
grid_sf <- st_make_grid(iraq_sf, n = c(150, 150), what = "centers")
grid_sf <- st_sf(geometry = grid_sf)

## Precise bounding box for zooming the map tightly to Iraq
bbox <- st_bbox(iraq_sf)

## Dummy off-map data, used only to force ggplot to draw the complete
## 5-class legend even if a given year's IDW surface happens not to
## realize every class after clipping to Iraq's borders
legend_df <- data.frame(
  x = rep(bbox$xmin - 100, length(humidity_labels)),
  y = rep(bbox$ymin - 100, length(humidity_labels)),
  Humidity_Class = factor(
    humidity_labels,
    levels = humidity_labels
  )
)

## ---------------------------------------------------------
## 5. Main loop: one map per year, 2015-2024
## ---------------------------------------------------------
for (yr in 2015:2024) {

  cat("Processing and interpolating year:", yr, "...\n")

  ## Governorate-level annual means for the current year
  city_data <- raw %>%
    filter(Year == yr) %>%
    group_by(City) %>%
    summarise(
      Humidity        = mean(Humidity, na.rm = TRUE),
      Temperature     = mean(Temp, na.rm = TRUE),
      Top_Atmosphere  = mean(TOA, na.rm = TRUE),
      Wind_Speed      = mean(Wind, na.rm = TRUE),
      Precipitation   = mean(Precip, na.rm = TRUE),
      Soil_Moisture   = mean(Soil, na.rm = TRUE),
      Latitude        = mean(Lat, na.rm = TRUE),
      Longitude       = mean(Lon, na.rm = TRUE)
    )

  city_sf <- st_as_sf(city_data, coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)

  ## Inverse-distance-weighted (IDW) interpolation of humidity across
  ## the 150x150 grid, using each governorate's own point location
  idw_model <- gstat(formula = Humidity ~ 1, data = city_sf, nmax = 7, set = list(idp = 2))
  idw_pred  <- predict(idw_model, grid_sf)
  idw_clipped <- st_intersection(idw_pred, iraq_outer)   # clip to Iraq's outer boundary

  coords <- st_coordinates(idw_clipped)
  idw_df <- data.frame(
    x = coords[, 1],
    y = coords[, 2],
    Humidity_Val = idw_clipped$var1.pred
  )

  ## Classify the interpolated surface using the FIXED breaks defined
  ## above (identical for every year -- do not recompute per year)
  idw_df$Humidity_Class <- cut(
    idw_df$Humidity_Val,
    breaks = humidity_breaks,
    labels = humidity_labels,
    include.lowest = TRUE,
    right = TRUE
  )

  idw_df$Humidity_Class <- factor(
    idw_df$Humidity_Class,
    levels = humidity_labels
  )

  ## Normalize the five pie-chart predictors to [0,1] for this year
  pie_data <- city_data %>%
    mutate(across(c(Temperature, Top_Atmosphere, Wind_Speed, Precipitation, Soil_Moisture), min_max_norm))

  ## Governorate label positions (shapefile centroids, then manually
  ## nudged so labels don't overlap pie charts or each other)
  lab <- st_centroid(iraq_sf)
  lab <- cbind(
    st_drop_geometry(lab),
    st_coordinates(lab)
  )

  lab$X <- lab$X + c(
    -0.20, 0.30, -0.20, 0.45, -0.10,
    0.15, 0.20, 0.20, 0.45, 0.35,
    0.15, -0.15, 0.15, -0.20, 0.20,
    -0.20, 0.15, 0.25
  )

  lab$Y <- lab$Y + c(
    -0.11, -0.11, -0.15, 0.18, 0.5,
    0.10, -0.18, -0.20, -0.10, 0.10,
    -0.21, 0.10, -0.15, 0.15, -0.20,
    0.10, -0.20, -0.16
  )
  ## ---------------------------------------------------------
  ## 6. Assemble the map
  ## ---------------------------------------------------------
  p <- ggplot() +
    geom_sf(data = ocean, fill = "grey85", colour = NA) +
    ## Layer 0: neighboring-country boundaries
    geom_sf(data = world, fill = "white", color = "grey60", linewidth = 0.4) +

    ## Invisible layer that exists ONLY to force the complete legend to
    ## render (alpha = 0), even if some humidity classes are absent
    ## from this particular year's clipped surface
    geom_tile(
      data = legend_df,
      aes(x = x, y = y, fill = Humidity_Class),
      alpha = 0
    ) +

    ## Layer 1: the interpolated humidity choropleth
    geom_tile(
      data = idw_df,
      aes(x = x, y = y, fill = Humidity_Class)
    ) +

    scale_fill_manual(
      values = humidity_colors,
      limits = humidity_labels,
      breaks = humidity_labels,
      drop = FALSE,
      name = "Humidity",
      guide = guide_legend(
        override.aes = list(alpha = 1)
      )
    ) +

    ## Layer 2: governorate boundary lines
    geom_sf(
      data = iraq_sf,
      fill = NA,
      color = "grey20",
      linewidth = 0.5
    ) +

    ## Start a new fill scale so the pie-chart colors don't clash with
    ## the humidity choropleth's fill scale above
    new_scale_fill() +

    ## Layer 3: per-governorate scatter-pie charts of the five
    ## (min-max-normalized) direct predictors
    geom_scatterpie(
      data = pie_data,
      aes(x = Longitude, y = Latitude, r = 0.25),
      cols = c("Temperature", "Top_Atmosphere", "Wind_Speed", "Precipitation", "Soil_Moisture"),
      alpha = .9
    ) +
    scale_fill_manual(
      values = c(
        "Temperature"    = "#FFC000",
        "Top_Atmosphere" = "#ED7D31",
        "Wind_Speed"     = "#C00000",
        "Precipitation"  = "#FF00FF",
        "Soil_Moisture"  = "#0000FF"
      ),
      labels = c(
        "Temperature"    = "Temperature",
        "Top_Atmosphere" = "Top of Atmosphere",
        "Wind_Speed"     = "Wind speed",
        "Precipitation"  = "Precipitation",
        "Soil_Moisture"  = "Soil Moisture"
      ),
      name = "Independent\nvariables"
    ) +

    ## Layer 4: governorate name labels
    geom_text(
      data = lab,
      aes(x = X, y = Y, label = NAME_1),
      size = 2.8,
      color = "black",
      fontface = "bold"
    ) +

    geom_text(
      data = country_labels,
      aes(X, Y, label = admin),
      size = 3.8,
      colour = "grey40"
    ) +
    ## Zoom tightly to Iraq's extent
    coord_sf(
      xlim = c(bbox$xmin - 0.5, bbox$xmax + 0.5),
      ylim = c(bbox$ymin - 0.5, bbox$ymax + 0.5),
      expand = FALSE
    ) +

    ## North arrow and scale bar, stacked in the top-right corner
    annotation_north_arrow(
      location = "tr",
      which_north = "true",
      pad_x = unit(0.8, "in"),
      pad_y = unit(0.8, "in"),
      style = north_arrow_orienteering,
      height = unit(0.8, "cm"),
      width  = unit(0.8, "cm")
    ) +
    annotation_scale(
      location = "tr",
      pad_x = unit(0.12, "in"),
      pad_y = unit(1.3, "in"),
      width_hint = 0.22
    ) +

    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      ## Legend in the top-left corner, inside the map frame
      legend.position = c(0.03, 0.98),
      legend.justification = c(0, 1),

      legend.direction = "vertical",
      legend.box = "vertical",

      legend.box.background = element_rect(
        fill = alpha("white", 0.85),
        colour = "black"
      ),

      legend.key.size = unit(0.45, "cm"),
      legend.spacing.y = unit(0.1, "cm"),
      legend.margin = margin(4, 4, 4, 4),

      ## Black frame around the whole map
      plot.background = element_rect(fill = NA, color = "black", linewidth = 1),
      plot.margin = margin(2, 2, 2, 2)
    )

  ## Save this year's map
  file_name <- paste0(output_dir, "/Iraq_Climate_Interpolated_", yr, ".png")
  ggsave(file_name, plot = p, width = 9, height = 8, dpi = 300)

}

cat("All maps generated and saved successfully!\n")
