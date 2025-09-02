###############################################################################
# TOPOGRAPHIC PREDICTORS → PLOT/AREA SUMMARIES
# -----------------------------------------------------------------------------
#
#  Functions: topographics
#  Author:  Konstantin Engelmayer
#  Date:  2025-06-12
#
# What this script does:
#   1) Load libraries (sf, terra, dplyr, tidyr, lubridate)
#   2) Define input paths for polygons (areas), sampling points, and predictor raster
#   3) Read polygons and points; attach area names
#   4) Build a point metadata table (pt_id ↔ plot)
#   5) Read 4-band predictor raster (slope, aspect, dem, ai), align CRS, extract
#   6) Summarise predictors:
#        6.1 per-plot (mean DEM, slope, circular-mean aspect, mean AI)
#        6.2 per-area  (same metrics; points are assigned to areas via spatial join)
#   7) Write AREA-level summary to CSV (for merging with your climate table)
#
# Output:
#   • data/analysis_ready_data/topo_ai_df.csv   (AREA-level summary)
#
# Notes:
#   • Aspect is averaged circularly (proper handling of angles).
#   • If your shapefile already has area names, keep them; otherwise we set them.
###############################################################################

# 1 ── Libraries ──────────────────────────────────────────────────────────────
library(sf)         # vector data handling
library(terra)      # raster data handling
library(dplyr)      # data manipulation
library(tidyr)      # reshaping data
library(lubridate)  # date handling

# 2 ── Paths (edit if needed) ─────────────────────────────────────────────────
shp_file    <- "data/vector_data/all_merged.shp"                # polygons (areas)
points_file <- "data/vector_data/trees_all_plots.gpkg"          # sampling points
rast_file   <- "data/satellite_data/predictors/predictors.tif"  # 4-band predictors

# 3 ── Input vectors ──────────────────────────────────────────────────────────

## 3.1 Areas: read polygons and ensure area names are present
areas <- sf::st_read(shp_file, quiet = TRUE)

# If your file already carries an "area" or "name" column, keep it.
# Otherwise assign names here (order must match row order in 'areas').
if (!("area" %in% names(areas))) {
  areas$area <- c("Eifel", "Kellerwald", "Koenigsforst",
                  "Calderner Wald", "Lindenberger Wald")
}
if (!("name" %in% names(areas))) areas$name <- areas$area

## 3.2 Sampling points: keep only the plot identifier
points <- sf::st_read(points_file, quiet = TRUE) |>
  dplyr::select(plot)

# 4 ── Point IDs & metadata (no geometry) ─────────────────────────────────────
#    Create a stable key (pt_id) to link extracted values back to plots/areas.
pt_meta <- points |>
  sf::st_drop_geometry() |>
  dplyr::mutate(pt_id = dplyr::row_number()) |>
  dplyr::select(pt_id, plot)

# Also prepare a pt_id ↔ area map via spatial join (used for AREA summaries)
pt_area_map <- points |>
  dplyr::mutate(pt_id = dplyr::row_number()) |>
  sf::st_join(areas[, "area"]) |>
  sf::st_drop_geometry() |>
  dplyr::select(pt_id, plot, area)

# 5 ── Extract raster predictors at points ────────────────────────────────────

## 5.1 Load predictor raster (4 bands) and name layers
predictors <- terra::rast(rast_file)
names(predictors) <- c("slope", "aspect", "dem", "ai")  # rename for clarity

## 5.2 Ensure CRS alignment and extract at point locations
pts_v <- terra::vect(points)  # sf → terra

if (!terra::compareGeom(pts_v, predictors, stopOnError = FALSE)) {
  # reproject points to raster CRS if needed
  pts_v <- terra::project(pts_v, terra::crs(predictors))
}

topo_df <- terra::extract(predictors, pts_v, ID = FALSE) |>
  tibble::as_tibble() |>
  dplyr::mutate(pt_id = dplyr::row_number())  # same key as pt_meta

# 6 ── Summaries ──────────────────────────────────────────────────────────────

# Helper: circular mean for aspect in degrees (0–360)
circ_mean_deg <- function(deg) {
  rad <- deg * pi / 180
  ang <- atan2(mean(sin(rad), na.rm = TRUE), mean(cos(rad), na.rm = TRUE)) * 180 / pi
  (ang %% 360 + 360) %% 360  # normalize to [0,360)
}

## 6.1 Plot-level summary (kept for reference/use elsewhere)
topo_plot <- topo_df |>
  dplyr::left_join(pt_meta, by = "pt_id") |>
  dplyr::group_by(plot) |>
  dplyr::summarise(
    dem_mean    = mean(dem,   na.rm = TRUE),   # m a.s.l.
    slope_mean  = mean(slope, na.rm = TRUE),   # degrees
    aspect_mean = circ_mean_deg(aspect),       # circular mean (0–360°)
    ai_mean     = mean(ai,    na.rm = TRUE),   # aridity index (if applicable)
    .groups = "drop"
  )

## 6.2 AREA-level summary (as per your comment “summarise to AREA level”)
topo_area <- topo_df |>
  dplyr::left_join(pt_area_map, by = "pt_id") |>
  dplyr::group_by(area) |>
  dplyr::summarise(
    dem_mean    = mean(dem,   na.rm = TRUE),
    slope_mean  = mean(slope, na.rm = TRUE),
    aspect_mean = circ_mean_deg(aspect),
    ai_mean     = mean(ai,    na.rm = TRUE),
    .groups = "drop"
  )

# 7 ── Save outputs ───────────────────────────────────────────────────────────
# AREA-level table (intended to be merged with your area-level climate table)
dir.create("data/analysis_ready_data", recursive = TRUE, showWarnings = FALSE)
write.csv(topo_area, "data/analysis_ready_data/topo_ai_df.csv", row.names = FALSE)
