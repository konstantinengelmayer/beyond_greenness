# ─────────────────────────────────────────────────────────────
#  POINT-BASED CLIMATE, PHENOLOGY & TOPOGRAPHY SUMMARY → PLOT & AREA LEVEL
#  Study area: Five regions in Germany (DE)
#
#  Functions: site_variables
#  Author:  Konstantin Engelmayer
#  Date:  2025-06-12
#
#  This script:
#    1. Extracts climate variables at tree sampling points
#    2. Aggregates monthly climate data to plot level
#    3. Filters climate to plot-year vegetation periods (SOS → EOS)
#    4. Computes climate, phenology, topography summaries per area
#    5. Builds a colour-shaded LaTeX table
#
#  Input data (all relative to project root):
#    ├─ data/vector_data/all_merged.shp              ← area polygons
#    ├─ data/vector_data/trees_all_plots.gpkg        ← sampling points (with plot id)
#    ├─ data/climate/dwd_4_variables/…               ← 4 stacked raster bricks
#    ├─ data/analysis_ready_data/phenology_df.csv    ← plot-year SOS / EOS
#    └─ data/analysis_ready_data/topo_ai_df.csv      ← plot topography (elev_m, slope_deg)
#
#  Output:
#    ├─ data/analysis_ready_data/climate_df.csv      ← plot-month climate (Step 8)
#    └─ veg_table.tex                                ← LaTeX table (Step 8)
#
#  Author:  Konstantin Engelmayer
#  Date:    2025-08-05
# ─────────────────────────────────────────────────────────────

# 1 ── libraries ──────────────────────────────────────────────
library(sf)         # vector data handling
library(terra)      # raster data handling
library(dplyr)      # data manipulation
library(tidyr)      # reshaping data
library(lubridate)  # date handling
library(kableExtra) # LaTeX tables
library(scales)     # rescale for colour shading

# 2 ── paths (edit if needed) ─────────────────────────────────
#    - Vector data
shp_file    <- "data/vector_data/all_merged.shp"        # polygons (area names)
points_file <- "data/vector_data/trees_all_plots.gpkg"  # sampling points

#    - Climate rasters
raster_dir  <- "data/climate/dwd_4_variables/"
temp_file   <- file.path(raster_dir, "air_temp_mean_1991_2024_merged.tif") # °C ×10
prec_file   <- file.path(raster_dir, "precipitation_1991_2024_merged.tif") # mm
evapo_file  <- file.path(raster_dir, "evapo_p_1991_2024_merged.tif")       # mm
soilmoist_file <- file.path(raster_dir, "soi_moist_1991_2024_merged.tif")  # %

#    - CSVs
phen_file <- "data/analysis_ready_data/phenology_df.csv"
topo_file <- "data/analysis_ready_data/topo_ai_df.csv"

#    - Outputs
clim_out  <- "data/analysis_ready_data/climate_df.csv"
tex_out   <- "veg_table.tex"

# 3 ── vectors & rasters ──────────────────────────────────────

## 3.1 Read area polygons and assign area names ---------------------------
areas <- st_read(shp_file, quiet = TRUE)
areas$area <- areas$name <- c("Eifel", "Kellerwald", "Koenigsforst",
                              "Calderner Wald", "Lindenberger Wald")

## 3.2 Read sampling points (keep only `plot` column) ---------------------
points <- st_read(points_file, quiet = TRUE) |>
  select(plot)

## 3.3 Load climate rasters ----------------------------------------------
Tmean  <- rast(temp_file)  / 10  # Convert 0.1 °C → °C
Pmm    <- rast(prec_file)
ETp    <- rast(evapo_file)
SoilMs <- rast(soilmoist_file)

## 3.4 Label raster layers with "YYYY-MM" ---------------------------------
start_date  <- ymd("1991-01-01")
layer_names <- format(seq(start_date, by = "1 month", length.out = nlyr(Tmean)),
                      "%Y-%m")

names(Tmean)  <- layer_names
names(Pmm)    <- layer_names
names(ETp)    <- layer_names
names(SoilMs) <- layer_names

# 4 ── harmonise CRS & attach area names to points ────────────────────────
crs_target <- crs(Tmean)

points <- st_transform(points, crs_target) |>
  st_join(st_transform(areas, crs_target), left = FALSE)

# 5 ── extract raster values at points ────────────────────────────────────
pts_v <- vect(points)

temp_df  <- terra::extract(Tmean,   pts_v, ID = FALSE) |> as_tibble()
prec_df  <- terra::extract(Pmm,     pts_v, ID = FALSE) |> as_tibble()
evapo_df <- terra::extract(ETp,     pts_v, ID = FALSE) |> as_tibble()
soil_df  <- terra::extract(SoilMs,  pts_v, ID = FALSE) |> as_tibble()

# 6 ── tidy to long format (one row per point × month) ────────────────────

## 6.1 Create point metadata table ----------------------------------------
pt_meta <- points |>
  st_drop_geometry() |>
  mutate(pt_id = row_number()) |>
  select(pt_id, area, plot)

## 6.2 Reshape extracted raster tables ------------------------------------
temp_long <- temp_df |>
  mutate(pt_id = row_number()) |>
  pivot_longer(-pt_id, names_to = "date", values_to = "temp")

prec_long <- prec_df |>
  mutate(pt_id = row_number()) |>
  pivot_longer(-pt_id, names_to = "date", values_to = "prec")

evapo_long <- evapo_df |>
  mutate(pt_id = row_number()) |>
  pivot_longer(-pt_id, names_to = "date", values_to = "evapo")

soil_long <- soil_df |>
  mutate(pt_id = row_number()) |>
  pivot_longer(-pt_id, names_to = "date", values_to = "soil_moist")

# 7 ── merge variables into single climate table ──────────────────────────
clim_long <- temp_long |>
  left_join(prec_long,  by = c("pt_id", "date")) |>
  left_join(evapo_long, by = c("pt_id", "date")) |>
  left_join(soil_long,  by = c("pt_id", "date")) |>
  mutate(date  = ymd(paste0(date, "-01")),
         month = month(date)) |>
  left_join(pt_meta, by = "pt_id")

# 8 ── aggregate to mean per plot × month & save --------------------------
clim_plot_month <- clim_long |>
  group_by(plot, area, date, month) |>
  summarise(temp       = mean(temp, na.rm = TRUE),
            prec       = mean(prec, na.rm = TRUE),
            evapo      = mean(evapo, na.rm = TRUE),
            soil_moist = mean(soil_moist, na.rm = TRUE),
            .groups    = "drop")

write.csv(clim_plot_month, clim_out, row.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────
#  PART B – Climate, Phenology & Topography Summary
# ─────────────────────────────────────────────────────────────────────────

# ---- (0) read phenology & helper tables ---------------------------------
phenology  <- read.csv(phen_file)
plot_area  <- clim_plot_month |>
  distinct(plot, area)            # 1 row per plot

# ---- (0b) read topography ----------------------------------------------
#  File must contain columns: plot, elev_m, slope_deg
topo_ai <- read.csv(topo_file)

# ---- (0c) median topography per area ------------------------------------
topo_area <- topo_ai |>
  left_join(plot_area, by = "plot") |>
  group_by(area) |>
  summarise(
    median_elev  = median(dem_mean,    na.rm = TRUE),  # m a.s.l.
    median_slope = median(slope_mean, na.rm = TRUE),  # degrees
    .groups      = "drop"
  )

# ---- (0d) SOS & EOS medians by area -------------------------------------
sos_eos_area <- phenology |>
  left_join(plot_area, by = "plot") |>
  dplyr::filter(!is.na(SOS_doy), !is.na(EOS_doy)) |>
  group_by(area) |>
  summarise(
    median_SOS = median(SOS_doy, na.rm = TRUE),
    median_EOS = median(EOS_doy, na.rm = TRUE),
    .groups    = "drop"
  )

# ---- (1) add year + mid-month DOY to climate data -----------------------
clim_tbl <- clim_plot_month |>
  mutate(
    year    = year(date),
    doy_mid = yday(date) + days_in_month(date) / 2
  )

# ---- (2) keep months inside veg period ----------------------------------
veg_clim <- clim_tbl |>
  inner_join(phenology, by = c("plot", "year")) |>
  dplyr::filter(!is.na(SOS_doy), !is.na(EOS_doy)) |>
  dplyr::filter(doy_mid >= SOS_doy & doy_mid <= EOS_doy)

# ---- (3) plot-year aggregates ------------------------------------------
plot_year <- veg_clim |>
  group_by(area, plot, year) |>
  summarise(
    temp_mean_plot_year = mean(temp, na.rm = TRUE),  # °C
    prec_sum_plot_year  = sum(prec,  na.rm = TRUE),  # mm
    .groups = "drop"
  )

# ---- (4) area-year averages of plot values -----------------------------
area_year <- plot_year |>
  group_by(area, year) |>
  summarise(
    mean_temp_yr     = mean(temp_mean_plot_year, na.rm = TRUE),
    mean_prec_sum_yr = mean(prec_sum_plot_year,  na.rm = TRUE),
    .groups = "drop"
  )

# ---- (5) long-term means per area --------------------------------------
area_means <- area_year |>
  group_by(area) |>
  summarise(
    mean_temp_veg     = mean(mean_temp_yr,     na.rm = TRUE), # °C
    mean_sum_prec_veg = mean(mean_prec_sum_yr, na.rm = TRUE), # mm
    .groups = "drop"
  ) |>
  left_join(sos_eos_area, by = "area") |>
  left_join(topo_area,    by = "area")        # add elevation & slope

# ---- (6) assemble final summary table (no colours) ---------------------
summary_tbl <- area_year |>
  group_by(area) |>
  summarise(
    mean_temp_veg     = mean(mean_temp_yr,     na.rm = TRUE),
    mean_sum_prec_veg = mean(mean_prec_sum_yr, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(sos_eos_area, by = "area") |>
  left_join(topo_area,    by = "area") |>
  select(
    Area                            = area,
    `Median Elevation (m)`          = median_elev,
    `Median Slope (\\textdegree)`   = median_slope,
    `Median SOS (DOY)`              = median_SOS,
    `Median EOS (DOY)`              = median_EOS,
    `Veg-period mean $T$ (°C)`      = mean_temp_veg,
    `Veg-period mean $\\Sigma P$ (mm)` = mean_sum_prec_veg
  ) |>
  mutate(across(where(is.numeric), ~ round(.x, 0)))   # <-- NEW: no decimals



# ---- (7) plain LaTeX table (no cell_spec, no shading) -------------------
veg_table_tex <- kableExtra::kbl(
  summary_tbl,
  format   = "latex",
  booktabs = TRUE,
  escape   = FALSE,
  caption  = "Vegetation-period climate, phenology and topography by study area."
) |>
  kableExtra::kable_styling(position = "center",
                            latex_options = "hold_position")

# ---- (8) output LaTeX code ---------------------------------------------
cat(veg_table_tex)
writeLines(veg_table_tex, "veg_table.tex")
