################################################################################
# Script: rolling_vi_tri_analysis.R
# Author: Konstantin Engelmayer
# Date:   27.08.2025
#
# Description:
#   Rolling-window analysis linking daily Landsat vegetation indices (VIs)
#   to per-plot tree-ring chronologies (TRI). The pipeline:
#     1) Load extracted Landsat surface reflectance pixels
#     2) Compute 8-daily VIs
#     3) Load (or build) detrended TRI chronologies
#     4) Compute rolling cumulative sums of each VI across windows (1–91 days)
#     5) Correlate cumulative VIs to TRI (same-year lag = 0 and next-year lag = 1)
#     6) Save a tidy correlation table (heat_all)
#
# Inputs:
#   - data/analysis_ready_data/SDC_extracted.csv
#   - data/analysis_ready_data/TRI_chronologies.csv
#   - scripts/r/functions/extract_TRW.R
#   - scripts/thesis_ready_scripts/functions/extract_SDC_grouped_by_plot.R
#
# Outputs:
#   - corr_table_vi_tri.rds  (tidy correlation table)
#
# Notes:
#   - forecast::tsclean() is applied with frequency = 46 (~8‑day composites)
#     to remove outliers while preserving gaps (replace.missing = FALSE).
################################################################################

# ─────────────────────────────────────────────────────────────────────────────
# 0. LIBRARIES & HELPER SOURCES
# ─────────────────────────────────────────────────────────────────────────────
# Core data wrangling & dates
library(dplyr)       # mutate(), across(), group_by(), summarise(), joins
library(tidyr)       # pivot_longer()/pivot_wider()
library(lubridate)   # year(), yday(), as_date()

# Time series & rolling ops
library(zoo)         # rollapply()
library(forecast)    # tsclean()

# Functional & plotting
library(purrr)       # map_dfc(), walk()

source("scripts/r/functions/extract_TRW.R")
source("scripts/thesis_ready_scripts/functions/extract_SDC_grouped_by_plot.R")

# ─────────────────────────────────────────────────────────────────────────────
# 1. SATELLITE DATA  – daily VIs per plot
# ─────────────────────────────────────────────────────────────────────────────

## 1A. (Optional) run extractor to create analysis-ready reflectance pixels ----
# SDC_extracted <- extract_SDC_points(
#   folder_paths = c(
#     "/data/satellite_data/SDC/kellerwald_lahntal",
#     "/data/satellite_data/SDC/lindenber_eifel_koenigsforst"
#   ),
#   points_path  = "data/vector_data/trees_all_plots.gpkg",
#   method       = "original"
# )
# write.csv(SDC_extracted, "data/analysis_ready_data/SDC_extracted.csv", row.names = FALSE)

## 1B. load reflectance & add calendar fields ---------------------------------
SDC_df <- read.csv("data/analysis_ready_data/SDC_extracted.csv")

names(SDC_df)[1] <- "plot"                 # ensure first column is named 'plot'

SDC_df <- SDC_df |>
  mutate(
    date = as.Date(date),                   # guarantee Date type (not POSIXct/IDate)
    year = year(date),                      # calendar year (e.g., 2019)
    doy  = yday(date)                       # day of year (1–365/366)
  )

## 1C. outlier-clean reflectance bands per plot -------------------------------
band_cols <- c("blue", "green", "red", "nir", "swir1", "swir2")

SDC_clean <- SDC_df %>%
  arrange(plot, date) %>%                   # keep dates ordered within plot
  group_by(plot) %>%
  mutate(across(
    all_of(band_cols),
    ~ {
      tsx <- ts(.x, frequency = 46)         # ≈ 8‑day composites across a year
      as.numeric(tsclean(tsx, iterate = 2, replace.missing = FALSE))
    }
  )) %>%
  ungroup()

## 1D. compute vegetation indices (per row) -----------------------------------
VI_raw <- SDC_clean %>%
  mutate(
    NDVI  = (nir - red)   / (nir + red),
    GNDVI = (nir - green) / (nir + green),
    EVI   = 2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1),
    GRVI  = (green - red) / (green + red),
    NIRv  = NDVI * nir,
    kNDVI = (1 - exp(- (nir - red)^2 / (2 * 0.2^2))) /
      (1 + exp(- (nir - red)^2 / (2 * 0.2^2))),
    SIPI  = ((nir - blue) / (nir - red)) * -1,
    NPCI  = (red - blue) / (red + blue) * -1,   # (R680 – R430)/(R680 + R430)
    NMDI  = (nir - (swir1 - swir2)) / (nir + (swir1 + swir2)),
    GVMI  = ((nir + 0.10) - (swir2 + 0.02)) / ((nir + 0.10) + (swir2 + 0.02)),
    CIG   = (nir / green) - 1,
    NDWI  = (nir - swir1) / (nir + swir1)
  ) %>%
  select(-all_of(band_cols))                # drop raw bands if not needed downstream

## 1E. outlier-clean VI columns per plot --------------------------------------
index_cols <- setdiff(names(VI_raw), c("plot", "date", "year", "doy", "area"))

VI_clean <- VI_raw %>%
  arrange(plot, date) %>%
  group_by(plot) %>%
  mutate(across(
    all_of(index_cols),
    ~ {
      tsx <- ts(.x, frequency = 46)
      as.numeric(tsclean(tsx, iterate = 2, replace.missing = FALSE))
    }
  )) %>%
  ungroup()

# ─────────────────────────────────────────────────────────────────────────────
# 2. TREE‑RING CHRONOLOGIES  – detrended TRI per plot
# ─────────────────────────────────────────────────────────────────────────────
## 2A. (Optional) build chronologies from raw TRW ------------------------------
# TRW_df <- extract_TRW(
#   xls_folder  = "data/dendro_data/MW3_Daten/Dendro-Data/d-XLS",
#   points_path = "data/vector_data/trees_all_plots.gpkg"
# )
# chronologies_by_plot <- TRW_df |>
#   group_by(plot) |>
#   group_modify(~{
#     wide  <- pivot_wider(.x, Year, names_from = tree_id, values_from = ring_width)
#     rwi   <- detrend(as.rwl(wide[, -1]), method = "Spline", nyrs = 13, f = 0.6)
#     chron <- chron(as.data.frame(rwi), prefix = paste0("CH_", unique(.x$plot)))
#     chron |>
#       tibble::rownames_to_column("year") |>
#       mutate(year = as.numeric(year) + 1649,  # calendar offset
#              plot = unique(.x$plot))
#   }) |>
#   ungroup() |>
#   rename(TRI = starts_with("std"))         # final TRI column
# chronologies_by_plot <- na.omit(chronologies_by_plot)
# write.csv(chronologies_by_plot,
#           "data/analysis_ready_data/TRI_chronologies.csv", row.names = FALSE)

## 2B. load precomputed chronologies ------------------------------------------
chronologies_by_plot <- read.csv("data/analysis_ready_data/TRI_chronologies.csv") %>%
  dplyr::filter(year < 2023)

# ─────────────────────────────────────────────────────────────────────────────
# 3. ROLLING‑SUM VIs  (1–91‑day windows, daily resolution)
# ─────────────────────────────────────────────────────────────────────────────
accum_windows <- 1:91                            # candidate window lengths

## 3A. long table: one VI per row ------------------------------------------------
SDC_long <- VI_clean |>
  pivot_longer(
    -c(plot, date, year, doy, area),           # keep date + calendar cols (+ area if present)
    names_to = "VI", values_to = "value"
  )

## 3B. rolling cumulative sums per window -------------------------------------
SDC_cum_continuous <- SDC_long |>
  group_by(plot, VI) |>
  group_modify(~{
    bind_cols(
      .x,
      map_dfc(
        accum_windows,
        \(w) tibble(!!paste0("cum", w) := zoo::rollapply(.x$value, w, sum,
                                                         align = "right",
                                                         fill = NA))
      )
    )
  }) |>
  ungroup()

## 3C. tidy: one row = one window length --------------------------------------
SDC_cum_long <- SDC_cum_continuous |>
  pivot_longer(
    starts_with("cum"),
    names_to = "window", names_prefix = "cum",
    values_to = "cum_value"
  ) |>
  mutate(window = as.integer(window))

# ─────────────────────────────────────────────────────────────────────────────
# 4. JOIN TRI and tag forest type (coniferous / deciduous)
# ─────────────────────────────────────────────────────────────────────────────

cum_TRI_filtered <- SDC_cum_long |>
  left_join(chronologies_by_plot, by = c("plot", "year")) |>
  dplyr::filter(year >= 2000) |>
  mutate(
    forest_type = case_when(
      plot %in% c("PF01","SF03","SF06","SF07","SF10","SF11",
                  "SF13","SF15","SF16","SF21") ~ "coniferous",
      plot %in% c("PF02","PF03","SF01","SF02","SF04","SF05","SF08",
                  "SF09","SF12","SF14","SF17","SF18","SF19","SF20") ~ "deciduous",
      TRUE ~ NA_character_
    )
  )

# ─────────────────────────────────────────────────────────────────────────────
# 5. CORRELATION HEAT‑TABLE  &  best (window, DOY)
# ─────────────────────────────────────────────────────────────────────────────
## 5A. same‑year (lag = 0) correlations ---------------------------------------
heat_curr <- cum_TRI_filtered |>
  group_by(plot, forest_type, VI, window, doy) |>
  summarise(
    n           = sum(!is.na(cum_value) & !is.na(TRI)),
    correlation = if (n >= 3) cor(cum_value, TRI, use = "complete.obs") else NA_real_,
    p_value     = if (n >= 3) cor.test(cum_value, TRI)$p.value else NA_real_,
    .groups = "drop"
  ) |>
  mutate(lag = 0)

## 5B. next‑year (lag = 1) correlations ---------------------------------------
# Align cumulative VI from year t with TRI from year t+1 within each (plot, VI, window, DOY)
heat_prev <- cum_TRI_filtered %>%
  select(plot, forest_type, VI, window, doy,
         cum_prev  = cum_value,
         year_prev = year) %>%
  mutate(year = year_prev + 1) %>%
  inner_join(
    chronologies_by_plot %>% select(plot, year, TRI_current = TRI),
    by = c("plot", "year")
  ) %>%
  group_by(plot, forest_type, VI, window, doy) %>%
  summarise(
    n           = sum(!is.na(cum_prev) & !is.na(TRI_current)),
    correlation = if (n >= 3) cor(cum_prev, TRI_current, use = "complete.obs") else NA_real_,
    p_value     = if (n >= 3) cor.test(cum_prev, TRI_current)$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(lag = 1)

## 5C. combine & save -----------------------------------------------------------
heat_all <- bind_rows(heat_curr, heat_prev)
saveRDS(heat_all, "data/analysis_ready_data/corr_table_vi_tri.rds")

# Optional: write CSV for quick inspection
# write.csv(heat_all, "corr_table_vi_tri.csv", row.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# End of file
# ─────────────────────────────────────────────────────────────────────────────
