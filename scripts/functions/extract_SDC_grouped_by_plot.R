################################################################################
# Function: extract_SDC_grouped_by_plot
# Author: [Konstantin Engelmayer]
# Date:   [23.07.2025]
#
# Description:
#   Reads multi-band SDC raster stacks (blue, green, red, nir, swir1, swir2)
#   stored in separate folders per study area, extracts pixel values at
#   tree-location points, averages them per plot and image date, optionally
#   applies temporal smoothing, and returns the result as a tidy data frame.
#
# Arguments:
#   folder_paths (character vector): One path per study-area folder.
#   points_path  (character): Vector layer path readable by {terra}; must
#                             contain attributes 'area' (study-area code) and
#                             'plot' (plot identifier).
#   method       (character): One of "original" | "sgolay" | "rolling7" | "loess".
#                             For Master-thesis data use "original" (other
#                             smoothing was applied later).
#   sg_n         (integer):   Window length for the Savitzky–Golay filter
#                             (method = "sgolay").
#   sg_p         (integer):   Polynomial order for the Savitzky–Golay filter
#                             (method = "sgolay").
#   roll_k       (integer):   Window size for the centred moving average
#                             (method = "rolling7").
#   loess_span   (numeric):   Span parameter for LOESS smoothing
#                             (method = "loess").
#
# Details:
#   - File naming: acquisition dates are read from the seven-digit Julian-day
#     code 'YYYYDDD' positioned 10–4 characters before the end of each filename
#     (e.g., "AreaA_2022265.tif" → 2022-09-22).
#   - Scaling: raw 16-bit integers are multiplied by 0.0001 to obtain
#     reflectance in the range 0–1.
#   - Temporal smoothing (applied separately for each plot_id × area):
#       • original  : raw values (no smoothing)
#       • sgolay    : Savitzky–Golay via signal::sgolayfilt()
#       • rolling7  : centred moving average via zoo::rollmean(k = roll_k)
#       • loess     : locally weighted regression via stats::loess()
#
# Returns:
#   A base-R data.frame in wide format with columns:
#     - plot_id, area, date, blue, green, red, nir, swir1, swir2
#
# Dependencies:
#   - terra
#   - sf
#   - dplyr
#   - zoo
#   - signal
#   - stats
#   - lubridate
#
# Example usage:
#   df <- extract_SDC_grouped_by_plot(
#           folder_paths = c("/data/AreaA", "/data/AreaB"),
#           points_path  = "plots.gpkg",
#           method       = "sgolay",
#           sg_n = 7, sg_p = 2
#        )
################################################################################

extract_SDC_points <- function(folder_paths,
                               points_path,
                               method       = c("original", "sgolay",
                                                "rolling7", "loess"),
                               sg_n         = 11,   # Savitzky–Golay window
                               sg_p         = 2,    # Savitzky–Golay poly order
                               roll_k       = 7,    # rolling‑mean window
                               loess_span   = 0.1)  # LOESS span
{
  # Match the chosen smoothing method ------------------------------------------------
  method <- match.arg(method)
  
  # ---- 0 · Load required packages & sanity checks -----------------------------------
  pkgs <- c("terra", "sf", "data.table", "signal", "zoo", "parallel", "dplyr")
  for (p in pkgs)
    if (!requireNamespace(p, quietly = TRUE))
      stop("Package '", p, "' is required but not installed.")
  
  if (!length(folder_paths))
    stop("'folder_paths' is empty.")
  if (!file.exists(points_path))
    stop("'points_path' does not exist: ", points_path)
  
  # ---- 1 · Reference raster & re‑project point layer --------------------------------
  first_tif <- list.files(folder_paths[1], "\\.tif$", full.names = TRUE)[1]
  if (is.na(first_tif))
    stop("No .tif files found in the first folder: ", folder_paths[1])
  
  ref_rast <- terra::rast(first_tif)                   # use first image as CRS/extent
  
  pts <- terra::vect(points_path) |>
    terra::project(ref_rast)                           # re‑project to raster CRS
  
  # Ensure mandatory attributes exist -------------------------------------------------
  if (!all(c("area", "plot") %in% names(pts)))
    stop("Points layer must contain both 'area' and 'plot' columns.")
  
  # Helper · Fast YYYYDDD → Date ------------------------------------------------------
  fast_dates <- function(tifs) {
    as.Date(substr(basename(tifs),
                   nchar(basename(tifs)) - 10,
                   nchar(basename(tifs)) - 4),
            format = "%Y%j")
  }
  
  # ---- 2 · Split points by area (processing unit) -----------------------------------
  pts_by_area <- split(pts, pts$area)
  
  # ---- 3 · Main loop over areas ------------------------------------------------------
  results <- vector("list", length(pts_by_area))
  names(results) <- names(pts_by_area)
  
  bands <- c("blue", "green", "red", "nir", "swir1", "swir2")
  
  for (ar in names(pts_by_area)) {
    
    # Find matching folder (case‑insensitive) -----------------------------------------
    folder <- folder_paths[grepl(ar, folder_paths, ignore.case = TRUE)][1]
    if (is.na(folder)) next                                  # skip if none found
    
    tif_paths <- list.files(folder, "\\.tif$", full.names = TRUE)
    if (!length(tif_paths)) next                              # skip if folder empty
    
    dates <- fast_dates(tif_paths)
    n_img <- length(tif_paths)
    
    # Read all GeoTIFFs into a single SpatRaster --------------------------------------
    ras <- terra::rast(tif_paths)
    names(ras) <- paste0(rep(bands, n_img), "_",
                         rep(format(dates, "%Y%j"), each = 6))
    
    # Extract pixel values for every point --------------------------------------------
    vals <- terra::extract(ras,
                           pts_by_area[[ar]],
                           ID   = TRUE,    # row index of pts
                           bind = FALSE)   # lean output (matrix)
    
    # Convert to data.table and attach plot IDs ---------------------------------------
    vals_dt <- data.table::as.data.table(vals)
    vals_dt[, plot_id := pts_by_area[[ar]]$plot[ID]]
    vals_dt[, ID      := NULL]                       # drop helper column
    
    # Reshape long → wide and aggregate per plot/date ---------------------------------
    long <- data.table::melt(vals_dt,
                             id.vars         = "plot_id",
                             variable.name   = "band_date",
                             value.name      = "value",
                             variable.factor = FALSE)
    
    long[, c("band", "date_str") := tstrsplit(band_date, "_", fixed = TRUE)]
    long[, date := as.Date(date_str, format = "%Y%j")]
    
    wide <- data.table::dcast(long,
                              plot_id + date ~ band,
                              value.var     = "value",
                              fun.aggregate = mean)  # average across points
    
    wide[, area := ar]
    results[[ar]] <- wide
  }
  
  # Combine areas into a single data.table --------------------------------------------
  df <- data.table::rbindlist(results, use.names = TRUE, fill = TRUE)
  
  # ---- 4 · Scale reflectances (16‑bit → SR) -----------------------------------------
  df[, (bands) := lapply(.SD, `*`, 0.0001), .SDcols = bands]
  
  # ---- 5 · Optional temporal smoothing ----------------------------------------------
  if (method != "original") {
    
    # Define smoothing helpers --------------------------------------------------------
    sgolay <- function(x) if (length(x) < sg_n) NA_real_ else
      signal::sgolayfilt(x, p = sg_p, n = sg_n)
    
    roll7  <- function(x) zoo::rollmean(x, roll_k, fill = NA)
    
    loess1 <- function(x, d) {
      if (all(is.na(x))) return(x)
      predict(loess(x ~ d, span = loess_span, na.action = na.exclude), d)
    }
    
    # Apply chosen smoother per plot × area -------------------------------------------
    data.table::setorder(df, plot_id, area, date)
    
    df[, (bands) :=
         lapply(.SD,
                function(v)
                  switch(method,
                         sgolay   = sgolay(v),
                         rolling7 = roll7(v),
                         loess    = loess1(v, as.numeric(date)))),
       by = .(plot_id, area),
       .SDcols = bands]
  }
  
  # Return as base data.frame ----------------------------------------------------------
  as.data.frame(df)
}
