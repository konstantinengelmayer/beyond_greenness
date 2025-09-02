# ─────────────────────────────────────────────────────────────
#  PLOT-WISE START-/END-OF-SEASON (SOS & EOS) ANALYSIS  •  FIVE AREAS (DE)
#
#  Functions: sos_eos
#  Author:  Konstantin Engelmayer
#  Date:  2025-06-12
#
#  This script:
#    ✓ reads daily Landsat surface-reflectance pixels per forest plot
#    ✓ computes EVI2 and smooths with Savitzky-Golay filter
#    ✓ fits one double-logistic curve per plot-year (fixed d & amplitudes)
#    ✓ uses grid search + nlsLM for robust fitting
#    ✓ applies filtering/snapping to realistic DOY windows:
#        - SOS: DOY 90–160
#        - EOS: DOY 240–330
#    ✓ outputs tidy table (plot, year, SOS_doy, EOS_doy, conv status)
#    ✓ optionally generates helper plots (time series, boxplots)
#
#  Input data:
#    - Landsat daily surface reflectance pixels
#    - Points assigned to plots
#
#  Output:
#    - CSV: one row per plot-year with SOS & EOS DOY
#    - Helper plots (optional): EVI2 curves + SOS/EOS, boxplots
#

# ─────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
#  0 • LIBRARIES
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dplyr)      # wrangling
  library(tidyr)      # pivots
  library(lubridate)  # date helpers
  library(zoo)        # rollmedian
  library(signal)     # Savitzky-Golay
  library(minpack.lm) # nlsLM
  library(ggplot2)    # plots
  library(purrr)      # map_* helpers
  library(tibble)
  library(terra); library(sf)  # (only needed by extractor)
})

# ─────────────────────────────────────────────────────────────────────────────
#  1 • CUSTOM EXTRACTORS  (adjust paths to your repo) 
# ─────────────────────────────────────────────────────────────────────────────
# source("scripts/r/functions/extract_SDC_buffered.R")  # Landsat reflectance

# ─────────────────────────────────────────────────────────────────────────────
#  2 • READ PIXELS  →  DAILY EVI2 PER PLOT
# ─────────────────────────────────────────────────────────────────────────────

SDC_df <- read.csv("data/analysis_ready_data/SDC_extracted.csv")
evi_ts <- SDC_df %>%
  rename(plot = 1L) %>%                             # make first col explicit
  mutate(EVI2 = 2.5 * (nir - red) / (nir + 2.4 * red + 1)) %>%
  select(plot, date, EVI2, area) %>%
  arrange(plot, date)

# ─────────────────────────────────────────────────────────────────────────────
#  3 • CORE HELPERS
# ─────────────────────────────────────────────────────────────────────────────
kv_part <- function(t, a, b, c) {                   # Zhang curvature (K_V)
  z  <- exp(a + b*t);  bc <- b*c
  num1 <- b^3*c*z*(3*z*(1-z)*(1+z)^3)
  num2 <- (1+z)^2*(1 + 2*z - 5*z^2)
  den  <- ((1+z)^4 + (bc*z)^2)^(5/2)
  z * (num1 - (bc*z)^3*num2) / den
}
extrema <- function(v) which(diff(sign(diff(v))) != 0) + 1

# ─────────────────────────────────────────────────────────────────────────────
#  4 • ONE PLOT-YEAR FIT  →  SOS / EOS
# ─────────────────────────────────────────────────────────────────────────────
one_fit <- function(df_one,
                    sg_p = 2, sg_n = 11,          # smoothing
                    amp_min   = 0.05,             # min total amplitude
                    sos_range = c(90, 160),       # plausible DOY windows
                    eos_range = c(240, 330)) {
  
  # --- initialise return row -----------------------------------------------
  out <- tibble(plot = first(df_one$plot),
                year = year(first(df_one$date)),
                SOS_doy = NA_real_, EOS_doy = NA_real_, conv = FALSE)
  
  if (nrow(df_one) < 15) return(out)
  
  # --- 1  Smooth EVI2 -------------------------------------------------------
  df_one <- df_one %>%
    arrange(date) %>%
    mutate(doy  = yday(date),
           EVI2 = signal::sgolayfilt(EVI2, p = sg_p, n = sg_n)) %>%
    drop_na()
  
  if (diff(range(df_one$EVI2)) < amp_min) return(out)
  
  # --- 2  Fixed background & amplitudes ------------------------------------
  d_fix  <- min(df_one$EVI2)
  c1_fix <- max(df_one$EVI2) - d_fix        # spring amplitude  (+)
  c2_fix <- -c1_fix                         # autumn amplitude (−)
  
  # --- 3  Scaled time axis  (t ≈ −6…+6) ------------------------------------
  df_one <- df_one %>% mutate(t = (doy - 182.5) / 30)  # 30 d ≈ 1 t unit
  
  # --- 4  Grid search for robust start values ------------------------------
  grid <- expand.grid(
    a1 = seq(-6,  0, 1),                    # spring inflection
    b1 = seq(-0.20, -0.03, 0.03),
    a2 = seq( 1,  6, 1),                    # autumn inflection
    b2 = seq( 0.03, 0.20, 0.03))
  
  y <- df_one$EVI2 ; t <- df_one$t
  rss <- apply(grid, 1, \(g){
    g <- as.numeric(g)
    y_hat <- d_fix +
      c1_fix/(1+exp(g[1] + g[2]*t)) +
      c2_fix/(1+exp(g[3] + g[4]*t))
    sum((y - y_hat)^2)
  })
  start_best <- as.list(grid[which.min(rss), ])
  
  # --- 5  Fine fit with nlsLM ----------------------------------------------
  form <- EVI2 ~ d_fix +
    c1_fix/(1+exp(a1 + b1*t)) +
    c2_fix/(1+exp(a2 + b2*t))
  
  fit <- try(
    nlsLM(form, data = df_one,
          start = start_best,
          control = nls.lm.control(maxiter = 800)),
    silent = TRUE)
  if (inherits(fit,"try-error")) return(out)
  out$conv <- TRUE
  pa <- coef(fit)
  
  # --- 6  SOS / EOS via curvature ------------------------------------------
  t_grid <- seq(-6, 6, 0.1)
  kv_s   <- kv_part(t_grid, pa["a1"], pa["b1"], c1_fix)
  kv_a   <- kv_part(t_grid, pa["a2"], pa["b2"], c2_fix)
  sp_e   <- extrema(kv_s)
  au_e   <- extrema(kv_a)
  if (!length(sp_e) || !length(au_e)) return(out)
  
  # translate to DOY
  doy_sp <- t_grid[sp_e] * 30 + 182.5
  doy_au <- t_grid[au_e] * 30 + 182.5
  
  # ensure spring < autumn by swapping both sets if necessary
  if (mean(doy_au) < mean(doy_sp)) {
    tmp <- doy_sp; doy_sp <- doy_au; doy_au <- tmp
    tmp <- sp_e;   sp_e   <- au_e;   au_e   <- tmp
  }
  
  # ---- choose SOS (earliest) ----------------------------------------------
  in_sp <- which(doy_sp >= sos_range[1] & doy_sp <= sos_range[2])
  SOS_t <- if (length(in_sp)) {
    t_grid[sp_e[in_sp[1]]]
  } else {
    t_grid[sp_e[which.min(abs(doy_sp - mean(sos_range)))]]
  }
  
  # ---- choose EOS (latest) -------------------------------------------------
  in_au <- which(doy_au >= eos_range[1] & doy_au <= eos_range[2])
  EOS_t <- if (length(in_au)) {
    t_grid[au_e[tail(in_au, 1)]]
  } else {
    t_grid[au_e[which.min(abs(doy_au - mean(eos_range)))]]
  }
  
  # safety swap & 60-d rule
  SOS <- SOS_t*30 + 182.5
  EOS <- EOS_t*30 + 182.5
  if (EOS < SOS)  { tmp <- SOS; SOS <- EOS; EOS <- tmp }
  if (EOS >= SOS + 60) {
    out$SOS_doy <- round(SOS)
    out$EOS_doy <- round(EOS)
  }
  out
}

# ─────────────────────────────────────────────────────────────────────────────
#  5 • BATCH RUN  →  TABLE WITH SOS / EOS
# ─────────────────────────────────────────────────────────────────────────────
phenology_df <- evi_ts %>%
  mutate(year = year(date)) %>%
  group_by(plot, year) %>%
  group_split() %>%
  map_dfr(one_fit)

# ─────────────────────────────────────────────────────────────────────────────
#  6 • QUICK INSPECTION
# ─────────────────────────────────────────────────────────────────────────────
print(head(phenology_df, 10))
cat("converged:", sum(phenology_df$conv), "/", nrow(phenology_df), "\n")

# ─────────────────────────────────────────────────────────────────────────────
#  7 • OPTIONAL BOXPLOTS
# ─────────────────────────────────────────────────────────────────────────────
phen_long <- phenology_df %>%
  dplyr::filter(conv) %>%
  pivot_longer(c(SOS_doy, EOS_doy), names_to = "phase", values_to = "doy")

ggplot(phen_long, aes(phase, doy, fill = phase)) +
  geom_boxplot(width = .5, alpha = .6, outlier.shape = 21) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Distribution of SOS / EOS (all plots & years)",
       x = NULL, y = "DOY") +
  theme_minimal() +
  theme(legend.position = "none")


###############################################################################
#  EXTRA:  Jahres-Zeitreihen mit eingezeichnetem SOS & EOS
###############################################################################
#  Voraussetzung:  • evi_ts           – tägliche EVI2-Serie
#                  • phenology_df     – Ergebnis-Tabelle (plot, year, SOS_doy, EOS_doy, conv)
###############################################################################

walk(unique(year(evi_ts$date)), function(yr) {
  
  ## 1)  SOS/EOS für dieses Jahr herausziehen  ───────────────────────────────
  anno <- phenology_df %>%
    dplyr::filter(year == yr, conv) %>%                       # nur konvergierte
    mutate(
      SOS_date = as.Date(SOS_doy - 1, origin = paste0(yr, "-01-01")),
      EOS_date = as.Date(EOS_doy - 1, origin = paste0(yr, "-01-01"))
    ) %>%
    select(plot, SOS_date, EOS_date) %>%
    pivot_longer(cols = c(SOS_date, EOS_date),
                 names_to  = "phase",
                 values_to = "date")
  
  ## 2)  Zeitreihe + vertikale Linien plotten  ───────────────────────────────
  p <- evi_ts %>%
    dplyr::filter(year(date) == yr) %>%
    ggplot(aes(date, EVI2)) +
    geom_line(colour = "grey40") +
    facet_wrap(~ plot, scales = "free_y") +
    geom_vline(data = anno,
               aes(xintercept = date, colour = phase),
               linetype = "dashed", linewidth = .6, show.legend = FALSE) +
    scale_colour_manual(values = c(SOS_date = "blue", EOS_date = "darkgreen")) +
    labs(title = paste("EVI2-Zeitreihe mit SOS / EOS  –  Jahr", yr),
         x = NULL, y = "EVI2") +
    theme_minimal()
  
  print(p)   # zeigt Plot im Plot-Pane
})

write.csv(phenology_df, "data/analysis_ready_data/phenology_df.csv", row.names = FALSE)
