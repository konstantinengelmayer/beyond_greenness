###############################################################################
# Function: final_figures
# Author: Konstantin Engelmayer
# Date:   27.08.02025
# -----------------------------------------------------------------------------
# Figures produced (counting paired panels as ONE figure, as requested):
#
#  1) Dendrograms (paired):
#        • "VIs Clustered by TRI Correlation Pattern"
#        • "Plots Clustered by EVI ↔ TRI Correlation Pattern"
#     Files: figures/final_figures/dendrogram_vi.png
#            figures/final_figures/dendrogram_plot.png
#
#  2) kNDVI – TRI Median Rolling Correlation of the FASY Cluster
#     File:  figures/final_figures/kNDVI_FASY_Cluster.png
#
#  3) Top-10 VI–TRI correlations by VI
#     File:  figures/final_figures/VI_top10_species_medians.png
#
#  4) kNDVI - TRI Rolling Correlation — Plot SF01
#     File:  figures/final_figures/kNDVI_SF01.png
#
#  5) kNDVI & GVMI Rolling Correlation with TRI — Plot SF18 (FASY-Kellerwald)
#     File:  figures/final_figures/SF18_kNDVI_GVMI_grid.png
#
#  6) SIPI–TRI rolling correlation — plots SF12 & SF21
#     File:  figures/final_figures/SIPI_TRI_SF12_SF21_grid.png
#
#  7) GNDVI–TRI rolling correlation — plots SF04 & SF08
#     File:  figures/final_figures/GNDVI_TRI_SF04_SF08_grid.png
#
#  8) NIRV–TRI rolling correlation — plots SF19 & SF20
#     File:  figures/final_figures/NIRV_TRI_SF19_SF20_grid.png
#
#  9) Correlation heatmaps (paired):
#        • "Plot × Plot (kNDVI)"
#        • "VI × VI (all plots)"
#     Files: figures/final_figures/plot_plot_kndvi_cor.png
#            figures/final_figures/Vi_VI_cor.png
#
# 10) Batch export: VI–TRI rolling correlation — one figure per plot (faceted by VI)
#     Folder: figures/VI_TRI_heatmaps_by_plot
#
# Notes:
#  • Packages are loaded ONCE below.
#  • Data (corr_table + plot metadata) are loaded ONCE and reused everywhere.
#  • Plot titles are preserved exactly as in your original code.
###############################################################################

# ───────────────────────── 0) PACKAGES (load once) ───────────────────────────
# install.packages(c(
#   "dplyr","tidyr","ggplot2","dendextend","pals","readxl","RColorBrewer",
#   "scales","kableExtra","ggcorrplot","patchwork","ggrepel","lubridate"
# ))
library(dplyr)
library(tidyr)
library(ggplot2)
library(dendextend)
library(pals)
library(readxl)
library(RColorBrewer)
library(scales)
library(kableExtra)
library(ggcorrplot)
library(patchwork)
library(ggrepel)
library(lubridate)     # used in some sections
library(grid)          # for unit()

# ───────────────────────── 1) PATHS & DATA (load once) ───────────────────────
# Try multiple candidate paths for the correlation RDS, then read the first that exists.
rds_candidates <- c(
  "data/analysis_ready_data/corr_table_vi_tri.rds"
)
rds_path <- rds_candidates[which(file.exists(rds_candidates))[1]]
if (is.na(rds_path)) stop("Could not find corr_table_vi_tri.rds in expected locations.")

plot_meta_path <- "data/vector_data/plot_metadata.xlsx"
if (!file.exists(plot_meta_path)) stop("plot_metadata.xlsx not found at: ", plot_meta_path)

# Load once
heat_all  <- readRDS(rds_path)
plot_meta <- read_xlsx(plot_meta_path)

# Attach metadata once (species, area, etc.) and reuse downstream
heat_all <- left_join(heat_all, plot_meta, by = "plot")

# Ensure output dirs exist
dir.create("figures/final_figures", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/VI_TRI_heatmaps_by_plot", recursive = TRUE, showWarnings = FALSE)

# Global VI order used in several sections (keep original choices/case)
vi_order <- c("NDVI","GNDVI","CIG","EVI","NIRv","kNDVI","NMDI","NDWI","GVMI","SIPI","NPCI","GRVI")

# ───────────────────────── 2) HELPERS (reused) ───────────────────────────────
## 2a) Make a coloured dendrogram from a correlation matrix
make_dend <- function(corr_mat, h_cut, palette_fun = glasbey) {
  d    <- as.dist(1 - corr_mat)              # dissimilarity = 1 − r
  hc   <- hclust(d, method = "average")
  dend <- as.dendrogram(hc)
  
  k    <- length(unique(cutree(hc, h = h_cut)))
  cols <- palette_fun(k)
  
  dend <- color_branches(dend, h = h_cut, col = cols)
  labels_colors(dend) <- get_leaves_branches_col(dend)
  labels_cex(dend)    <- 0.9
  hang.dendrogram(dend, hang = 0.1)
}

## 2b) Right margin lines for dendrogram labels (dynamic)
right_lines <- function(dend) {
  w_in   <- max(strwidth(labels(dend), units = "inches"))
  char_w <- par("cin")[1]                     # width of one margin line (inches)
  ceiling(w_in / char_w) + 1
}

## 2c) Common contour binning (keeps your original bin logic)
make_bins <- function(x, step = 0.1, palette = "RdYlBu") {
  rng    <- range(x, na.rm = TRUE)
  breaks <- seq(floor(rng[1]/step)*step, ceiling(rng[2]/step)*step, by = step)
  nbins  <- length(breaks) - 1
  pal    <- colorRampPalette(rev(brewer.pal(11, palette)))(nbins)
  mids   <- (head(breaks,-1) + tail(breaks,-1)) / 2
  labs   <- sprintf("%.2f", mids)
  labs[labs == "-0.00"] <- "0.00"
  list(breaks = breaks, nbins = nbins, pal = pal, labels = labs)
}

## 2d) Grid lines (every 2 months on x, every 2 units on y)
make_grids <- function(date_vec, win_vec) {
  xg <- seq(min(date_vec), max(date_vec), by = "2 month")
  yg <- seq(0, ceiling(max(win_vec, na.rm = TRUE)), by = 2)
  list(x_grid = xg, y_grid = yg)
}

# ───────────────────────── 3) FIGURE 1 — DENDROGRAMS (paired) ────────────────
## 3a) Dendrogram 1: plots × plots (kNDVI only)
kndvi_mat <- heat_all %>%
  dplyr::filter(tolower(VI) == "kndvi") %>%
  mutate(time_id = paste0("w", window, "_d", doy, "_lag", lag)) %>%
  group_by(time_id, plot) %>%
  summarise(corr = mean(correlation, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = plot, values_from = corr, values_fill = NA) %>%
  select(-time_id) %>%
  cor(use = "pairwise.complete.obs")

dend_plot <- make_dend(kndvi_mat, h_cut = 0.20)

# Replace labels with "Plot – Species" (preserve branch colours)
labs       <- labels(dend_plot)
sp_vec     <- plot_meta$species[match(labs, plot_meta$plot)]
new_labels <- paste(labs, sp_vec, sep = " – ")

labels(dend_plot)        <- new_labels
labels_colors(dend_plot) <- get_leaves_branches_col(dend_plot)

## 3b) Dendrogram 2: VI × VI (all plots)
vi_mat <- heat_all %>%
  mutate(rec_id = paste(plot, "w", window, "d", doy, "lag", lag, sep = "_")) %>%
  group_by(rec_id, VI) %>%
  summarise(corr = mean(correlation, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = VI, values_from = corr, values_fill = NA) %>%
  select(-rec_id) %>%
  cor(use = "pairwise.complete.obs")

dend_vi <- make_dend(vi_mat, h_cut = 0.20)

## 3c) Export the paired dendrograms
fig_w  <- 2000  # px
fig_h  <- 1750  # px
fig_res <- 300  # dpi
op <- par(no.readonly = TRUE)

png("figures/final_figures/dendrogram_vi.png", width = fig_w, height = fig_h, res = fig_res)
par(mfrow = c(1,1), mar = c(4, 2, 2, right_lines(dend_vi)))
plot(dend_vi, horiz = TRUE,
     main = "VIs Clustered by TRI Correlation Pattern",
     xlab = "1 − Pearson r")
abline(v = 0.20, lty = 2)
dev.off()

png("figures/final_figures/dendrogram_plot.png", width = fig_w, height = fig_h, res = fig_res)
par(mfrow = c(1,1), mar = c(4, 2, 2, right_lines(dend_plot)))
plot(dend_plot, horiz = TRUE,
     main = "Plots Clustered by EVI ↔ TRI Correlation Pattern",
     xlab = "1 − Pearson r")
abline(v = 0.20, lty = 2)
dev.off()
par(op)

# ───────────────────────── 4) FIGURE 2 — FASY MEDIAN HEATMAP ─────────────────
vi_target <- "kNDVI"
plot_set  <- c("SF19","SF18","SF05","SF20","SF14","SF02","PF03","SF15","SF21")

df_med <- heat_all %>%
  dplyr::filter(plot %in% plot_set,
         tolower(VI) == tolower(vi_target),
         !is.na(correlation)) %>%
  group_by(window, doy, lag) %>%
  summarise(med_r = median(correlation), .groups = "drop") %>%
  mutate(
    ext_doy      = if_else(lag == 1, doy, doy + 365),
    month_date   = as.Date(ext_doy - 1, origin = "2000-01-01"),
    window_month = window * 8 / 30.437
  )

bins <- make_bins(df_med$med_r)
grds <- make_grids(df_med$month_date, df_med$window_month)

p1 <- ggplot(df_med, aes(month_date, window_month, z = med_r,
                         fill = after_stat(as.numeric(level)))) +
  geom_contour_filled(breaks = bins$breaks, colour = NA) +
  geom_contour(breaks = bins$breaks, colour = "grey70", linewidth = 0.15) +
  scale_x_date(date_breaks = "2 month",
               labels = date_format("%b", locale = "en"),
               expand = c(0,0)) +
  scale_y_continuous(breaks = grds$y_grid, expand = c(0,0)) +
  geom_vline(xintercept = grds$x_grid, colour = "grey50", linetype = "dotted") +
  geom_hline(yintercept = grds$y_grid, colour = "grey50", linetype = "dotted") +
  scale_fill_stepsn(
    colours = bins$pal, limits = c(1, bins$nbins), breaks = seq_len(bins$nbins),
    labels = bins$labels, name = expression("Median "~italic(r)),
    guide = guide_colourbar(barheight = unit(5,"cm"), barwidth = unit(0.5,"cm"),
                            ticks.colour = "black")
  ) +
  labs(title = paste(vi_target, "– TRI Median Rolling Correlation of the FASY Cluster"),
       x = "Month (previous → current year)",
       y = "Window length (months)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("figures/final_figures/kNDVI_FASY_Cluster.png", p1,
       width = 6, height = 4, dpi = 300, bg = "white")

# ───────────────────────── TABLE (not numbered as a figure) ──────────────────
# Median of top-10 |r| per plot × VI, shaded by significance counts
top10 <- heat_all %>%
  dplyr::filter(!is.na(correlation)) %>%
  group_by(plot, VI) %>%
  slice_max(abs(correlation), n = 10, with_ties = FALSE) %>%
  ungroup()

medians <- top10 %>%
  group_by(plot, VI, species) %>%
  summarise(median_r = median(abs(correlation)), .groups = "drop")

table_data <- medians %>%
  pivot_wider(names_from = VI, values_from = median_r) %>%
  select(plot, species, all_of(vi_order))

# order rows like the dendrogram (plots) and color plot/species labels
leaf_labels <- labels(dend_plot)
leaf_cols   <- labels_colors(dend_plot)
plot_ids    <- sub("^([^ ]+).*", "\\1", leaf_labels)
plot_col_vec <- setNames(leaf_cols, plot_ids)

table_data <- table_data %>%
  mutate(plot = factor(plot, levels = rev(plot_ids))) %>%
  arrange(plot)

sig_counts <- heat_all %>%
  dplyr::filter(p_value < 0.05) %>%
  group_by(plot, VI) %>%
  summarise(n_sig = n(), .groups = "drop") %>%
  pivot_wider(names_from = VI, values_from = n_sig) %>%
  select(plot, all_of(vi_order))

gp_red        <- "#F8766D"
all_vi_values <- unlist(table_data[vi_order], use.names = FALSE)
shade_fun     <- col_numeric(c("#FFFFFF", gp_red),
                             domain = range(all_vi_values, na.rm = TRUE))

for (vi in vi_order) {
  cnt <- sig_counts[[vi]][match(table_data$plot, sig_counts$plot)]
  cnt[is.na(cnt)] <- 0
  table_data[[vi]] <- mapply(function(val, n_sig){
    if (n_sig < 10){
      cell_spec(sprintf("%.2f", val), format = "latex",
                background = "#D9D9D9", color = "black")
    } else {
      cell_spec(sprintf("%.2f", val), format = "latex",
                background = shade_fun(val), color = "black")
    }
  }, table_data[[vi]], cnt, SIMPLIFY = FALSE)
}

table_data <- table_data %>%
  mutate(row_hex = plot_col_vec[as.character(plot)],
         plot    = cell_spec(as.character(plot),   format = "latex", color = row_hex),
         species = cell_spec(as.character(species),format = "latex", color = row_hex)) %>%
  select(-row_hex)

old_fmt <- getOption("knitr.table.format"); options(knitr.table.format = "latex")
kbl(table_data,
    format = "latex", booktabs = TRUE, escape = FALSE,
    caption = paste(
      "Median of the ten strongest VI–TRI correlations per plot.",
      "Cells with fewer than 10 significant correlations are greyed out;",
      "deeper", gp_red, "= stronger correlation."
    ),
    col.names = c("Plot","Species", vi_order)) %>%
  kable_styling(latex_options = "hold_position")
options(knitr.table.format = old_fmt)

# ───────────────────────── 5) FIGURE 3 — SPECIES LINES/BOXES ─────────────────
best_combo <- heat_all %>%
  dplyr::filter(!is.na(correlation)) %>%
  group_by(plot, VI) %>%
  slice_max(abs(correlation), n = 10, with_ties = FALSE) %>%
  ungroup()

med_tbl <- best_combo %>%
  group_by(VI) %>%
  summarise(median_r = median(abs(correlation)), .groups = "drop")

# ordered by descending median |r| (local order to this figure only)
vi_order_sp <- med_tbl %>% arrange((median_r)) %>% pull(VI)
best_combo  <- mutate(best_combo, VI = factor(VI, levels = vi_order_sp))
med_tbl     <- mutate(med_tbl,    VI = factor(VI, levels = vi_order_sp))

species_level <- best_combo %>%
  group_by(species, VI) %>%
  summarise(species_median = median(abs(correlation)), .groups = "drop")

species_counts <- best_combo %>%
  distinct(plot, species) %>% count(species, name = "n")

species_labs <- setNames(
  paste0(species_counts$species, " (n = ", species_counts$n, ")"),
  species_counts$species
)

p_species <- ggplot() +
  geom_boxplot(
    data    = best_combo,
    aes(VI, abs(correlation)),
    outlier.shape = NA, width = .6,
    colour = "grey40",  fill = "grey90"
  ) +
  geom_line(
    data = species_level,
    aes(VI, species_median, group = species, colour = species),
    linewidth = .7, alpha = .2
  ) +
  geom_point(
    data = species_level,
    aes(VI, species_median, colour = species),
    size = 2
  ) +
  geom_text(
    data = med_tbl,
    aes(VI, median_r, label = sprintf("%.2f", median_r)),
    vjust = -0.6, size = 3
  ) +
  scale_colour_brewer(palette = "Dark2", name = "Species", labels = species_labs) +
  labs(
    title    = "Top-10 VI–TRI correlations by VI",
    subtitle = "Lines connect each species’ median |r| across vegetation indices",
    x        = "Vegetation index (ordered by overall median |r|)",
    y        = "Pearson correlation (|r|)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold"),
    panel.grid.major.x = element_line(colour = "grey85"),
    panel.grid.major.y = element_blank(),
    legend.position    = "right"
  )

ggsave("figures/final_figures/VI_top10_species_medians.png", p_species,
       width = 11, height = 5, dpi = 300, bg = "white")

# ───────────────────────── 6) FIGURE 4 — SF01 kNDVI ──────────────────────────
sig_thr <- 0.05
vi_set  <- c("kNDVI")
plot_id <- "SF01"

df <- heat_all %>%
  dplyr::filter(plot == plot_id, VI %in% vi_set, !is.na(correlation)) %>%
  mutate(
    sig          = p_value < sig_thr,
    ext_doy      = if_else(lag == 1, doy, doy + 365),
    month_date   = as.Date(ext_doy - 1, origin = "2000-01-01"),
    window_month = window * 8 / 30.437,
    VI           = factor(VI, levels = vi_set),
    panel_lab    = paste(VI)
  )
df_nsig <- dplyr::filter(df, !sig)
df_sig  <- dplyr::filter(df,  sig)
best_pts <- df_sig %>%
  group_by(panel_lab) %>% slice_max(abs(correlation), n = 1, with_ties = FALSE) %>% ungroup()

bins <- make_bins(df$correlation); grds <- make_grids(df$month_date, df$window_month)

p <- ggplot() +
  geom_contour_filled(
    data = df,
    aes(month_date, window_month, z = correlation, fill = after_stat(as.numeric(level))),
    breaks = bins$breaks, colour = NA
  ) +
  geom_contour_filled(
    data = df_nsig,
    aes(month_date, window_month, z = correlation, fill = after_stat(as.numeric(level))),
    breaks = bins$breaks, fill = "white", alpha = 0.5, colour = NA
  ) +
  geom_contour(
    data = df, aes(month_date, window_month, z = correlation),
    colour = "grey70", size = 0.15, breaks = bins$breaks
  ) +
  geom_point(
    data = best_pts, aes(month_date, window_month),
    size = 2, shape = 21, stroke = .8, fill = "yellow", colour = "black"
  ) +
  geom_label_repel(
    data = best_pts, aes(month_date, window_month, label = sprintf("r = %.2f", correlation)),
    size = 3, fill = scales::alpha("white", 0.3), colour = "black", label.size = 0,
    label.r = unit(0.1, "lines"), segment.color = NA, box.padding = 0.3, point.padding = 0.2
  ) +
  scale_x_date(date_breaks = "2 month", labels = date_format("%b", locale = "en"), expand = c(0, 0)) +
  scale_y_continuous(breaks = grds$y_grid, expand = c(0,0)) +
  geom_vline(xintercept = grds$x_grid, colour = "grey50", linetype = "dotted") +
  geom_hline(yintercept = grds$y_grid, colour = "grey50", linetype = "dotted") +
  scale_fill_stepsn(
    colours = bins$pal, limits = c(1, bins$nbins), breaks = seq_len(bins$nbins),
    labels = bins$labels, name = expression("Pearson "~italic(r)),
    guide = guide_colourbar(barheight = unit(5,"cm"), barwidth = unit(0.5,"cm"),
                            ticks.colour = "black")
  ) +
  labs(
    title = "kNDVI - TRI Rolling Correlation — Plot SF01",
    x     = "Month (prev year → current year)",
    y     = "Window length (months)"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        strip.text = element_text(face = "bold"))

ggsave("figures/final_figures/kNDVI_SF01.png", p, width = 6, height = 4, dpi = 300, bg = "white")

# ───────────────────────── 7) FIGURE 5 — SF18 kNDVI & GVMI ───────────────────
sig_thr <- 0.05
vi_set  <- c("kNDVI", "GVMI")
plot_id <- "SF18"

df <- heat_all %>%
  dplyr::filter(plot == plot_id, VI %in% vi_set, !is.na(correlation)) %>%
  mutate(
    sig          = p_value < sig_thr,
    ext_doy      = if_else(lag == 1, doy, doy + 365),
    month_date   = as.Date(ext_doy - 1, origin = "2000-01-01"),
    window_month = window * 8 / 30.437,
    VI           = factor(VI, levels = vi_set),
    panel_lab    = paste(VI)
  )
df_nsig <- dplyr::filter(df, !sig); df_sig <- dplyr::filter(df, sig)
best_pts <- df_sig %>% group_by(panel_lab) %>% slice_max(abs(correlation), n = 1, with_ties = FALSE) %>% ungroup()
bins <- make_bins(df$correlation); grds <- make_grids(df$month_date, df$window_month)

p <- ggplot() +
  geom_contour_filled(
    data = df,
    aes(month_date, window_month, z = correlation, fill = after_stat(as.numeric(level))),
    breaks = bins$breaks, colour = NA
  ) +
  geom_contour_filled(
    data = df_nsig,
    aes(month_date, window_month, z = correlation, fill = after_stat(as.numeric(level))),
    breaks = bins$breaks, fill = "white", alpha = 0.5, colour = NA
  ) +
  geom_contour(data = df, aes(month_date, window_month, z = correlation),
               colour = "grey70", size = 0.15, breaks = bins$breaks) +
  geom_point(data = best_pts, aes(month_date, window_month),
             size = 2, shape = 21, stroke = .8, fill = "yellow", colour = "black") +
  geom_label_repel(
    data = best_pts, aes(month_date, window_month, label = sprintf("r = %.2f", correlation)),
    size = 3, fill = scales::alpha("white", 0.3), colour = "black", label.size = 0,
    label.r = unit(0.1, "lines"), box.padding = 0.3, point.padding = 0.2, segment.color = NA
  ) +
  facet_grid(~ VI) +
  scale_x_date(date_breaks = "2 month", labels = date_format("%b", locale = "en"), expand = c(0, 0)) +
  scale_y_continuous(breaks = grds$y_grid, expand = c(0,0)) +
  geom_vline(xintercept = grds$x_grid, colour = "grey50", linetype = "dotted") +
  geom_hline(yintercept = grds$y_grid, colour = "grey50", linetype = "dotted") +
  scale_fill_stepsn(
    colours = bins$pal, limits = c(1, bins$nbins), breaks = seq_len(bins$nbins),
    labels = bins$labels, name = expression("Pearson "~italic(r)),
    guide = guide_colourbar(barheight = unit(5,"cm"), barwidth = unit(0.5,"cm"),
                            ticks.colour = "black")
  ) +
  labs(
    title = "kNDVI & GVMI Rolling Correlation with TRI — Plot SF18 (FASY-Kellerwald)",
    x     = "Month (prev year → current year)",
    y     = "Window length (months)"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        strip.text = element_text(face = "bold"))

ggsave("figures/final_figures/SF18_kNDVI_GVMI_grid.png", p, width = 12, height = 5, dpi = 300, bg = "white")

# ───────────────────────── 8) FIGURE 6 — SIPI (SF12 & SF21) ──────────────────
sig_thr <- 0.05
df <- heat_all %>%
  dplyr::filter(plot %in% c("SF12","SF21"), tolower(VI) == "sipi", !is.na(correlation)) %>%
  mutate(
    sig          = p_value < sig_thr,
    ext_doy      = if_else(lag == 1, doy, doy + 365),
    month_date   = as.Date(ext_doy - 1, origin = "2000-01-01"),
    window_month = window * 8 / 30.437,
    plot_lab     = sprintf("%s – %s – %s", plot, species, area)
  )
df_nsig  <- dplyr::filter(df, !sig); df_sig <- dplyr::filter(df, sig)
best_pts <- df_sig %>% group_by(plot_lab) %>% slice_max(abs(correlation), n = 1, with_ties = FALSE) %>% ungroup()
bins <- make_bins(df$correlation); grds <- make_grids(df$month_date, df$window_month)

p <- ggplot() +
  geom_contour_filled(
    data = df,
    aes(month_date, window_month, z = correlation, fill = after_stat(as.numeric(level))),
    breaks = bins$breaks, colour = NA, linewidth = 0
  ) +
  geom_contour_filled(
    data = df_nsig,
    aes(month_date, window_month, z = correlation, fill = after_stat(as.numeric(level))),
    breaks = bins$breaks, fill = "white", alpha = 0.5, colour = NA, linewidth = 0
  ) +
  geom_contour(data = df, aes(month_date, window_month, z = correlation),
               colour = "grey70", size = 0.15, breaks = bins$breaks) +
  { if (nrow(best_pts) > 0)
    list(
      geom_point(data = best_pts, aes(month_date, window_month),
                 size = 2, shape = 21, stroke = .8, fill = "yellow", colour = "black"),
      geom_label_repel(
        data = best_pts, aes(month_date, window_month, label = sprintf("r = %.2f", correlation)),
        size = 3, fill = scales::alpha("white", 0.3), colour = "black",
        label.size = 0, label.r = unit(0.1, "lines"), box.padding = 0.3,
        point.padding = 0.2, segment.color = NA
      )
    )
  } +
  facet_wrap(~ plot_lab, ncol = 2) +
  scale_x_date(date_breaks = "2 month", labels = date_format("%b", locale = "en"), expand = c(0,0)) +
  scale_y_continuous(breaks = grds$y_grid, expand = c(0,0)) +
  geom_vline(xintercept = grds$x_grid, colour = "grey50", linetype = "dotted", linewidth = 0.3) +
  geom_hline(yintercept = grds$y_grid, colour = "grey50", linetype = "dotted", linewidth = 0.3) +
  scale_fill_stepsn(
    colours = bins$pal, limits = c(1, bins$nbins), breaks = seq_len(bins$nbins),
    labels = bins$labels, name = expression("Pearson "~italic(r)),
    guide = guide_colourbar(barheight = unit(5, "cm"), barwidth = unit(0.5, "cm"),
                            ticks.colour = "black")
  ) +
  labs(
    title = "SIPI–TRI rolling correlation — plots SF12 & SF21",
    x     = "Month (prev year → current year)",
    y     = "Window length (months)"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        strip.text = element_text(face = "bold"))

ggsave("figures/final_figures/SIPI_TRI_SF12_SF21_grid.png", p, width = 12, height = 5, dpi = 300, bg = "white")

# ───────────────────────── 9) FIGURE 7 — GNDVI (SF04 & SF08) ─────────────────
sig_thr <- 0.05
df <- heat_all %>%
  dplyr::filter(plot %in% c("SF04","SF08"), tolower(VI) == "gndvi", !is.na(correlation)) %>%
  mutate(
    sig          = p_value < sig_thr,
    ext_doy      = if_else(lag == 1, doy, doy + 365),
    month_date   = as.Date(ext_doy - 1, origin = "2000-01-01"),
    window_month = window * 8 / 30.437,
    plot_lab     = sprintf("%s – %s – %s", plot, species, area)
  )
df_nsig  <- dplyr::filter(df, !sig); df_sig <- dplyr::filter(df, sig)
best_pts <- df_sig %>% group_by(plot_lab) %>% slice_max(abs(correlation), n = 1, with_ties = FALSE) %>% ungroup()
bins <- make_bins(df$correlation); grds <- make_grids(df$month_date, df$window_month)

p <- ggplot() +
  geom_contour_filled(
    data = df,
    aes(month_date, window_month, z = correlation, fill = after_stat(as.numeric(level))),
    breaks = bins$breaks, colour = NA, linewidth = 0
  ) +
  geom_contour_filled(
    data = df_nsig,
    aes(month_date, window_month, z = correlation, fill = after_stat(as.numeric(level))),
    breaks = bins$breaks, fill = "white", alpha = 0.5, colour = NA, linewidth = 0
  ) +
  geom_contour(data = df, aes(month_date, window_month, z = correlation),
               colour = "grey70", size = 0.15, breaks = bins$breaks) +
  { if (nrow(best_pts) > 0)
    list(
      geom_point(data = best_pts, aes(month_date, window_month),
                 size = 2, shape = 21, stroke = .8, fill = "yellow", colour = "black"),
      geom_label_repel(
        data = best_pts, aes(month_date, window_month, label = sprintf("r = %.2f", correlation)),
        size = 3, fill = scales::alpha("white", 0.3), colour = "black",
        label.size = 0, label.r = unit(0.1, "lines"), box.padding = 0.3,
        point.padding = 0.2, segment.color = NA
      )
    )
  } +
  facet_wrap(~ plot_lab, ncol = 2) +
  scale_x_date(date_breaks = "2 month", labels = date_format("%b", locale = "en"), expand = c(0,0)) +
  scale_y_continuous(breaks = grds$y_grid, expand = c(0,0)) +
  geom_vline(xintercept = grds$x_grid, colour = "grey50", linetype = "dotted", linewidth = 0.3) +
  geom_hline(yintercept = grds$y_grid, colour = "grey50", linetype = "dotted", linewidth = 0.3) +
  scale_fill_stepsn(
    colours = bins$pal, limits = c(1, bins$nbins), breaks = seq_len(bins$nbins),
    labels = bins$labels, name = expression("Pearson "~italic(r)),
    guide = guide_colourbar(barheight = unit(5, "cm"), barwidth = unit(0.5, "cm"),
                            ticks.colour = "black")
  ) +
  labs(
    title = "GNDVI–TRI rolling correlation — plots SF04 & SF08",
    x     = "Month (prev year → current year)",
    y     = "Window length (months)"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        strip.text = element_text(face = "bold"))

ggsave("figures/final_figures/GNDVI_TRI_SF04_SF08_grid.png", p, width = 12, height = 5, dpi = 300, bg = "white")

# ───────────────────────── 10) FIGURE 8 — NIRV (SF19 & SF20) ─────────────────
sig_thr <- 0.05
df <- heat_all %>%
  dplyr::filter(plot %in% c("SF19","SF20"), tolower(VI) == "nirv", !is.na(correlation)) %>%
  mutate(
    sig          = p_value < sig_thr,
    ext_doy      = if_else(lag == 1, doy, doy + 365),
    month_date   = as.Date(ext_doy - 1, origin = "2000-01-01"),
    window_month = window * 8 / 30.437,
    plot_lab     = sprintf("%s – %s – %s", plot, species, area)
  )
df_nsig  <- dplyr::filter(df, !sig); df_sig <- dplyr::filter(df, sig)
best_pts <- df_sig %>% group_by(plot_lab) %>% slice_max(abs(correlation), n = 1, with_ties = FALSE) %>% ungroup()
bins <- make_bins(df$correlation); grds <- make_grids(df$month_date, df$window_month)

p <- ggplot() +
  geom_contour_filled(
    data = df,
    aes(month_date, window_month, z = correlation, fill = after_stat(as.numeric(level))),
    breaks = bins$breaks, colour = NA, linewidth = 0
  ) +
  geom_contour_filled(
    data = df_nsig,
    aes(month_date, window_month, z = correlation, fill = after_stat(as.numeric(level))),
    breaks = bins$breaks, fill = "white", alpha = 0.5, colour = NA, linewidth = 0
  ) +
  geom_contour(data = df, aes(month_date, window_month, z = correlation),
               colour = "grey70", size = 0.15, breaks = bins$breaks) +
  { if (nrow(best_pts) > 0)
    list(
      geom_point(data = best_pts, aes(month_date, window_month),
                 size = 2, shape = 21, stroke = .8, fill = "yellow", colour = "black"),
      geom_label_repel(
        data = best_pts, aes(month_date, window_month, label = sprintf("r = %.2f", correlation)),
        size = 3, fill = scales::alpha("white", 0.3), colour = "black",
        label.size = 0, label.r = unit(0.1, "lines"), box.padding = 0.3,
        point.padding = 0.2, segment.color = NA
      )
    )
  } +
  facet_wrap(~ plot_lab, ncol = 2) +
  scale_x_date(date_breaks = "2 month", labels = date_format("%b", locale = "en"), expand = c(0,0)) +
  scale_y_continuous(breaks = grds$y_grid, expand = c(0,0)) +
  geom_vline(xintercept = grds$x_grid, colour = "grey50", linetype = "dotted", linewidth = 0.3) +
  geom_hline(yintercept = grds$y_grid, colour = "grey50", linetype = "dotted", linewidth = 0.3) +
  scale_fill_stepsn(
    colours = bins$pal, limits = c(1, bins$nbins), breaks = seq_len(bins$nbins),
    labels = bins$labels, name = expression("Pearson "~italic(r)),
    guide = guide_colourbar(barheight = unit(5, "cm"), barwidth = unit(0.5, "cm"),
                            ticks.colour = "black")
  ) +
  labs(
    title = "NIRV–TRI rolling correlation — plots SF19 & SF20",
    x     = "Month (prev year → current year)",
    y     = "Window length (months)"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        strip.text = element_text(face = "bold"))

ggsave("figures/final_figures/NIRV_TRI_SF19_SF20_grid.png", p, width = 12, height = 5, dpi = 300, bg = "white")

# ───────────────────────── 11) FIGURE 9 — CORRELATION HEATMAPS (paired) ─────
corr_plot <- heat_all %>%
  dplyr::filter(tolower(VI) == "kndvi") %>%
  mutate(id = paste0("w", window, "_d", doy, "_lag", lag)) %>%
  group_by(id, plot) %>%
  summarise(corr = mean(correlation, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = plot, values_from = corr, values_fill = NA) %>%
  select(-id) %>%
  cor(use = "pairwise.complete.obs")

p1_cor <- ggcorrplot(corr_plot, type = "upper", lab = TRUE, title = "Plot × Plot (kNDVI)")

corr_vi <- heat_all %>%
  mutate(id = paste(plot, "w", window, "d", doy, "lag", lag, sep = "_")) %>%
  group_by(id, VI) %>%
  summarise(corr = mean(correlation, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = VI, values_from = corr, values_fill = NA) %>%
  select(-id) %>%
  cor(use = "pairwise.complete.obs")

p2_cor <- ggcorrplot(corr_vi, type = "upper", lab = TRUE, title = "VI × VI (all plots)")

ggsave("figures/final_figures/plot_plot_kndvi_cor.png", p1_cor, bg = "white", height = 11, width = 12)
ggsave("figures/final_figures/Vi_VI_cor.png",         p2_cor, bg = "white", height = 8,  width = 9)

# ───────────────────────── 12) FIGURE 10 — BATCH: ALL PLOTS (faceted by VI) ─
alpha_thr <- 0.05   # significance threshold for shading

for (pl in unique(heat_all$plot)) {
  
  df <- heat_all %>%
    dplyr::filter(plot == pl, !is.na(correlation)) %>%
    mutate(
      VI           = factor(VI, levels = vi_order),
      ext_doy      = if_else(lag == 1, doy, doy + 365),
      month_date   = as.Date(ext_doy - 1, origin = "2000-01-01"),
      window_month = window * 8 / 30.437,
      sig          = p_value < alpha_thr
    )
  if (nrow(df) == 0) next
  
  df_nonsig <- dplyr::filter(df, !sig)
  best_pts  <- df %>% dplyr::filter(sig) %>%
    group_by(VI) %>% slice_max(abs(correlation), n = 1, with_ties = FALSE) %>% ungroup()
  
  bins <- make_bins(df$correlation)
  grds <- make_grids(df$month_date, df$window_month)
  
  p <- ggplot() +
    geom_contour_filled(
      data = df,
      aes(month_date, window_month, z = correlation, fill = after_stat(as.numeric(level))),
      breaks = bins$breaks, colour = NA, linewidth = 0, na.rm = TRUE
    ) +
    geom_contour_filled(
      data = df_nonsig,
      aes(month_date, window_month, z = correlation, fill = after_stat(as.numeric(level))),
      breaks = bins$breaks, fill = "white", alpha = 0.5, colour = NA, linewidth = 0, na.rm = TRUE
    ) +
    geom_contour(data = df, aes(month_date, window_month, z = correlation),
                 colour = "grey70", size = 0.15, breaks = bins$breaks) +
    scale_fill_stepsn(
      colours = bins$pal, limits = c(1, bins$nbins), breaks = seq_len(bins$nbins),
      labels = bins$labels, name = expression("Pearson "~italic(r)),
      guide = guide_colourbar(barheight = unit(5, "cm"), barwidth = unit(0.5, "cm"),
                              ticks.colour = "black")
    ) +
    { if (nrow(best_pts) > 0)
      list(
        geom_point(data = best_pts, aes(month_date, window_month),
                   size = 2, shape = 21, stroke = .8, fill = "yellow", colour = "black"),
        geom_label_repel(
          data = best_pts, aes(month_date, window_month, label = sprintf("r = %.2f", correlation)),
          size = 3, fill = scales::alpha("white", 0.3), colour = "black",
          label.size = 0, label.r = unit(0.1, "lines"), box.padding = 0.3,
          point.padding = 0.2, segment.color = NA
        )
      )
    } +
    facet_wrap(~ VI, ncol = 4) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b", expand = c(0,0)) +
    scale_y_continuous(breaks = grds$y_grid, expand = c(0,0)) +
    geom_vline(xintercept = grds$x_grid, colour = "grey70", linetype = "dotted") +
    geom_hline(yintercept = grds$y_grid, colour = "grey70", linetype = "dotted") +
    labs(
      title = sprintf("VI–TRI rolling correlation — plot %s  (%s, %s)",
                      pl, unique(df$species), unique(df$area)),
      x     = "Month (prev year → current year)",
      y     = "Window length (months)"
    ) +
    theme_minimal(base_size = 11) +
    theme(strip.text = element_text(face = "bold"))
  
  ggsave(
    file.path("figures/VI_TRI_heatmaps_by_plot",
              sprintf("%s_%s_VI_TRI_contour_%s.png",
                      unique(df$species), unique(df$area), pl)),
    plot = p, width = 16, height = 8, dpi = 300, bg = "white"
  )
}

