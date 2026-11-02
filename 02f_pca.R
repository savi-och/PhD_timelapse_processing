# =============================================================================
# CELL TRACKING ANALYSIS - STEP 02f: PCA & QC DIAGNOSTICS
# =============================================================================
# Purpose:
#   1. Add inv_outreach_ratio (TrackLength / MaxDis) to turning_msd_final.csv
#      and save the result as tracks_clean_v1.csv (versioned master track file).
#   2. Run Principal Component Analysis (PCA) on the combined motility metric
#      set for two purposes:
#      a. QC/noise filtering: identify tracks with a noise-like multivariate
#         signature (low alpha_30pct, high ratio, low MaxDis, low speed,
#         low straightness).
#      b. Population clustering: assess whether biologically distinct motility
#         groups exist in multivariate space.
#
# PCA variable set (z-score standardised before PCA):
#   log_Speed, log_MaxDis, Straightness, mean_turning_angle, alpha_30pct,
#   inv_outreach_ratio, TrackDuration [Variant A: raw | Variant B: log-transformed]
#
# Two PCA variants are run:
#   A — raw TrackDuration
#   B — log_TrackDuration
# This lets you assess whether TrackDuration's right-skew changes the PC structure.
#
# Analysis sequence:
#   Section 1 — Configuration
#   Section 2 — Load data and compute inv_outreach_ratio
#   Section 3 — Save tracks_clean_v1.csv
#   Section 4 — PCA setup: variable selection and complete-case filtering
#   Section 5 — Run PCA (both variants)
#   Section 6 — Scree / variance explained plots
#   Section 7 — Loading plots (heatmap + bar chart)
#   Section 8 — Score biplots (4 grouping variables × 2 variants = 8 pages)
#   Section 9 — [Gated] Mahalanobis distance outlier flagging
#   Section 10 — Save all outputs
#
# Statistical references:
#   - PCA:  Jolliffe IT (2002) Principal Component Analysis, 2nd ed. Springer.
#   - Kaiser criterion for PC retention:
#             Kaiser HF (1960) The application of electronic computers to factor
#             analysis. Educational and Psychological Measurement 20:141-151.
#   - Mahalanobis distance for outlier detection:
#             Mahalanobis PC (1936) On the generalised distance in statistics.
#             Proc National Institute of Sciences of India 2:49-55.
#
# Inputs:
#   turning_msd_final.csv  — from Script 02e (59/77 ROIs matched to tracks_clean)
#
# Outputs:
#   tracks_clean_v1.csv    — turning_msd_final + inv_outreach_ratio + flags
#   pca_scores_v1.csv      — PC scores (PCA A) + metadata + QC flag (if enabled)
#   pca_loadings_v1.csv    — variable loadings per PC for both variants (stacked)
#   pca_diagnostics.pdf    — scree, loading, and biplot pages for both variants
#
# Note: turning_msd_final.csv is NOT overwritten.
#       tracks_clean.csv   is NOT touched.
# =============================================================================

rm(list = ls())
gc()


# --- 0. PACKAGES --------------------------------------------------------------

# install.packages(c("tidyverse", "patchwork"))
#install.packages("GGally") #if not already installed.
#install.packages("ggtext")
library(GGally)
library(ggtext)
library(tidyverse)
library(patchwork)


# =============================================================================
# ----SECTION 1 — CONFIGURATION----
# =============================================================================
# All user-adjustable settings are collected here so you never have to hunt
# through the rest of the script to change a threshold or file name.

setwd("Z:\\lightmicro-inst\\scharalambous\\Resident_cell_motility_tr_v3\\ROI files")

input_file  <- "tracks_clean_v3_alldata.csv" #"turning_msd_final.csv"


gmm_biplot_file <- "Z:\\lightmicro-inst\\scharalambous\\Resident_cell_motility_tr_v3\\ROI files\\Complete data -retrospective fitting G=5 VEE (from v2)\\gmm_component_assignments_VEE_G5_SpaceB_newROIs_alldata.csv"

# Pixel size and non-translocating threshold
# Tracks with MaxDis < 1.1 µm (= 2 pixels at 0.55 µm/px) are flagged.
# Their inv_outreach_ratio is set to NA and they are excluded from PCA
# (ratio is unreliable when the cell barely moved).
px_size_um       <- 0.55
maxdis_threshold <- 2 * px_size_um    # 1.1 µm

# PCA variable sets — differ only in TrackDuration transformation --- !!!CHANGE FOR MATRIX PLOT!!!

### ADD TRACK_LENGTH IF YOU WANT

pca_vars_A <- c(
  "log_Speed", "log_MaxDis", "log_TrackLength","Track_Straightness",
  "mean_turning_angle", "alpha_30pct",
  "log_inv_outreach_ratio", "TrackDuration"
)

pca_vars_B <- c(
  "log_Speed", "log_MaxDis","Track_Straightness",
  "mean_turning_angle", "alpha_30pct",
  "log_inv_outreach_ratio", "log_TrackDuration"
)

# Grouping variables used to colour biplots
colour_vars <- c("Treatment", "Day", "Location", "MouseID","GMM_component")

# Consistent colour palettes (same as rest of pipeline)
pal_treatment <- c("Mock" = "#457B9D", "RP"       = "#E63946")
pal_day       <- c("D1"   = "#2A9D8F", "D2"       = "#F4A261")
pal_location  <- c("Core" = "#264653", "CoreEdge" = "#E9C46A", "Edge" = "#E63946")
pal_mouse     <- c("M1"   = "#457B9D", "M2"       = "#E63946", "M3"   = "#F4A261","M4" = "#2A9D8F","M5" = "#6A0572" )
pal_gmm       <- c("1" = "#457B9D", "2" = "#E63946","3" = "#2A9D8F", "4" = "#F4A261", "5" = "#6A0572" )

palettes <- list(
  Treatment = pal_treatment,
  Day       = pal_day,
  Location  = pal_location,
  MouseID   = pal_mouse,
  GMM_component = pal_gmm
)

# Mahalanobis outlier flagging — gated flag
# Set to TRUE only after reviewing pca_diagnostics.pdf
flag_outliers  <- FALSE
mahal_pcs      <- 2      # number of PCs used for distance calculation
mahal_alpha    <- 0.001  # chi-squared tail probability (lower = stricter)
# PC axis pairs to plot in biplots.
# Each element is a length-2 integer vector: c(x_pc, y_pc).
# Add or remove pairs freely — all combinations will be plotted per grouping variable.
# Example: list(c(1,2), c(1,3), c(2,3)) gives three biplot pages per colour variable.
pc_pairs <- list(c(1, 2), c(1, 3), c(2, 3))

# =============================================================================
# ----SECTION 2 — LOAD DATA----
# =============================================================================
# Loads tracks_clean_v1.csv — the versioned master track file produced by
# Scripts 02e and 02f (previous run). The following columns are expected to
# already be present:
#   inv_outreach_ratio  — computed in 02f Section 2
#   Translocating       — flag: MaxDis >= 1.1 µm
#   log_TrackDuration   — log(TrackDuration)
#   log_Speed, log_MaxDis, log_TrackLength — log-transformed motility metrics
#   mean_turning_angle, alpha_30pct        — from Script 02e batch processing
#   GMM_component, GMM_max_prob            — from Script 02d (if flagging done)
#
# This section only applies factor encoding for grouping variables.
# No new columns are computed here.

tracks <- read_csv(input_file, show_col_types = FALSE) %>%
  mutate(
    log_inv_outreach_ratio = log(inv_outreach_ratio),
    log_TrackDuration  = log(TrackDuration),
    Location       = factor(Location,  levels = c("Core", "CoreEdge", "Edge")),
    Treatment      = factor(Treatment, levels = c("Mock", "RP")),
    Day            = factor(Day,       levels = c("D1", "D2")),
    MouseID        = factor(MouseID),
    TumourID_global  = factor(TumourID_global),
    SliceID_global   = factor(SliceID_global),
    ROIID_global     = factor(ROIID_global),
   # GMM_component    = factor(GMM_component)   # if already joined from 02d
  )

cat("Tracks loaded:", nrow(tracks), "\n")
cat("Columns available:", paste(names(tracks), collapse = ", "), "\n\n")


# =============================================================================
# ----SECTION 4 — PCA SETUP: VARIABLE SELECTION AND COMPLETE-CASE FILTERING----
# =============================================================================
# PCA cannot handle missing values — any track with NA on any PCA variable
# is excluded (complete-case analysis).
#
# Main reasons a track may be excluded:
#   - inv_outreach_ratio = NA  : MaxDis = 0 (cell never moved at all)
#   - alpha_30pct = NA         : track too short for MSD regression (< 3 lags)
#   - mean_turning_angle = NA  : single-step track (angle not computable)
#
# The helper function below selects only the PCA variables and ID columns,
# counts NAs per variable so you can see what is being dropped and why,
# and then removes incomplete rows.
#
# Non-translocating tracks (MaxDis < 1.1 µm) have NA for inv_outreach_ratio
# and will be excluded automatically — no separate filter needed.

id_cols <- c("ROIID_global", "ID", "MouseID", "Treatment", "Day", "Location",
             "TumourID_global", "SliceID_global", "Translocating")

make_pca_data <- function(df, vars, label) {
  cat("=== PCA VARIANT", label, "— complete-case filter ===\n")

  # Per-variable NA counts before filtering
  na_counts <- df %>% summarise(across(all_of(vars), ~ sum(is.na(.))))
  cat("NA count per PCA variable:\n")
  walk2(names(na_counts), as.numeric(na_counts[1, ]),
        ~ cat("  ", formatC(.x, width = 22, flag = "-"), ":", .y, "\n"))

  df_sel  <- df %>% select(all_of(c(id_cols, vars)))
  n_before <- nrow(df_sel)
  df_cc    <- df_sel %>% filter(complete.cases(across(all_of(vars))))
  n_after  <- nrow(df_cc)

  cat("Tracks before filter:", n_before, "\n")
  cat("Tracks after  filter:", n_after,  "\n")
  cat("Tracks excluded:     ", n_before - n_after, "\n\n")

  df_cc
}

pca_data_A <- make_pca_data(tracks, pca_vars_A, "A (raw TrackDuration)")
pca_data_B <- make_pca_data(tracks, pca_vars_B, "B (log TrackDuration)")


# Join GMM component assignments to PCA data for biplot colouring
#gmm_file <- "gmm_component_assignments_VEE_SpaceB.csv"
if (file.exists(gmm_biplot_file)) {
  gmm_comps <- read_csv(gmm_biplot_file, show_col_types = FALSE) %>%
    select(ROIID_global, ID, GMM_component, GMM_max_prob) %>%
    mutate(GMM_component = factor(GMM_component))
  pca_data_A <- pca_data_A %>%
    left_join(gmm_comps, by = c("ROIID_global", "ID"))
  pca_data_B <- pca_data_B %>%
    left_join(gmm_comps, by = c("ROIID_global", "ID"))
  cat("GMM joined:", sum(!is.na(pca_data_A$GMM_component)), "/",
      nrow(pca_data_A), "tracks\n")
} else {
  cat("GMM file not found — GMM_component biplot will be skipped\n")
}

# =============================================================================
# ----SECTION 5 — RUN PCA----
# =============================================================================
# prcomp() performs PCA with centering and scaling (z-score standardisation).
#
# Why z-score standardise?
#   The variables are on very different scales:
#     mean_turning_angle : 0–180 degrees
#     Straightness       : 0–1   (dimensionless)
#     alpha_30pct        : ~0–2  (dimensionless)
#     TrackDuration      : tens to hundreds of seconds
#     log_Speed, etc.    : log µm/min units
#   Without scaling, variables with larger numeric ranges would dominate the PCs
#   simply because of their units, not because they vary more biologically.
#   After z-scoring, each variable has mean = 0 and SD = 1 before PCA.
#
# What PCA produces:
#   fit$x        : score matrix — one row per track, one column per PC.
#                  The score tells you where that track sits on each PC axis.
#   fit$rotation : loading matrix — one row per variable, one column per PC.
#                  The loading tells you how much each original variable
#                  contributes to each PC (range −1 to +1).
#   fit$sdev     : standard deviation of each PC.
#                  sdev^2 = eigenvalue = variance explained by that PC.
#
# The first PC (PC1) captures the most variance in the data, PC2 the second
# most, and so on. PCs are orthogonal (uncorrelated with each other).
#
# For noise QC: a noise track signature (low speed, low MaxDis, high ratio,
# low alpha, low straightness) should project to a specific corner of PC space
# if these variables co-vary consistently. If they do, they will cluster there
# and can be separated from biological tracks by a PC threshold.

run_pca <- function(pca_data, vars, label) {
  mat <- pca_data %>% select(all_of(vars)) %>% as.matrix()
  fit <- prcomp(mat, center = TRUE, scale. = TRUE)

  var_exp <- fit$sdev^2 / sum(fit$sdev^2) * 100
  cumvar  <- cumsum(var_exp)

  cat("--- PCA", label, "---\n")
  cat("Variance explained per PC (%):\n")
  walk2(seq_along(var_exp), round(var_exp, 1),
        ~ cat("  PC", .x, ":", .y, "%\n"))
  cat("Cumulative variance:\n")
  walk2(seq_along(cumvar), round(cumvar, 1),
        ~ cat("  PC1 to PC", .x, ":", .y, "%\n"))
  cat("\n")

  list(fit = fit, var_exp = var_exp, cumvar = cumvar, label = label)
}

pca_A <- run_pca(pca_data_A, pca_vars_A, "A (raw TrackDuration)")
pca_B <- run_pca(pca_data_B, pca_vars_B, "B (log TrackDuration)")


# =============================================================================
# ----SECTION 6 — SCREE / VARIANCE EXPLAINED PLOTS----
# =============================================================================
# Two panels per PCA variant:
#
# Left panel — Scree plot (bar chart):
#   Each bar = % variance explained by that PC.
#   The dashed red line = Kaiser criterion: 100 / p % where p = number of
#   variables (here p = 7, so the line is at ~14.3%).
#   PCs above this line explain more variance than any single original variable
#   would on average, and are generally worth retaining for interpretation.
#   The "elbow" (where bars drop off sharply) is a visual guide to the
#   natural dimensionality of the data.
#   Reference: Kaiser HF (1960) Educational and Psychological Measurement 20:141.
#
# Right panel — Cumulative variance curve:
#   Shows how many PCs together explain 80% of total variance (dotted line).
#   This is a common practical threshold for "enough PCs to summarise the data".

make_scree <- function(pca_result) {
  p  <- length(pca_result$var_exp)
  df <- data.frame(
    PC     = factor(paste0("PC", seq_len(p)), levels = paste0("PC", seq_len(p))),
    VarExp = pca_result$var_exp,
    CumVar = pca_result$cumvar
  )
  kaiser_line <- 100 / p

  p_bar <- ggplot(df, aes(x = PC, y = VarExp)) +
    geom_col(fill = "steelblue4", alpha = 0.75, width = 0.6) +
    geom_hline(yintercept = kaiser_line,
               linetype = "dashed", colour = "firebrick", linewidth = 0.7) +
    annotate("text", x = Inf, y = kaiser_line,
             label = paste0("Kaiser criterion (", round(kaiser_line, 1), "%)"),
             hjust = 1.05, vjust = -0.5, size = 3, colour = "firebrick") +
    labs(title = paste("Scree plot — PCA", pca_result$label),
         subtitle = "Bars above dashed line explain more than average per-variable variance",
         x = NULL, y = "Variance explained (%)") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title    = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9, colour = "grey40"))

  p_cum <- ggplot(df, aes(x = PC, y = CumVar, group = 1)) +
    geom_line(colour = "steelblue4", linewidth = 0.9) +
    geom_point(colour = "steelblue4", size = 2.5) +
    geom_hline(yintercept = 80, linetype = "dotted",
               colour = "grey40", linewidth = 0.7) +
    annotate("text", x = Inf, y = 80,
             label = "80% threshold", hjust = 1.05, vjust = -0.5,
             size = 3, colour = "grey40") +
    ylim(0, 100) +
    labs(title = "Cumulative variance explained",
         subtitle = "Dotted line at 80%",
         x = NULL, y = "Cumulative variance (%)") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title    = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9, colour = "grey40"))

  p_bar | p_cum
}

scree_A <- make_scree(pca_A)
scree_B <- make_scree(pca_B)


# =============================================================================
# ----SECTION 7 — LOADING PLOTS----
# =============================================================================
# Loadings (stored in fit$rotation) describe how each original variable
# contributes to each PC.
#
#   Loading = +1 : the variable perfectly aligns with that PC (positive direction)
#   Loading = -1 : the variable perfectly aligns but is inverted
#   Loading =  0 : the variable does not contribute to that PC
#
# Two displays per PCA variant:
#
# Top — Heatmap across all PCs (up to PC4):
#   Shows the full loading matrix at a glance. Blue = positive, red = negative.
#   Useful for seeing which variables "own" which PCs.
#
# Bottom — Bar chart for PC1 and PC2 only:
#   More precise reading of the two most important PCs.
#   Blue bars = positive loading (variable increases → PC score increases).
#   Red bars  = negative loading (variable increases → PC score decreases).
#   The variable ordering (by loading magnitude) helps rank contributions.

make_loading_plot <- function(pca_result) {
  n_pcs <- min(4, ncol(pca_result$fit$rotation))

  load_df <- as.data.frame(pca_result$fit$rotation[, 1:n_pcs]) %>%
    rownames_to_column("Variable") %>%
    pivot_longer(-Variable, names_to = "PC", values_to = "Loading") %>%
    mutate(PC = factor(PC, levels = paste0("PC", 1:n_pcs)))

  p_heat <- ggplot(load_df, aes(x = PC, y = Variable, fill = Loading)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = round(Loading, 2)), size = 2.8) +
    scale_fill_gradient2(low = "#E63946", mid = "white", high = "#457B9D",
                         midpoint = 0, limits = c(-1, 1), name = "Loading") +
    labs(title    = paste("Variable loadings — PCA", pca_result$label),
         subtitle = "Blue = positive contribution, Red = negative contribution",
         x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(panel.grid    = element_blank(),
          plot.title    = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9, colour = "grey40"))

  p_bar <- load_df %>%
    filter(PC %in% c("PC1", "PC2")) %>%
    mutate(Direction = if_else(Loading > 0, "Positive", "Negative")) %>%
    ggplot(aes(x = reorder(Variable, Loading), y = Loading, fill = Direction)) +
    geom_col(alpha = 0.8, width = 0.6) +
    geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.4) +
    coord_flip() +
    facet_wrap(~ PC) +
    scale_fill_manual(values = c("Positive" = "#457B9D", "Negative" = "#E63946"),
                      guide = "none") +
    labs(title = "PC1 and PC2 loadings (bar chart)",
         x = NULL, y = "Loading") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 11))

  p_heat / p_bar + plot_layout(heights = c(1, 1.2))
}

load_A <- make_loading_plot(pca_A)
load_B <- make_loading_plot(pca_B)

# =============================================================================
# ----SECTION 7b — FEATURE PAIRS MATRIX----
# =============================================================================
# A pairs matrix shows every pairwise combination of PCA input variables in
# one grid, with univariate histograms on the diagonal.
#
# Layout (p variables → p×p grid):
#   Diagonal   : histogram of that variable's distribution
#   Off-diagonal: scatter of variable[row] vs variable[col], coloured by a
#                 grouping variable (default: MouseID — most informative for
#                 checking animal-level clustering before PCA)
#
# This is the multivariate equivalent of the univariate distribution checks
# already in Script 02c. It helps you see:
#   - Which pairs of variables are strongly correlated (linear scatter)
#   - Which are independent (circular/diffuse cloud)
#   - Whether any variable has extreme outliers pulling the scatter
#   - Whether MouseID (or Treatment) separates in raw feature space,
#     before PCA compresses everything into PCs
#
# Implementation uses GGally::ggpairs() which handles the grid construction.
# The colour_var_pairs setting below controls which grouping variable is used.
# Re-run with a different variable to get Treatment or Day version.
#


# Which grouping variable to colour the pairs matrix by
colour_var_pairs <- "GMM_component"   #####        ---->! change to "Treatment", "Day", or "Location", or "MouseID"

# Use the PCA Variant A complete-case data (raw TrackDuration version)
# Rename log_TrackDuration absent here — pca_data_A uses TrackDuration
pairs_plot_data <- pca_data_A %>%     # tracks
  select(all_of(c(pca_vars_A, colour_var_pairs))) %>%
  mutate(
    log_TrackDuration = log(TrackDuration),   # needed since it exists in pca_vars_B not tracks yet
    across(all_of(colour_var_pairs), as.factor)
  ) %>%
  filter(complete.cases(.))

# Variable display labels (same order as pca_vars_A)
var_labels <- c(
  log_Speed          = "Log Speed\n(µm/min)",
  log_MaxDis         = "Log_MaxDis\n(µm)",
  log_TrackLength    = "Log_Track Length\n(μm)",
  Track_Straightness = "Straightness\n(0-1)",
  mean_turning_angle = "Turning\nAngle (°)",
  alpha_30pct        = "MSD α\n(30%)",
  log_inv_outreach_ratio = "log_Inv \nOutreach Ratio",
  TrackDuration      = "Track\nDuration (s)"
)

# Active palette for chosen grouping variable
pairs_palette <- palettes[[colour_var_pairs]]

# Build coloured HTML legend string — must be outside the ggpairs chain
colour_legend <- paste(
  sapply(names(pairs_palette), function(g) {
    sprintf("<span style='color:%s'>**%s**</span>", pairs_palette[[g]], g)
  }),
  collapse = " &nbsp;|&nbsp; "
)

pairs_matrix <- ggpairs(
  pairs_plot_data,
  columns      = pca_vars_A,
  columnLabels = var_labels[pca_vars_A],
  aes(colour = .data[[colour_var_pairs]], alpha = 0.3),
  diag  = list(continuous = wrap("barDiag", bins = 40, alpha = 0.6)),
  lower = list(continuous = wrap("points",  size = 0.3, alpha = 0.2)),
  upper = list(continuous = wrap("cor",     size = 2.8))
) +
  labs(
    title    = paste("Feature pairs matrix — coloured by", colour_var_pairs),
    subtitle = paste0(
      "Lower: scatter &nbsp;|&nbsp; Diagonal: histogram &nbsp;|&nbsp; Upper: Pearson r",
      " &nbsp;|&nbsp; n = ", nrow(pairs_plot_data), " tracks<br>",
      colour_legend
    )                          # <- paste0 closes here
  ) +
  theme_bw(base_size = 9) +
  theme(
    plot.title       = element_text(face = "bold", size = 11),
    plot.subtitle    = element_markdown(size = 9, colour = "grey40"),
    strip.text       = element_text(size = 7),
    panel.grid.minor = element_blank()
  )
for (i in 1:pairs_matrix$nrow) {
  for (j in 1:pairs_matrix$ncol) {
    if (i != j) {   # skip diagonal histograms — fill handled separately
      pairs_matrix[i, j] <- pairs_matrix[i, j] +
        scale_colour_manual(values = pairs_palette)
    }
    if (i == j) {   # diagonal histograms use fill not colour
      pairs_matrix[i, j] <- pairs_matrix[i, j] +
        scale_fill_manual(values = pairs_palette)
    }
  }
}
  
print(pairs_matrix)

# =============================================================================
# ----SECTION 8 — SCORE BIPLOTS----
# =============================================================================
# A biplot places each track as a point in PC space (PC1 on x-axis, PC2 on y),
# and overlays arrows showing each variable's loading direction.
#
# How to read a biplot:
#   - Points: each dot is one track. Tracks that are similar across all
#     variables cluster together. Tracks in different corners of the plot
#     have different motility profiles.
#   - Arrows: each arrow shows which direction in PC space corresponds to
#     higher values of that variable. Arrow length reflects how strongly
#     the variable contributes to the two displayed PCs.
#   - Angle between arrows: small angle = variables are correlated;
#     opposite directions = anti-correlated; 90° = uncorrelated.
#
# Because there are thousands of tracks, points are made small and semi-
# transparent so you see density rather than individual points.
# Each page colours points by one grouping variable (Treatment, Day,
# Location, MouseID) — if a condition clusters in a distinct region,
# that PC axis may be capturing a biological or technical effect related
# to that variable.
#
# Arrow scaling:
#   Arrows are scaled so the longest arrow reaches ~40% of the score range.
#   This is a visual convenience only — arrows do not share the same axis
#   units as the points. Labels are placed at 112% of arrow length to avoid
#   overlapping the arrowhead.

make_biplot <- function(pca_result, pca_data, colour_var, palette,
                        pc_x = 1, pc_y = 2) {
  
  # Extract requested PC columns by index
  scores <- as.data.frame(pca_result$fit$x[, c(pc_x, pc_y)]) %>%
    setNames(c("PCx", "PCy")) %>%
    bind_cols(pca_data %>% select(all_of(c(colour_vars, "ROIID_global", "ID",
                                           "Translocating"))))
  
  loadings <- as.data.frame(pca_result$fit$rotation[, c(pc_x, pc_y)]) %>%
    setNames(c("PCx", "PCy")) %>%
    rownames_to_column("Variable")
  
  scale_arr <- min(
    diff(range(scores$PCx)) * 0.40,
    diff(range(scores$PCy)) * 0.40
  ) / max(abs(c(loadings$PCx, loadings$PCy)))
  
  var_exp_x <- round(pca_result$var_exp[pc_x], 1)
  var_exp_y <- round(pca_result$var_exp[pc_y], 1)
  
  ggplot() +
    geom_point(data = scores,
               aes(x = PCx, y = PCy, colour = .data[[colour_var]]),
               size = 1.5, alpha = 0.4) +
    geom_segment(data = loadings,
                 aes(x = 0, y = 0,
                     xend = PCx * scale_arr,
                     yend = PCy * scale_arr),
                 arrow     = arrow(length = unit(0.15, "cm"), type = "closed"),
                 colour    = "grey20",
                 linewidth = 0.5) +
    geom_text(data = loadings,
              aes(x     = PCx * scale_arr * 1.12,
                  y     = PCy * scale_arr * 1.12,
                  label = Variable),
              size = 2.8, colour = "grey20") +
    scale_colour_manual(values = palette, name = colour_var) +
    labs(
      title    = paste0("Biplot — PCA ", pca_result$label,
                        " | PC", pc_x, " vs PC", pc_y),
      subtitle = paste0("Coloured by: ", colour_var,
                        " | PC", pc_x, ": ", var_exp_x,
                        "% | PC", pc_y, ": ", var_exp_y, "%"),
      x = paste0("PC", pc_x, " (", var_exp_x, "%)"),
      y = paste0("PC", pc_y, " (", var_exp_y, "%)")
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title       = element_text(face = "bold", size = 11),
          plot.subtitle    = element_text(size = 9, colour = "grey40"),
          panel.grid.minor = element_blank(),
          legend.position  = "right") +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
}

# Build biplots: all PC pairs × all grouping variables × both PCA variants
# Stored as named lists: biplots_A[["PC1v2_Treatment"]] etc.
build_biplots <- function(pca_result, pca_data, pairs, cvars, pals) {
  out <- list()
  for (pair in pairs) {
    pair_label <- paste0("PC", pair[1], "v", pair[2])
    for (cv in cvars) {
      key <- paste0(pair_label, "_", cv)
      out[[key]] <- make_biplot(pca_result, pca_data, cv, pals[[cv]],
                                pc_x = pair[1], pc_y = pair[2])
    }
  }
  out
}


biplots_A <- build_biplots(pca_A, pca_data_A, pc_pairs, colour_vars, palettes)
biplots_B <- build_biplots(pca_B, pca_data_B, pc_pairs, colour_vars, palettes)

# =============================================================================
# ----SECTION 9 — [GATED] MAHALANOBIS DISTANCE OUTLIER FLAGGING----
# =============================================================================
# This section runs only when flag_outliers <- TRUE (see Section 1).
# Set it to TRUE only after reviewing the biplots in pca_diagnostics.pdf.
#
# What Mahalanobis distance measures:
#   Standard Euclidean distance from the centre of a cloud treats all axes
#   equally. Mahalanobis distance adjusts for the fact that PCs may have
#   different spreads — a track that is 3 SD out on a wide axis is less
#   extreme than one that is 3 SD out on a narrow axis.
#   Practically: it flags tracks that are genuinely unusual in multivariate
#   PC space, not just outliers on one variable.
#
# Threshold:
#   Under multivariate normality, Mahalanobis distance squared follows a
#   chi-squared distribution with df = mahal_pcs (number of PCs used).
#   A track whose distance exceeds the (1 - mahal_alpha) quantile of this
#   distribution is flagged as a multivariate outlier.
#   Default: alpha = 0.001 (0.1% tail) — conservative to avoid over-flagging.
#
# IMPORTANT: This is a geometric criterion. A flagged track is unusual in
# PC space — it is not automatically "noise". Always inspect flagged tracks'
# raw variable values before deciding to remove them.
#
# When enabled, outputs:
#   - Mahal_dist_A and QC_outlier_A columns appended to pca_data_A
#   - An additional biplot page highlighting outlier tracks in red

if (flag_outliers) {

  flag_mahal <- function(pca_result, pca_data, n_pcs, alpha, suffix) {

    scores_mat <- pca_result$fit$x[, 1:n_pcs, drop = FALSE]
    cov_mat    <- cov(scores_mat)
    centre     <- colMeans(scores_mat)
    dists      <- mahalanobis(scores_mat, center = centre, cov = cov_mat)
    threshold  <- qchisq(1 - alpha, df = n_pcs)

    cat("PCA", pca_result$label, "— Mahalanobis outlier flagging\n")
    cat("  Chi-sq threshold (df =", n_pcs, ", alpha =", alpha, "):",
        round(threshold, 3), "\n")
    cat("  Tracks flagged:", sum(dists > threshold), "/", length(dists), "\n\n")

    df_flag <- data.frame(
      Mahal_dist  = round(dists, 4),
      QC_outlier  = dists > threshold
    )
    names(df_flag) <- paste0(names(df_flag), "_", suffix)
    df_flag
  }

  cat("=== MAHALANOBIS OUTLIER FLAGGING ===\n")
  mahal_A <- flag_mahal(pca_A, pca_data_A, mahal_pcs, mahal_alpha, "A")
  mahal_B <- flag_mahal(pca_B, pca_data_B, mahal_pcs, mahal_alpha, "B")

  pca_data_A <- bind_cols(pca_data_A, mahal_A)
  pca_data_B <- bind_cols(pca_data_B, mahal_B)

  make_outlier_biplot <- function(pca_result, pca_data, outlier_col, label) {
    scores <- as.data.frame(pca_result$fit$x[, 1:2]) %>%
      setNames(c("PC1", "PC2")) %>%
      mutate(QC_outlier = pca_data[[outlier_col]])

    var_exp <- round(pca_result$var_exp[1:2], 1)

    ggplot(scores, aes(x = PC1, y = PC2)) +
      geom_point(data = scores %>% filter(!QC_outlier),
                 colour = "steelblue4", size = 0.4, alpha = 0.20) +
      geom_point(data = scores %>% filter(QC_outlier),
                 colour = "#E63946", size = 0.7, alpha = 0.65) +
      labs(title    = paste("Outlier overlay — PCA", label),
           subtitle = paste0("Red = Mahalanobis outlier | chi-sq threshold,",
                             " alpha = ", mahal_alpha,
                             ", df = ", mahal_pcs,
                             " | n flagged = ", sum(scores$QC_outlier)),
           x = paste0("PC1 (", var_exp[1], "%)"),
           y = paste0("PC2 (", var_exp[2], "%)")) +
      theme_bw(base_size = 11) +
      theme(plot.title       = element_text(face = "bold", size = 11),
            plot.subtitle    = element_text(size = 9, colour = "grey40"),
            panel.grid.minor = element_blank())
  }

  outlier_biplot_A <- make_outlier_biplot(
    pca_A, pca_data_A, "QC_outlier_A", pca_A$label)
  outlier_biplot_B <- make_outlier_biplot(
    pca_B, pca_data_B, "QC_outlier_B", pca_B$label)
}

# =============================================================================
# ----SECTION 10 — SAVE OUTPUTS----
# =============================================================================

# -- Plots ---------------------------------------------------------------------

pdf("pca_diagnostics.pdf", width = 12, height = 6)

suppressWarnings(print(pairs_matrix))

# PCA A
print(scree_A)
print(load_A)
invisible(lapply(biplots_A, print))
if (flag_outliers) print(outlier_biplot_A)

# PCA B
print(scree_B)
print(load_B)
invisible(lapply(biplots_B, print))
if (flag_outliers) print(outlier_biplot_B)

dev.off()
cat("Saved: pca_diagnostics.pdf\n")

# -- Loadings CSV --------------------------------------------------------------

loadings_A_df <- as.data.frame(pca_A$fit$rotation) %>%
  rownames_to_column("Variable") %>%
  mutate(PCA_variant = "A_raw_TrackDuration")

loadings_B_df <- as.data.frame(pca_B$fit$rotation) %>%
  rownames_to_column("Variable") %>%
  mutate(PCA_variant = "B_log_TrackDuration")

write_csv(bind_rows(loadings_A_df, loadings_B_df), "pca_loadings_v1.csv")
cat("Saved: pca_loadings_v1.csv\n")

# -- Scores CSV (PCA A) --------------------------------------------------------

scores_A_df <- as.data.frame(pca_A$fit$x) %>%
  rename_with(~ paste0(.x, "_A"))

pca_scores_out <- bind_cols(
  pca_data_A %>% select(ROIID_global, ID, MouseID, Treatment, Day, Location,
                        TumourID_global, SliceID_global, Translocating,
                        all_of(pca_vars_A)),
  scores_A_df
)

if (flag_outliers) {
  pca_scores_out <- pca_scores_out %>%
    bind_cols(pca_data_A %>% select(Mahal_dist_A, QC_outlier_A))
}

# Append GMM component if joined
if ("GMM_component" %in% names(pca_data_A)) {
  pca_scores_out <- pca_scores_out %>%
    left_join(
      pca_data_A %>% select(ROIID_global, ID, GMM_component, GMM_max_prob),
      by = c("ROIID_global", "ID")
    )
}

write_csv(pca_scores_out, "pca_scores_v1.csv")
cat("Saved: pca_scores_v1.csv\n")
cat("Tracks in scores file:", nrow(pca_scores_out), "\n\n")

# --- Final summary ------------------------------------------------------------

cat("=== COMPLETE ===\n\n")
cat("Files written:\n")
cat("   pca_diagnostics.pdf    — scree, loadings, biplots\n")
cat("   pca_loadings_v1.csv    — loading matrix (both PCA variants)\n")
cat("   pca_scores_v1.csv      — PC scores + metadata (PCA A)\n\n")

