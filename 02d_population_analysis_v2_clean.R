# =============================================================================
# CELL TRACKING ANALYSIS - STEP 2d: POPULATION IDENTIFICATION
# =============================================================================
# Purpose: Determine whether the track-level speed distribution supports
#          distinct cell subpopulations using a sequential approach.
#
# Analysis sequence:
#   Step 1 — Hartigan's dip test: formal test of unimodality on log_Speed
#   Step 2 — Univariate GMM (Section 4): BIC-driven component selection on
#             log_Speed; fits 1–G_max components, selects optimal G
#   Step 3 — Univariate component profiling (Section 5): do components differ
#             on other motility features?
#   Step 5 — Bivariate GMM (Section 6): noise detection in 2D feature spaces
#   Step 6 — 3D GMM (Section 6d): noise detection in 3D feature spaces
#             (primary analysis, three variables + mean_turning_angle)
#   Section 6e — Component profiling and export for 3D GMM
#   Section 6f — [OPTIONAL] Retrospective prediction for new ROI batches
#
# All analyses pooled across conditions; stratification deferred to Scripts 03/04.
#
# Statistical notes:
#   - Univariate GMM: mclust::Mclust(), model "V" (variable variance per
#     component). Allows each component its own mean and variance.
#   - 3D GMM: EEE preferred for biological interpretability.
#       EEE = equal volume, shape, orientation across components
#       VEE = variable volume, equal shape and orientation
#       VVV = fully variable (biologically implausible; rejected despite BIC)
#   - mclust BIC convention: HIGHER = better (opposite sign to standard BIC).
#   - Dip test null: distribution is unimodal. p < 0.05 rejects unimodality.
#     Reference: Hartigan JA & Hartigan PM (1985) Ann Stat 13(1):70-84.
#   - mclust reference: Scrucca L et al. (2016) mclust 5. R J 8(1):289-317.
#
# Match key throughout pipeline: ROIID_global + ID
#
# Input:  tracks_clean_v1.csv
#         Required columns: log_Speed, inv_outreach_ratio, mean_turning_angle,
#                           alpha_30pct, TrackDuration, ROIID_global, ID
#
# Outputs:
#   population_analysis_TrLen_MaxDis_ratio_vsMaxDis.pdf
#   population_analysis_bivariate_GMM.pdf
#   population_analysis_trivariate_GMM_SpaceAB_<model>_G<n>.pdf
#   component_profiles_<model>.pdf              [saved in Section 6e]
#   gmm_component_assignments_<model>_Space<A/B>.csv
#   gmm_3D_<A/B>_<model>_fit.rds               [saved in Section 6e]
#   population_gmm_summary.csv                  [univariate GMM, Section 9]
#   population_profile_summary.csv              [univariate GMM, Section 9]
# =============================================================================

rm(list = ls())
gc()

setwd("Z:\\lightmicro-inst\\scharalambous\\Resident_cell_motility_tr_v3\\ROI files")


# --- 0. PACKAGES --------------------------------------------------------------
# install.packages(c("diptest", "mclust", "tidyverse", "patchwork", "cowplot"))

library(tidyverse)
library(diptest)    # Hartigan's dip test
library(mclust)     # Gaussian mixture models
library(patchwork)  # plot assembly
library(cowplot)


# --- 1. CONFIGURATION ---------------------------------------------------------

# ---- Univariate GMM ----------------------------------------------------------
G_max                <- 5     # max components to evaluate for univariate GMM
assignment_threshold <- 0.80  # posterior probability threshold for hard assignment - ! Important for assignment of components

# ---- Input file --------------------------------------------------------------
input_file <- "tracks_clean_v3_alldata.csv"

# ---- Bivariate GMM (Section 6) -----------------------------------------------
flag_noise   <- FALSE   # set TRUE after reviewing bivariate GMM plots
G_max_biv    <- 5     # max components for bivariate and 3D GMM
                        # biologically motivated: noise / static / intermediate /
                        # migrating. Increase to 5 to test an additional subtype.
noise_comp_A <- NA      # component number = noise in Space A (review plots first)
noise_comp_B <- NA      # component number = noise in Space B (review plots first)

# ---- 3D GMM (Section 6d) -----------------------------------------------------
# gmm_model_types: evaluated when force_model = NULL (BIC selects best combination).
#   All combinations of model × G = 1:G_max_biv are evaluated; highest BIC wins.
#
#   To test additional model geometries, extend this vector. Common candidates:
#     "EEE" — equal volume, shape, orientation (current default; biologically preferred)
#     "VEE" — variable volume, equal shape/orientation
#     "VVV" — fully variable (biologically implausible; rejected in initial analysis)
#     "EEV" — equal volume/shape, variable orientation
#     "VEV" — variable volume/orientation, equal shape
#     "EVE" — equal volume/orientation, variable shape
#     "EVI" — equal volume, variable shape, axes aligned with coordinates
#     "VVI" — variable volume/shape, axes aligned (diagonal covariance)
#     "VVE" — variable volume/shape, equal orientation
#   See mclust documentation (?mclustModelNames) for the full list.
#
#   To test G = 3–5: set G_max_biv <- 5 below.
#   To force a single model for a dedicated run: set force_model <- "EEE" (or other).
#   Output filenames embed the model label — run once per model to keep traceability.
gmm_model_types <- c("VEE")

gmm_seed        <- 42     # fixed seed for reproducibility
flag_noise_3D   <- FALSE  # set TRUE after reviewing 3D GMM plots
noise_comp_3D_A <- NA     # noise component in 3D Space A (review plots first)
noise_comp_3D_B <- NA     # noise component in 3D Space B (review plots first)

# force_model: override BIC model selection.
#   NULL   — mclust evaluates all models in gmm_model_types across G = 1:G_max_biv
#             and selects the highest BIC combination automatically.
#             BIC plot will show all model lines for comparison.
#   "EEE"  — forces EEE model (current default; biologically preferred).
#   "VEE"  — forces VEE model.
#   To test multiple models: change force_model, re-run Sections 6d–6e each time.
#   Output filenames include the model label for traceability.
force_model <- "VEE" #! IMPORTANT TO DEFINE MODEL

# ---- Retrospective prediction (Section 6f) -----------------------------------
# READ ME!!
# Set run_retro_predict <- TRUE once new ROI tracks are available.
# new_tracks_file : path to CSV with new ROI tracks (same format as input_file).
# retro_rds_B     : full path (or bare filename if in working directory) to the
#                   Space B .rds model file saved by Section 6e.
#                   Model name and G are read directly from the file — no need
#                   to specify them separately.
# append_to_existing : if TRUE, also writes a combined CSV (original + new).
#   Inspect new-ROI output first; set append_to_existing <- TRUE only once
#   satisfied, to produce the combined file for downstream scripts.
run_retro_predict  <- TRUE
new_tracks_file    <- "Z:\\lightmicro-inst\\scharalambous\\Resident_cell_motility_tr_v3\\ROI files\\tracks_clean_v3_alldata.csv"
retro_rds_B        <- "Z:\\lightmicro-inst\\scharalambous\\Resident_cell_motility_tr_v2\\ROI files\\GMM modelling_noise\\G=5\\VEE\\gmm_3D_B_VEE_fit.rds"   # update path/filename as needed
retro_vars_B <- c("log_inv_outreach_ratio", "log_MaxDis", "mean_turning_angle")
append_to_existing <- FALSE

# ---- GMM component join (Section 10) ----------------------------------------- !!!!
# Set these to the assignment CSV(s) you want joined onto tracks_clean.
# Use the _combined.csv if new ROIs have been retrospectively classified,
# otherwise use the standard Space A/B output from Section 6e.
# Set to NA to skip that space.
gmm_join_file <- "Z:\\lightmicro-inst\\scharalambous\\Resident_cell_motility_tr_v3\\ROI files\\GMM modelling - complete dataset - new fitting G=5 VEE\\gmm_component_assignments_VEE_SpaceB.csv"

  #"Z:\\lightmicro-inst\\scharalambous\\Resident_cell_motility_tr_v2\\ROI files\\GMM modelling_noise\\G=5\\06JUN26_Outputs from use of VEE G=5 Space B on 60ROIs\\gmm_component_assignments_VEE_G5_SpaceB_newROIs.csv"



# READ ME : CHANGE LABELS IN 6f TO REFLECT THE SPACE YOU USED FOR FITTING - SHOULD BE OKAY AS SPACE BE WAS DECIDED
# BUT JUST IN CASE



# ----  Pipeline colour palettes (consistent across all scripts in this project) ----
comp_colours  <- c("1" = "#457B9D", "2" = "#E63946",
                   "3" = "#2A9D8F", "4" = "#F4A261", "5" = "#6A0572","6"="black","7"="grey","8"="blue")
pal_treatment <- c("Mock" = "#457B9D", "RP"  = "#E63946")
pal_day       <- c("D1"   = "#2A9D8F", "D2"  = "#F4A261")
pal_mouse     <- c("M1"   = "#457B9D", "M2"  = "#E63946", "M3" = "#F4A261")


# --- 2. LOAD DATA -------------------------------------------------------------

tracks <- read_csv(input_file, show_col_types = FALSE) %>%
  mutate(
    Location        = factor(Location,  levels = c("Core", "CoreEdge", "Edge")),
    Treatment       = factor(Treatment, levels = c("Mock", "RP")),
    Day             = factor(Day,       levels = c("D1", "D2")),
    MouseID         = factor(MouseID),
    TumourID_global = factor(TumourID_global),
    SliceID_global  = factor(SliceID_global),
    ROIID_global    = factor(ROIID_global)
  )

cat("Tracks loaded:", nrow(tracks), "tracks\n")
cat("Columns available:", paste(names(tracks), collapse = ", "), "\n\n")

speed_vec  <- tracks$log_Speed
n_complete <- sum(!is.na(speed_vec))
speed_vec  <- speed_vec[!is.na(speed_vec)]

cat("Tracks with non-missing log_Speed:", n_complete, "\n")
cat("Speed range (log scale): [", round(min(speed_vec), 3), ",",
    round(max(speed_vec), 3), "]\n")
cat("Speed range (µm/min):    [", round(exp(min(speed_vec)), 2), ",",
    round(exp(max(speed_vec)), 2), "]\n\n")


# --- 3. STEP 1: HARTIGAN'S DIP TEST -------------------------------------------
# Tests H0: distribution is unimodal. p < 0.05 rejects unimodality.
# D < 0.01 with large n: statistically detectable but biologically negligible.
# Reference: Hartigan JA & Hartigan PM (1985) Ann Stat 13(1):70-84.

cat("=== STEP 1: HARTIGAN'S DIP TEST ===\n")

dip_result <- dip.test(speed_vec, simulate.p.value = FALSE)

cat("Dip statistic D =", round(dip_result$statistic, 5), "\n")
cat("p-value          =", round(dip_result$p.value,   5), "\n")
cat("n tracks         =", length(speed_vec), "\n")
cat("Interpretation: D < 0.01 = negligible | 0.01-0.05 = modest | > 0.05 = clear\n")
cat("Observed D:", round(dip_result$statistic, 5), "->",
    ifelse(dip_result$statistic < 0.01, "negligible",
           ifelse(dip_result$statistic < 0.05, "modest", "clear")),
    "departure from unimodality\n")

dip_label <- paste0(
  "Dip test: D = ", round(dip_result$statistic, 4),
  ", p = ", ifelse(dip_result$p.value < 0.0001, "< 0.0001",
                   round(dip_result$p.value, 4)),
  "  (n = ", length(speed_vec), ")"
)


# --- 4. STEP 2: UNIVARIATE GAUSSIAN MIXTURE MODEL ----------------------------
# Model "V": variable variance per component (each component has own mean + SD).
# mclust BIC: HIGHER = better (opposite sign convention to standard BIC).
# Reference: Scrucca L et al. (2016) mclust 5. R J 8(1):289-317.

cat("\n=== STEP 2: UNIVARIATE GMM (log_Speed) ===\n")
cat("Fitting 1 to", G_max, "components (model: V)\n")

gmm_fit   <- Mclust(speed_vec, G = 1:G_max, modelNames = "V", verbose = FALSE)
G_optimal <- gmm_fit$G
cat("Optimal components (BIC):", G_optimal, "\n\n")

bic_df <- data.frame(
  G   = 1:G_max,
  BIC = as.numeric(gmm_fit$BIC[, "V"])
) %>%
  mutate(
    Delta_BIC = BIC - max(BIC, na.rm = TRUE),
    Selected  = G == G_optimal
  )

cat("BIC table (higher = better):\n")
print(bic_df, digits = 4)

comp_params <- data.frame(
  Component     = seq_len(G_optimal),
  Mean_logSpeed = round(gmm_fit$parameters$mean, 4),
  SD_logSpeed   = round(sqrt(gmm_fit$parameters$variance$sigmasq), 4),
  Proportion    = round(gmm_fit$parameters$pro, 4)
) %>%
  mutate(
    Mean_Speed_ummin = round(exp(Mean_logSpeed), 2),
    GSD_Speed        = round(exp(SD_logSpeed), 3)
  ) %>%
  arrange(Mean_logSpeed)

cat("\nGMM component parameters (slowest to fastest):\n")
print(comp_params)
cat("Note: Mean_Speed_ummin = geometric mean; GSD_Speed = geometric SD factor.\n")

prob_matrix               <- as.data.frame(gmm_fit$z)
names(prob_matrix)        <- paste0("P_comp", seq_len(G_optimal))
prob_matrix$GMM_component <- apply(gmm_fit$z, 1, which.max)
prob_matrix$GMM_max_prob  <- apply(gmm_fit$z, 1, max)

tracks_complete <- tracks %>% filter(!is.na(log_Speed))
stopifnot(nrow(tracks_complete) == nrow(prob_matrix))
tracks_complete <- bind_cols(tracks_complete, prob_matrix)

cat("\nComponent membership counts:\n")
print(table(tracks_complete$GMM_component))
cat("Tracks with max_prob >=", assignment_threshold, ":",
    sum(tracks_complete$GMM_max_prob >= assignment_threshold), "\n")
cat("Tracks with max_prob < ", assignment_threshold, "(uncertain):",
    sum(tracks_complete$GMM_max_prob <  assignment_threshold), "\n")


# --- 5. STEP 3: UNIVARIATE COMPONENT PROFILING --------------------------------
# Descriptive only: mean ± SD per univariate GMM component on other features.
# Informs whether speed alone distinguishes populations or multivariate
# analysis (PCA / 3D GMM) is warranted.

cat("\n=== STEP 3: UNIVARIATE COMPONENT PROFILING ===\n")

profile_summary <- tracks_complete %>%
  group_by(GMM_component) %>%
  summarise(
    n_tracks          = n(),
    Speed_geomean     = round(exp(mean(log_Speed,  na.rm = TRUE)), 3),
    Speed_GSD         = round(exp(sd(log_Speed,    na.rm = TRUE)), 3),
    Straightness_mean = round(mean(Straightness,   na.rm = TRUE), 4),
    Straightness_sd   = round(sd(Straightness,     na.rm = TRUE), 4),
    MaxDis_mean       = round(mean(MaxDis,         na.rm = TRUE), 2),
    MaxDis_sd         = round(sd(MaxDis,           na.rm = TRUE), 2),
    TrackLength_mean  = round(mean(TrackLength,    na.rm = TRUE), 2),
    TrackLength_sd    = round(sd(TrackLength,      na.rm = TRUE), 2),
    .groups = "drop"
  )

cat("Profile summary per univariate GMM component:\n")
print(profile_summary, n = Inf)


# --- 5b. PAIRWISE FEATURE SCATTER PLOTS ---------------------------------------
# Six orthogonal pairwise scatter plots across:
#   log_Speed, log_TrackLength, log_MaxDis, Straightness
#
# Multiple colouring options are defined below (Options A–E).
# Only ONE definition of make_scatter should be active at a time.
# To switch colouring: comment out the active block and uncomment the desired one.
# Pipeline colour palettes — defined in Section 1.

# Scatter data: compute log transforms not already in tracks_complete
scatter_data <- tracks_complete %>%
  mutate(
    log_MaxDis      = log(MaxDis      + 1e-6),
    log_TrackLength = log(TrackLength + 1e-6),
    Component       = factor(GMM_component)
  ) %>%
  filter(
    !is.na(log_Speed), !is.na(log_MaxDis),
    !is.na(log_TrackLength), !is.na(Straightness)
  )

# [OPTION A — coloured by GMM component]
# make_scatter <- function(xvar, yvar, xlabel, ylabel) {
#   ggplot(scatter_data, aes(x = .data[[xvar]], y = .data[[yvar]],
#                            colour = Component)) +
#     geom_point(size = 0.3, alpha = 0.3) +
#     geom_density2d(linewidth = 0.5, alpha = 0.8, adjust = 1.5) +
#     scale_colour_manual(values = comp_colours, name = "Component") +
#     labs(x = xlabel, y = ylabel) +
#     theme_bw(base_size = 10) +
#     theme(legend.position = "right", panel.grid.minor = element_blank(),
#           axis.title = element_text(size = 9))
# }

# [OPTION B — no colour, points only]
# make_scatter <- function(xvar, yvar, xlabel, ylabel) {
#   ggplot(scatter_data, aes(x = .data[[xvar]], y = .data[[yvar]])) +
#     geom_point(size = 0.3, alpha = 0.12, colour = "steelblue4") +
#     labs(x = xlabel, y = ylabel) +
#     theme_bw(base_size = 10) +
#     theme(panel.grid.minor = element_blank(),
#           axis.title = element_text(size = 9))
# }

# [OPTION C — coloured by Treatment]
# make_scatter <- function(xvar, yvar, xlabel, ylabel) {
#   ggplot(scatter_data, aes(x = .data[[xvar]], y = .data[[yvar]],
#                            colour = Treatment)) +
#     geom_point(size = 0.3, alpha = 0.3) +
#     geom_density2d(linewidth = 0.5, alpha = 0.8, adjust = 1.5) +
#     scale_colour_manual(values = pal_treatment, name = "Treatment") +
#     labs(x = xlabel, y = ylabel) +
#     theme_bw(base_size = 10) +
#     theme(legend.position = "right", panel.grid.minor = element_blank(),
#           axis.title = element_text(size = 9))
# }

# [OPTION D — coloured by Day]
# make_scatter <- function(xvar, yvar, xlabel, ylabel) {
#   ggplot(scatter_data, aes(x = .data[[xvar]], y = .data[[yvar]],
#                            colour = Day)) +
#     geom_point(size = 0.3, alpha = 0.3) +
#     geom_density2d(linewidth = 0.5, alpha = 0.8, adjust = 1.5) +
#     scale_colour_manual(values = pal_day, name = "Day") +
#     labs(x = xlabel, y = ylabel) +
#     theme_bw(base_size = 10) +
#     theme(legend.position = "right", panel.grid.minor = element_blank(),
#           axis.title = element_text(size = 9))
# }

# [OPTION E — coloured by Mouse] *** ACTIVE ***
make_scatter <- function(xvar, yvar, xlabel, ylabel) {
  ggplot(scatter_data, aes(x = .data[[xvar]], y = .data[[yvar]],
                           colour = MouseID)) +
    geom_point(size = 0.3, alpha = 0.3) +
    geom_density2d(linewidth = 0.5, alpha = 0.8, adjust = 1.5) +
    scale_colour_manual(values = pal_mouse, name = "Mouse") +
    labs(x = xlabel, y = ylabel) +
    theme_bw(base_size = 10) +
    theme(legend.position  = "right",
          panel.grid.minor = element_blank(),
          axis.title       = element_text(size = 9))
}

p_s_l  <- make_scatter("log_Speed",       "log_TrackLength",
                       "Log Speed (µm/min)",    "Log Track Length (µm)")
p_s_d  <- make_scatter("log_Speed",       "log_MaxDis",
                       "Log Speed (µm/min)",    "Log Max Displacement (µm)")
p_s_st <- make_scatter("log_Speed",       "Straightness",
                       "Log Speed (µm/min)",    "Straightness (0-1)")
p_l_d  <- make_scatter("log_TrackLength", "log_MaxDis",
                       "Log Track Length (µm)", "Log Max Displacement (µm)")
p_l_st <- make_scatter("log_TrackLength", "Straightness",
                       "Log Track Length (µm)", "Straightness (0-1)")
p_d_st <- make_scatter("log_MaxDis",      "Straightness",
                       "Log Max Displacement (µm)", "Straightness (0-1)")


# --- 5c. TRACKLENGTH/MAXDIS RATIO PLOTS ---------------------------------------
# Ratio = TrackLength / MaxDis: 1 = linear path, > 1 = tortuous path.
# MaxDis filter: removes tracks with MaxDis < 1.1 µm (< 2 pixels at 0.55 µm/px).
# Two plot pages: raw MaxDis x-axis (p6) and log MaxDis x-axis (p7).

maxdis_threshold <- 1.1   # µm — 2 pixels at 0.55 µm/px

ratio_data_all <- tracks_complete %>%
  filter(!is.na(MaxDis), !is.na(TrackLength), !is.na(SpeedMean_ummin),
         MaxDis > 0, TrackLength > 0) %>%
  mutate(
    Ratio            = TrackLength / MaxDis,
    Translocating    = MaxDis >= maxdis_threshold,
    Inv_Straightness = 1 / Straightness,
    log_duration     = log(TrackDuration),
    log_MaxDis       = log(MaxDis + 1e-6)   # computed here for p7 x-axis
  )

ratio_data_filtered <- ratio_data_all %>% filter(Translocating)

cat("\n--- Ratio data summary ---\n")
cat("All tracks (valid MaxDis):", nrow(ratio_data_all), "\n")
cat("Non-translocating removed (MaxDis <", maxdis_threshold, "µm):",
    sum(!ratio_data_all$Translocating), "\n")
cat("Retained for filtered plot:", nrow(ratio_data_filtered), "\n")
cat("Ratio range (all):     [", round(min(ratio_data_all$Ratio), 2),
    ",", round(max(ratio_data_all$Ratio), 2), "]\n")
cat("Ratio range (filtered):[", round(min(ratio_data_filtered$Ratio), 2),
    ",", round(max(ratio_data_filtered$Ratio), 2), "]\n")

make_ratio_panel <- function(data, xvar, xlabel, title, xlim = NULL, ylim = NULL) {
  
  p <- ggplot(data, aes(x = .data[[xvar]], y = Ratio,
                        colour = factor(MouseID))) +
    geom_point(size = 1.5, alpha = 0.5) +
    # Note: mouse colours here are intentionally hardcoded to preserve original
    # plot appearance. To align with pipeline pal_mouse, replace with:
    #   scale_colour_manual(values = pal_mouse, name = "Mouse")
    scale_colour_manual(values = c("M1" = "#E63946",
                                   "M2" = "#457B9D",
                                   "M3" = "#FFAA00"),
                        name = "Mouse") +
    # [OPTION: filled density contours — comment out geom_point block above,
    #  then uncomment below to switch from scatter to density display]
    # geom_density2d_filled(alpha = 0.85, contour_var = "ndensity", bins = 12) +
    # geom_density2d(colour = "grey30", linewidth = 0.25, alpha = 0.6,
    #                contour_var = "ndensity", bins = 12) +
    # scale_fill_viridis_d(option = "viridis", direction = -1,
    #                      name = "Density", guide = guide_legend(reverse = TRUE)) +
    labs(title = title, x = xlabel, y = "TrackLength / MaxDis") +
    theme_bw(base_size = 10) +
    theme(plot.title       = element_text(face = "bold", size = 10),
          panel.grid.minor = element_blank(),
          legend.position  = "right",
          legend.key.size  = unit(0.4, "cm"),
          legend.text      = element_text(size = 7),
          legend.title     = element_text(size = 8))
  
  if (!is.null(xlim)) p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
  if (is.null(xlim) && !is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}


# --- 6a. STEP 5: BIVARIATE GMM — NOISE DETECTION --------------------------------
# Fits 2D GMM (EEE) in two feature spaces to detect noise tracks from
# mis-segmentation (inflated TrackLength for static cells).
#
# Space A: log_MaxDis + log_TrackLength
#   Noise = high TrackLength relative to MaxDis (off the main diagonal).
# Space B: log_inv_outreach_ratio + log_MaxDis
#   Noise = high ratio + low MaxDis. Requires tracks_clean_v1.csv.
#
# Decision boundary shown as 200×200 grid overlay.
# Reference: Scrucca L et al. (2016) mclust 5. R J 8(1):289-317.

cat("\n=== STEP 5: BIVARIATE GMM — NOISE DETECTION ===\n")

run_bivar_gmm <- function(df, var1, var2, label) {
  
  cat("\n--- Space", label, ":", var1, "+", var2, "---\n")
  
  biv_data <- df %>%
    filter(!is.na(.data[[var1]]), !is.na(.data[[var2]])) %>%
    select(all_of(c("ROIID_global", "ID", "MouseID", var1, var2)))
  
  mat <- biv_data %>% select(all_of(c(var1, var2))) %>% as.matrix()
  
  models_to_fit <- if (!is.null(force_model)) force_model else gmm_model_types
  set.seed(gmm_seed)
  fit        <- Mclust(mat, G = 1:G_max_biv, modelNames = models_to_fit, verbose = FALSE)
  G_opt      <- fit$G
  best_model <- fit$modelName
  cat("Force model:", ifelse(is.null(force_model), "None (BIC selects)", force_model), "\n")
  cat("Selected: G =", G_opt, "| Model =", best_model, "\n")
  
  bic_vals <- fit$BIC[, best_model]
  bic_tbl  <- data.frame(
    G        = seq_along(bic_vals),
    BIC      = as.numeric(bic_vals),
    Selected = seq_along(bic_vals) == G_opt
  ) %>%
    mutate(Delta_BIC = BIC - max(BIC, na.rm = TRUE))
  cat("BIC table (model:", best_model, "):\n"); print(bic_tbl, digits = 4)
  
  biv_data$component <- apply(fit$z, 1, which.max)
  biv_data$max_prob  <- apply(fit$z, 1, max)
  
  comp_sum <- biv_data %>%
    group_by(component) %>%
    summarise(
      n         = n(),
      mean_var1 = round(mean(.data[[var1]], na.rm = TRUE), 3),
      sd_var1   = round(sd(.data[[var1]],   na.rm = TRUE), 3),
      mean_var2 = round(mean(.data[[var2]], na.rm = TRUE), 3),
      sd_var2   = round(sd(.data[[var2]],   na.rm = TRUE), 3),
      .groups   = "drop"
    )
  names(comp_sum)[3:6] <- c(paste0("mean_", var1), paste0("sd_", var1),
                             paste0("mean_", var2), paste0("sd_", var2))
  cat("\nComponent summary:\n"); print(comp_sum)
  cat("\nSet noise_comp_", label, " in Section 1 after reviewing plots.\n", sep = "")
  
  x_seq          <- seq(min(mat[, 1]) - 0.2, max(mat[, 1]) + 0.2, length.out = 200)
  y_seq          <- seq(min(mat[, 2]) - 0.2, max(mat[, 2]) + 0.2, length.out = 200)
  grid_df        <- expand.grid(x = x_seq, y = y_seq)
  names(grid_df) <- c(var1, var2)
  grid_pred      <- predict(fit, newdata = as.matrix(grid_df))
  grid_df$component <- factor(apply(grid_pred$z, 1, which.max))
  
  list(fit = fit, G_opt = G_opt, best_model = best_model, bic_tbl = bic_tbl,
       biv_data = biv_data, comp_sum = comp_sum,
       grid_df = grid_df, var1 = var1, var2 = var2, label = label)
}

# Compute log transforms (checks if already present in file before recomputing)
tracks_biv <- tracks %>%
  mutate(
    log_TrackLength        = if ("log_TrackLength"        %in% names(tracks))
      log_TrackLength        else log(TrackLength        + 1e-6),
    log_MaxDis             = if ("log_MaxDis"             %in% names(tracks))
      log_MaxDis             else log(MaxDis             + 1e-6),
    log_inv_outreach_ratio = if ("log_inv_outreach_ratio" %in% names(tracks))
      log_inv_outreach_ratio else log(inv_outreach_ratio + 1e-6)
  )

biv_A <- run_bivar_gmm(tracks_biv, "log_MaxDis", "log_TrackLength", "A")

space_B_available <- "inv_outreach_ratio" %in% names(tracks_biv) &&
  any(!is.na(tracks_biv$inv_outreach_ratio))

if (space_B_available) {
  biv_B <- run_bivar_gmm(tracks_biv, "log_inv_outreach_ratio", "log_MaxDis", "B")
} else {
  cat("\nSpace B skipped: inv_outreach_ratio not found or all NA.\n")
  cat("Set input_file <- 'tracks_clean_v1.csv' in Section 1 to enable Space B.\n")
}


# --- 6b. BIVARIATE GMM PLOTS --------------------------------------------------

make_bivar_plots <- function(biv_result) {
  
  var1     <- biv_result$var1
  var2     <- biv_result$var2
  label    <- biv_result$label
  G_opt      <- biv_result$G_opt
  best_model <- biv_result$best_model
  bic_tbl    <- biv_result$bic_tbl
  biv_data <- biv_result$biv_data %>% mutate(component = factor(component))
  grid_df  <- biv_result$grid_df
  comp_pal <- comp_colours[as.character(seq_len(G_opt))]
  
  p_bic <- ggplot(bic_tbl, aes(x = G, y = BIC)) +
    geom_line(colour = "grey50", linewidth = 0.8) +
    geom_point(aes(colour = Selected, size = Selected)) +
    scale_colour_manual(values = c("FALSE" = "grey50", "TRUE" = "#E63946"),
                        guide = "none") +
    scale_size_manual(values = c("FALSE" = 2.5, "TRUE" = 4.5), guide = "none") +
    geom_vline(xintercept = G_opt, linetype = "dashed",
               colour = "#E63946", linewidth = 0.7) +
    scale_x_continuous(breaks = 1:G_max_biv) +
    labs(title    = paste0("Space ", label, " — BIC (", best_model, ")"),
         subtitle = "Higher = better | Red = selected G",
         x = "Components (G)", y = "mclust BIC") +
    theme_bw(base_size = 11) +
    theme(plot.title       = element_text(face = "bold", size = 11),
          plot.subtitle    = element_text(size = 9, colour = "grey40"),
          panel.grid.minor = element_blank())
  
  p_scatter <- ggplot() +
    geom_raster(data = grid_df,
                aes(x = .data[[var1]], y = .data[[var2]], fill = component),
                alpha = 0.15, interpolate = TRUE) +
    geom_point(data = biv_data,
               aes(x = .data[[var1]], y = .data[[var2]], colour = component),
               size = 0.85, alpha = 0.3) +
    scale_fill_manual(values   = comp_pal, name = "Region") +
    scale_colour_manual(values = comp_pal, name = "Component") +
    labs(title    = paste0("Space ", label, " — decision boundary"),
         subtitle = paste0("Shaded = decision regions | G = ", G_opt, "(", best_model, ")"),
         x = var1, y = var2) +
    theme_bw(base_size = 11) +
    theme(plot.title       = element_text(face = "bold", size = 11),
          plot.subtitle    = element_text(size = 9, colour = "grey40"),
          panel.grid.minor = element_blank())
  
  mouse_tbl <- biv_data %>%
    group_by(MouseID, component) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(MouseID) %>%
    mutate(pct = n / sum(n) * 100) %>%
    ungroup()
  
  p_mouse <- ggplot(mouse_tbl, aes(x = MouseID, y = pct, fill = component)) +
    geom_col(width = 0.6, colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = comp_pal, name = "Component") +
    labs(title    = paste0("Space ", label, " — membership by Mouse"),
         subtitle = "Noise component should be M2-enriched if mis-segmentation hypothesis holds",
         x = "Mouse", y = "Tracks (%)") +
    theme_bw(base_size = 11) +
    theme(plot.title       = element_text(face = "bold", size = 11),
          plot.subtitle    = element_text(size = 9, colour = "grey40"),
          panel.grid.minor = element_blank())
  
  (p_bic | p_scatter | p_mouse) +
    plot_annotation(
      title = paste0("Bivariate GMM — Space ", label,
                     " (", var1, " + ", var2, ")"),
      theme = theme(plot.title = element_text(face = "bold", size = 13))
    )
}

p_biv_A <- make_bivar_plots(biv_A)
if (space_B_available) p_biv_B <- make_bivar_plots(biv_B)


# # --- 6c. BIVARIATE NOISE FLAGGING [CONDITIONAL] --------------------------------
# 
# if (flag_noise) {
#   
#   if (is.na(noise_comp_A) && is.na(noise_comp_B)) {
#     cat("\n[Step 5] flag_noise = TRUE but noise_comp_A and noise_comp_B are NA.\n")
#     cat("         Set at least one in Section 1 and re-run.\n")
#     
#   } else {
#     
#     cat("\n=== BIVARIATE NOISE FLAGGING ===\n")
#     
#     flag_A <- biv_A$biv_data %>%
#       select(ROIID_global, ID) %>%
#       mutate(noise_A = biv_A$biv_data$component == noise_comp_A)
#     
#     if (space_B_available && !is.na(noise_comp_B)) {
#       flag_B <- biv_B$biv_data %>%
#         select(ROIID_global, ID) %>%
#         mutate(noise_B = biv_B$biv_data$component == noise_comp_B)
#       flags <- full_join(flag_A, flag_B, by = c("ROIID_global", "ID")) %>%
#         mutate(bivariate_GMM_noise = replace_na(noise_A, FALSE) |
#                  replace_na(noise_B, FALSE)) %>%
#         select(ROIID_global, ID, bivariate_GMM_noise)
#     } else {
#       flags <- flag_A %>%
#         mutate(bivariate_GMM_noise = replace_na(noise_A, FALSE)) %>%
#         select(ROIID_global, ID, bivariate_GMM_noise)
#     }
#     
#     cat("Flag counts:\n")
#     print(table(flags$bivariate_GMM_noise, useNA = "ifany"))
#     cat("Proportion flagged:", round(mean(flags$bivariate_GMM_noise) * 100, 1), "%\n")
#     
#     flag_mouse <- tracks %>%
#       select(ROIID_global, ID, MouseID) %>%
#       left_join(flags, by = c("ROIID_global", "ID")) %>%
#       group_by(MouseID) %>%
#       summarise(n_total   = n(),
#                 n_noise   = sum(bivariate_GMM_noise, na.rm = TRUE),
#                 pct_noise = round(n_noise / n_total * 100, 1),
#                 .groups   = "drop")
#     cat("\nNoise flags by Mouse:\n"); print(flag_mouse)
#     
#     tracks_out <- tracks %>%
#       left_join(flags, by = c("ROIID_global", "ID")) %>%
#       mutate(bivariate_GMM_noise = replace_na(bivariate_GMM_noise, FALSE))
#     write_csv(tracks_out, "tracks_clean_v1.csv")
#     cat("\ntracks_clean_v1.csv overwritten with bivariate_GMM_noise column.\n")
#   }
#   
# } else {
#   cat("\n[Step 5] flag_noise = FALSE — bivariate flagging skipped.\n")
#   cat("  1. Review bivariate GMM plots\n")
#   cat("  2. Identify noise component in Space A (and B if available)\n")
#   cat("  3. Set noise_comp_A / noise_comp_B in Section 1\n")
#   cat("  4. Set flag_noise <- TRUE and re-run\n")
# }


# # ---- 6d. STEP 6: 3D GMM — NOISE DETECTION ------------------------------------- <----- NEED TO RUN 6a !!!
# Extends bivariate GMM to three dimensions by adding mean_turning_angle.
# Rationale: mis-segmented static cells have erratic turning angles (spurious
# detections jump randomly) independent of the TrackLength/MaxDis signal.
#
# Space A: log_TrackLength + log_MaxDis + mean_turning_angle
# Space B: log_inv_outreach_ratio + log_MaxDis + mean_turning_angle
#
# Model preference: EEE (biological interpretability). VVV rejected despite
# potentially higher BIC (biologically implausible: all components fully free).
# G_max = 4: noise / static / intermediate / migrating.
#
# To test more models: add to gmm_model_types in Section 1.
# To test more G: increase G_max_biv in Section 1.
# force_model = NULL: evaluates all gmm_model_types, selects by BIC.
#
# Visualisation: three 2D projection panels per space; decision boundary in each
# panel computed at the mean of the fixed (third) variable.
#
# Reference: Scrucca L et al. (2016) mclust 5. R J 8(1):289-317.

cat("\n=== STEP 6: 3D GMM — NOISE DETECTION ===\n")

run_gmm_2D_3D <- function(df, vars, label) {
  
  n_vars <- length(vars)
  cat("\n--- Space", label, ":", paste(vars, collapse = " + "), "---\n")
  cat("Dimensionality:", n_vars, "D\n")
  
  biv_data <- df %>%
    filter(complete.cases(select(., all_of(vars)))) %>%
    select(all_of(c("ROIID_global", "ID", "MouseID", vars)))
  
  mat <- biv_data %>% select(all_of(vars)) %>% as.matrix()
  cat("Tracks entering GMM:", nrow(mat), "\n")
  
  models_to_fit <- if (!is.null(force_model)) force_model else gmm_model_types
  
  set.seed(gmm_seed)
  fit        <- Mclust(mat, G = 1:G_max_biv,
                       modelNames = models_to_fit, verbose = FALSE)
  G_opt      <- fit$G
  best_model <- fit$modelName
  cat("Force model:", ifelse(is.null(force_model), "None (BIC selects)", force_model), "\n")
  cat("Selected: G =", G_opt, "| Model =", best_model, "\n")
  
  bic_mat  <- fit$BIC
  bic_long <- expand.grid(
    G     = as.integer(rownames(bic_mat)),
    Model = colnames(bic_mat)
  ) %>%
    mutate(
      BIC = mapply(function(g, m) bic_mat[as.character(g), m], G, Model)
    ) %>%
    filter(!is.na(BIC)) %>%
    mutate(
      Selected  = G == G_opt & Model == best_model,
      Delta_BIC = BIC - max(BIC, na.rm = TRUE)
    )
  
  cat("\nBIC table (all models, higher = better):\n")
  print(pivot_wider(select(bic_long, G, Model, BIC),
                    names_from = Model, values_from = BIC), digits = 5)
  
  biv_data$component <- apply(fit$z, 1, which.max)
  biv_data$max_prob  <- apply(fit$z, 1, max)
  
  comp_sum <- biv_data %>%
    group_by(component) %>%
    summarise(
      n = n(),
      across(all_of(vars),
             list(mean = ~ round(mean(.x, na.rm = TRUE), 3),
                  sd   = ~ round(sd(.x,   na.rm = TRUE), 3))),
      .groups = "drop"
    )
  cat("\nComponent summary:\n"); print(comp_sum)
  cat("\nSet noise_comp_3D_", label, " in Section 1 after reviewing plots.\n", sep = "")
  
  # Decision grids — one per variable pair; third variable fixed at its mean
  var_pairs <- combn(vars, 2, simplify = FALSE)
  
  grids <- lapply(var_pairs, function(pair) {
    v1         <- pair[1]
    v2         <- pair[2]
    other_vars <- setdiff(vars, pair)
    
    x_seq   <- seq(min(mat[, v1]) - 0.2, max(mat[, v1]) + 0.2, length.out = 150)
    y_seq   <- seq(min(mat[, v2]) - 0.2, max(mat[, v2]) + 0.2, length.out = 150)
    grid_df <- expand.grid(x = x_seq, y = y_seq)
    names(grid_df) <- c(v1, v2)
    
    for (ov in other_vars) grid_df[[ov]] <- mean(mat[, ov], na.rm = TRUE)
    
    grid_mat          <- as.matrix(grid_df[, vars, drop = FALSE])
    grid_pred         <- predict(fit, newdata = grid_mat)
    grid_df$component <- factor(apply(grid_pred$z, 1, which.max))
    
    list(grid_df = grid_df, v1 = v1, v2 = v2, other_vars = other_vars)
  })
  
  list(fit        = fit,
       G_opt      = G_opt,
       best_model = best_model,
       bic_long   = bic_long,
       biv_data   = biv_data,
       comp_sum   = comp_sum,
       grids      = grids,
       vars       = vars,
       label      = label)
}

if (!"mean_turning_angle" %in% names(tracks_biv))
  stop("mean_turning_angle not found. Ensure input_file is tracks_clean_v1.csv.")

gmm_3D_A <- run_gmm_2D_3D(
  tracks_biv,
  vars  = c("log_TrackLength", "log_MaxDis", "mean_turning_angle"),
  label = "A"
)

if (space_B_available) {
  gmm_3D_B <- run_gmm_2D_3D(
    tracks_biv,
    vars  = c("log_inv_outreach_ratio", "log_MaxDis", "mean_turning_angle"),
    label = "B"
  )
} else {
  cat("\n3D Space B skipped: inv_outreach_ratio not available.\n")
}


make_3D_plots <- function(gmm_result) {
  
  vars       <- gmm_result$vars
  label      <- gmm_result$label
  G_opt      <- gmm_result$G_opt
  best_model <- gmm_result$best_model
  bic_long   <- gmm_result$bic_long
  biv_data   <- gmm_result$biv_data %>% mutate(component = factor(component))
  grids      <- gmm_result$grids
  comp_pal   <- comp_colours[as.character(seq_len(G_opt))]
  # model_pal: add a colour for every model you add to gmm_model_types.
  # mclust will error on the BIC plot if a model name is absent from this palette.
  model_pal  <- c(
    "EEE" = "#457B9D", "VEE" = "#2A9D8F", "VVV" = "#E63946",
    "EEV" = "#F4A261", "VEV" = "#6A0572", "EVE" = "#264653",
    "EVI" = "#E9C46A", "VVI" = "#606C38", "VVE" = "#BC6C25","EVV" = "gray"
  )
  
  p_bic <- ggplot(bic_long, aes(x = G, y = BIC, colour = Model, group = Model)) +
    geom_line(linewidth = 0.8) +
    geom_point(aes(size = Selected, shape = Selected)) +
    geom_point(data   = bic_long %>% filter(Selected),
               colour = "black", size = 5, shape = 1, stroke = 1.2) +
    scale_colour_manual(values = model_pal, name = "Model") +
    scale_size_manual(values   = c("FALSE" = 2,  "TRUE" = 4),  guide = "none") +
    scale_shape_manual(values  = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
    scale_x_continuous(breaks  = 1:G_max_biv) +
    labs(title    = paste0("3D Space ", label, " — BIC across models"),
         subtitle = paste0("Selected: G = ", G_opt, ", Model = ", best_model,
                           " | Higher = better"),
         x = "Components (G)", y = "mclust BIC") +
    theme_bw(base_size = 11) +
    theme(plot.title       = element_text(face = "bold", size = 11),
          plot.subtitle    = element_text(size = 9, colour = "grey40"),
          panel.grid.minor = element_blank())
  
  scatter_panels <- lapply(grids, function(g) {
    v1      <- g$v1
    v2      <- g$v2
    fix_str <- if (length(g$other_vars) > 0) paste0(g$other_vars, " @ mean") else ""
    
    ggplot() +
      geom_raster(data = g$grid_df,
                  aes(x = .data[[v1]], y = .data[[v2]], fill = component),
                  alpha = 0.15, interpolate = TRUE) +
      geom_point(data = biv_data,
                 aes(x = .data[[v1]], y = .data[[v2]], colour = component),
                 size = 0.8, alpha = 0.25) +
      scale_fill_manual(values   = comp_pal, name = "Region") +
      scale_colour_manual(values = comp_pal, name = "Component") +
      labs(title = paste0(v1, " vs ", v2), subtitle = fix_str, x = v1, y = v2) +
      theme_bw(base_size = 10) +
      theme(plot.title       = element_text(face = "bold", size = 10),
            plot.subtitle    = element_text(size = 8, colour = "grey40"),
            panel.grid.minor = element_blank())
  })
  
  mouse_tbl <- biv_data %>%
    group_by(MouseID, component) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(MouseID) %>%
    mutate(pct = n / sum(n) * 100) %>%
    ungroup()
  
  p_mouse <- ggplot(mouse_tbl, aes(x = MouseID, y = pct, fill = component)) +
    geom_col(width = 0.6, colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = comp_pal, name = "Component") +
    labs(title    = paste0("3D Space ", label, " — membership by Mouse"),
         subtitle = "Noise component should be M2-enriched if mis-segmentation holds",
         x = "Mouse", y = "Tracks (%)") +
    theme_bw(base_size = 11) +
    theme(plot.title       = element_text(face = "bold", size = 11),
          plot.subtitle    = element_text(size = 9, colour = "grey40"),
          panel.grid.minor = element_blank())
  
  top_row    <- wrap_plots(scatter_panels, nrow = 1)
  bottom_row <- (p_bic | p_mouse) + plot_layout(widths = c(1.2, 0.8))
  
  (top_row / bottom_row) +
    plot_layout(heights = c(1.6, 1)) +
    plot_annotation(
      title = paste0("3D GMM Space ", label, " (",
                     paste(vars, collapse = " + "),
                     ") — G = ", G_opt, ", Model = ", best_model),
      theme = theme(plot.title = element_text(face = "bold", size = 13))
    )
}

p_3D_A <- make_3D_plots(gmm_3D_A)
if (space_B_available) p_3D_B <- make_3D_plots(gmm_3D_B)


# ---- 3D Noise flagging [gated] -----------------------------------------------

if (flag_noise_3D) {
  
  if (is.na(noise_comp_3D_A) && is.na(noise_comp_3D_B)) {
    cat("\n[Step 6] flag_noise_3D = TRUE but noise_comp_3D_A and _B are both NA.\n")
    cat("         Set at least one in Section 1 and re-run.\n")
    
  } else {
    
    cat("\n=== 3D NOISE FLAGGING ===\n")
    
    flag_3D_A <- gmm_3D_A$biv_data %>%
      select(ROIID_global, ID) %>%
      mutate(noise_3D_A = gmm_3D_A$biv_data$component == noise_comp_3D_A)
    
    if (space_B_available && !is.na(noise_comp_3D_B)) {
      flag_3D_B <- gmm_3D_B$biv_data %>%
        select(ROIID_global, ID) %>%
        mutate(noise_3D_B = gmm_3D_B$biv_data$component == noise_comp_3D_B)
      flags_3D <- full_join(flag_3D_A, flag_3D_B, by = c("ROIID_global", "ID")) %>%
        mutate(bivariate_GMM_noise_3D = replace_na(noise_3D_A, FALSE) |
                 replace_na(noise_3D_B, FALSE)) %>%
        select(ROIID_global, ID, bivariate_GMM_noise_3D)
    } else {
      flags_3D <- flag_3D_A %>%
        mutate(bivariate_GMM_noise_3D = replace_na(noise_3D_A, FALSE)) %>%
        select(ROIID_global, ID, bivariate_GMM_noise_3D)
    }
    
    cat("3D flag counts:\n")
    print(table(flags_3D$bivariate_GMM_noise_3D, useNA = "ifany"))
    cat("Proportion flagged:",
        round(mean(flags_3D$bivariate_GMM_noise_3D) * 100, 1), "%\n")
    
    flag_mouse_3D <- tracks %>%
      select(ROIID_global, ID, MouseID) %>%
      left_join(flags_3D, by = c("ROIID_global", "ID")) %>%
      group_by(MouseID) %>%
      summarise(n_total   = n(),
                n_noise   = sum(bivariate_GMM_noise_3D, na.rm = TRUE),
                pct_noise = round(n_noise / n_total * 100, 1),
                .groups   = "drop")
    cat("\n3D noise flags by Mouse:\n"); print(flag_mouse_3D)
    
    tracks_out <- tracks %>%
      left_join(flags_3D, by = c("ROIID_global", "ID")) %>%
      mutate(bivariate_GMM_noise_3D = replace_na(bivariate_GMM_noise_3D, FALSE))
    write_csv(tracks_out, "tracks_clean_v1.csv")
    cat("\ntracks_clean_v1.csv overwritten with bivariate_GMM_noise_3D column.\n")
  }
  
} else {
  cat("\n[Step 6] flag_noise_3D = FALSE — 3D flagging skipped.\n")
  cat("  1. Review 3D GMM plots\n")
  cat("  2. Identify noise component in 3D Space A (and B if available)\n")
  cat("  3. Set noise_comp_3D_A / noise_comp_3D_B in Section 1\n")
  cat("  4. Set flag_noise_3D <- TRUE and re-run\n")
}


# --- 6e. COMPONENT PROFILING AND EXPORT ---------------------------------------
# Profiles the 3D GMM components on variables NOT used in GMM fitting.
# Separation here is independent external validation of the component structure.
# To profile a different model: change force_model in Section 1 and re-run 6d–6e.

cat("\n=== COMPONENT PROFILING AND EXPORT ===\n")
cat("Model:", gmm_3D_A$best_model, "| G =", gmm_3D_A$G_opt, "\n")

make_profile_plots <- function(gmm_result, tracks_full) {
  
  comp_assign <- gmm_result$biv_data %>%
    select(ROIID_global, ID, component) %>%
    mutate(component = factor(component))
  
  profile_data <- tracks_full %>%
    inner_join(comp_assign, by = c("ROIID_global", "ID"))
  
  comp_pal <- comp_colours[as.character(seq_len(gmm_result$G_opt))]
  
  # 6 profiling variables — 2×3 grid
  # log_Speed, Track_Straightness, alpha_30pct, TrackDuration: independent of GMM
  # MaxDis (raw µm) and TrackLength (raw µm): NOT log-transformed here so units
  #   are directly interpretable; these ARE used (log-transformed) in GMM fitting,
  #   so separation here reflects the raw-scale distribution across components
  profile_vars <- list(
    list(col = "log_Speed",          ylabel = "Log Speed (µm/min)",    title = "Log Speed"),
    list(col = "Track_Straightness", ylabel = "Straightness (0-1)",    title = "Straightness"),
    list(col = "alpha_30pct",        ylabel = "MSD α (30%)",           title = "MSD α (30%)"),
    list(col = "TrackDuration",      ylabel = "Duration (s)",          title = "Track Duration"),
    list(col = "MaxDis",             ylabel = "Max Displacement (µm)", title = "Max Displacement"),
    list(col = "log_TrackLength",        ylabel = "Log Track Length (µm)",     title = "Track Length")    #log_TrackLength, Log Track Length
  )
  
  panels <- lapply(profile_vars, function(v) {
    plot_data <- profile_data %>% filter(!is.na(.data[[v$col]]))
    
    ggplot(plot_data, aes(x = component, y = .data[[v$col]], fill = component)) +
      geom_violin(alpha = 0.4, colour = NA, trim = FALSE) +
      geom_boxplot(width = 0.2, outlier.size = 0.4,
                   outlier.alpha = 0.3, colour = "grey30") +
      scale_fill_manual(values = comp_pal, guide = "none") +
      labs(title = v$title, x = "Component", y = v$ylabel) +
      theme_bw(base_size = 11) +
      theme(plot.title       = element_text(face = "bold", size = 11),
            panel.grid.minor = element_blank())
  })
  
  wrap_plots(panels, nrow = 2, ncol = 3) +
    plot_annotation(
      title    = paste0("Component profiling — Space ", gmm_result$label,
                        " | G = ", gmm_result$G_opt,
                        ", Model = ", gmm_result$best_model),
      subtitle = "Variables NOT used in GMM fitting — separation is independent external validation",
      theme    = theme(plot.title    = element_text(face = "bold", size = 13),
                       plot.subtitle = element_text(size = 10, colour = "grey40"))
    )
}

p_profile_A <- make_profile_plots(gmm_3D_A, tracks_biv)
if (space_B_available)
  p_profile_B <- make_profile_plots(gmm_3D_B, tracks_biv)

print(p_profile_A)
if (space_B_available) print(p_profile_B)

# -- Component assignment export -----------------------------------------------
# Writes: gmm_component_assignments_<model>_Space<A/B>.csv
# Columns: ROIID_global, ID, MouseID, Day, GMM_component, GMM_max_prob,
#          P_comp1 ... P_comp<G>

export_assignments <- function(gmm_result, space_label) {
  
  prob_mat        <- as.data.frame(gmm_result$fit$z)
  names(prob_mat) <- paste0("P_comp", seq_len(gmm_result$G_opt))
  
  day_lookup <- tracks_biv %>%
    select(ROIID_global, ID, Day) %>%
    distinct()
  
  out <- gmm_result$biv_data %>%
    select(ROIID_global, ID, MouseID, component, max_prob) %>%
    rename(GMM_component = component,
           GMM_max_prob  = max_prob) %>%
    left_join(day_lookup, by = c("ROIID_global", "ID")) %>%
    bind_cols(prob_mat)
  
  filename <- paste0("gmm_component_assignments_",
                     gmm_result$best_model,
                     "_Space", space_label, ".csv")
  write_csv(out, filename)
  cat("Saved:", filename, "\n")
  cat("  Tracks:", nrow(out), "| G =", gmm_result$G_opt,
      "| Model:", gmm_result$best_model, "\n")
  cat("  Component counts:\n")
  print(table(out$GMM_component))
  cat("  Tracks with max_prob >= 0.80:", sum(out$GMM_max_prob >= 0.80), "\n\n")
}

export_assignments(gmm_3D_A, "A")
if (space_B_available) export_assignments(gmm_3D_B, "B")

# ---- Export cluster descriptives from fitted GMM parameters --------------------
# Means and SDs come from fit$parameters (the actual fitted Gaussians), not
# from empirical summaries of assigned tracks. For EEE the covariance matrix
# is shared across components (SDs will be equal); for VEE/VVV it differs.
# Columns: Space, Model, G_opt, Component, n_tracks, proportion,
#          mean_<var>, sd_<var> for each variable in the GMM space.

export_cluster_descriptives <- function(gmm_result, space_label) {
  fit   <- gmm_result$fit
  vars  <- gmm_result$vars
  G     <- gmm_result$G_opt
  model <- gmm_result$best_model
  
  # Covariance: EEE gives a p×p matrix; VEE/VVV give a p×p×G array
  Sigma_all  <- fit$parameters$variance$sigma
  shared_cov <- (length(dim(Sigma_all)) == 2)
  
  rows <- lapply(seq_len(G), function(k) {
    means_k <- fit$parameters$mean[, k]
    sds_k   <- if (shared_cov) sqrt(diag(Sigma_all)) else
      sqrt(diag(Sigma_all[, , k]))
    n_k     <- round(sum(fit$z[, k]))
    prop_k  <- round(fit$parameters$pro[k], 4)
    
    row <- data.frame(Space = space_label, Model = model,
                      G_opt = G, Component = k,
                      n_tracks = n_k, proportion = prop_k)
    for (v in vars) {
      row[[paste0("mean_", v)]] <- round(means_k[v], 4)
      row[[paste0("sd_",   v)]] <- round(sds_k[v],   4)
    }
    row
  })
  
  out   <- do.call(rbind, rows)
  
  ### filename
  fname <- paste0("gmm_cluster_descriptives_", model,
                  "_G", G, "_Space", space_label, ".csv")
  write_csv(out, fname)
  cat("Saved:", fname, "\n")
  invisible(out)
}

# ---- Export profile summary of tracks per component ----------------------------
# Empirical summaries (mean, SD, median, IQR) of profiling variables for tracks
# hard-assigned to each component. Variables are independent of GMM fitting
# except MaxDis and TrackLength which were used log-transformed in fitting —
# raw µm values are reported here for interpretability.
# Columns: Space, Model, G_opt, Component, n_tracks,
#          <var>_mean, <var>_sd, <var>_median, <var>_IQR for each profiling var.

export_profile_summary <- function(gmm_result, tracks_full, space_label) {
  model <- gmm_result$best_model
  G     <- gmm_result$G_opt
  
  comp_assign  <- gmm_result$biv_data %>%
    select(ROIID_global, ID, component)
  
  profile_data <- tracks_full %>%
    inner_join(comp_assign, by = c("ROIID_global", "ID"))
  
  profile_vars <- c("log_Speed", "Track_Straightness", "alpha_30pct",
                    "TrackDuration", "MaxDis", "TrackLength", "log_TrackLength")
  
  rows <- lapply(seq_len(G), function(k) {
    d   <- profile_data %>% filter(component == k)
    row <- data.frame(Space = space_label, Model = model,
                      G_opt = G, Component = k, n_tracks = nrow(d))
    for (v in profile_vars) {
      if (v %in% names(d)) {
        vals <- d[[v]][!is.na(d[[v]])]
        row[[paste0(v, "_mean")]]   <- round(mean(vals),   4)
        row[[paste0(v, "_sd")]]     <- round(sd(vals),     4)
        row[[paste0(v, "_median")]] <- round(median(vals), 4)
        row[[paste0(v, "_IQR")]]    <- round(IQR(vals),    4)
      }
    }
    row
  })
  
  out   <- do.call(rbind, rows)
  
  ## filename
  fname <- paste0("gmm_profile_summary_", model,
                  "_G", G, "_Space", space_label, ".csv")
  write_csv(out, fname)
  cat("Saved:", fname, "\n")
  invisible(out)
}

export_cluster_descriptives(gmm_3D_A, "A")
if (space_B_available) export_cluster_descriptives(gmm_3D_B, "B")

export_profile_summary(gmm_3D_A, tracks_biv, "A")
if (space_B_available) export_profile_summary(gmm_3D_B, tracks_biv, "B")


# ---- Save profile plots to PDF -------------------------------------------------

model_label <- gmm_3D_A$best_model
pdf(paste0("component_profiles_", model_label, ".pdf"), width = 14, height = 8)
print(p_profile_A)
if (space_B_available) print(p_profile_B)
dev.off()
cat("Saved: component_profiles_", model_label, ".pdf\n", sep = "")

# ---- Save fitted model objects for retrospective prediction (Section 6f) -------

saveRDS(gmm_3D_A$fit,
        paste0("gmm_3D_A_", gmm_3D_A$best_model, "_fit.rds"))
if (space_B_available)
  saveRDS(gmm_3D_B$fit,
          paste0("gmm_3D_B_", gmm_3D_B$best_model, "_fit.rds"))
cat("Saved: .rds model objects (used by Section 6f for retrospective prediction)\n")


# --- 6f. RETROSPECTIVE PREDICTION [CONDITIONAL] --------------------------------
# Loads a previously saved mclust model fit (.rds set in Section 1) and applies
# predict() to classify new ROI tracks into the existing components.
# The GMM is NOT re-fitted — all component boundaries come from the original model.
# Model name and G are read directly from the .rds file — no need to specify them.
#
# Gate: run_retro_predict <- TRUE in Section 1 (default FALSE).
# Required settings in Section 1:
#   new_tracks_file    : CSV with new ROI tracks (same format as input_file).
#   retro_rds_B        : path to the Space A .rds file saved by Section 6e.
#   append_to_existing : if TRUE, also writes a combined CSV (original + new)
#                        for downstream scripts (track visualisation, 03/04).
#
# Output format matches export_assignments() exactly:
#   ROIID_global, ID, MouseID, Day, GMM_component, GMM_max_prob, P_comp1..P_comp4
#
# Suggested workflow:
#   1. Generate tracks_clean_vX.csv for the new ROI batch (Script 01)
#   2. Set run_retro_predict <- TRUE, update new_tracks_file and retro_rds_B
#   3. Run Section 6f only — inspect the newROIs CSV and profile PDF
#   4. Once satisfied, set append_to_existing <- TRUE and re-run to produce
#      the combined file ready for Scripts 03/04 and track_visualisation.R

if (run_retro_predict) {
  
  cat("\n=== SECTION 6f: RETROSPECTIVE PREDICTION ===\n")
  
  # Load and prepare new tracks
  if (!file.exists(new_tracks_file))
    stop("[6f] File not found: ", new_tracks_file,
         "\nSet new_tracks_file in Section 1.")
  
  cat("Loading:", new_tracks_file, "\n")
  tracks_new <- read_csv(new_tracks_file, show_col_types = FALSE) %>%
    mutate(
      Location        = factor(Location,  levels = c("Core", "CoreEdge", "Edge")),
      Treatment       = factor(Treatment, levels = c("Mock", "RP")),
      Day             = factor(Day,       levels = c("D1", "D2")),
      MouseID         = factor(MouseID),
      TumourID_global = factor(TumourID_global),
      SliceID_global  = factor(SliceID_global),
      ROIID_global    = factor(ROIID_global)
    )
  cat("New tracks loaded:", nrow(tracks_new), "\n")
  
  # Compute same log transforms as tracks_biv (checks file first)
  tracks_new_biv <- tracks_new %>%
    mutate(
      log_TrackLength = if ("log_TrackLength" %in% names(tracks_new))
        log_TrackLength else log(TrackLength + 1e-6),
      log_MaxDis      = if ("log_MaxDis" %in% names(tracks_new))
        log_MaxDis      else log(MaxDis      + 1e-6),
      log_inv_outreach_ratio = if ("log_inv_outreach_ratio" %in% names(tracks_new))
        log_inv_outreach_ratio else if ("inv_outreach_ratio" %in% names(tracks_new))
          log(inv_outreach_ratio + 1e-6) else NA_real_
    )
  
  # Helper: load RDS, predict, export
  # Model name and G are read from the RDS — no external label needed.
  predict_space <- function(space_label, vars, tracks_df, rds_path) {
    
    if (is.na(rds_path) || !file.exists(rds_path)) {
      cat("[6f] RDS not found:", rds_path, "— Space", space_label, "skipped.\n")
      return(invisible(NULL))
    }
    
    cat("\n--- Space", space_label, ":", paste(vars, collapse = " + "), "---\n")
    fit         <- readRDS(rds_path)
    G_fit       <- fit$G
    model_label <- fit$modelName
    cat("Loaded model: G =", G_fit, "| Model =", model_label, "\n")
    
    pred_data <- tracks_df %>%
      filter(complete.cases(select(., all_of(vars)))) %>%
      select(all_of(c("ROIID_global", "ID", "MouseID", vars)))
    cat("Tracks entering prediction:", nrow(pred_data), "\n")
    
    mat         <- pred_data %>% select(all_of(vars)) %>% as.matrix()
    pred_result <- predict(fit, newdata = mat)
    
    prob_mat        <- as.data.frame(pred_result$z)
    names(prob_mat) <- paste0("P_comp", seq_len(G_fit))
    # Pad P_comp columns to G_fit so rbind works across models with different G
    for (i in seq_len(G_fit)) {
      col_nm <- paste0("P_comp", i)
      if (!col_nm %in% names(prob_mat)) prob_mat[[col_nm]] <- NA_real_
    }
    
    day_lookup <- tracks_df %>% select(ROIID_global, ID, Day) %>% distinct()
    
    out <- pred_data %>%
      select(ROIID_global, ID, MouseID) %>%
      mutate(
        GMM_component = apply(pred_result$z, 1, which.max),
        GMM_max_prob  = apply(pred_result$z, 1, max)
      ) %>%
      left_join(day_lookup, by = c("ROIID_global", "ID")) %>%
      bind_cols(prob_mat)
    
    new_file <- paste0("gmm_component_assignments_", model_label,
                       "_G", G_fit, "_Space", space_label, "_newROIs.csv")
    write_csv(out, new_file)
    cat("Saved:", new_file, "\n")
    cat("  Component counts:\n"); print(table(out$GMM_component))
    cat("  Tracks with max_prob >= 0.80:", sum(out$GMM_max_prob >= 0.80), "\n")
    
    if (append_to_existing) {
      orig_file <- paste0("gmm_component_assignments_", model_label,
                          "_Space", space_label, ".csv")
      if (file.exists(orig_file)) {
        orig      <- read_csv(orig_file, show_col_types = FALSE)
        combined  <- bind_rows(orig, out)
        comb_file <- paste0("gmm_component_assignments_", model_label,
                            "_G", G_fit, "_Space", space_label, "_combined.csv")
        write_csv(combined, comb_file)
        cat("Saved combined:", comb_file, "(", nrow(orig), "original +",
            nrow(out), "new =", nrow(combined), "tracks)\n")
      } else {
        cat("[6f] Original file not found:", orig_file,
            "— combined output skipped.\n")
      }
    }
    
    invisible(out)
  }
  
  # Space B                                                             #<--------------- #cHANGE LABEL HERE
  retro_space <- predict_space(
    space_label = "B",
    vars        = c("log_inv_outreach_ratio", "log_MaxDis", "mean_turning_angle"),
    tracks_df   = tracks_new_biv,
    rds_path    = retro_rds_B
  )
  
  # -- Profile plots for new ROIs -----------------------------------------------
  # Reuses make_profile_plots from Section 6e.
  # Model name and G are derived from the same RDS — no extra settings needed.
  # Output: component_profiles_<model>_G<n>_SpaceA_newROIs.pdf
  
  make_retro_profiles <- function(retro_out, space_label, tracks_df, rds_path) {
    if (is.null(retro_out)) return(invisible(NULL))
    fit         <- readRDS(rds_path)
    model_label <- fit$modelName
    G_fit       <- fit$G
    retro_gmm <- list(
      biv_data   = retro_out %>%
        rename(component = GMM_component, max_prob = GMM_max_prob),
      G_opt      = G_fit,
      best_model = model_label,
      label      = space_label
    )
    p     <- make_profile_plots(retro_gmm, tracks_df)                   #<--- run make profile plot helper in section 6e - change to log_TrackLength
    fname <- paste0("component_profiles_", model_label,
                    "_G", G_fit, "_Space", space_label, "_newROIs.pdf")
    pdf(fname, width = 14, height = 8)
    print(p)
    dev.off()
    cat("Saved:", fname, "\n")
    invisible(p)
  }
  
  make_retro_profiles(retro_space, "B", tracks_new_biv, retro_rds_B)  #<--------------- #CHANGE LABLE HERE
  
  # -- Scatter and mouse membership plots for new ROIs --------------------------
  # Equivalent to the scatter_panels + p_mouse from make_3D_plots.
  # Decision boundary grids are omitted — not available for retrospective fits.
  # Output: scatter_<model>_G<n>_Space<label>_newROIs.pdf
  
  make_retro_scatter <- function(retro_out, space_label, vars, tracks_df, rds_path) {
    if (is.null(retro_out)) return(invisible(NULL))
    
    fit         <- readRDS(rds_path)
    model_label <- fit$modelName
    G_fit       <- fit$G
    comp_pal    <- comp_colours[as.character(seq_len(G_fit))]
    
    # Join component assignments onto feature values
    plot_data <- tracks_df %>%
      inner_join(
        retro_out %>% select(ROIID_global, ID, GMM_component),
        by = c("ROIID_global", "ID")
      ) %>%
      mutate(component = factor(GMM_component))
    
    # One scatter panel per variable pair
    var_pairs <- combn(vars, 2, simplify = FALSE)
    
    scatter_panels <- lapply(var_pairs, function(pair) {
      v1 <- pair[1]
      v2 <- pair[2]
      ggplot(plot_data,
             aes(x = .data[[v1]], y = .data[[v2]], colour = component)) +
        geom_point(size = 0.8, alpha = 0.3) +
        scale_colour_manual(values = comp_pal, name = "Component") +
        labs(title = paste0(v1, " vs ", v2), x = v1, y = v2) +
        theme_bw(base_size = 10) +
        theme(plot.title       = element_text(face = "bold", size = 10),
              panel.grid.minor = element_blank())
    })
    
    # Mouse membership bar chart
    mouse_tbl <- plot_data %>%
      group_by(MouseID, component) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(MouseID) %>%
      mutate(pct = n / sum(n) * 100) %>%
      ungroup()
    
    p_mouse <- ggplot(mouse_tbl, aes(x = MouseID, y = pct, fill = component)) +
      geom_col(width = 0.6, colour = "white", linewidth = 0.3) +
      scale_fill_manual(values = comp_pal, name = "Component") +
      labs(title    = paste0("Space ", space_label, " — membership by Mouse (new ROIs)"),
           x = "Mouse", y = "Tracks (%)") +
      theme_bw(base_size = 11) +
      theme(plot.title       = element_text(face = "bold", size = 11),
            panel.grid.minor = element_blank())
    
    combined <- wrap_plots(scatter_panels, nrow = 1) / p_mouse +
      plot_layout(heights = c(1.6, 1)) +
      plot_annotation(
        title = paste0("3D GMM Space ", space_label, " — new ROIs | G = ",
                       G_fit, ", Model = ", model_label),
        theme = theme(plot.title = element_text(face = "bold", size = 13))
      )
    
    fname <- paste0("scatter_", model_label, "_G", G_fit,
                    "_Space", space_label, "_newROIs.pdf")
    pdf(fname, width = 14, height = 10)
    print(combined)
    dev.off()
    cat("Saved:", fname, "\n")
    invisible(combined)
  }
  
  make_retro_scatter(retro_space, "B", retro_vars_B, tracks_new_biv, retro_rds_B)
  
  } else {
  cat("\n[6f] run_retro_predict = FALSE — retrospective prediction skipped.\n")
  cat("  Set run_retro_predict <- TRUE in Section 1 once new ROI tracks are ready.\n")
  cat("  Set new_tracks_file and retro_rds_B in Section 1.\n")
  cat("  Set append_to_existing <- TRUE to also write a combined output file.\n")
}

# -- Descriptive stats and profile summary for new ROIs -----------------------
# Mirrors the 6e exports but for the retrospectively classified batch.
# Cluster descriptives: means/SDs from the loaded RDS fit parameters.
# Profile summary: empirical summaries per component.
#   tracks_df used = tracks_new_biv (new ROIs only) when append_to_existing = FALSE.
#   When append_to_existing = TRUE and tracks_biv is in session, both original
#   and new ROI tracks are combined for a comprehensive profile.
# Output filenames include _newROIs suffix to distinguish from 6e outputs.

export_retro_descriptives <- function(retro_out, rds_path, space_label, tracks_df) {
  if (is.null(retro_out)) return(invisible(NULL))
  
  fit         <- readRDS(rds_path)
  model_label <- fit$modelName
  G_fit       <- fit$G
  vars        <- rownames(fit$parameters$mean)
  Sigma_all   <- fit$parameters$variance$sigma
  shared_cov  <- (length(dim(Sigma_all)) == 2)
  
  
  
  # Cluster descriptives from fitted parameters
  desc_rows <- lapply(seq_len(G_fit), function(k) {
    means_k <- fit$parameters$mean[, k]
    sds_k   <- if (shared_cov) sqrt(diag(Sigma_all)) else
      sqrt(diag(Sigma_all[, , k]))
    n_k     <- sum(retro_out$GMM_component == k, na.rm = TRUE)
    prop_k  <- round(n_k / nrow(retro_out), 4)
    
    row <- data.frame(Space = space_label, Model = model_label,
                      G_opt = G_fit, Component = k,
                      n_tracks = n_k, proportion = prop_k)
    for (v in vars) {
      row[[paste0("mean_", v)]] <- round(means_k[v], 4)
      row[[paste0("sd_",   v)]] <- round(sds_k[v],   4)
    }
    row
  })
  
  desc_out  <- do.call(rbind, desc_rows)
  fname_desc <- paste0("gmm_cluster_descriptives_", model_label,
                       "_G", G_fit, "_Space", space_label, "_newROIs.csv")
  write_csv(desc_out, fname_desc)
  cat("Saved:", fname_desc, "\n")
  
  # Profile summary — empirical per component
  comp_assign  <- retro_out %>%
    select(ROIID_global, ID, GMM_component) %>%
    rename(component = GMM_component)
  
  profile_data <- tracks_df %>%
    inner_join(comp_assign, by = c("ROIID_global", "ID"))
  
  profile_vars <- c("log_Speed", "Track_Straightness", "alpha_30pct",
                    "TrackDuration", "MaxDis", "TrackLength", "log_TrackLength")
  
  prof_rows <- lapply(seq_len(G_fit), function(k) {
    d   <- profile_data %>% filter(component == k)
    row <- data.frame(Space = space_label, Model = model_label,
                      G_opt = G_fit, Component = k, n_tracks = nrow(d))
    for (v in profile_vars) {
      if (v %in% names(d)) {
        vals <- d[[v]][!is.na(d[[v]])]
        row[[paste0(v, "_mean")]]   <- round(mean(vals),   4)
        row[[paste0(v, "_sd")]]     <- round(sd(vals),     4)
        row[[paste0(v, "_median")]] <- round(median(vals), 4)
        row[[paste0(v, "_IQR")]]    <- round(IQR(vals),    4)
      }
    }
    row
  })
  
  prof_out   <- do.call(rbind, prof_rows)
  fname_prof <- paste0("gmm_profile_summary_", model_label,
                       "_G", G_fit, "_Space", space_label, "_newROIs.csv")
  write_csv(prof_out, fname_prof)
  cat("Saved:", fname_prof, "\n")
  
  invisible(list(descriptives = desc_out, profile = prof_out))
}

# Determine tracks to profile: new ROIs only, or combined if available
if (append_to_existing && exists("tracks_biv")) {
  cat("[6f] append_to_existing = TRUE and tracks_biv in session — profiling combined dataset.\n")
  tracks_for_retro_profile <- bind_rows(tracks_biv, tracks_new_biv) %>%
    distinct(ROIID_global, ID, .keep_all = TRUE)
} else {
  if (append_to_existing)
    cat("[6f] append_to_existing = TRUE but tracks_biv not in session — profiling new ROIs only.\n")
  tracks_for_retro_profile <- tracks_new_biv
}

export_retro_descriptives(retro_space, retro_rds_B, "B", tracks_for_retro_profile)

} else {
  cat("\n[6f] run_retro_predict = FALSE — retrospective prediction skipped.\n")
  cat("  Set run_retro_predict <- TRUE in Section 1 once new ROI tracks are ready.\n")
  cat("  Set new_tracks_file and retro_rds_A in Section 1.\n")
  cat("  Set append_to_existing <- TRUE to also write a combined output file.\n")
}



# =============================================================================
# [OPTIONAL] UNIVARIATE STEP 4: SUB-CLASSIFICATION WITHIN COMPONENTS
# =============================================================================
# Purpose: Apply hard population labels (Pop1, Pop2, ...) derived from the
#          UNIVARIATE GMM on log_Speed (Section 4). Intended for sub-classifying
#          tracks WITHIN an identified 3D GMM component — e.g. splitting a
#          "motile" component into slow-motile vs fast-motile using speed alone.
#
# To activate:
#   1. Define  assign_populations <- TRUE  in Section 1
#   2. Review population_analysis.pdf (Section 8) and population_gmm_summary.csv
#      to confirm G_optimal >= 2 and components separate biologically
#   3. Uncomment the code block below and re-run
#
# Dependencies (computed in Sections 4–5, always available without re-running):
#   G_optimal, comp_params, tracks_complete, gmm_fit, assignment_threshold
#
# if (assign_populations) {
#
#   if (G_optimal < 2) {
#     cat("\n[Step 4] G_optimal = 1: only one component found.\n")
#     cat("         Population assignment skipped.\n")
#
#   } else {
#
#     cat("\n=== STEP 4: POPULATION ASSIGNMENT ===\n")
#     cat("G_optimal =", G_optimal, "| Threshold:", assignment_threshold, "\n\n")
#
#     # Map components to population labels ordered by mean log_Speed (slowest = Pop1)
#     comp_order <- comp_params %>%
#       arrange(Mean_logSpeed) %>%
#       mutate(PopLabel = paste0("Pop", seq_len(n())))
#
#     comp_map <- setNames(comp_order$PopLabel, comp_order$Component)
#
#     tracks_complete <- tracks_complete %>%
#       mutate(
#         Population = ifelse(
#           GMM_max_prob >= assignment_threshold,
#           comp_map[as.character(GMM_component)],
#           "Uncertain"
#         ),
#         Population = factor(Population,
#                             levels = c(comp_order$PopLabel, "Uncertain"))
#       )
#
#     cat("Population assignment counts:\n")
#     print(table(tracks_complete$Population, useNA = "ifany"))
#     cat("Proportion uncertain:",
#         round(mean(tracks_complete$Population == "Uncertain") * 100, 1), "%\n")
#
#     pop_join <- tracks_complete %>%
#       select(Filepath, Filename, ID, Population, GMM_component,
#              GMM_max_prob, starts_with("P_comp"))
#
#     tracks_out <- tracks %>%
#       left_join(pop_join, by = c("Filepath", "Filename", "ID"))
#
#     write_csv(tracks_out, "tracks_clean.csv")
#     cat("\ntracks_clean.csv overwritten with Population column.\n")
#     cat("Columns added: Population, GMM_component, GMM_max_prob, P_comp*\n")
#   }
#
# } else {
#   cat("\n[Step 4] assign_populations = FALSE — population assignment skipped.\n")
# }


# --- 7. PLOTS -----------------------------------------------------------------
# ACTIVE plots (built and saved in Section 8):
#   p5b_scatter : pairwise feature scatter plots (Mouse-coloured)
#   p6_ratio_raw / p7_ratio_log : TrackLength/MaxDis ratio plots
#
# COMMENTED OUT — univariate GMM plots (built from Sections 3–5).
# These relate to the log_Speed GMM and were used in early-stage exploration.
# They are no longer part of the routine workflow.
# To reactivate all four:
#   1. Uncomment the four blocks below (p1_rug through p4_profile)
#   2. Define  assign_populations <- FALSE  in Section 1
#   3. Uncomment the p5_uncertainty block further below
#   4. Uncomment the population_analysis.pdf block in Section 8


# -- [INACTIVE] Page 1: log_Speed distribution + dip test annotation -----------
# Requires: speed_vec, dip_label (Sections 3–4)

# p1_rug <- ggplot(data.frame(x = speed_vec), aes(x = x)) +
#   geom_histogram(aes(y = after_stat(density)),
#                  bins = 60, fill = "steelblue", alpha = 0.45, colour = NA) +
#   geom_density(colour = "steelblue4", linewidth = 0.9, adjust = 1.2) +
#   geom_rug(aes(x = x), alpha = 0.03, colour = "steelblue4", linewidth = 0.3) +
#   annotate("text", x = Inf, y = Inf, label = dip_label,
#            hjust = 1.05, vjust = 1.5, size = 3.5, colour = "grey30") +
#   labs(title    = "Track-level log Speed distribution (pooled, all conditions)",
#        subtitle = "Histogram + KDE + rug | Dip test annotation",
#        x = "Log Speed (µm/min)", y = "Density") +
#   theme_bw(base_size = 11) +
#   theme(plot.title       = element_text(face = "bold", size = 12),
#         plot.subtitle    = element_text(size = 10, colour = "grey40"),
#         panel.grid.minor = element_blank())


# -- [INACTIVE] Page 2: Univariate GMM BIC selection curve --------------------
# Requires: bic_df, G_optimal, G_max (Section 4)

# p2_bic <- ggplot(bic_df, aes(x = G, y = BIC)) +
#   geom_line(colour = "grey50", linewidth = 0.8) +
#   geom_point(aes(colour = Selected, size = Selected)) +
#   scale_colour_manual(values = c("FALSE" = "grey50", "TRUE" = "#E63946"),
#                       guide = "none") +
#   scale_size_manual(values = c("FALSE" = 2.5, "TRUE" = 4.5), guide = "none") +
#   geom_vline(xintercept = G_optimal, linetype = "dashed",
#              colour = "#E63946", linewidth = 0.7) +
#   annotate("text", x = G_optimal + 0.15, y = min(bic_df$BIC, na.rm = TRUE),
#            label = paste("G =", G_optimal, "(optimal)"),
#            hjust = 0, size = 3.5, colour = "#E63946") +
#   scale_x_continuous(breaks = 1:G_max) +
#   labs(title    = "GMM BIC by number of components",
#        subtitle = "Model type V (variable variance) | Higher BIC = better",
#        x = "Number of components (G)", y = "mclust BIC") +
#   theme_bw(base_size = 11) +
#   theme(plot.title       = element_text(face = "bold", size = 12),
#         plot.subtitle    = element_text(size = 10, colour = "grey40"),
#         panel.grid.minor = element_blank())


# -- [INACTIVE] Page 3: GMM component density overlay on histogram -------------
# Requires: speed_vec, gmm_fit, G_optimal, comp_colours (Sections 4–5)
# Uses lapply + do.call(rbind) per pipeline convention (no map_dfr).

# x_grid <- seq(min(speed_vec) - 0.5, max(speed_vec) + 0.5, length.out = 500)
#
# comp_densities <- do.call(rbind, lapply(seq_len(G_optimal), function(g) {
#   mu  <- gmm_fit$parameters$mean[g]
#   sig <- sqrt(gmm_fit$parameters$variance$sigmasq[g])
#   pro <- gmm_fit$parameters$pro[g]
#   data.frame(
#     x         = x_grid,
#     density   = pro * dnorm(x_grid, mean = mu, sd = sig),
#     Component = as.character(g)
#   )
# }))
#
# mixture_density <- comp_densities %>%
#   group_by(x) %>%
#   summarise(density = sum(density), .groups = "drop")
#
# p3_gmm <- ggplot() +
#   geom_histogram(data = data.frame(x = speed_vec),
#                  aes(x = x, y = after_stat(density)),
#                  bins = 60, fill = "grey80", alpha = 0.5, colour = NA) +
#   geom_line(data = comp_densities,
#             aes(x = x, y = density, colour = Component), linewidth = 0.9) +
#   geom_line(data = mixture_density,
#             aes(x = x, y = density),
#             colour = "black", linewidth = 1.1, linetype = "dashed") +
#   scale_colour_manual(values = comp_colours,
#                       name   = "Component",
#                       labels = paste0("Comp ", seq_len(G_optimal),
#                                       " (µ=", round(exp(gmm_fit$parameters$mean), 1),
#                                       " µm/min)")) +
#   annotate("text", x = Inf, y = Inf,
#            label = paste("G =", G_optimal, "components (BIC-optimal)"),
#            hjust = 1.05, vjust = 1.5, size = 3.5) +
#   labs(title    = "GMM component decomposition of log Speed",
#        subtitle = "Coloured lines = individual components | Dashed = mixture density",
#        x = "Log Speed (µm/min)", y = "Density") +
#   theme_bw(base_size = 11) +
#   theme(plot.title       = element_text(face = "bold", size = 12),
#         plot.subtitle    = element_text(size = 10, colour = "grey40"),
#         legend.position  = "right",
#         panel.grid.minor = element_blank())


# -- [INACTIVE] Page 4: Univariate GMM component profile boxplots -------------
# Requires: tracks_complete (Section 4), comp_colours

# profile_plot_data <- tracks_complete %>% mutate(Component = factor(GMM_component))
#
# make_profile_box <- function(col, ylabel, title) {
#   ggplot(profile_plot_data,
#          aes(x = Component, y = .data[[col]], fill = Component)) +
#     geom_violin(alpha = 0.4, colour = NA, trim = FALSE) +
#     geom_boxplot(width = 0.25, outlier.size = 0.4,
#                  outlier.alpha = 0.3, colour = "grey30") +
#     scale_fill_manual(values = comp_colours, guide = "none") +
#     labs(title = title, x = "GMM Component", y = ylabel) +
#     theme_bw(base_size = 11) +
#     theme(plot.title       = element_text(face = "bold", size = 11),
#           panel.grid.minor = element_blank())
# }
#
# p4a <- make_profile_box("Straightness", "Straightness (0-1)",    "Straightness by component")
# p4b <- make_profile_box("MaxDis",       "Max Displacement (µm)", "Max Displacement by component")
# p4c <- make_profile_box("TrackLength",  "Track Length (µm)",     "Track Length by component")
#
# p4_profile <- (p4a | p4b | p4c) +
#   plot_annotation(
#     title    = "Step 3: Feature profiling by univariate GMM component",
#     subtitle = "Does speed-based component membership separate on other features?",
#     theme    = theme(plot.title    = element_text(face = "bold", size = 13),
#                      plot.subtitle = element_text(size = 10, colour = "grey40"))
#   )


# -- Page 5 [CONDITIONAL]: Assignment uncertainty histogram --------------------
# Requires assign_populations to be defined and the Step 4 block above to have run.
# To reactivate:
#   1. Define  assign_populations <- FALSE  (or TRUE) in Section 1
#   2. Run the optional Step 4 block
#   3. Uncomment the block below
#   4. Uncomment the corresponding print line in the population_analysis.pdf block

# if (assign_populations && G_optimal >= 2) {
#   p5_uncertainty <- ggplot(tracks_complete,
#                            aes(x = GMM_max_prob, fill = factor(GMM_component))) +
#     geom_histogram(bins = 50, colour = NA, alpha = 0.7, position = "identity") +
#     geom_vline(xintercept = assignment_threshold,
#                linetype = "dashed", colour = "black", linewidth = 0.8) +
#     annotate("text", x = assignment_threshold + 0.01, y = Inf,
#              label = paste("Threshold:", assignment_threshold),
#              hjust = 0, vjust = 1.5, size = 3.5) +
#     scale_fill_manual(values = comp_colours, name = "Component") +
#     labs(title    = "Step 4: Assignment uncertainty (max posterior probability)",
#          subtitle = "Tracks left of threshold are labelled Uncertain",
#          x = "Max posterior probability", y = "Count") +
#     theme_bw(base_size = 11) +
#     theme(plot.title       = element_text(face = "bold", size = 12),
#           plot.subtitle    = element_text(size = 10, colour = "grey40"),
#           panel.grid.minor = element_blank())
# }


# -- Page 5: Pairwise feature scatter plots (Mouse-coloured) -------------------

p5b_scatter <- (p_s_l | p_s_d | p_s_st) / (p_l_d | p_l_st | p_d_st) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Pairwise feature scatter plots",
    subtitle = "Points coloured by Mouse | Contours show density | Log scale for Speed, Length, MaxDis",
    theme    = theme(plot.title    = element_text(face = "bold", size = 13),
                     plot.subtitle = element_text(size = 10, colour = "grey40"))
  )


# -- Pages 6–7: TrackLength/MaxDis ratio plots ---------------------------------

ratio_99 <- quantile(ratio_data_all$Ratio,          0.99, na.rm = TRUE)
speed_99 <- quantile(ratio_data_all$SpeedMean_ummin, 0.99, na.rm = TRUE)
cat("Ratio 99th percentile:", round(ratio_99, 1),
    "| SpeedMean 99th percentile:", round(speed_99, 1), "µm/min\n")

p6a_raw <- make_ratio_panel(ratio_data_all,      xvar = "MaxDis",
                            xlabel = "MaxDis (µm)", title = "All tracks (unfiltered)")
p6b_raw <- make_ratio_panel(ratio_data_filtered, xvar = "MaxDis",
                            xlabel = "MaxDis (µm)",
                            title  = paste0("Filtered: MaxDis ≥ ", maxdis_threshold, " µm"))

p6_ratio_raw <- (p6a_raw | p6b_raw) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "TrackLength / MaxDis ratio vs Max Displacement (raw scale)",
    subtitle = paste0("Ratio = 1: linear | Ratio > 1: tortuous",
                      "  |  Non-translocating tracks removed in right panel",
                      " (MaxDis < ", maxdis_threshold, " µm)"),
    theme = theme(plot.title    = element_text(face = "bold", size = 12),
                  plot.subtitle = element_text(size = 8, colour = "grey40"))
  )

p7a_log <- make_ratio_panel(ratio_data_all,      xvar = "log_MaxDis",
                            xlabel = "Log Max Displacement", title = "All tracks (unfiltered)")
p7b_log <- make_ratio_panel(ratio_data_filtered, xvar = "log_MaxDis",
                            xlabel = "Log Max Displacement",
                            title  = paste0("Filtered: MaxDis ≥ ", maxdis_threshold, " µm"))

p7_ratio_log <- (p7a_log | p7b_log) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "TrackLength / MaxDis ratio vs Max Displacement (log scale)",
    subtitle = paste0("Ratio = 1: linear | Ratio > 1: tortuous",
                      "  |  Non-translocating tracks removed in right panel",
                      " (MaxDis < ", maxdis_threshold, " µm)"),
    theme = theme(plot.title    = element_text(face = "bold", size = 12),
                  plot.subtitle = element_text(size = 8, colour = "grey40"))
  )


# --- 8. SAVE PLOTS ------------------------------------------------------------

# -- [COMMENTED OUT] population_analysis.pdf -----------------------------------
# Contains univariate GMM plots: p1_rug, p2_bic, p3_gmm (Sections 3–4),
# p4_profile (Section 5), p5b_scatter, p6_ratio_raw, p7_ratio_log.
# These are not part of the routine 3D GMM workflow.
#
# To reactivate:
#   1. Define  assign_populations <- FALSE  in Section 1
#   2. Uncomment the p1_rug through p4_profile blocks in Section 7
#   3. Uncomment the p5_uncertainty block in Section 7
#   4. Uncomment this entire block

# pdf("population_analysis.pdf", width = 11, height = 6)
# print(p1_rug)
# print(p2_bic)
# print(p3_gmm)
# print(p4_profile)
# print(p5b_scatter)
# print(p6_ratio_raw)
# print(p7_ratio_log)
# if (assign_populations && G_optimal >= 2) print(p5_uncertainty)
# dev.off()
# cat("\nSaved: population_analysis.pdf\n")


# -- TrackLength/MaxDis ratio plots --------------------------------------------

pdf("population_analysis_TrLen_MaxDis_ratio_vsMaxDis.pdf", width = 11, height = 6)
print(p6_ratio_raw)
print(p7_ratio_log)
dev.off()
cat("\nSaved: population_analysis_TrLen_MaxDis_ratio_vsMaxDis.pdf\n")


# ---- Bivariate GMM -------------------------------------------------------------

pdf("population_analysis_bivariate_GMM.pdf", width = 11, height = 6)
print(p_biv_A)
if (space_B_available) print(p_biv_B)
dev.off()
cat("\nSaved: population_analysis_bivariate_GMM.pdf\n")


# ----Trivariate (3D) GMM -------------------------------------------------------
# Filename is dynamic: encodes model and G for traceability across force_model runs.

pdf(paste0("population_analysis_trivariate_GMM_SpaceAB_",
           gmm_3D_A$best_model, "_G", gmm_3D_A$G_opt, ".pdf"),
    width = 14, height = 10)
print(p_3D_A)
if (space_B_available) print(p_3D_B)
dev.off()
cat("Saved: population_analysis_trivariate_GMM_SpaceAB_",
    gmm_3D_A$best_model, "_G", gmm_3D_A$G_opt, ".pdf\n", sep = "")


# --- 9. SAVE NUMERICAL OUTPUTS ------------------------------------------------
# Saves univariate GMM component parameters and feature profile summary.
#
# IMPORTANT: comp_params and profile_summary describe the UNIVARIATE log_Speed
# GMM (Sections 4–5) — they are NOT from the 3D GMM. The 3D GMM component
# assignments are already exported to CSV in Section 6e.
#
# population_gmm_summary.csv    : per-component mean, SD, proportion, geometric
#                                  speed (µm/min) from the univariate GMM
# population_profile_summary.csv: per-component means of Straightness, MaxDis,
#                                  TrackLength stratified by univariate component

write_csv(comp_params,     "population_gmm_summary.csv")
write_csv(profile_summary, "population_profile_summary.csv")
cat("Saved: population_gmm_summary.csv (univariate GMM component parameters)\n")
cat("Saved: population_profile_summary.csv (univariate GMM feature profiles)\n")

cat("\n--- Run summary ---\n")
cat("Dip test D =", round(dip_result$statistic, 5),
    "| p =", round(dip_result$p.value, 5), "\n")
cat("Univariate GMM components:", G_optimal, "\n")
cat("3D GMM Space A: G =", gmm_3D_A$G_opt, "| Model =", gmm_3D_A$best_model, "\n")
if (space_B_available)
  cat("3D GMM Space B: G =", gmm_3D_B$G_opt, "| Model =", gmm_3D_B$best_model, "\n")
cat("\nNext steps:\n")
cat("  1. Review population_analysis_bivariate_GMM.pdf\n")
cat("  2. Review population_analysis_trivariate_GMM_SpaceAB_",
    gmm_3D_A$best_model, "_G", gmm_3D_A$G_opt, ".pdf\n", sep = "")
cat("  3. Review component_profiles_", model_label, ".pdf\n", sep = "")
cat("  4. Set noise_comp_3D_A / _B and flag_noise_3D <- TRUE to flag noise\n")
cat("  5. Once new ROI tracks ready: set run_retro_predict <- TRUE (Section 6f)\n")
cat("  6. Proceed to Scripts 03/04 for Treatment x Day mixed models\n")


# ---- 10. JOIN GMM COMPONENTS TO TRACKS ----------------------------------------
# Left-joins GMM component assignments (Space A and/or B) onto the most recent
# tracks_clean CSV and writes a new output file for use in Scripts 03/04 and
# track_visualisation.R.
# Files to join are set in Section 1 (gmm_join_file_A, gmm_join_file_B).
# Set either to NA to skip that space.
# Match key: ROIID_global + ID (consistent with pipeline convention).
#
# Output columns added per space:
#   GMM_component_A, GMM_max_prob_A, P_comp1_A ... P_comp<G>_A
#   GMM_component_B, GMM_max_prob_B, P_comp1_B ... P_comp<G>_B
# Tracks with no match get NA in those columns.

cat("\n=== SECTION 10: JOIN GMM COMPONENTS TO TRACKS ===\n")

tracks_join <- tracks   # start from the loaded tracks

join_gmm <- function(df, gmm_file) {
  if (is.na(gmm_file)) {
    cat("gmm_join_file set to NA — skipped.\n")
    return(df)
  }
  if (!file.exists(gmm_file)) {
    cat("File not found:", gmm_file, "— skipped.\n")
    return(df)
  }
  gmm <- read_csv(gmm_file, show_col_types = FALSE) %>%
    select(-any_of(c("MouseID", "Day")))
  df_out <- df %>%
    left_join(gmm, by = c("ROIID_global", "ID"))
  n_matched <- sum(!is.na(df_out$GMM_component))
  unmatched <- df_out %>% filter(is.na(GMM_component))
  if (nrow(unmatched) > 0) {
    cat("Unmatched tracks (", nrow(unmatched), "):\n", sep = "")
    print(unmatched %>% select(ROIID_global, ID, MouseID, Day) %>% as.data.frame())
  }
  cat("Joined:", n_matched, "/", nrow(df_out), "tracks |", gmm_file, "\n")
  df_out
}

tracks_join <- join_gmm(tracks_join, gmm_join_file)


out_file <- paste0(tools::file_path_sans_ext(input_file), "_GMM.csv")
write_csv(tracks_join, out_file)
cat("Saved:", out_file, "\n")
cat("Columns added:",
    paste(setdiff(names(tracks_join), names(tracks)), collapse = ", "), "\n")

