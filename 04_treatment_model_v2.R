# =============================================================================
# CELL TRACKING ANALYSIS - STEP 4: TREATMENT AND DAY MODEL
# =============================================================================
# Purpose: Test whether motility parameters differ by Treatment (Mock/RP),
#          Day (D1/D2), and their interaction, with Location as a covariate.
#
# Key distinction from Script 03:
#   - Treatment * Day is the primary fixed effect of interest
#   - Location is a nuisance covariate (adjusted for, not tested)
#   - EMMs focus on Treatment x Day combinations
#   - Interaction test answers: does the RP effect differ between D1 and D2?
#
# Models fitted for each variable:
#   ROI level:   variable ~ Treatment * Day + Location +
#                           (1|MouseID/TumourID_global/SliceID_global)
#   Track level: variable ~ Treatment * Day + Location +
#                           (1|MouseID/TumourID_global/SliceID_global/ROIID_global)
#
# Random effects note:
#   Primary:     random intercept only (1|...)
#   Sensitivity: random slope for Treatment commented out below —
#                asks whether the RP response varies across tumours.
#                Likely singular with n=5; enable and check isSingular().
#
# DF approximation:
#   Default: Kenward-Roger (more accurate for small n ROI level).
#   Reference: Kenward & Roger (1997) Biometrics 53:983-997
#
# Variables: Speed (raw + log), TrackLength (log only), Straightness (raw + log)
#
# Outputs (all suffixed with output_suffix):
#   treatment_model_summaries_INTERACTION_<suffix>.txt : full model summaries + R2
#   treatment_emms_INTERACTION_<suffix>.csv            : EMMs per Treatment x Day
#   treatment_pairwise_INTERACTION_<suffix>.csv        : pairwise contrasts (Tukey)
#   treatment_r2_INTERACTION_<suffix>.csv              : R2 marginal/conditional
#   treatment_assumption_plots_<suffix>.pdf            : residual + QQ plots
#   treatment_emm_plots_interaction_<suffix>.pdf       : EMM plots + forest plots
#   treatment_effect_sizes_INTERACTION_<suffix>.csv    : Cohen's d per effect
#   treatment_geometric_means_INTERACTION_<suffix>.csv : back-transformed group means
#   treatment_fold_changes_INTERACTION_<suffix>.csv    : fold changes from log contrasts
#   treatment_model_comparisons_<suffix>.csv           : interaction + slope LRT/BIC
#   treatment_exploratory_plots_<suffix>.pdf           : violin + jitter plots
#   ks_track_results_<suffix>.csv                      : KS tests (track level only)
#   per_mouse_speed_summary_<suffix>.csv               : per-mouse mean summary
#   per_mouse_direction_check_<suffix>.csv             : RP vs Mock fold changes
#   per_mouse_trajectory_plots_<suffix>.pdf            : D1→D2 trajectories
#   per_mouse_facet_trajectories_<suffix>.pdf          : per-mouse 1x3 facet plots
#   per_mouse_facet_trajectories_normalised_<suffix>.pdf : normalised trajectories
#
# Inputs: tracks_clean_v2.csv, roi_summary_v2.csv (from 02e_batch_revised.R)
#         GMM_component column required if gmm_keep != c() (from 02d)
# =============================================================================


# --- PARAMETERS ---------------------------------------------------------------

data_dir          <- "Z:\\lightmicro-inst\\scharalambous\\Resident_cell_motility_tr_v3\\ROI files"  # <-- UPDATE PATH
tracks_clean_path <- "Z:\\lightmicro-inst\\scharalambous\\Resident_cell_motility_tr_v3\\ROI files\\Complete data - retrospective fitting G=5 VEE (from v2)\\tracks_clean_v3_alldata_GMM_retrospective.csv"
roi_summary_path  <- "roi_summary_v3_GMM_1_3_4_5.csv"


# GMM filter
# Set to c() for full dataset; e.g. c(1, 2) to keep components 1 and 2.
# Requires GMM_component column (run 02d first).
gmm_keep <- c(1,3,4,5)


# Output suffix — auto-constructed
output_suffix <- if (length(gmm_keep) == 0) "all" else
                 paste0("GMM", paste(gmm_keep, collapse = "_"))


# ROI-level analysis gate
# Set TRUE only when gmm_keep is finalised AND roi_summary matches it.
# When FALSE: track-level models, KS tests run; ROI models, per-mouse and
# trajectory plots are skipped.
roi_analysis_enabled <- TRUE

# --- Variable toggle ----------------------------------------------------------
# Mirrors var_include in scripts 02 and 03. Controls which variables enter
# models, comparisons, exploratory plots and KS tests.
var_include <- c(
  "Speed"        = TRUE,
  "TrackLength"  = TRUE,
  "Straightness" = TRUE,
  "MaxDis"       = FALSE,
  "TurningAngle" = TRUE,
  "Alpha"        = TRUE
)

# MaxDis displacement filter
# Applied after GMM filter within retained components.
# Removes near-stationary tracks that never meaningfully displaced
# from their starting position regardless of GMM component assignment.
# Threshold: MaxDis > 1.1 µm (equivalent to log_MaxDis > 0.095)
maxdis_threshold <- 0  # µm

# Alpha R2 threshold — must match script 02
alpha_r2_threshold <- 0.5

# --- Variable-specific random formula overrides -------------------------------
# Variables listed here use their own random formula instead of the level
# default below. Use when ICC results justify a richer structure for a
# specific variable only (e.g. log Speed shows non-negligible Tumour-level
# variance while others approach singularity with TumourID included).
#
# ROI default:   MouseID/SliceID_global          (TumourID dropped)
# Track default: MouseID/SliceID_global/ROIID_global
#
# Speed overrides to full Mouse/Tumour/Slice[/ROI] to capture tumour-level
# variance identified in ICC (script 02) and support tumour-size analyses.
random_overrides_roi <- list(
  "Speed" = "(1 | MouseID/TumourID_global/SliceID_global)"
)
random_overrides_track <- list(
  "Speed" = "(1 | MouseID/TumourID_global/SliceID_global/ROIID_global)"
)


# --- 0. PACKAGES --------------------------------------------------------------
# install.packages(c("lme4", "lmerTest", "emmeans", "effectsize",
#                    "performance", "tidyverse", "patchwork", "pbkrtest"))
library(pbkrtest)
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(effectsize)
library(performance)
library(patchwork)

# # Allow Kenward-Roger emmeans with large n (set above your track n)
# emm_options(pbkrtest.limit = 7000,
#             lmerTest.limit = 7000)


# --- 1. LOAD DATA -------------------------------------------------------------

setwd(data_dir)

tracks <- read_csv(tracks_clean_path, show_col_types = FALSE) %>%
  mutate(
    MouseID         = factor(MouseID),
    TumourID_global = factor(TumourID_global),
    SliceID_global  = factor(SliceID_global),
    ROIID_global    = factor(ROIID_global),
    Location        = factor(Location,  levels = c("Core", "CoreEdge", "Edge")),
    Treatment       = factor(Treatment, levels = c("Mock", "RP")),
    Day             = factor(Day,       levels = c("D1", "D2"))
  )

roi <- read_csv(roi_summary_path, show_col_types = FALSE) %>%
  mutate(
    MouseID         = factor(MouseID),
    TumourID_global = factor(TumourID_global),
    SliceID_global  = factor(SliceID_global),
    ROIID_global    = factor(ROIID_global),
    Location        = factor(Location,  levels = c("Core", "CoreEdge", "Edge")),
    Treatment       = factor(Treatment, levels = c("Mock", "RP")),
    Day             = factor(Day,       levels = c("D1", "D2"))
  )

cat("Track level loaded:", nrow(tracks), "tracks\n")
cat("ROI level loaded:  ", nrow(roi),    "ROIs\n")
cat("Output suffix:     ", output_suffix, "\n\n")


# --- 1b. GMM FILTER -----------------------------------------------------------

if (length(gmm_keep) > 0) {
  if (!"GMM_component" %in% names(tracks)) {
    stop("GMM_component column not found in tracks — run script 02d first.")
  }
  n_before <- nrow(tracks)
  tracks   <- tracks %>% filter(GMM_component %in% gmm_keep)
  cat("GMM filter applied: kept components", paste(gmm_keep, collapse = ", "),
      "->", nrow(tracks), "tracks retained (", n_before - nrow(tracks), "removed)\n\n")
} else {
  cat("No GMM filter applied (gmm_keep = c()). Using full dataset.\n\n")
}

# --- 1b-ii. MAXDIS FILTER -----------------------------------------------------
# Secondary filter applied within retained GMM components.
# Removes tracks where maximum displacement from origin <= 1.1 µm,
# indicating near-stationary behaviour regardless of GMM assignment.
# No track duration limit is set, so this acts as the displacement-based
# quality gate to remove tracks with no meaningful productive displacement.
# Full exclusion log saved for transparency and methods reporting.

if (!"MaxDis" %in% names(tracks)) {
  warning("MaxDis column not found — MaxDis filter skipped.")
} else {
  n_before_maxdis        <- nrow(tracks)
  tracks_excluded_maxdis <- tracks %>% filter(MaxDis <= maxdis_threshold)
  tracks                 <- tracks %>% filter(MaxDis >  maxdis_threshold)
  
  cat("MaxDis filter: MaxDis >", maxdis_threshold, "µm\n")
  cat("  Removed: ", nrow(tracks_excluded_maxdis),
      " (", round(100 * nrow(tracks_excluded_maxdis) / n_before_maxdis, 1), "%)\n",
      sep = "")
  cat("  Retained:", nrow(tracks), "\n\n")
  
  # Exclusion log — key columns for reporting and auditing
  write_csv(
    tracks_excluded_maxdis %>%
      select(any_of(c(
        "ROIID_global", "SliceID_global", "TumourID_global", "MouseID",
        "Location", "Treatment", "Day", "GMM_component",
        "MaxDis", "log_MaxDis", "SpeedMean_ummin",
        "Track_Length", "Track_Straightness", "TrackDuration"
      ))),
    paste0("excluded_tracks_maxdis_", output_suffix, ".csv")
  )
  cat("Saved: excluded_tracks_maxdis_", output_suffix, ".csv\n\n", sep = "")
}


# Complete-case filter
# Removes tracks with NA in any core variable column.
# Alpha is excluded from this filter — it runs on its own subset via filter_col.
# The 20 excluded tracks have no turning angle due to complete immobility.
core_cols <- c("SpeedMean_ummin", "log_Speed", "log_TrackLength",
               "Track_Straightness", "log_Straightness",
               "MaxDis", "log_MaxDis", "mean_turning_angle")

# --- 1c. COMPLETE-CASE FILTER -------------------------------------------------
n_before <- nrow(tracks)
tracks   <- tracks %>% filter(if_all(any_of(core_cols), ~ !is.na(.)))
cat("Complete-case filter: ", n_before - nrow(tracks),
    "tracks removed (NA in core columns)\n")
cat("Tracks retained:", nrow(tracks), "\n\n")

# Alpha bounds filter — remove physically impossible values
# Applied before model fitting; Alpha must be in (0, 2.25]
n_alpha_before <- sum(!is.na(tracks$alpha_30pct))
tracks <- tracks %>%
  mutate(alpha_30pct = ifelse(alpha_30pct <= 0 | alpha_30pct > 2.25,
                              NA_real_, alpha_30pct))
cat("Alpha bounds filter: ",
    n_alpha_before - sum(!is.na(tracks$alpha_30pct)),
    "values set to NA (outside 0-2 range)\n\n")


# --- 1d. GMM COMPOSITION BY LOCATION ------------------------------------------
# Quantifies the proportion of each GMM component within each
# Location x Treatment x Day combination.
# Interpretation guide:
#   Compare Mock vs RP within each Location and Day:
#   - If imbalance present on D1: likely pre-existing state confound
#   - If imbalance emerges by D2: consistent with treatment-driven composition shift
#   - If spatially uniform: Location covariate in mixed model absorbs confounding
# Note: requires GMM_component column — if gmm_keep = c() this reflects
#       full unfiltered composition; if gmm_keep != c() reflects retained components only.
# =============================================================================
#  GMM COMPOSITION WITH MOTILE SUBCLASSIFICATION
# =============================================================================
# Components 1 and 3 are split into:
#   "1"        — all tracks in component 1 (MaxDis <= 7 µm, low displacement)
#   "1 motile" — component 1 tracks with MaxDis > 7 µm (productive displacement)
#   "3"        — all tracks in component 3 (MaxDis <= 7 µm, low displacement)
#   "3 motile" — component 3 tracks with MaxDis > 7 µm
# Components 2, 4, 5 are unchanged.
# Low displacement sits below high displacement in the stacked bar.
# Stacking order: 1 | 1 motile | 2 | 3 | 3 motile | 4 | 5
# =============================================================================

if ("GMM_component" %in% names(tracks)) {
  
  # Reclassify components 1 and 3 by MaxDis threshold
  maxdis_motile_threshold <- 7  # µm — tracks above this in comp 1/3 are "motile"
  
  tracks_comp_plot <- tracks %>%
    mutate(GMM_component_plot = case_when(
      GMM_component == 1 & MaxDis >  maxdis_motile_threshold ~ "1 motile",
      GMM_component == 1 & MaxDis <= maxdis_motile_threshold ~ "1",
      GMM_component == 3 & MaxDis >  maxdis_motile_threshold ~ "3 motile",
      GMM_component == 3 & MaxDis <= maxdis_motile_threshold ~ "3",
      TRUE ~ as.character(GMM_component)
    ))
  
  # Report split
  cat("\n=== GMM COMPOSITION — MOTILE SUBCLASSIFICATION ===\n")
  cat("MaxDis threshold for motile subclass: >", maxdis_motile_threshold, "µm\n\n")
  
  tracks_comp_plot %>%
    count(GMM_component, GMM_component_plot) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
    print(n = Inf)
  cat("\n")
  
  # Compute composition
  # Stacking order: low-displacement below, motile above within each component
  component_levels <- c("1", "1 motile", "2", "3", "3 motile", "4", "5")
  
  gmm_composition <- tracks_comp_plot %>%
    group_by(ROIID_global,Location,Treatment, GMM_component_plot) %>%    # Location,Treatment, Day, GMM_component_plot| ROIID_global,Location,Treatment, GMM_component_plot
    summarise(n = n(), .groups = "drop") %>%
    group_by(ROIID_global,Location, Treatment) %>%          #Location, Treatment, Day
    mutate(
      pct         = round(100 * n / sum(n), 1),
      n_condition = sum(n),
      motile_n    = sum(n[GMM_component_plot %in% c("1 motile", "3 motile")]),
      motile_pct  = round(100 * motile_n / sum(n), 1)
    ) %>%
    ungroup() %>%
    rename(GMM_component = GMM_component_plot) %>%
    mutate(GMM_component = factor(GMM_component,
                                  levels = component_levels)) %>%
    arrange(ROIID_global,Location,Treatment, GMM_component)
  
  cat("=== GMM COMPOSITION BY TREATMENT AND DAY ===\n")
  cat("Rows show % of each GMM subclass within Treatment x Day\n\n")
  print(gmm_composition, n = Inf)
  
  # # write_csv(gmm_composition,
  #           paste0("gmm_composition_by_Loc_trt_ROI_", output_suffix, "_MaxDis",
  #                  maxdis_motile_threshold, ".csv"))
  cat("\nSaved: gmm_composition_by_Loc_trt_ROI", output_suffix,
      "_MaxDis", maxdis_motile_threshold, ".csv\n\n", sep = "")
  
} else {
  cat("\nGMM_component column not found — skipping composition table.\n")
  cat("Run script 02d first if GMM analysis has been completed.\n\n")
}


  # 
  # gmm_composition <- tracks_comp_filtered %>%
  #   group_by(Location,Treatment, Day, GMM_component) %>%   #Location,
  #   summarise(n = n(), .groups = "drop") %>%
  #   group_by(Location, Treatment, Day) %>%                 #Location,
  #   mutate(
  #     pct         = round(100 * n / sum(n), 1),
  #     n_condition = sum(n)
  #   ) %>%
  #   ungroup() %>%
  #   arrange(Location,Treatment, Day, GMM_component)         #Location,
  # 
  # cat("\n=== GMM COMPOSITION BY LOCATION ===\n")
  # cat("Rows show % of each GMM component within Location x Treatment x Day\n")
  # cat("Focus: Mock vs RP within same Location and Day\n\n")
  # print(gmm_composition, n = Inf)
  # 
  # write_csv(gmm_composition,
  #           paste0("gmm_composition_by_location_", output_suffix, "_MaxDis7.csv"))    #_location_
  # cat("Saved: gmm_composition_by_location_", output_suffix, "_MaxDis7.csv\n\n", sep = "")  #_location_
  
# } else {
#   cat("\nGMM_component column not found — skipping composition table.\n")
#   cat("Run script 02d first if GMM analysis has been completed.\n\n")
# }


# --- 2. DEFINE VARIABLES ------------------------------------------------------

# Track level column names (from tracks_clean_v2.csv)
# Track_Straightness is the canonical raw column name in v2 files.
track_vars <- list(
  list(name = "Speed",        raw = "SpeedMean_ummin",   log = "log_Speed",
       raw_label = "Speed (µm/min)",        log_label = "Log Speed (µm/min)",
       run_raw = TRUE),
  list(name = "TrackLength",  raw = NULL,                log = "log_TrackLength",
       raw_label = NULL,                   log_label = "Log Track Length (µm)",
       run_raw = FALSE),
  list(name = "Straightness", raw = "Track_Straightness", log = "log_Straightness",
       raw_label = "Straightness (-)",      log_label = "Log Straightness",
       run_raw = TRUE),
  list(name = "MaxDis",       raw = "MaxDis",              log = "log_MaxDis",
       raw_label = "MaxDis (um)",            log_label = "Log MaxDis (um)",
       run_raw = TRUE),
  list(name = "TurningAngle", raw = "mean_turning_angle",  log = NULL,
       raw_label = "Mean Turning Angle (deg)", log_label = NULL,
       run_raw = TRUE),
  list(name = "Alpha",        raw = "alpha_30pct",         log = NULL,
       raw_label = "Alpha (MSD)",             log_label = NULL,
       run_raw = TRUE,
       filter_col = "msd_r2_30pct",           filter_val = alpha_r2_threshold)
)

# ROI level column names (from roi_summary_v2.csv)
roi_vars <- list(
  list(name = "Speed",        raw = "Speed_mean",        log = "log_Speed_mean",
       raw_label = "Speed (µm/min)",        log_label = "Log Speed (µm/min)",
       run_raw = TRUE),
  list(name = "TrackLength",  raw = NULL,                log = "log_Length_mean",
       raw_label = NULL,                   log_label = "Log Track Length (µm)",
       run_raw = FALSE),
  list(name = "Straightness", raw = "Straightness_mean", log = "log_Straight_mean",
       raw_label = "Straightness (-)",      log_label = "Log Straightness",
       run_raw = TRUE),
  list(name = "MaxDis",       raw = "MaxDis_mean",         log = "log_MaxDis_mean",
       raw_label = "MaxDis (um)",            log_label = "Log MaxDis (um)",
       run_raw = TRUE),
  list(name = "TurningAngle", raw = "TurningAngle_mean",   log = NULL,
       raw_label = "Mean Turning Angle (deg)", log_label = NULL,
       run_raw = TRUE),
  list(name = "Alpha",        raw = "Alpha30_filt_mean",   log = NULL,
       raw_label = "Alpha (MSD, filtered)",   log_label = NULL,
       run_raw = TRUE)
)

# Colours and shapes
treatment_colours <- c("Mock" = "#457B9D", "RP" = "#E63946")
day_shapes        <- c("D1" = 16, "D2" = 17)

location_shapes <- c(
  "Core"     = 21,
  "CoreEdge" = 22,
  "Edge"     = 24
)

location_colours <- c(
  "Core"     = "#E63946",
  "CoreEdge" = "#457B9D",
  "Edge"     = "#2A9D8F"
)


# --- 3. HELPER: FIT MODEL -----------------------------------------------------
# Fixed effects: Treatment * Day + Location
# Random effects:
#   Default ROI:   MouseID/SliceID_global (TumourID dropped — singular for most variables)
#   Default Track: MouseID/SliceID_global/ROIID_global
#   Overrides:     declared in random_overrides_roi / random_overrides_track above
#                  applied per variable inside run_models_for_level

random_formula_roi_default   <- "(1 | MouseID/SliceID_global)"
random_formula_track_default <- "(1 | MouseID/SliceID_global/ROIID_global)"

# Random slope sensitivity — uncomment to test; likely singular with n=5 tumours
# random_formula_roi_default   <- "(1 + Treatment | TumourID_global) + (1 | SliceID_global)"
# random_formula_track_default <- "(1 + Treatment | TumourID_global) + (1 | SliceID_global/ROIID_global)"

# Convenience aliases used in run_model_comparisons and Section C
random_formula_roi   <- random_formula_roi_default
random_formula_track <- random_formula_track_default


fit_treatment_model <- function(data, col, random_formula) {
  formula <- as.formula(
    paste(col, "~ Treatment * Day + Location +", random_formula)
  )
  tryCatch(
    lmer(formula, data = data, REML = TRUE),
    error   = function(e) { cat("  Model error:", e$message, "\n"); NULL },
    warning = function(w) {
      m <- lmer(formula, data = data, REML = TRUE)
      cat("  Model warning:", w$message, "\n")
      m
    }
  )
}


# --- 4. HELPER: EXTRACT EMMs AND PAIRWISE -------------------------------------
# EMMs for the Treatment x Day grid (2 x 2 = 4 combinations), averaged over Location.
# Pairwise contrasts:
#   1. Treatment within each Day (Mock vs RP | D1, Mock vs RP | D2)
#   2. Day within each Treatment (D1 vs D2 | Mock, D1 vs D2 | RP)

extract_emm_trt <- function(model, var_name, level_name, transform) {
  
  # KR for ROI level (small n); Satterthwaite for Track level (large n)
  # Reference: Kenward & Roger (1997) Biometrics 53:983-997
  df_method <- if (level_name == "ROI") "kenward-roger" else "satterthwaite"

  emm <- emmeans(model, ~ Treatment * Day, lmer.df = df_method)

  pw_trt_by_day <- pairs(emm, simple = "Treatment",
                         adjust = "tukey", infer = TRUE)
  pw_day_by_trt <- pairs(emm, simple = "Day",
                         adjust = "tukey", infer = TRUE)

  emm_raw <- as.data.frame(emm)
  names(emm_raw)[names(emm_raw) == "emmean"] <- "emmean_val"
  names(emm_raw)[grepl("lower|LCL", names(emm_raw), ignore.case = TRUE)] <- "lower_CL"
  names(emm_raw)[grepl("upper|UCL", names(emm_raw), ignore.case = TRUE)] <- "upper_CL"
  emm_df <- emm_raw %>%
    mutate(Variable = var_name, Level = level_name, Transform = transform)

  clean_pw <- function(pw_obj, contrast_type) {
    df <- as.data.frame(pw_obj)
    names(df)[grepl("lower|LCL", names(df), ignore.case = TRUE)] <- "lower_CL"
    names(df)[grepl("upper|UCL", names(df), ignore.case = TRUE)] <- "upper_CL"
    df %>% mutate(
      Variable      = var_name,
      Level         = level_name,
      Transform     = transform,
      Contrast_type = contrast_type,
      Significant   = p.value < 0.05
    )
  }

  pw_df <- bind_rows(
    clean_pw(pw_trt_by_day, "Treatment within Day"),
    clean_pw(pw_day_by_trt, "Day within Treatment")
  )

  list(
    emm      = emm_df,
    pairwise = pw_df,
    emm_obj  = emm,
    pw_trt   = pw_trt_by_day,
    pw_day   = pw_day_by_trt
  )
}


# --- 5. HELPER: ASSUMPTION PLOTS ----------------------------------------------

make_assumption_plots <- function(model, title) {

  resid_df <- data.frame(
    fitted   = fitted(model),
    residual = residuals(model)
  )

  p_resid <- ggplot(resid_df, aes(x = fitted, y = residual)) +
    geom_point(alpha = 0.3, size = 0.8) +
    geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
    geom_smooth(method = "loess", se = FALSE, colour = "blue", linewidth = 0.8) +
    labs(title = paste(title, "— Residuals vs Fitted"),
         x = "Fitted values", y = "Residuals") +
    theme_bw(base_size = 10)

  p_qq <- ggplot(resid_df, aes(sample = residual)) +
    stat_qq(alpha = 0.3, size = 0.8) +
    stat_qq_line(colour = "red") +
    labs(title = paste(title, "— QQ Plot"),
         x = "Theoretical quantiles", y = "Sample quantiles") +
    theme_bw(base_size = 10)

  p_resid | p_qq
}


# --- 6. HELPER: EMM PLOT ------------------------------------------------------
# Layout: X-axis = Day, facets = Treatment (Mock | RP)
# Points = EMM, error bars = 95% CI, coloured by Treatment.
# Forest plots: Treatment effect within each Day; Day effect within each Treatment.

make_emm_plot_trt <- function(emm_obj, pw_trt, pw_day, title, ylabel) {

  emm_df <- as.data.frame(emm_obj)
  names(emm_df)[names(emm_df) == "emmean"] <- "mean"
  names(emm_df)[grepl("lower|LCL", names(emm_df), ignore.case = TRUE)] <- "lower"
  names(emm_df)[grepl("upper|UCL", names(emm_df), ignore.case = TRUE)] <- "upper"

  p_emm <- ggplot(emm_df, aes(x = Day, y = mean, colour = Treatment)) +
    geom_point(size = 4, position = position_dodge(0.3)) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  width = 0.15, linewidth = 0.8,
                  position = position_dodge(0.3)) +
    geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE),
                  width = 0.0, linewidth = 3, alpha = 0.5,
                  position = position_dodge(0.3)) +
    facet_wrap(~ Treatment) +
    scale_colour_manual(values = treatment_colours) +
    labs(title = title, x = "Day", y = ylabel, caption = "Point = EMM | Thick bar = ±SE | Thin bar = 95% CI") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none",
          strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold", size = 11))

  make_forest <- function(pw_obj, plot_title, x_label, y_var) {
    df <- as.data.frame(pw_obj)
    names(df)[grepl("lower|LCL", names(df), ignore.case = TRUE)] <- "lower"
    names(df)[grepl("upper|UCL", names(df), ignore.case = TRUE)] <- "upper"
    df <- df %>%
      mutate(stars = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE            ~ "ns"
      ))

    ggplot(df, aes(x = estimate, y = as.character(.data[[y_var]]))) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_point(size = 3) +
      geom_errorbar(aes(xmin = lower, xmax = upper),
                    width = 0.2, linewidth = 0.8,
                    orientation = "y") +
      geom_text(aes(label = stars, x = upper),
                hjust = -0.3, size = 4) +
      labs(title = plot_title, x = x_label, y = NULL) +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 10))
  }

  p_forest_day <- make_forest(
    pw_obj     = pw_trt,
    plot_title = "Mock vs RP | by Day",
    x_label    = "Estimate (Mock - RP)",
    y_var      = "Day"
  )

  p_forest_trt <- make_forest(
    pw_obj     = pw_day,
    plot_title = "D1 vs D2 | by Treatment",
    x_label    = "Estimate (D1 - D2)",
    y_var      = "Treatment"
  )

  p_emm | p_forest_day | p_forest_trt
}


# --- 7. RUN MODELS ------------------------------------------------------------

cat("\n=== FITTING TREATMENT MODELS ===\n")

all_emms_trt        <- list()
all_pairwise_trt    <- list()
all_r2_trt          <- list()
all_models          <- list()
assumption_plots    <- list()
emm_plots           <- list()

model_summary_file <- paste0("treatment_model_summaries_INTERACTION_", output_suffix, ".txt")

sink(model_summary_file)
cat("TREATMENT MODEL SUMMARIES\n")
cat("=========================\n")
cat("Fixed effects: Treatment * Day + Location\n")
cat("Random effects: MouseID/TumourID_global/SliceID_global (ROI level)\n")
cat("              MouseID/TumourID_global/SliceID_global/ROIID_global (Track level)\n")
cat("Output suffix:", output_suffix, "\n\n")
sink()

run_models_for_level <- function(var_list, data, random_formula,
                                 level_name, rf_overrides = NULL,
                                 point_size = 1, alpha = 0.4) {

  for (vp in var_list) {
    for (type in c("raw", "log")) {

      if (type == "raw" && !isTRUE(vp$run_raw)) next
      if (type == "log" && is.null(vp$log)) next

      col       <- if (type == "raw") vp$raw else vp$log
      ylabel    <- if (type == "raw") vp$raw_label else vp$log_label
      transform <- if (type == "raw") "Raw" else "Log"
      key       <- paste(level_name, vp$name, transform, sep = "_")
      title     <- paste(level_name, "—", vp$name, paste0("(", transform, ")"))

      cat("  Fitting:", title, "\n")

      # Apply variable-specific random formula override if declared
      rf <- if (!is.null(rf_overrides) && !is.null(rf_overrides[[vp$name]])) {
        cat("  Random formula override for", vp$name, ":", rf_overrides[[vp$name]], "\n")
        rf_overrides[[vp$name]]
      } else {
        random_formula
      }

      # Apply variable-specific data filter (e.g. alpha_30pct)
      data_v <- if (!is.null(vp$filter_col)) {
        n_all  <- nrow(data)
        d      <- data %>% filter(.data[[vp$filter_col]] >= vp$filter_val)
        cat("  Processing:", ylabel,    
            "| filtered:", vp$filter_col, ">=", vp$filter_val,
            "| n (post-R2 filter):", nrow(d), "/", n_all,
            "| n (model-ready, excl. NA):", sum(!is.na(d[[col]])), "\n")
        d
      } else {
        data
      }

      model <- fit_treatment_model(data_v, col, rf)
      if (is.null(model)) next

      all_models[[key]] <<- model

      singular <- isSingular(model)
      if (singular) cat("  NOTE: Singular fit for", title, "\n")

      # R² — Nakagawa & Schielzeth (2013) Methods Ecol Evol 4:133-142
      r2_vals <- tryCatch(
        r2(model),
        error = function(e) { cat("  R2 failed:", e$message, "\n"); NULL }
      )

      if (!is.null(r2_vals)) {
        r2m <- round(r2_vals$R2_marginal,    3)
        r2c <- round(r2_vals$R2_conditional, 3)
        cat("  R2 marginal (fixed effects):  ", r2m, "\n")
        cat("  R2 conditional (fixed+random):", r2c, "\n")
      } else {
        r2m <- NA; r2c <- NA
      }

      all_r2_trt[[key]] <<- tibble(
        Model          = key,
        Variable       = vp$name,
        Level          = level_name,
        Transform      = transform,
        Singular_fit   = singular,
        R2_marginal    = r2m,
        R2_conditional = r2c,
        R2_random      = round(r2c - r2m, 3)
      )

      sink(model_summary_file, append = TRUE)
      cat("\n", paste(rep("=", 60), collapse = ""), "\n")
      cat(title, "\n")
      cat(paste(rep("=", 60), collapse = ""), "\n")
      cat("Singular fit:", singular, "\n")
      cat("R2 marginal (fixed effects):  ", r2m, "\n")
      cat("R2 conditional (fixed+random):", r2c, "\n")
      cat("R2 random effects only:       ", round(r2c - r2m, 3), "\n\n")
      print(summary(model))
      cat("\nAnova (Type III) — interaction test:\n")
      print(anova(model))
      sink()

      emm_result <- tryCatch(
        extract_emm_trt(model, vp$name, level_name, transform),
        error = function(e) { cat("  EMM failed:", e$message, "\n"); NULL }
      )

      if (!is.null(emm_result)) {
        all_emms_trt[[key]]     <<- emm_result$emm
        all_pairwise_trt[[key]] <<- emm_result$pairwise

        emm_plots[[key]] <<- make_emm_plot_trt(
          emm_obj = emm_result$emm_obj,
          pw_trt  = emm_result$pw_trt,
          pw_day  = emm_result$pw_day,
          title   = title,
          ylabel  = ylabel
        )
      }

      assumption_plots[[key]] <<- make_assumption_plots(model, title)
    }
  }
}

# ROI level — gated
if (!roi_analysis_enabled) {
  cat("\nROI-level treatment models skipped (roi_analysis_enabled = FALSE).\n")
  cat("Set roi_analysis_enabled <- TRUE once final roi_summary is confirmed.\n\n")
} else {
  cat("\n--- ROI level ---\n")
  run_models_for_level(
    var_list       = roi_vars[map_lgl(roi_vars, ~ isTRUE(var_include[.x$name]))],
    data           = roi,
    random_formula = random_formula_roi_default,
    rf_overrides   = random_overrides_roi,
    level_name     = "ROI"
  )
}

# Track level — always runs
cat("\n--- Track level ---\n")
run_models_for_level(
  var_list       = track_vars[map_lgl(track_vars, ~ isTRUE(var_include[.x$name]))],
  data           = tracks,
  random_formula = random_formula_track_default,
  rf_overrides   = random_overrides_track,
  level_name     = "Track"
)


# --- 8. SAVE NUMERICAL RESULTS ------------------------------------------------

emm_table_trt      <- bind_rows(all_emms_trt)
pairwise_table_trt <- bind_rows(all_pairwise_trt) %>%
  mutate(
    p_display = case_when(
      is.na(p.value)   ~ "NA",
      p.value < 0.0001 ~ "< 0.0001",
      TRUE             ~ as.character(round(p.value, 4))
    ),
    stars = case_when(
      is.na(p.value)  ~ "NA",
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )

r2_table_trt <- bind_rows(all_r2_trt)

write_csv(emm_table_trt,      paste0("treatment_emms_INTERACTION_",      output_suffix, ".csv"))
write_csv(pairwise_table_trt, paste0("treatment_pairwise_INTERACTION_",  output_suffix, ".csv"))
write_csv(r2_table_trt,       paste0("treatment_r2_INTERACTION_",        output_suffix, ".csv"))

cat("Saved: treatment_model_summaries_INTERACTION_", output_suffix, ".txt\n", sep = "")
cat("Saved: treatment_emms_INTERACTION_",            output_suffix, ".csv\n", sep = "")
cat("Saved: treatment_pairwise_INTERACTION_",        output_suffix, ".csv\n", sep = "")
cat("Saved: treatment_r2_INTERACTION_",              output_suffix, ".csv\n", sep = "")

cat("\n--- R2 Summary ---\n")
print(r2_table_trt, n = Inf)


# --- 9. SAVE ASSUMPTION PLOTS -------------------------------------------------

pdf(paste0("treatment_assumption_plots_", output_suffix, ".pdf"), width = 10, height = 4)
for (nm in names(assumption_plots)) print(assumption_plots[[nm]])
dev.off()
cat("Saved: treatment_assumption_plots_", output_suffix, ".pdf\n", sep = "")


# --- 10. SAVE EMM PLOTS -------------------------------------------------------

pdf(paste0("treatment_emm_plots_interaction_", output_suffix, ".pdf"), width = 10, height = 5)
for (vp in roi_vars) {
  for (type in c("Raw", "Log")) {
    key_roi   <- paste("ROI",   vp$name, type, sep = "_")
    key_track <- paste("Track", vp$name, type, sep = "_")
    if (!is.null(emm_plots[[key_roi]]))   print(emm_plots[[key_roi]])
    if (!is.null(emm_plots[[key_track]])) print(emm_plots[[key_track]])
  }
}
dev.off()
cat("Saved: treatment_emm_plots_interaction_", output_suffix, ".pdf\n", sep = "")
cat("\nAll outputs saved.\n")


# =============================================================================
# SECTION A: EFFECT SIZES
# =============================================================================
# Cohen's d from EMM contrasts using residual SD (sigma)
# Satterthwaite or Kenward-Roger df matched to each fixed effect
# Fold change from back-transforming log-scale estimates
# Runs for all stored models (track always; ROI only if gate was open)
# Reference: Cohen (1988) Statistical Power Analysis for the Behavioral Sciences
# =============================================================================

cat("\n=== EFFECT SIZES ===\n")

effect_size_results <- list()
geomean_results     <- list()
foldchange_results  <- list()

for (key in names(all_models)) {

  model <- all_models[[key]]

  parts     <- strsplit(key, "_")[[1]]
  level_nm  <- parts[1]
  var_nm    <- parts[2]
  transform <- parts[3]

  cat("\n  Effect sizes for:", key, "\n")

  aov_tbl <- tryCatch(
    anova(model),
    error = function(e) { cat("  ANOVA failed:", e$message, "\n"); NULL }
  )
  if (is.null(aov_tbl)) next

  sigma_resid <- sigma(model)
  sigma_total <- sqrt(sum(as.data.frame(VarCorr(model))$vcov))

  # --- EMMs for each fixed effect ---
  for (effect in c("Treatment", "Day", "Location")) {

    emm_obj <- tryCatch(
      emmeans(model, as.formula(paste("~", effect)),
              lmer.df = if (level_nm == "ROI") "kenward-roger" else "satterthwaite"),
      error = function(e) { cat("  EMM failed for", effect, ":", e$message, "\n"); NULL }
    )
    if (is.null(emm_obj)) next

    effect_df <- tryCatch(
      aov_tbl[effect, "DenDF"],
      error = function(e) df.residual(model)
    )
    if (is.null(effect_df) || is.na(effect_df)) {
      effect_df <- df.residual(model)
    }

    eff_resid <- tryCatch(
      eff_size(emm_obj, sigma = sigma_resid, edf = effect_df),
      error = function(e) { cat("  eff_size failed:", e$message, "\n"); NULL }
    )
    eff_total <- tryCatch(
      eff_size(emm_obj, sigma = sigma_total, edf = effect_df),
      error = function(e) NULL
    )

    if (is.null(eff_resid)) next

    emm_raw <- as.data.frame(emm_obj)
    names(emm_raw)[names(emm_raw) == "emmean"] <- "emmean_val"
    names(emm_raw)[grepl("lower|LCL", names(emm_raw), ignore.case = TRUE)] <- "lower_CL"
    names(emm_raw)[grepl("upper|UCL", names(emm_raw), ignore.case = TRUE)] <- "upper_CL"

    if (transform == "Log") {
      emm_raw <- emm_raw %>%
        mutate(
          geometric_mean     = round(exp(emmean_val), 3),
          geometric_mean_LCL = round(exp(lower_CL),   3),
          geometric_mean_UCL = round(exp(upper_CL),   3)
        )
    }

    pw_raw <- as.data.frame(pairs(emm_obj, infer = TRUE))
    names(pw_raw)[grepl("lower|LCL", names(pw_raw), ignore.case = TRUE)] <- "lower_CL"
    names(pw_raw)[grepl("upper|UCL", names(pw_raw), ignore.case = TRUE)] <- "upper_CL"

    if (transform == "Log") {
      pw_raw <- pw_raw %>%
        mutate(
          fold_change     = round(exp(estimate), 3),
          fold_change_LCL = round(exp(lower_CL), 3),
          fold_change_UCL = round(exp(upper_CL), 3),
          direction = case_when(
            fold_change > 1 ~ "first group faster",
            fold_change < 1 ~ "second group faster",
            TRUE            ~ "no difference"
          )
        )
    }

    eff_df <- as.data.frame(eff_resid) %>%
      mutate(
        Model            = key,
        Variable         = var_nm,
        Level            = level_nm,
        Transform        = transform,
        Fixed_effect     = effect,
        Den_df           = round(effect_df, 1),
        sigma_residual   = round(sigma_resid, 4),
        sigma_total      = round(sigma_total, 4),
        d_residual       = round(effect.size, 3),
        d_total          = if (!is.null(eff_total)) {
          round(as.data.frame(eff_total)$effect.size, 3)
        } else NA,
        d_interpretation = case_when(
          abs(effect.size) >= 0.8 ~ "Large",
          abs(effect.size) >= 0.5 ~ "Medium",
          abs(effect.size) >= 0.2 ~ "Small",
          TRUE                    ~ "Negligible"
        )
      )

    effect_size_results[[paste(key, effect, "effectsize", sep = "_")]] <- eff_df

    if (transform == "Log") {
      geomean_results[[paste(key, effect, "geomean", sep = "_")]] <- emm_raw %>%
        mutate(Model = key, Variable = var_nm, Level = level_nm,
               Transform = transform, Fixed_effect = effect)

      foldchange_results[[paste(key, effect, "foldchange", sep = "_")]] <- pw_raw %>%
        mutate(Model = key, Variable = var_nm, Level = level_nm,
               Transform = transform, Fixed_effect = effect)
    }
  }

  # --- Conditional EMMs: Treatment within each Day ---
  emm_trt_by_day <- tryCatch(
    emmeans(model, ~ Treatment | Day,
            lmer.df = if (level_nm == "ROI") "kenward-roger" else "satterthwaite"),
    error = function(e) { cat("  EMM failed Treatment|Day:", e$message, "\n"); NULL }
  )

  if (!is.null(emm_trt_by_day)) {

    eff_trt_by_day <- tryCatch(
      eff_size(emm_trt_by_day, sigma = sigma_resid,
               edf = aov_tbl["Treatment", "DenDF"]),
      error = function(e) NULL
    )
    if (!is.null(eff_trt_by_day)) {
      effect_size_results[[paste(key, "Treatment_by_Day", sep = "_")]] <-
        as.data.frame(eff_trt_by_day) %>%
        mutate(Model = key, Variable = var_nm, Level = level_nm,
               Transform = transform, Fixed_effect = "Treatment within Day",
               d_residual = round(effect.size, 3),
               sigma_residual = round(sigma_resid, 4),
               d_interpretation = case_when(
                 abs(effect.size) >= 0.8 ~ "Large",
                 abs(effect.size) >= 0.5 ~ "Medium",
                 abs(effect.size) >= 0.2 ~ "Small",
                 TRUE                    ~ "Negligible"
               ))
    }

    emm_trt_day_raw <- as.data.frame(emm_trt_by_day)
    names(emm_trt_day_raw)[names(emm_trt_day_raw) == "emmean"] <- "emmean_val"
    names(emm_trt_day_raw)[grepl("lower|LCL", names(emm_trt_day_raw), ignore.case = TRUE)] <- "lower_CL"
    names(emm_trt_day_raw)[grepl("upper|UCL", names(emm_trt_day_raw), ignore.case = TRUE)] <- "upper_CL"

    if (transform == "Log") {
      emm_trt_day_raw <- emm_trt_day_raw %>%
        mutate(
          geometric_mean     = round(exp(emmean_val), 3),
          geometric_mean_LCL = round(exp(lower_CL),   3),
          geometric_mean_UCL = round(exp(upper_CL),   3)
        )
      geomean_results[[paste(key, "Treatment_by_Day", "geomean", sep = "_")]] <-
        emm_trt_day_raw %>%
        mutate(Model = key, Variable = var_nm, Level = level_nm,
               Transform = transform, Fixed_effect = "Treatment within Day")
    }

    pw_trt_by_day <- tryCatch(
      as.data.frame(pairs(emm_trt_by_day, infer = TRUE)),
      error = function(e) NULL
    )
    if (!is.null(pw_trt_by_day) && transform == "Log") {
      names(pw_trt_by_day)[grepl("lower|LCL", names(pw_trt_by_day), ignore.case = TRUE)] <- "lower_CL"
      names(pw_trt_by_day)[grepl("upper|UCL", names(pw_trt_by_day), ignore.case = TRUE)] <- "upper_CL"
      pw_trt_by_day <- pw_trt_by_day %>%
        mutate(
          fold_change     = round(exp(estimate), 3),
          fold_change_LCL = round(exp(lower_CL), 3),
          fold_change_UCL = round(exp(upper_CL), 3),
          direction = case_when(
            fold_change > 1 ~ "first group faster",
            fold_change < 1 ~ "second group faster",
            TRUE            ~ "no difference"
          )
        )
      foldchange_results[[paste(key, "Treatment_by_Day", "foldchange", sep = "_")]] <-
        pw_trt_by_day %>%
        mutate(Model = key, Variable = var_nm, Level = level_nm,
               Transform = transform, Fixed_effect = "Treatment within Day")
    }
  }

  # --- Conditional EMMs: Day within each Treatment ---
  emm_day_by_trt <- tryCatch(
    emmeans(model, ~ Day | Treatment,
            lmer.df = if (level_nm == "ROI") "kenward-roger" else "satterthwaite"),
    error = function(e) { cat("  EMM failed Day|Treatment:", e$message, "\n"); NULL }
  )

  if (!is.null(emm_day_by_trt)) {

    eff_day_by_trt <- tryCatch(
      eff_size(emm_day_by_trt, sigma = sigma_resid,
               edf = aov_tbl["Day", "DenDF"]),
      error = function(e) NULL
    )
    if (!is.null(eff_day_by_trt)) {
      effect_size_results[[paste(key, "Day_by_Treatment", sep = "_")]] <-
        as.data.frame(eff_day_by_trt) %>%
        mutate(Model = key, Variable = var_nm, Level = level_nm,
               Transform = transform, Fixed_effect = "Day within Treatment",
               d_residual = round(effect.size, 3),
               sigma_residual = round(sigma_resid, 4),
               d_interpretation = case_when(
                 abs(effect.size) >= 0.8 ~ "Large",
                 abs(effect.size) >= 0.5 ~ "Medium",
                 abs(effect.size) >= 0.2 ~ "Small",
                 TRUE                    ~ "Negligible"
               ))
    }

    emm_day_trt_raw <- as.data.frame(emm_day_by_trt)
    names(emm_day_trt_raw)[names(emm_day_trt_raw) == "emmean"] <- "emmean_val"
    names(emm_day_trt_raw)[grepl("lower|LCL", names(emm_day_trt_raw), ignore.case = TRUE)] <- "lower_CL"
    names(emm_day_trt_raw)[grepl("upper|UCL", names(emm_day_trt_raw), ignore.case = TRUE)] <- "upper_CL"

    if (transform == "Log") {
      emm_day_trt_raw <- emm_day_trt_raw %>%
        mutate(
          geometric_mean     = round(exp(emmean_val), 3),
          geometric_mean_LCL = round(exp(lower_CL),   3),
          geometric_mean_UCL = round(exp(upper_CL),   3)
        )
      geomean_results[[paste(key, "Day_by_Treatment", "geomean", sep = "_")]] <-
        emm_day_trt_raw %>%
        mutate(Model = key, Variable = var_nm, Level = level_nm,
               Transform = transform, Fixed_effect = "Day within Treatment")
    }

    pw_day_by_trt <- tryCatch(
      as.data.frame(pairs(emm_day_by_trt, infer = TRUE)),
      error = function(e) NULL
    )
    if (!is.null(pw_day_by_trt) && transform == "Log") {
      names(pw_day_by_trt)[grepl("lower|LCL", names(pw_day_by_trt), ignore.case = TRUE)] <- "lower_CL"
      names(pw_day_by_trt)[grepl("upper|UCL", names(pw_day_by_trt), ignore.case = TRUE)] <- "upper_CL"
      pw_day_by_trt <- pw_day_by_trt %>%
        mutate(
          fold_change     = round(exp(estimate), 3),
          fold_change_LCL = round(exp(lower_CL), 3),
          fold_change_UCL = round(exp(upper_CL), 3),
          direction = case_when(
            fold_change > 1 ~ "first group faster",
            fold_change < 1 ~ "second group faster",
            TRUE            ~ "no difference"
          )
        )
      foldchange_results[[paste(key, "Day_by_Treatment", "foldchange", sep = "_")]] <-
        pw_day_by_trt %>%
        mutate(Model = key, Variable = var_nm, Level = level_nm,
               Transform = transform, Fixed_effect = "Day within Treatment")
    }
  }

  # --- Conditional EMMs: Treatment within each Location ---
  emm_trt_by_loc <- tryCatch(
    emmeans(model, ~ Treatment | Location,
            lmer.df = if (level_nm == "ROI") "kenward-roger" else "satterthwaite"),
    error = function(e) { cat("  EMM failed Treatment|Location:", e$message, "\n"); NULL }
  )

  if (!is.null(emm_trt_by_loc)) {

    eff_trt_by_loc <- tryCatch(
      eff_size(emm_trt_by_loc, sigma = sigma_resid,
               edf = aov_tbl["Treatment", "DenDF"]),
      error = function(e) NULL
    )
    if (!is.null(eff_trt_by_loc)) {
      effect_size_results[[paste(key, "Treatment_by_Location", sep = "_")]] <-
        as.data.frame(eff_trt_by_loc) %>%
        mutate(Model = key, Variable = var_nm, Level = level_nm,
               Transform = transform, Fixed_effect = "Treatment within Location",
               d_residual = round(effect.size, 3),
               sigma_residual = round(sigma_resid, 4),
               d_interpretation = case_when(
                 abs(effect.size) >= 0.8 ~ "Large",
                 abs(effect.size) >= 0.5 ~ "Medium",
                 abs(effect.size) >= 0.2 ~ "Small",
                 TRUE                    ~ "Negligible"
               ))
    }

    emm_trt_loc_raw <- as.data.frame(emm_trt_by_loc)
    names(emm_trt_loc_raw)[names(emm_trt_loc_raw) == "emmean"] <- "emmean_val"
    names(emm_trt_loc_raw)[grepl("lower|LCL", names(emm_trt_loc_raw), ignore.case = TRUE)] <- "lower_CL"
    names(emm_trt_loc_raw)[grepl("upper|UCL", names(emm_trt_loc_raw), ignore.case = TRUE)] <- "upper_CL"

    if (transform == "Log") {
      emm_trt_loc_raw <- emm_trt_loc_raw %>%
        mutate(
          geometric_mean     = round(exp(emmean_val), 3),
          geometric_mean_LCL = round(exp(lower_CL),   3),
          geometric_mean_UCL = round(exp(upper_CL),   3)
        )
      geomean_results[[paste(key, "Treatment_by_Location", "geomean", sep = "_")]] <-
        emm_trt_loc_raw %>%
        mutate(Model = key, Variable = var_nm, Level = level_nm,
               Transform = transform, Fixed_effect = "Treatment within Location")
    }

    pw_trt_by_loc <- tryCatch(
      as.data.frame(pairs(emm_trt_by_loc, infer = TRUE)),
      error = function(e) NULL
    )
    if (!is.null(pw_trt_by_loc) && transform == "Log") {
      names(pw_trt_by_loc)[grepl("lower|LCL", names(pw_trt_by_loc), ignore.case = TRUE)] <- "lower_CL"
      names(pw_trt_by_loc)[grepl("upper|UCL", names(pw_trt_by_loc), ignore.case = TRUE)] <- "upper_CL"
      pw_trt_by_loc <- pw_trt_by_loc %>%
        mutate(
          fold_change     = round(exp(estimate), 3),
          fold_change_LCL = round(exp(lower_CL), 3),
          fold_change_UCL = round(exp(upper_CL), 3),
          direction = case_when(
            fold_change > 1 ~ "first group faster",
            fold_change < 1 ~ "second group faster",
            TRUE            ~ "no difference"
          )
        )
      foldchange_results[[paste(key, "Treatment_by_Location", "foldchange", sep = "_")]] <-
        pw_trt_by_loc %>%
        mutate(Model = key, Variable = var_nm, Level = level_nm,
               Transform = transform, Fixed_effect = "Treatment within Location")
    }
  }
}

write_csv(bind_rows(effect_size_results), paste0("treatment_effect_sizes_INTERACTION_",    output_suffix, ".csv"))
write_csv(bind_rows(geomean_results),     paste0("treatment_geometric_means_INTERACTION_", output_suffix, ".csv"))
write_csv(bind_rows(foldchange_results),  paste0("treatment_fold_changes_INTERACTION_",    output_suffix, ".csv"))

cat("Saved: treatment_effect_sizes_INTERACTION_",    output_suffix, ".csv\n", sep = "")
cat("Saved: treatment_geometric_means_INTERACTION_", output_suffix, ".csv\n", sep = "")
cat("Saved: treatment_fold_changes_INTERACTION_",    output_suffix, ".csv\n", sep = "")


# --- 11. EXPLORATORY PLOTS ----------------------------------------------------
# X-axis = Treatment, facets = Day (treatment focus)
# Track-level always runs; ROI-level gated.
# Point shapes encode Location for spatial context.

make_exploratory_trt <- function(data, col, ylabel, title,
                                 point_size = 1.5, alpha = 0.2) {

  plot_data <- data %>%
    rename(y = all_of(col)) %>%
    filter(!is.na(y), !is.na(Location), !is.na(Treatment), !is.na(Day))

  summary_data <- plot_data %>%
    group_by(Day, Treatment) %>%
    summarise(
      mean_y = mean(y, na.rm = TRUE),
      sd_y   = sd(y,   na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(ymin = mean_y - sd_y,
           ymax = mean_y + sd_y)

  summary_data <<- summary_data

  ggplot() +
    geom_violin(
      data   = plot_data,
      aes(x = Treatment, y = y),
      fill   = "grey85",
      colour = "grey50",
      linewidth = 0.4,
      alpha  = 0.6,
      trim   = FALSE
    ) +
    geom_jitter(
      data   = plot_data,
      aes(x = Treatment, y = y, fill = Location, colour = Location, shape = Location),
      width  = 0.12,
      size   = point_size,
      alpha  = alpha,
      stroke = 0.4
    ) +
    geom_point(
      data   = summary_data,
      aes(x = Treatment, y = mean_y),
      shape  = 21, size = 2,
      fill   = "white", colour = "black", stroke = 0.9
    ) +
    geom_errorbar(
      data      = summary_data,
      aes(x = Treatment, ymin = ymin, ymax = ymax),
      width     = 0.15, colour = "black", linewidth = 0.6
    ) +
    facet_wrap(~ Day) +
    scale_fill_manual(values   = location_colours, name = "Location") +
    scale_colour_manual(values = location_colours, name = "Location") +
    scale_shape_manual(values  = location_shapes,  name = "Location") +
    labs(title = title, x = "Day", y = ylabel) +
    theme_bw(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 12),
      strip.background = element_rect(fill = "grey90", colour = "grey60"),
      strip.text       = element_text(face = "bold", size = 10),
      legend.position  = "right",
      panel.grid.minor = element_blank()
    )
}

# Exploratory plot variables derived from toggled track_vars / roi_vars.
# Avoids duplicating the variable list — inherits the same toggle.
build_expl_vars <- function(var_list, level_name) {
  out <- list()
  for (vp in var_list[map_lgl(var_list, ~ isTRUE(var_include[.x$name]))]) {
    if (isTRUE(vp$run_raw) && !is.null(vp$raw)) {
      out <- c(out, list(list(
        col   = vp$raw,
        label = vp$raw_label,
        title = paste(level_name, "—", vp$raw_label),
        filter_col = if (!is.null(vp$filter_col)) vp$filter_col else NULL,
        filter_val = if (!is.null(vp$filter_col)) vp$filter_val else NULL
      )))
    }
    if (!isFALSE(vp$run_log) && !is.null(vp$log)) {
      out <- c(out, list(list(
        col   = vp$log,
        label = vp$log_label,
        title = paste(level_name, "—", vp$log_label),
        filter_col = NULL, filter_val = NULL
      )))
    }
  }
  out
}

track_expl_vars <- build_expl_vars(track_vars, "Track")
roi_expl_vars   <- build_expl_vars(roi_vars,   "ROI")

pdf(paste0("treatment_exploratory_plots_", output_suffix, ".pdf"), width = 8, height = 5)

# Track level — always runs
for (vp in track_expl_vars) {
  data_v <- if (!is.null(vp$filter_col)) {
    tracks %>% filter(.data[[vp$filter_col]] >= vp$filter_val)
  } else tracks
  p <- make_exploratory_trt(
    data       = data_v,
    col        = vp$col,
    ylabel     = vp$label,
    title      = vp$title,
    point_size = 0.8,
    alpha      = 0.2
  )
  print(p)
}

# ROI level — gated
if (roi_analysis_enabled) {
  for (vp in roi_expl_vars) {
    p <- make_exploratory_trt(
      data       = roi,
      col        = vp$col,
      ylabel     = vp$label,
      title      = vp$title,
      point_size = 1.8,
      alpha      = 0.6
    )
    print(p)
  }
} else {
  cat("ROI-level exploratory plots skipped (roi_analysis_enabled = FALSE).\n")
}

dev.off()
cat("Saved: treatment_exploratory_plots_", output_suffix, ".pdf\n", sep = "")


# --- 12. MODEL COMPARISONS (LRT + BIC) ----------------------------------------
# Comparison 1: Treatment x Day interaction necessary?
#   Full:    variable ~ Treatment * Day + Location + (random)
#   Reduced: variable ~ Treatment + Day + Location + (random)
#
# Comparison 2: Random slope vs random intercept only
#   Full:    (1 + Treatment | TumourID_global) + (1 | Slice[/ROI])
#   Reduced: (1 | Mouse/Tumour/Slice[/ROI])
#   NOTE: likely singular with n=5 tumours — check isSingular().
#
# Both comparisons use REML=FALSE (required for LRT between nested models).
# BIC_diff = BIC(reduced) - BIC(full)
#   > 10  : strong evidence to keep full model
#   0-10  : weak evidence
#   < 0   : simpler model preferred

run_model_comparisons <- function(var_list, data, random_intercept,
                                  random_slope, level_name) {

  results <- list()

  for (vp in var_list) {
    for (type in c("raw", "log")) {

      if (type == "raw" && !isTRUE(vp$run_raw)) next

      col       <- if (type == "raw") vp$raw else vp$log
      transform <- if (type == "raw") "Raw" else "Log"
      key       <- paste(level_name, vp$name, transform, sep = "_")

      cat("  Comparing:", key, "\n")

      fit_ml <- function(formula_str) {
        tryCatch(
          lmer(as.formula(formula_str), data = data, REML = FALSE),
          error   = function(e) { cat("    Model error:", e$message, "\n"); NULL },
          warning = function(w) {
            suppressWarnings(lmer(as.formula(formula_str), data = data, REML = FALSE))
          }
        )
      }

      f_full        <- paste(col, "~ Treatment * Day + Location +", random_intercept)
      f_no_interact <- paste(col, "~ Treatment + Day + Location +", random_intercept)
      f_slope       <- paste(col, "~ Treatment * Day + Location +", random_slope)

      m_full        <- fit_ml(f_full)
      m_no_interact <- fit_ml(f_no_interact)
      m_slope       <- fit_ml(f_slope)

      if (is.null(m_full)) next

      # Comparison 1: interaction test
      int_row <- if (!is.null(m_no_interact)) {
        lrt_int <- tryCatch(anova(m_no_interact, m_full), error = function(e) NULL)
        if (!is.null(lrt_int)) {
          bic_diff <- round(BIC(m_no_interact) - BIC(m_full), 2)
          tibble(
            Model         = key,
            Variable      = vp$name,
            Level         = level_name,
            Transform     = transform,
            Comparison    = "Interaction (Trt*Day vs Trt+Day)",
            LRT_Chi2      = round(lrt_int[2, "Chisq"],       3),
            LRT_df        = lrt_int[2, "Df"],
            LRT_p         = round(lrt_int[2, "Pr(>Chisq)"],  4),
            LRT_p_display = ifelse(lrt_int[2, "Pr(>Chisq)"] < 0.0001,
                                   "< 0.0001",
                                   as.character(round(lrt_int[2, "Pr(>Chisq)"], 4))),
            BIC_full      = round(BIC(m_full),        2),
            BIC_reduced   = round(BIC(m_no_interact), 2),
            BIC_diff      = bic_diff,
            BIC_decision  = case_when(
              bic_diff > 10  ~ "Keep interaction",
              bic_diff > 0   ~ "Weak evidence for interaction",
              TRUE           ~ "Prefer additive model"
            ),
            Singular_full    = isSingular(m_full),
            Singular_reduced = isSingular(m_no_interact)
          )
        }
      } else NULL

      # Comparison 2: random slope test
      slope_row <- if (!is.null(m_slope)) {
        lrt_slope <- tryCatch(anova(m_full, m_slope), error = function(e) NULL)
        if (!is.null(lrt_slope)) {
          bic_diff_s <- round(BIC(m_full) - BIC(m_slope), 2)
          tibble(
            Model         = key,
            Variable      = vp$name,
            Level         = level_name,
            Transform     = transform,
            Comparison    = "Random slope (intercept+slope vs intercept only)",
            LRT_Chi2      = round(lrt_slope[2, "Chisq"],      3),
            LRT_df        = lrt_slope[2, "Df"],
            LRT_p         = round(lrt_slope[2, "Pr(>Chisq)"], 4),
            LRT_p_display = ifelse(lrt_slope[2, "Pr(>Chisq)"] < 0.0001,
                                   "< 0.0001",
                                   as.character(round(lrt_slope[2, "Pr(>Chisq)"], 4))),
            BIC_full      = round(BIC(m_slope), 2),
            BIC_reduced   = round(BIC(m_full),  2),
            BIC_diff      = bic_diff_s,
            BIC_decision  = case_when(
              bic_diff_s > 10  ~ "Keep random slope",
              bic_diff_s > 0   ~ "Weak evidence for slope",
              TRUE             ~ "Prefer intercept only"
            ),
            Singular_full    = isSingular(m_slope),
            Singular_reduced = isSingular(m_full)
          )
        }
      } else NULL

      results[[paste0(key, "_interaction")]] <- int_row
      results[[paste0(key, "_slope")]]       <- slope_row
    }
  }

  bind_rows(results)
}

# Random slope formulas for comparison (REML=FALSE versions)
random_slope_roi   <- "(1 + Treatment | TumourID_global) + (1 | SliceID_global)"
random_slope_track <- "(1 + Treatment | TumourID_global) + (1 | SliceID_global/ROIID_global)"

cat("\n=== MODEL COMPARISONS ===\n")
cat("Using REML=FALSE for all LRT comparisons\n\n")

comp_results <- list()

# ROI level — gated
if (roi_analysis_enabled) {
  cat("--- ROI level ---\n")
  comp_results[["ROI"]] <- run_model_comparisons(
    var_list         = roi_vars[map_lgl(roi_vars, ~ isTRUE(var_include[.x$name]))],
    data             = roi,
    random_intercept = random_formula_roi_default,
    random_slope     = random_slope_roi,
    level_name       = "ROI"
  )
} else {
  cat("ROI-level model comparisons skipped (roi_analysis_enabled = FALSE).\n\n")
}

# Track level — always runs
cat("\n--- Track level ---\n")
comp_results[["Track"]] <- run_model_comparisons(
  var_list         = track_vars[map_lgl(track_vars, ~ isTRUE(var_include[.x$name]))],
  data             = tracks,
  random_intercept = random_formula_track_default,
  random_slope     = random_slope_track,
  level_name       = "Track"
)

comparison_table <- bind_rows(comp_results)

write_csv(comparison_table, paste0("treatment_model_comparisons_", output_suffix, ".csv"))
cat("Saved: treatment_model_comparisons_", output_suffix, ".csv\n", sep = "")

cat("\n--- Model Comparison Summary ---\n")
print(comparison_table %>%
        select(Model, Comparison, LRT_p_display, BIC_diff, BIC_decision,
               Singular_full, Singular_reduced),
      n = Inf)


# =============================================================================
# ---- KS TEST — TRACK LEVEL DISTRIBUTIONAL COMPARISONS ----
# =============================================================================
# NOTE: Standard KS test used as descriptive complement to mixed model.
# Tracks are not independent (nested within ROI/Slice/Tumour).
# Results should be interpreted as exploratory — the mixed model remains
# the primary inferential analysis.
# Benjamini-Hochberg correction applied across all pairwise comparisons
# per variable.
# Reference: Benjamini & Hochberg (1995) J R Stat Soc B 57:289-300
#
# Uses the already-loaded and GMM-filtered `tracks` object for consistency
# with the mixed model analysis above.
# =============================================================================

cat("\n=== KS TESTS (TRACK LEVEL — EXPLORATORY) ===\n")
cat("Note: tracks are not independent; treat as descriptive only.\n\n")

# KS test variables — derived from toggled track_vars (same toggle as models)
ks_vars <- track_expl_vars

run_ks <- function(x, y, group1, group2, comparison_type, variable) {

  x <- x[!is.na(x)]
  y <- y[!is.na(y)]

  if (length(x) < 2 | length(y) < 2) {
    return(tibble(
      Variable = variable, Comparison_type = comparison_type,
      Group1 = group1, Group2 = group2,
      n1 = length(x), n2 = length(y),
      D_statistic = NA, p_value = NA
    ))
  }

  ks <- ks.test(x, y)

  tibble(
    Variable        = variable,
    Comparison_type = comparison_type,
    Group1          = group1,
    Group2          = group2,
    n1              = length(x),
    n2              = length(y),
    D_statistic     = round(ks$statistic, 4),
    p_value         = ks$p.value
  )
}

run_all_ks <- function(data, col, label) {

  results <- list()
  get_vals <- function(df) df[[col]]

  loc_pairs <- list(c("Core", "CoreEdge"), c("Core", "Edge"), c("CoreEdge", "Edge"))

  results[["trt_overall"]] <- run_ks(get_vals(filter(data, Treatment == "Mock")),
                                     get_vals(filter(data, Treatment == "RP")),
                                     "Mock", "RP", "Treatment overall", label)

  results[["day_overall"]] <- run_ks(get_vals(filter(data, Day == "D1")),
                                     get_vals(filter(data, Day == "D2")),
                                     "D1", "D2", "Day overall", label)

  for (pair in loc_pairs) {
    key <- paste("loc", pair[1], pair[2], sep = "_")
    results[[key]] <- run_ks(get_vals(filter(data, Location == pair[1])),
                             get_vals(filter(data, Location == pair[2])),
                             pair[1], pair[2], "Location overall", label)
  }

  for (d in c("D1", "D2")) {
    sub <- filter(data, Day == d)
    results[[paste("trt_day", d, sep = "_")]] <- run_ks(
      get_vals(filter(sub, Treatment == "Mock")),
      get_vals(filter(sub, Treatment == "RP")),
      "Mock", "RP", paste("Treatment within", d), label)
  }

  for (loc in c("Core", "CoreEdge", "Edge")) {
    sub <- filter(data, Location == loc)
    results[[paste("trt_loc", loc, sep = "_")]] <- run_ks(
      get_vals(filter(sub, Treatment == "Mock")),
      get_vals(filter(sub, Treatment == "RP")),
      "Mock", "RP", paste("Treatment within", loc), label)
  }

  for (trt in c("Mock", "RP")) {
    sub <- filter(data, Treatment == trt)
    results[[paste("day_trt", trt, sep = "_")]] <- run_ks(
      get_vals(filter(sub, Day == "D1")),
      get_vals(filter(sub, Day == "D2")),
      "D1", "D2", paste("Day within", trt), label)
  }

  for (trt in c("Mock", "RP")) {
    for (pair in loc_pairs) {
      sub <- filter(data, Treatment == trt)
      key <- paste("loc_trt", trt, pair[1], pair[2], sep = "_")
      results[[key]] <- run_ks(get_vals(filter(sub, Location == pair[1])),
                               get_vals(filter(sub, Location == pair[2])),
                               pair[1], pair[2], paste("Location within", trt), label)
    }
  }

  for (d in c("D1", "D2")) {
    for (pair in loc_pairs) {
      sub <- filter(data, Day == d)
      key <- paste("loc_day", d, pair[1], pair[2], sep = "_")
      results[[key]] <- run_ks(get_vals(filter(sub, Location == pair[1])),
                               get_vals(filter(sub, Location == pair[2])),
                               pair[1], pair[2], paste("Location within", d), label)
    }
  }

  bind_rows(results) %>%
    mutate(
      p_adjusted  = p.adjust(p_value, method = "BH"),
      Significant = p_adjusted < 0.05,
      p_display   = case_when(
        is.na(p_adjusted)    ~ "NA",
        p_adjusted < 0.0001  ~ "< 0.0001",
        TRUE                 ~ as.character(round(p_adjusted, 4))
      ),
      stars = case_when(
        is.na(p_adjusted) ~ "NA",
        p_adjusted < 0.001 ~ "***",
        p_adjusted < 0.01  ~ "**",
        p_adjusted < 0.05  ~ "*",
        TRUE               ~ "ns"
      )
    )
}

# Duration balance check — confirms track duration distributions are comparable
# across Treatment arms after GMM filtering, validating TrackLength comparisons.
ks_duration <- bind_rows(
  run_ks(
    x = tracks %>% filter(Treatment == "Mock", Day == "D1") %>% pull(TrackDuration),
    y = tracks %>% filter(Treatment == "RP",   Day == "D1") %>% pull(TrackDuration),
    group1 = "Mock", group2 = "RP",
    comparison_type = "Treatment within D1",
    variable = "TrackDuration"
  ),
  run_ks(
    x = tracks %>% filter(Treatment == "Mock", Day == "D2") %>% pull(TrackDuration),
    y = tracks %>% filter(Treatment == "RP",   Day == "D2") %>% pull(TrackDuration),
    group1 = "Mock", group2 = "RP",
    comparison_type = "Treatment within D2",
    variable = "TrackDuration"
  ),
  run_ks(
    x = tracks %>% filter(Day == "D1") %>% pull(TrackDuration),
    y = tracks %>% filter(Day == "D2") %>% pull(TrackDuration),
    group1 = "D1", group2 = "D2",
    comparison_type = "Day overall",
    variable = "TrackDuration"
  )
) %>%
  mutate(
    p_adjusted  = p.adjust(p_value, method = "BH"),
    Significant = p_adjusted < 0.05,
    p_display   = ifelse(p_adjusted < 0.0001, "< 0.0001",
                         as.character(round(p_adjusted, 4)))
  )

write_csv(ks_duration, paste0("ks_duration_balance_", output_suffix, ".csv"))
cat("\n--- TrackDuration balance check ---\n")
print(ks_duration, n = Inf)
cat("Saved: ks_duration_balance_", output_suffix, ".csv\n", sep = "")

ks_results <- map_dfr(ks_vars, function(vp) {
  cat("Running KS tests for:", vp$label, "\n")
  data_v <- if (!is.null(vp$filter_col)) {
    tracks %>% filter(.data[[vp$filter_col]] >= vp$filter_val)
  } else tracks
  run_all_ks(data_v, vp$col, vp$label)
})

write_csv(ks_results, paste0("ks_track_results_", output_suffix, ".csv"))
cat("\nSaved: ks_track_results_", output_suffix, ".csv\n", sep = "")
cat("Total tests run:", nrow(ks_results), "\n")
cat("Significant after BH correction:", sum(ks_results$Significant, na.rm = TRUE), "\n")

cat("\n--- Significant results ---\n")
ks_results %>%
  filter(Significant) %>%
  select(Variable, Comparison_type, Group1, Group2,
         n1, n2, D_statistic, p_display, stars) %>%
  print(n = Inf)

# --- Duration distribution plots ----------------------------------------------
pdf(paste0("track_duration_distributions_", output_suffix, ".pdf"),
    width = 10, height = 8)

# Panel 1: ECDF by Treatment within each Day
# D statistic is the maximum vertical distance between the two ECDF lines
p_ecdf <- tracks %>%
  ggplot(aes(x = TrackDuration, colour = Treatment, linetype = Treatment)) +
  stat_ecdf(linewidth = 0.9) +
  facet_wrap(~ Day) +
  scale_colour_manual(values = treatment_colours) +
  labs(
    title    = "Track Duration — ECDF by Treatment",
    subtitle = paste0("D statistic shown as maximum vertical separation between curves",
                      "  |  suffix: ", output_suffix),
    x        = "Track Duration (s)",
    y        = "Cumulative proportion",
    caption  = "KS D ~ 0.14 = maximum 14 percentage point separation between CDFs"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 10, colour = "grey40"),
    strip.background = element_rect(fill = "grey90", colour = "grey60"),
    strip.text       = element_text(face = "bold", size = 10),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

# Panel 2: Density by Treatment within each Day
# Shows the bimodal structure and where the distributions separate
p_dens <- tracks %>%
  ggplot(aes(x = TrackDuration, fill = Treatment, colour = Treatment)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  facet_wrap(~ Day) +
  scale_fill_manual(values  = treatment_colours) +
  scale_colour_manual(values = treatment_colours) +
  labs(
    title    = "Track Duration — Density by Treatment",
    subtitle = "Shaded area shows distributional overlap",
    x        = "Track Duration (s)",
    y        = "Density"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 10, colour = "grey40"),
    strip.background = element_rect(fill = "grey90", colour = "grey60"),
    strip.text       = element_text(face = "bold", size = 10),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

# Panel 3: Density by Day overall
# Shows D1 vs D2 separation
p_day <- tracks %>%
  ggplot(aes(x = TrackDuration, fill = Day, colour = Day)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  scale_fill_manual(values   = c("D1" = "#457B9D", "D2" = "#E63946")) +
  scale_colour_manual(values = c("D1" = "#457B9D", "D2" = "#E63946")) +
  labs(
    title    = "Track Duration — Density by Day (overall)",
    subtitle = "D1 vs D2 regardless of treatment",
    x        = "Track Duration (s)",
    y        = "Density"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 10, colour = "grey40"),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

# Panel 4: Box/violin by Treatment x Day x Location
# Checks whether duration imbalance is location-specific
p_loc <- tracks %>%
  ggplot(aes(x = Location, y = TrackDuration,
             fill = Treatment, colour = Treatment)) +
  geom_violin(alpha = 0.35, linewidth = 0.5,
              position = position_dodge(0.8), scale = "width") +
  geom_boxplot(alpha = 0.6, linewidth = 0.5, width = 0.15, outlier.size = 0.3,
               position = position_dodge(0.8)) +
  facet_wrap(~ Day) +
  scale_fill_manual(values   = treatment_colours) +
  scale_colour_manual(values = treatment_colours) +
  labs(
    title    = "Track Duration — by Location, Treatment and Day",
    subtitle = "Checks whether duration imbalance is spatially uniform",
    x        = "Location",
    y        = "Track Duration (s)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 10, colour = "grey40"),
    strip.background = element_rect(fill = "grey90", colour = "grey60"),
    strip.text       = element_text(face = "bold", size = 10),
    axis.text.x      = element_text(angle = 30, hjust = 1),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

print(p_ecdf)
print(p_dens)
print(p_day)
print(p_loc)

dev.off()
cat("Saved: track_duration_distributions_", output_suffix, ".pdf\n", sep = "")
# =============================================================================
# ---- SECTION C: PER-MOUSE DIRECTION ANALYSIS ----
# =============================================================================
# Purpose: Assess whether treatment effects are consistent in direction
#          across individual animals — a biological validity check
#          independent of formal statistical testing.
#
# With n=3 mice, formal tests of between-mouse consistency are underpowered.
# Direction consistency is the appropriate validity criterion at this n.
#
# All outputs use roi data — gated by roi_analysis_enabled.
# =============================================================================

if (!roi_analysis_enabled) {
  cat("\nSection C (per-mouse direction analysis) skipped (roi_analysis_enabled = FALSE).\n")
  cat("Set roi_analysis_enabled <- TRUE once final roi_summary is confirmed.\n")
} else {

  cat("\n=== SECTION C: PER-MOUSE DIRECTION ANALYSIS ===\n")

  # Log-transformed variables — direction assessed via geometric mean + fold change
  per_mouse_vars <- map(
    roi_vars[map_lgl(roi_vars, ~ isTRUE(var_include[.x$name]) &
                       !is.null(.x$log))],
    function(vp) list(
      col       = vp$log,
      label     = vp$log_label,
      raw_col   = vp$raw,          # may be NULL for TrackLength
      raw_label = if (!is.null(vp$raw)) vp$raw_label else vp$log_label
    )
  )
  
  # Raw variables — direction assessed via arithmetic mean + raw difference
  # TurningAngle: raw difference interpretable (degrees shift)
  # Alpha: raw difference interpretable (shift in diffusion exponent)
  # Note: Alpha30_filt_mean in roi_summary already filtered to R2 >= 0.5
  per_mouse_vars_raw <- map(
    roi_vars[map_lgl(roi_vars, ~ isTRUE(var_include[.x$name]) & is.null(.x$log))],
    function(vp) list(
      col       = vp$raw,
      label     = vp$raw_label
    )
  )

  # --- 1. Per-mouse summary table ---
  per_mouse_long <- map_dfr(per_mouse_vars, function(vp) {
    roi %>%
      group_by(MouseID, TumourID_global, Treatment, Day) %>%
      summarise(
        mean_log   = mean(.data[[vp$col]],    na.rm = TRUE),
        mean_raw   = if (!is.null(vp$raw_col))
          mean(.data[[vp$raw_col]], na.rm = TRUE) else NA_real_,
        geom_mean  = exp(mean(.data[[vp$col]], na.rm = TRUE)),
        sd_log     = sd(.data[[vp$col]],      na.rm = TRUE),
        n_rois     = n(),
        .groups    = "drop"
      ) %>%
      mutate(
        Variable   = vp$label,
        mean_log   = round(mean_log,  3),
        mean_raw   = round(mean_raw,  3),
        geom_mean  = round(geom_mean, 3),
        sd_log     = round(sd_log,    3)
      )
  })

  write_csv(per_mouse_long, paste0("per_mouse_speed_summary_", output_suffix, ".csv"))
  cat("Saved: per_mouse_speed_summary_", output_suffix, ".csv\n", sep = "")

  # ---- 2. Direction check table ----
  direction_results <- list()

  for (vp in per_mouse_vars) {

    mock_vals <- roi %>%
      filter(Treatment == "Mock") %>%
      group_by(MouseID, Day) %>%
      summarise(mock_log  = mean(.data[[vp$col]], na.rm = TRUE),
                mock_geom = exp(mean(.data[[vp$col]], na.rm = TRUE)),
                mock_n    = n(),
                .groups   = "drop")

    rp_vals <- roi %>%
      filter(Treatment == "RP") %>%
      group_by(MouseID, Day) %>%
      summarise(rp_log  = mean(.data[[vp$col]], na.rm = TRUE),
                rp_geom = exp(mean(.data[[vp$col]], na.rm = TRUE)),
                rp_n    = n(),
                .groups = "drop")

    direction_results[[vp$label]] <- left_join(
      mock_vals, rp_vals, by = c("MouseID", "Day")
    ) %>%
      mutate(
        Variable    = vp$label,
        log_diff    = round(rp_log  - mock_log,  3),
        fold_change = round(exp(rp_log - mock_log), 3),
        direction   = case_when(
          fold_change > 1 ~ "RP higher",
          fold_change < 1 ~ "Mock higher",
          TRUE            ~ "equal"
        ),
        mock_log  = round(mock_log,  3),
        mock_geom = round(mock_geom, 3),
        rp_log    = round(rp_log,    3),
        rp_geom   = round(rp_geom,   3)
      ) %>%
      arrange(Day, MouseID) %>%
      select(Variable, MouseID, Day,
             mock_log, mock_geom, mock_n,
             rp_log,   rp_geom,   rp_n,
             log_diff, fold_change, direction)
  }

  direction_table <- bind_rows(direction_results)

  write_csv(direction_table, paste0("per_mouse_direction_check_", output_suffix, ".csv"))
  cat("Saved: per_mouse_direction_check_", output_suffix, ".csv\n", sep = "")

  for (vp in per_mouse_vars) {
    cat("\n---", vp$label, "---\n")
    direction_table %>% filter(Variable == vp$label) %>% print(n = Inf)
  }

  cat("\n--- Direction consistency summary ---\n")
  direction_table %>%
    group_by(Variable, Day, direction) %>%
    summarise(n_mice = n(), .groups = "drop") %>%
    pivot_wider(names_from = direction, values_from = n_mice, values_fill = 0) %>%
    mutate(
      n_consistent = if ("RP higher" %in% names(.)) `RP higher` else 0,
      verdict = case_when(
        n_consistent == 3 ~ "Consistent — RP higher in all 3 mice",
        n_consistent == 2 ~ "Mostly consistent — RP higher in 2/3 mice",
        n_consistent == 1 ~ "Inconsistent — RP higher in 1/3 mice",
        n_consistent == 0 ~ "Consistent — Mock higher in all 3 mice"
      )
    ) %>%
    select(Variable, Day, verdict) %>%
    print(n = Inf)
  
  # ----  2b. Raw variable direction check (TurningAngle, Alpha) ----
  if (length(per_mouse_vars_raw) > 0) {
    
    cat("\n--- Raw variable direction check (arithmetic mean) ---\n")
    cat("Note: fold change not computed — raw difference used instead\n\n")
    
    direction_results_raw <- list()
    
    for (vp in per_mouse_vars_raw) {
      
      mock_vals <- roi %>%
        filter(Treatment == "Mock") %>%
        group_by(MouseID, Day) %>%
        summarise(mock_mean = mean(.data[[vp$col]], na.rm = TRUE),
                  mock_n    = sum(!is.na(.data[[vp$col]])),
                  .groups   = "drop")
      
      rp_vals <- roi %>%
        filter(Treatment == "RP") %>%
        group_by(MouseID, Day) %>%
        summarise(rp_mean = mean(.data[[vp$col]], na.rm = TRUE),
                  rp_n    = sum(!is.na(.data[[vp$col]])),
                  .groups = "drop")
      
      direction_results_raw[[vp$label]] <- left_join(
        mock_vals, rp_vals, by = c("MouseID", "Day")
      ) %>%
        mutate(
          Variable  = vp$label,
          raw_diff  = round(rp_mean - mock_mean, 4),
          ratio     = round(rp_mean / mock_mean, 3),
          direction = case_when(
            raw_diff > 0  ~ "RP higher",
            raw_diff < 0  ~ "Mock higher",
            TRUE          ~ "equal"
          ),
          mock_mean = round(mock_mean, 4),
          rp_mean   = round(rp_mean,   4)
        ) %>%
        arrange(Day, MouseID) %>%
        select(Variable, MouseID, Day,
               mock_mean, mock_n,
               rp_mean,   rp_n,
               raw_diff,  ratio, direction)
    }
    
    direction_table_raw <- bind_rows(direction_results_raw)
    write_csv(direction_table_raw,
              paste0("per_mouse_direction_check_raw_", output_suffix, ".csv"))
    cat("Saved: per_mouse_direction_check_raw_", output_suffix, ".csv\n\n", sep = "")
    
    for (vp in per_mouse_vars_raw) {
      cat("\n---", vp$label, "---\n")
      direction_table_raw %>% filter(Variable == vp$label) %>% print(n = Inf)
    }
    
    cat("\n--- Raw direction consistency summary ---\n")
    direction_table_raw %>%
      group_by(Variable, Day, direction) %>%
      summarise(n_mice = n(), .groups = "drop") %>%
      print(n = Inf)
  }

  # ---- 3. Per-mouse trajectory plots ----

  mouse_colours <- c("M1" = "#457B9D", "M2" = "#E63946", "M3" = "#F4A261", "M4" = "#2A9D8F","M5" = "#6A4C93")
  treatment_linetypes <- c("Mock" = "dashed", "RP" = "solid")

  pdf(paste0("per_mouse_trajectory_plots_", output_suffix, ".pdf"), width = 9, height = 5)

  for (vp in per_mouse_vars) {

    # Mouse-level (averaged across tumours)
    plot_data_mouse <- roi %>%
      group_by(MouseID, Treatment, Day) %>%
      summarise(geom_mean = exp(mean(.data[[vp$col]], na.rm = TRUE)),
                n_rois    = n(), .groups = "drop")

    summary_data_mouse <- plot_data_mouse %>%
      group_by(Treatment, Day) %>%
      summarise(mean_geom = mean(geom_mean, na.rm = TRUE),
                sd_geom   = sd(geom_mean,   na.rm = TRUE), .groups = "drop") %>%
      mutate(ymin = mean_geom - sd_geom, ymax = mean_geom + sd_geom)

    p_mouse <- ggplot() +
      geom_line(data = plot_data_mouse,
                aes(x = Day, y = geom_mean,
                    group    = interaction(MouseID, Treatment),
                    colour   = MouseID,
                    linetype = Treatment),
                linewidth = 0.9, alpha = 0.8) +
      geom_point(data = plot_data_mouse,
                 aes(x = Day, y = geom_mean, colour = MouseID, shape = Treatment),
                 size = 3) +
      geom_point(data = summary_data_mouse,
                 aes(x = Day, y = mean_geom, group = Treatment),
                 shape = 21, size = 4, fill = "white", colour = "black", stroke = 1,
                 position = position_dodge(0.15)) +
      geom_errorbar(data = summary_data_mouse,
                    aes(x = Day, ymin = ymin, ymax = ymax, group = Treatment),
                    width = 0.1, linewidth = 0.7, colour = "black",
                    position = position_dodge(0.15)) +
      facet_wrap(~ Treatment) +
      scale_colour_manual(values = mouse_colours,         name = "Mouse") +
      scale_linetype_manual(values = treatment_linetypes, name = "Treatment") +
      scale_shape_manual(values = c("Mock" = 1, "RP" = 16), name = "Treatment") +
      labs(title    = paste("Mouse-level trajectories —", vp$label),
           subtitle = "Each line = one mouse | Open circle ± SD = overall mean",
           x = "Day", y = vp$raw_label) +
      theme_bw(base_size = 11) +
      theme(plot.title       = element_text(face = "bold", size = 12),
            plot.subtitle    = element_text(size = 10, colour = "grey40"),
            strip.background = element_rect(fill = "grey90", colour = "grey60"),
            strip.text       = element_text(face = "bold", size = 10),
            legend.position  = "right",
            panel.grid.minor = element_blank())

    # Tumour-level (separate lines per tumour)
    plot_data_tumour <- roi %>%
      group_by(MouseID, TumourID_global, Treatment, Day) %>%
      summarise(geom_mean = exp(mean(.data[[vp$col]], na.rm = TRUE)),
                n_rois    = n(), .groups = "drop")

    summary_data <- roi %>%
      group_by(Treatment, Day) %>%
      summarise(mean_geom = exp(mean(.data[[vp$col]], na.rm = TRUE)),
                sd_geom   = sd(exp(.data[[vp$col]]),  na.rm = TRUE), .groups = "drop") %>%
      mutate(ymin = mean_geom - sd_geom, ymax = mean_geom + sd_geom)

    p_tumour <- ggplot() +
      geom_line(data = plot_data_tumour,
                aes(x = Day, y = geom_mean,
                    group    = interaction(TumourID_global, Treatment),
                    colour   = MouseID,
                    linetype = Treatment),
                linewidth = 0.9, alpha = 0.8) +
      geom_point(data = plot_data_tumour,
                 aes(x = Day, y = geom_mean, colour = MouseID, shape = Treatment),
                 size = 3) +
      geom_point(data = summary_data,
                 aes(x = Day, y = mean_geom, group = Treatment),
                 shape = 21, size = 4, fill = "white", colour = "black", stroke = 1,
                 position = position_dodge(0.15)) +
      geom_errorbar(data = summary_data,
                    aes(x = Day, ymin = ymin, ymax = ymax, group = Treatment),
                    width = 0.1, linewidth = 0.7, colour = "black",
                    position = position_dodge(0.15)) +
      facet_wrap(~ Treatment) +
      scale_colour_manual(values = mouse_colours,         name = "Mouse") +
      scale_linetype_manual(values = treatment_linetypes, name = "Treatment") +
      scale_shape_manual(values = c("Mock" = 1, "RP" = 16), name = "Treatment") +
      labs(title    = paste("Tumour-level trajectories —", vp$label),
           subtitle = "Each line = one tumour | Open circle ± SD = overall mean",
           x = "Day", y = vp$raw_label) +
      theme_bw(base_size = 11) +
      theme(plot.title       = element_text(face = "bold", size = 12),
            plot.subtitle    = element_text(size = 10, colour = "grey40"),
            strip.background = element_rect(fill = "grey90", colour = "grey60"),
            strip.text       = element_text(face = "bold", size = 10),
            legend.position  = "right",
            panel.grid.minor = element_blank())

    print(p_mouse)
    print(p_tumour)
    
    # --- Raw variable trajectory plots (arithmetic mean) ---
    if (length(per_mouse_vars_raw) > 0) {
      for (vp in per_mouse_vars_raw) {
        
        # Note in subtitle distinguishes arithmetic mean from geometric mean pages
        plot_data_raw <- roi %>%
          group_by(MouseID, Treatment, Day) %>%
          summarise(arith_mean = mean(.data[[vp$col]], na.rm = TRUE),
                    n_rois     = sum(!is.na(.data[[vp$col]])),
                    .groups    = "drop")
        
        summary_data_raw <- plot_data_raw %>%
          group_by(Treatment, Day) %>%
          summarise(mean_val = mean(arith_mean, na.rm = TRUE),
                    sd_val   = sd(arith_mean,   na.rm = TRUE),
                    .groups  = "drop") %>%
          mutate(ymin = mean_val - sd_val, ymax = mean_val + sd_val)
        
        p_raw <- ggplot() +
          geom_line(data = plot_data_raw,
                    aes(x = Day, y = arith_mean,
                        group    = interaction(MouseID, Treatment),
                        colour   = MouseID,
                        linetype = Treatment),
                    linewidth = 0.9, alpha = 0.8) +
          geom_point(data = plot_data_raw,
                     aes(x = Day, y = arith_mean,
                         colour = MouseID, shape = Treatment),
                     size = 3) +
          geom_point(data = summary_data_raw,
                     aes(x = Day, y = mean_val, group = Treatment),
                     shape = 21, size = 4, fill = "white",
                     colour = "black", stroke = 1,
                     position = position_dodge(0.15)) +
          geom_errorbar(data = summary_data_raw,
                        aes(x = Day, ymin = ymin, ymax = ymax,
                            group = Treatment),
                        width = 0.1, linewidth = 0.7, colour = "black",
                        position = position_dodge(0.15)) +
          facet_wrap(~ Treatment) +
          scale_colour_manual(values = mouse_colours,          name = "Mouse") +
          scale_linetype_manual(values = treatment_linetypes,  name = "Treatment") +
          scale_shape_manual(values = c("Mock" = 1, "RP" = 16), name = "Treatment") +
          labs(title    = paste("Mouse-level trajectories —", vp$label),
               subtitle = "Arithmetic mean (no log transform) | Open circle ± SD = overall mean",
               x = "Day", y = vp$label) +
          theme_bw(base_size = 11) +
          theme(plot.title       = element_text(face = "bold", size = 12),
                plot.subtitle    = element_text(size = 10, colour = "grey40"),
                strip.background = element_rect(fill = "grey90", colour = "grey60"),
                strip.text       = element_text(face = "bold", size = 10),
                legend.position  = "right",
                panel.grid.minor = element_blank())
        
        print(p_raw)
      }
    }
    
  }

  dev.off()
  cat("Saved: per_mouse_trajectory_plots_", output_suffix, ".pdf\n", sep = "")
  cat("Pages 1-2: Speed | Pages 3-4: TrackLength | Pages 5-6: Straightness\n")
  cat("Odd pages = Mouse level | Even pages = Tumour level\n")

  # --- 4. Per-mouse 1x3 facet trajectory plots ---

  pdf(paste0("per_mouse_facet_trajectories_", output_suffix, ".pdf"), width = 10, height = 4)

  for (vp in per_mouse_vars) {

    plot_data_facet <- roi %>%
      group_by(MouseID, Treatment, Day) %>%
      summarise(geom_mean = exp(mean(.data[[vp$col]], na.rm = TRUE)),
                n_rois    = n(), .groups = "drop")

    p_facet <- ggplot(plot_data_facet,
                      aes(x = Day, y = geom_mean, colour = Treatment, group = Treatment)) +
      geom_line(linewidth = 0.9) +
      geom_point(aes(shape = Treatment), size = 3) +
      facet_wrap(~ MouseID, nrow = 1) +
      scale_colour_manual(values = treatment_colours, name = "Treatment") +
      scale_shape_manual(values = c("Mock" = 1, "RP" = 16), name = "Treatment") +
      labs(title    = paste("Per-mouse trajectories —", vp$raw_label),
           subtitle = "Each line = treatment group averaged across tumours within mouse",
           x = "Day", y = vp$raw_label) +
      theme_bw(base_size = 11) +
      theme(plot.title       = element_text(face = "bold", size = 12),
            plot.subtitle    = element_text(size = 10, colour = "grey40"),
            strip.background = element_rect(fill = "grey90", colour = "grey60"),
            strip.text       = element_text(face = "bold", size = 11),
            legend.position  = "right",
            panel.grid.minor = element_blank())

    print(p_facet)
    
    # --- Raw variable trajectory plots (arithmetic mean) ---
    if (length(per_mouse_vars_raw) > 0) {
      for (vp in per_mouse_vars_raw) {
        
        # Note in subtitle distinguishes arithmetic mean from geometric mean pages
        plot_data_raw <- roi %>%
          group_by(MouseID, Treatment, Day) %>%
          summarise(arith_mean = mean(.data[[vp$col]], na.rm = TRUE),
                    n_rois     = sum(!is.na(.data[[vp$col]])),
                    .groups    = "drop")
        
        summary_data_raw <- plot_data_raw %>%
          group_by(Treatment, Day) %>%
          summarise(mean_val = mean(arith_mean, na.rm = TRUE),
                    sd_val   = sd(arith_mean,   na.rm = TRUE),
                    .groups  = "drop") %>%
          mutate(ymin = mean_val - sd_val, ymax = mean_val + sd_val)
        
        p_raw <- ggplot() +
          geom_line(data = plot_data_raw,
                    aes(x = Day, y = arith_mean,
                        group    = interaction(MouseID, Treatment),
                        colour   = MouseID,
                        linetype = Treatment),
                    linewidth = 0.9, alpha = 0.8) +
          geom_point(data = plot_data_raw,
                     aes(x = Day, y = arith_mean,
                         colour = MouseID, shape = Treatment),
                     size = 3) +
          geom_point(data = summary_data_raw,
                     aes(x = Day, y = mean_val, group = Treatment),
                     shape = 21, size = 4, fill = "white",
                     colour = "black", stroke = 1,
                     position = position_dodge(0.15)) +
          geom_errorbar(data = summary_data_raw,
                        aes(x = Day, ymin = ymin, ymax = ymax,
                            group = Treatment),
                        width = 0.1, linewidth = 0.7, colour = "black",
                        position = position_dodge(0.15)) +
          facet_wrap(~ Treatment) +
          scale_colour_manual(values = mouse_colours,          name = "Mouse") +
          scale_linetype_manual(values = treatment_linetypes,  name = "Treatment") +
          scale_shape_manual(values = c("Mock" = 1, "RP" = 16), name = "Treatment") +
          labs(title    = paste("Mouse-level trajectories —", vp$label),
               subtitle = "Arithmetic mean (no log transform) | Open circle ± SD = overall mean",
               x = "Day", y = vp$label) +
          theme_bw(base_size = 11) +
          theme(plot.title       = element_text(face = "bold", size = 12),
                plot.subtitle    = element_text(size = 10, colour = "grey40"),
                strip.background = element_rect(fill = "grey90", colour = "grey60"),
                strip.text       = element_text(face = "bold", size = 10),
                legend.position  = "right",
                panel.grid.minor = element_blank())
        
        print(p_raw)
      }
    }
    
  }

  dev.off()
  cat("Saved: per_mouse_facet_trajectories_", output_suffix, ".pdf\n", sep = "")

  # --- 5. Normalised per-mouse facet trajectory plots ---

  pdf(paste0("per_mouse_facet_trajectories_normalised_", output_suffix, ".pdf"),
      width = 10, height = 4)

  for (vp in per_mouse_vars) {

    plot_data_facet <- roi %>%
      group_by(MouseID, Treatment, Day) %>%
      summarise(geom_mean = exp(mean(.data[[vp$col]], na.rm = TRUE)),
                n_rois    = n(), .groups = "drop")

    # Option 1: normalise per Mouse x Treatment (both D1 = 1)
    d1_by_trt <- plot_data_facet %>%
      filter(Day == "D1") %>%
      select(MouseID, Treatment, d1_ref = geom_mean)

    plot_data_opt1 <- plot_data_facet %>%
      left_join(d1_by_trt, by = c("MouseID", "Treatment")) %>%
      mutate(norm_geom = geom_mean / d1_ref)

    # Option 2: normalise per Mouse (Mock D1 as reference)
    d1_mock <- plot_data_facet %>%
      filter(Day == "D1", Treatment == "Mock") %>%
      select(MouseID, d1_mock_ref = geom_mean)

    plot_data_opt2 <- plot_data_facet %>%
      left_join(d1_mock, by = "MouseID") %>%
      mutate(norm_geom = geom_mean / d1_mock_ref)

    make_facet_plot <- function(plot_data, subtitle_text) {
      ggplot(plot_data, aes(x = Day, y = norm_geom,
                            colour = Treatment, group = Treatment)) +
        geom_hline(yintercept = 1, linetype = "dotted",
                   colour = "grey50", linewidth = 0.5) +
        geom_line(linewidth = 0.9) +
        geom_point(aes(shape = Treatment), size = 3) +
        facet_wrap(~ MouseID, nrow = 1) +
        scale_colour_manual(values = treatment_colours, name = "Treatment") +
        scale_shape_manual(values = c("Mock" = 1, "RP" = 16), name = "Treatment") +
        labs(title    = paste("Per-mouse trajectories —", vp$raw_label),
             subtitle = subtitle_text,
             x = "Day", y = paste0(vp$raw_label, " (normalised to D1)")) +
        theme_bw(base_size = 11) +
        theme(plot.title       = element_text(face = "bold", size = 12),
              plot.subtitle    = element_text(size = 10, colour = "grey40"),
              strip.background = element_rect(fill = "grey90", colour = "grey60"),
              strip.text       = element_text(face = "bold", size = 11),
              legend.position  = "right",
              panel.grid.minor = element_blank())
    }

    print(make_facet_plot(
      plot_data_opt1,
      "Option 1: each treatment normalised to its own D1 (both start at 1)"
    ))

    print(make_facet_plot(
      plot_data_opt2,
      "Option 2: both treatments normalised to Mock D1 (shared reference)"
    ))
    
    # --- Raw variable trajectory plots (arithmetic mean) ---
    if (length(per_mouse_vars_raw) > 0) {
      for (vp in per_mouse_vars_raw) {
        
        # Note in subtitle distinguishes arithmetic mean from geometric mean pages
        plot_data_raw <- roi %>%
          group_by(MouseID, Treatment, Day) %>%
          summarise(arith_mean = mean(.data[[vp$col]], na.rm = TRUE),
                    n_rois     = sum(!is.na(.data[[vp$col]])),
                    .groups    = "drop")
        
        summary_data_raw <- plot_data_raw %>%
          group_by(Treatment, Day) %>%
          summarise(mean_val = mean(arith_mean, na.rm = TRUE),
                    sd_val   = sd(arith_mean,   na.rm = TRUE),
                    .groups  = "drop") %>%
          mutate(ymin = mean_val - sd_val, ymax = mean_val + sd_val)
        
        p_raw <- ggplot() +
          geom_line(data = plot_data_raw,
                    aes(x = Day, y = arith_mean,
                        group    = interaction(MouseID, Treatment),
                        colour   = MouseID,
                        linetype = Treatment),
                    linewidth = 0.9, alpha = 0.8) +
          geom_point(data = plot_data_raw,
                     aes(x = Day, y = arith_mean,
                         colour = MouseID, shape = Treatment),
                     size = 3) +
          geom_point(data = summary_data_raw,
                     aes(x = Day, y = mean_val, group = Treatment),
                     shape = 21, size = 4, fill = "white",
                     colour = "black", stroke = 1,
                     position = position_dodge(0.15)) +
          geom_errorbar(data = summary_data_raw,
                        aes(x = Day, ymin = ymin, ymax = ymax,
                            group = Treatment),
                        width = 0.1, linewidth = 0.7, colour = "black",
                        position = position_dodge(0.15)) +
          facet_wrap(~ Treatment) +
          scale_colour_manual(values = mouse_colours,          name = "Mouse") +
          scale_linetype_manual(values = treatment_linetypes,  name = "Treatment") +
          scale_shape_manual(values = c("Mock" = 1, "RP" = 16), name = "Treatment") +
          labs(title    = paste("Mouse-level trajectories —", vp$label),
               subtitle = "Arithmetic mean (no log transform) | Open circle ± SD = overall mean",
               x = "Day", y = vp$label) +
          theme_bw(base_size = 11) +
          theme(plot.title       = element_text(face = "bold", size = 12),
                plot.subtitle    = element_text(size = 10, colour = "grey40"),
                strip.background = element_rect(fill = "grey90", colour = "grey60"),
                strip.text       = element_text(face = "bold", size = 10),
                legend.position  = "right",
                panel.grid.minor = element_blank())
        
        print(p_raw)
      }
    }
  }

  dev.off()
  cat("Saved: per_mouse_facet_trajectories_normalised_", output_suffix, ".pdf\n", sep = "")
  cat("2 normalisation options per variable (Options 1 and 2)\n")

}  # end roi_analysis_enabled gate for Section C

# =============================================================================
# ---- SECTION D: GMM COMPOSITIONAL ANALYSIS ----
# =============================================================================
# Purpose: Visualise how GMM component proportions shift across Treatment
#          and Day within each tumour Location.
#
# Two visualisations:
#
#   1. Stacked bar chart (1x3 Location facets)
#      X-axis: Mock D1 | RP D1 | Mock D2 | RP D2
#      Fill:   GMM component
#      Shows within-location compositional shifts across Treatment and Day
#
#   
# Outputs:
#   gmm_composition_barchart_<suffix>.pdf

#
# Requires: gmm_composition object from Section 1d
#           GMM_component column in tracks (script 02d must have been run)
# =============================================================================

# --- SECTION D PARAMETERS ----------------------------------------------------
# Toggle to use a pre-saved composition CSV instead of recomputing from
# the current tracks object. Useful for plotting compositional shifts
# from a different filtered population without rerunning the full pipeline.
#
# Set use_saved_composition = TRUE and point composition_path to the
# desired gmm_composition_by_location_<suffix>.csv file.
# Set use_saved_composition = FALSE to use the gmm_composition object
# already computed from the current tracks object in Section 1d.

use_saved_composition  <- TRUE
composition_path       <- "Z:\\lightmicro-inst\\scharalambous\\Resident_cell_motility_tr_v3\\ROI files\\gmm_composition_by_Loc_trt_GMM1_3_4_5_MaxDis7.csv"  # <-- update path
maxdis_motile_threshold <- 7   # must match the value used when CSV was generated

# Load or use existing
if (use_saved_composition) {
  if (!file.exists(composition_path)) {
    stop("Composition file not found: ", composition_path,
         "\nCheck composition_path in Section D parameters.")
  }
  gmm_composition <- read_csv(composition_path, show_col_types = FALSE) %>%
    mutate(
      Location      = factor(Location,
                             levels = c("Core", "CoreEdge", "Edge")),
      Treatment     = factor(Treatment, levels = c("Mock", "RP")),
      Day           = factor(Day,       levels = c("D1", "D2")),
      GMM_component = factor(GMM_component,
                             levels = c("1", "1 motile", "2",
                                        "3", "3 motile", "4", "5"))
    )
  cat("Section D: loaded composition from", composition_path, "\n")
  cat("  Rows:", nrow(gmm_composition), "\n\n")
} else {
  if (!exists("gmm_composition")) {
    stop("gmm_composition object not found.\n",
         "Run Section 1d first or set use_saved_composition = TRUE.")
  }
  cat("Section D: using gmm_composition from current session.\n\n")
}

if (!exists("gmm_composition") || nrow(gmm_composition) == 0) {
  cat("\nSection D skipped — gmm_composition not available.\n")
} else {
  
  
  cat("\n=== SECTION D: GMM COMPOSITIONAL ANALYSIS ===\n")
  
  # # --- Colour palette for GMM components ------------------------------------
  # # Extend if more than 5 components exist in your GMM solution
  # n_components  <- length(unique(gmm_composition$GMM_component))
  # gmm_colours <- c(
  #   `1` = "#457B9D",   
  #   `2` = "#E63946",   
  #   `3` = "#2A9D8F",   
  #   `4` = "#F4A261",
  #   `5` = "#6A4C93"
  # )
  # 
  # --- Colour palette -------------------------------------------------------
  # Low-displacement subclass: lighter shade of parent component colour
  # Motile subclass: full colour of parent component
  # Stacking order mirrors factor levels: low below, motile above
  gmm_colours <- c(
    "1"        = "#A8C8E0",   # comp 1 low displacement — light blue
    "1 motile" = "#457B9D",   # comp 1 motile           — full blue
    "2"        = "#E63946",   # comp 2
    "3"        = "#8ED0C7",   # comp 3 low displacement — light teal
    "3 motile" = "#2A9D8F",   # comp 3 motile           — full teal
    "4"        = "#F4A261",   # comp 4
    "5"        = "#6A4C93"    # comp 5
  )
  
  # --- Prepare data ---------------------------------------------------------
  # Create ordered Treatment_Day factor for x-axis
  # Order: Mock D1 | RP D1 | Mock D2 | RP D2
  # This groups by Day so within-Day Treatment comparison is immediate
  gmm_plot_data <- gmm_composition %>%
    mutate(
      Treatment_Day = factor(
        paste(Treatment, Day),
        levels = c("Mock D1", "RP D1", "Mock D2", "RP D2")
      ),
      Location = factor(Location, levels = c("Core", "CoreEdge", "Edge")),
      GMM_component = factor(GMM_component,
                             levels = c("1", "1 motile", "2",
                                        "3", "3 motile", "4", "5"))
    )
  
  # --- 1. STACKED BAR CHART -------------------------------------------------
  
  pdf(paste0("gmm_composition_barchart_", output_suffix,
             "_MaxDis", maxdis_motile_threshold, ".pdf"),
      width = 10, height = 5)
  p_bar <- ggplot(gmm_plot_data,
                  aes(x = Treatment_Day, y = pct,
                      fill = GMM_component)) +
    geom_bar(stat = "identity", position = "stack",
             width = 0.7, colour = "white", linewidth = 0.3) +
    geom_text(
      aes(label = ifelse(pct >= 8, paste0(pct, "%"), "")),
      position = position_stack(vjust = 0.5),
      size = 3, colour = "white", fontface = "bold"
    ) +
    # Vertical separator between D1 and D2 groups
    geom_vline(xintercept = 2.5,
               linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    annotate("text", x = 1.5, y = 105, label = "Day 1",
             size = 3.5, colour = "grey30", fontface = "bold") +
    annotate("text", x = 3.5, y = 105, label = "Day 2",
             size = 3.5, colour = "grey30", fontface = "bold") +
    facet_wrap(~ Location, nrow = 1) +
    scale_fill_manual(
      values = gmm_colours,
      name   = "GMM component",
      breaks = c("1", "1 motile", "2", "3", "3 motile", "4", "5"),
      labels = c(
        "1"        = "1 (<7 µm)",
        "1 motile" = "1 motile (>7 µm)",
        "2"        = "2",
        "3"        = "3 (<7 µm)",
        "3 motile" = "3 motile (>7 µm)",
        "4"        = "4",
        "5"        = "5"
      )
      ) +
    scale_y_continuous(limits = c(0, 107),
                       breaks = seq(0, 100, 25),
                       labels = function(x) paste0(x, "%")) +
    scale_x_discrete(labels = c("Mock D1" = "Mock",
                                "RP D1"   = "RP",
                                "Mock D2" = "Mock",
                                "RP D2"   = "RP")) +
    labs(
      title    = "GMM component composition by Treatment and Day",
      subtitle = paste0(
        "Components 1 and 3 split by MaxDis > ", maxdis_motile_threshold,
        " µm  |  Light = low displacement  |  Full = motile  |  suffix: ",
        output_suffix
      ),
      x       = NULL,
      y       = "% of tracks",
      caption = paste0(
        "Dashed line separates Day 1 (left) from Day 2 (right)\n",
        "Light shading = MaxDis ≤ ", maxdis_motile_threshold,
        " µm (low productive displacement within component)"
      )
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title         = element_text(face = "bold", size = 12),
      plot.subtitle      = element_text(size = 10, colour = "grey40"),
      plot.caption       = element_text(size = 9,  colour = "grey50"),
      axis.text.x        = element_text(size = 10),
      legend.position    = "right",
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank()
    )
  
  print(p_bar)
  dev.off()
  cat("Saved: gmm_composition_barchart_", output_suffix,
      "_MaxDis", maxdis_motile_threshold, ".pdf\n", sep = "")
  
}  # end GMM_component gate

# =============================================================================
# ---- SECTION D (ADDITION): MOUSE-AVERAGED MOTILE-FRACTION BARCHART ----
# =============================================================================
# APPEND this INSIDE the Section D gate, AFTER the existing stacked-bar
# dev.off() line and BEFORE the closing "}  # end GMM_component gate".
# The existing stacked chart is NOT modified.
#
# WHY A SEPARATE PLOT (not error bars on the stack):
#   Stacked bars cannot carry interpretable error bars — a segment's y-position
#   depends on the segments beneath it, and proportions summing to 100% makes
#   segment uncertainties mechanically dependent. So the motile fraction is
#   shown here as a SINGLE non-stacked bar per cell, with SEM + per-mouse dots.
#   (Weissgerber et al. 2015, PLoS Biol 13:e1002128 — show the data, not just
#   bar+error.)
#
# POOLING RULE (A): each mouse's motile% is computed at whatever grain
#   `group_by_vars` specifies (pool that mouse's OWN tracks within the grain),
#   THEN average across mice. Mouse = unit of replication. Changing
#   group_by_vars automatically re-pools. NOTE: pooled means will NOT equal the
#   arithmetic mean of finer sub-cell means when track counts are unequal
#   (mean-of-pooled != mean-of-means; Jensen's inequality) — this is correct.
#
# INPUT: tracks_comp_plot (must contain MouseID, Location, Treatment, Day,
#   GMM_component_plot). The pre-saved gmm_composition CSV canNOT be used here
#   because it has no MouseID.
# =============================================================================

#####!!!! MAKE SURE tracks_Comp_plot is present in 01c/d!!! #####
##### IF YOU NEED TO FACET OTHERWISE CHANGE D2D AND X ELEMENTS WITHING GGPLOT2

if (exists("tracks_comp_plot")) {
  
  # --- USER INPUT: which grouping to plot ------------------------------------
  # Edit this vector. "Treatment" should stay (it is the x-axis); any other
  # variables become facets (2-simple layout).
  #   c("Location","Treatment","Day") -> x=Treatment, facets = Location x Day
  #   c("Treatment","Day")            -> location-pooled, facet = Day
  #   c("Treatment","Location")       -> day-pooled,      facet = Location
  #   c("Treatment")                  -> fully pooled,    single panel
  #group_by_vars <- c("Treatment")
  group_by_vars <- c("Location")
  
  # --- D2a. Derive per-track motile flag -------------------------------------
  ma_dat <- tracks_comp_plot
  ma_dat$motile <- as.integer(
    ma_dat$GMM_component_plot %in% c("1 motile", "3 motile")
  )
  
  # --- D2b. POOLING RULE (A): per-mouse % at the chosen grain -----------------
  # Group by (group_by_vars + MouseID), pool that mouse's tracks, compute its %.
  per_mouse <- ma_dat %>%
    group_by(across(all_of(c(group_by_vars, "MouseID")))) %>%
    summarise(n_tracks = n(),
              motile_pct = 100 * mean(motile), .groups = "drop")
  
  # --- D2c. Average across mice + SEM -----------------------------------------
  # SEM = SD / sqrt(n_mice); undefined (NA) when n_mice == 1.
  summ <- per_mouse %>%
    group_by(across(all_of(group_by_vars))) %>%
    summarise(n_mice   = n(),
              mean_pct = mean(motile_pct),
              sd_pct   = sd(motile_pct),
              sem_pct  = ifelse(n_mice > 1, sd_pct / sqrt(n_mice), NA_real_),
              .groups  = "drop")
  
  # Ensure factor ordering where those columns are present
  if ("Location"  %in% names(summ))
    summ$Location  <- factor(summ$Location,  levels = c("Core","CoreEdge","Edge"))
  if ("Treatment" %in% names(summ))
    summ$Treatment <- factor(summ$Treatment, levels = c("Mock","RP"))
  if ("Day"       %in% names(summ))
    summ$Day       <- factor(summ$Day,       levels = c("D1","D2"))
  if ("Location"  %in% names(per_mouse))
    per_mouse$Location  <- factor(per_mouse$Location,  levels = c("Core","CoreEdge","Edge"))
  if ("Treatment" %in% names(per_mouse))
    per_mouse$Treatment <- factor(per_mouse$Treatment, levels = c("Mock","RP"))
  if ("Day"       %in% names(per_mouse))
    per_mouse$Day       <- factor(per_mouse$Day,       levels = c("D1","D2"))
  
  # --- D2d. Build facet spec from leftover vars (2-simple) --------------------
  facet_vars <- setdiff(group_by_vars, "Location")
  facet_layer <- if (length(facet_vars) == 0) {
    NULL
  } else {
    facet_wrap(reformulate(facet_vars), nrow = 1)
  }
  
  trt_colours <- c("Mock" = "#4C72B0", "RP" = "#C44E52")
  # loc_colours <- c("Core" = "#4C72B0", "CoreEdge" = "#C44E52","Edge" = "#F4A261")
  
  # --- D2e. Plot: bars (mean) + SEM + per-mouse dots --------------------------
  p_motile <- ggplot(summ, aes(x = Location, y = mean_pct, fill = Location)) +
    geom_col(width = 0.65, colour = "grey30", linewidth = 0.3, alpha = 0.85) +
    geom_errorbar(aes(ymin = pmax(mean_pct - sem_pct, 0),
                      ymax = mean_pct + sem_pct),
                  width = 0.2, linewidth = 0.5, na.rm = TRUE) +
    geom_jitter(data = per_mouse,
                aes(x = Location, y = motile_pct),
                inherit.aes = FALSE,
                width = 0.12, height = 0, size = 1.8,
                shape = 21, fill = "white", colour = "black", stroke = 0.4) +
    { if (!is.null(facet_layer)) facet_layer } +
    scale_fill_manual(values = trt_colours, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(
      title    = "Mouse-averaged motile fraction",
      subtitle = paste0("Grouping: ", paste(group_by_vars, collapse = " x "),
                        "  |  bar = mean across mice, error = +/- SEM, ",
                        "dots = per-mouse %  |  suffix: ", output_suffix),
      x = NULL, y = "Motile tracks (%)",
      caption = "Motile = GMM components '1 motile' + '3 motile' (MaxDis > 7 um). Mouse = replication unit. n=1 cells have no SEM."
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9,  colour = "grey40"),
          plot.caption  = element_text(size = 8,  colour = "grey50"),
          strip.text    = element_text(face = "bold"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor   = element_blank())
  
  ggsave(paste0("motile_fraction_mouseavg_barchart_", output_suffix,
                "_MaxDis", maxdis_motile_threshold, ".jpg"),
         plot = p_motile, width = 4.5, height = 4.5, dpi = 300)

  # # Save the underlying numbers too (audit trail)
  # write_csv(summ,
  #           paste0("motile_fraction_mouseavg_summary_", output_suffix, ".csv"))
  # write_csv(per_mouse,
  #           paste0("motile_fraction_per_mouse_", output_suffix, ".csv"))
  # 
  # cat("Saved: motile_fraction_mouseavg_barchart_", output_suffix,
  #     "_MaxDis", maxdis_motile_threshold, ".jpeg\n", sep = "")
  # cat("Saved: motile_fraction_mouseavg_summary_", output_suffix, ".csv\n", sep = "")
  # cat("Saved: motile_fraction_per_mouse_", output_suffix, ".csv\n", sep = "")

} else {
  cat("\nMouse-averaged motile barchart skipped — tracks_comp_plot not found.\n")
}

######## PER MOUSE PLOTS #############

if (exists("per_mouse")) {
  
  # Treatment must be present to dodge Mock vs RP within a mouse.
  if (!("Treatment" %in% names(per_mouse))) {
    cat("Per-mouse barchart skipped — 'Treatment' not in group_by_vars.\n")
  } else {
    
    # Facet by everything EXCEPT Treatment and MouseID (2-simple).
    pm_facet_vars <- setdiff(group_by_vars, "Treatment")
    pm_facet_layer <- if (length(pm_facet_vars) == 0) NULL else
      facet_wrap(reformulate(pm_facet_vars), nrow = 1)
    
    trt_colours <- c("Mock" = "#4C72B0", "RP" = "#C44E52")
    
    p_permouse <- ggplot(per_mouse,
                         aes(x = MouseID, y = motile_pct, fill = Treatment)) +
      geom_col(position = position_dodge(width = 0.8, preserve = "single"),
               width = 0.7, colour = "grey30", linewidth = 0.3) +
      { if (!is.null(pm_facet_layer)) pm_facet_layer } +
      scale_fill_manual(values = trt_colours, name = "Treatment") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
      labs(
        title    = "Per-mouse motile fraction",
        subtitle = paste0("Grouping: ", paste(group_by_vars, collapse = " x "),
                          "  |  each bar = one mouse's motile % (tracks pooled ",
                          "at this grain)  |  suffix: ", output_suffix),
        x = "Mouse", y = "Motile tracks (%)",
        caption = "Motile = '1 motile' + '3 motile' (MaxDis > 7 um). Missing Mock/RP bar = that mouse had no tracks in that cell."
      ) +
      theme_bw(base_size = 11) +
      theme(plot.title    = element_text(face = "bold", size = 12),
            plot.subtitle = element_text(size = 9,  colour = "grey40"),
            plot.caption  = element_text(size = 8,  colour = "grey50"),
            strip.text    = element_text(face = "bold"),
            axis.text.x   = element_text(size = 9),
            panel.grid.major.x = element_blank(),
            panel.grid.minor   = element_blank())
    
    ggsave(paste0("motile_fraction_per_mouse_barchart_", output_suffix,
                  "_MaxDis", maxdis_motile_threshold, ".jpg"),
           plot = p_permouse, width = 10, height = 4.5, dpi = 300)
    
    cat("Saved: motile_fraction_per_mouse_barchart_", output_suffix,
        "_MaxDis", maxdis_motile_threshold, ".jpg\n", sep = "")
  }
  
} else {
  cat("\nPer-mouse barchart skipped — per_mouse object not found ",
      "(run the D2 mouse-averaged block first).\n", sep = "")
}

# =============================================================================
# STANDALONE: PER-MOUSE D1->D2 MOTILE-FRACTION TRAJECTORY
# =============================================================================
# Self-contained. Append anywhere after tracks_comp_plot exists (e.g. after the
# Section D bar-chart additions). Does NOT depend on `per_mouse` or
# `group_by_vars` — computes its own per-mouse Treatment x Day means, so it
# always shows the same D1->D2 trajectory regardless of bar-chart settings.
#
# PURPOSE
#   Show whether the D1->D2 numerical increase in motile fraction is CONSISTENT
#   across mice or driven by a few. This directly illustrates WHY the Day effect
#   was non-significant (Day p ~ 0.12 in the Treatment*Day GLMM): inconsistent
#   per-mouse direction is the visual explanation of a null main effect.
#   Faceted by Treatment so all mice's trajectories are visible within each arm.
#
# HONESTY CAVEAT (in caption): D1 and D2 are SEPARATE slices from each mouse,
#   NOT the same cells tracked longitudinally. Dashed connectors are a visual
#   aid only, not a paired trajectory. Framed descriptively — the underlying
#   Day increase is NOT statistically supported.
# =============================================================================

if (exists("tracks_comp_plot")) {
  
  # --- T1. Per-mouse motile % at Treatment x Day grain (pooling rule A) -------
  traj_dat <- tracks_comp_plot
  traj_dat$motile <- as.integer(
    traj_dat$GMM_component_plot %in% c("1 motile", "3 motile"))
  
  per_mouse_td <- traj_dat %>%
    group_by(Treatment, Day, MouseID) %>%
    summarise(n_tracks = n(),
              motile_pct = 100 * mean(motile), .groups = "drop") %>%
    mutate(Treatment = factor(Treatment, levels = c("Mock", "RP")),
           Day       = factor(Day,       levels = c("D1", "D2")))
  
  # --- T1b. Mean across mice + SEM per Treatment x Day (for the bars) ---------
  summ_td <- per_mouse_td %>%
    group_by(Treatment, Day) %>%
    summarise(n_mice   = n(),
              mean_pct = mean(motile_pct),
              sem_pct  = ifelse(n_mice > 1, sd(motile_pct) / sqrt(n_mice), NA_real_),
              .groups  = "drop")
  
  # --- T2. Plot: bars + SEM (bottom layer) then mouse trajectories (top) ------
  p_traj <- ggplot() +
    # bars = group mean per Day, light grey fill so coloured dots read on top
    geom_col(data = summ_td,
             aes(x = Day, y = mean_pct),
             width = 0.6, fill = "grey85", colour = "grey40", linewidth = 0.3) +
    geom_errorbar(data = summ_td,
                  aes(x = Day,
                      ymin = pmax(mean_pct - sem_pct, 0),
                      ymax = mean_pct + sem_pct),
                  width = 0.15, linewidth = 0.5, colour = "grey30", na.rm = TRUE) +
    # per-mouse dashed trajectories + coloured points on top
    geom_line(data = per_mouse_td,
              aes(x = Day, y = motile_pct, group = MouseID, colour = MouseID),
              linetype = "dashed", linewidth = 0.5, alpha = 0.85) +
    geom_point(data = per_mouse_td,
               aes(x = Day, y = motile_pct, colour = MouseID),
               size = 2.4) +
    facet_wrap(~ Treatment, nrow = 1) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(
      title    = "Motile fraction D1 -> D2: group mean (bars) + per-mouse trajectories",
      subtitle = paste0("Bars = mean across mice +/- SEM; dots/lines = per-mouse ",
                        "motile %  |  suffix: ", output_suffix),
      x = "Day", y = "Motile tracks (%)", colour = "Mouse",
      caption = paste0("D1 and D2 are SEPARATE slices per mouse, not ",
                       "longitudinally tracked cells; dashed lines are a visual ",
                       "aid only.\nDescriptive: the day-associated increase was ",
                       "not statistically supported (Day p ~ 0.12). n=1 cells have no SEM.")
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", size = 11.5),
          plot.subtitle = element_text(size = 9,  colour = "grey40"),
          plot.caption  = element_text(size = 8,  colour = "grey50"),
          strip.text    = element_text(face = "bold"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor   = element_blank())
  
  ggsave(paste0("motile_fraction_trajectory_bars_", output_suffix,
                "_MaxDis", maxdis_motile_threshold, ".jpg"),
         plot = p_traj, width = 8, height = 4.5, dpi = 300)
  
  write_csv(per_mouse_td,
            paste0("motile_fraction_permouse_trajectory_", output_suffix, ".csv"))
  
  cat("Saved: motile_fraction_trajectory_bars_", output_suffix,
      "_MaxDis", maxdis_motile_threshold, ".jpg\n", sep = "")
  cat("Saved: motile_fraction_permouse_trajectory_", output_suffix, ".csv\n", sep = "")
  
} else {
  cat("\nTrajectory plot skipped — tracks_comp_plot not found.\n")
}

# =============================================================================
# ---- SECTION E: BINOMIAL GLMM — COMPOSITION (MOTILE vs STATIC PROPORTION) ----
# =============================================================================
# APPEND to the END of 04_treatment_model_v2.R, AFTER the composition plots
# (Section D). Nothing above this section is modified. Self-contained.
#
# PURPOSE
#   Inferential, properly-clustered counterpart to the descriptive
#   mouse-averaged composition percentages. Tests whether RP shifts the
#   PROPORTION of motile tracks (and whether that shift differs by Location),
#   accounting for the Mouse/Tumour/Slice/ROI nesting that a raw percentage
#   ignores. This is the correct instrument for the composition claim: it
#   models a binary motile/static outcome rather than inferring composition
#   indirectly from a speed model.
#
# MODEL
#   motile ~ Treatment * Location + Day + (1 | MouseID/TumourID_global/
#            SliceID_global/ROIID_global),  family = binomial
#   Outcome (definition a): motile = 1 for GMM_component_plot in
#   {"1 motile","3 motile"}; 0 for static components {1,3,4,5}.
#   Fitted on tracks_comp_plot (full/unfiltered GMM1_3_4_5 population).
#
# METHOD REFERENCES
#   Bolker BM et al. (2009) Generalized linear mixed models: a practical
#     guide for ecology and evolution. Trends Ecol Evol 24:127-135.
#   Agresti A (2013) Categorical Data Analysis, 3rd ed., ch.13 (GLMMs for
#     nested binary data).
#   Harrison XA et al. (2018) A brief introduction to mixed effects modelling
#     and multi-model inference in ecology. PeerJ 6:e4794 (convergence /
#     singularity diagnostics for binomial GLMMs).
#   Odds ratios via emmeans(..., type="response"): Lenth R, emmeans vignette.
# =============================================================================

cat("\n=== SECTION E: COMPOSITION BINOMIAL GLMM ===\n")

library(lme4)      # glmer
library(emmeans)
# (broom.mixed optional — not required here)

# --- E1. DERIVE BINARY OUTCOME ----------------------------------------------
# Use the EXISTING GMM_component_plot column (do NOT re-derive from MaxDis,
# to avoid any threshold mismatch). motile = 1 only for the displacement-
# positive subclasses of components 1 and 3.

glmm_data <- tracks_comp_plot   #READ ME: RUN THE FIRST PORTION TO HAVE THE DATA
glmm_data$motile <- as.integer(
  glmm_data$GMM_component_plot %in% c("1 motile", "3 motile")
)

cat("  Motile tracks:", sum(glmm_data$motile),
    "/", nrow(glmm_data),
    sprintf("(%.1f%%)\n", 100 * mean(glmm_data$motile)))
cat("  Motile % by Treatment x Location (raw, unclustered — sanity check):\n")
print(round(100 * tapply(glmm_data$motile,
                         list(glmm_data$Treatment, glmm_data$Location),
                         mean), 1))

# --- E2. FIXED-EFFECTS FORMULA (EDIT THIS LINE to change the model) ---------
# To run Location*Day later, replace the RHS below, e.g.:
#   "motile ~ Location * Day + Treatment + (random)"
# or the three-way "motile ~ Treatment * Location * Day + (random)"
# (three-way likely singular with ~5 mice — see note at end).
glmm_random  <- "(1 | MouseID/TumourID_global/SliceID_global/ROIID_global)"
glmm_formula <- as.formula(
  paste("motile ~ Treatment * Location + Day +", glmm_random)
)

# --- E3. FIT ----------------------------------------------------------------
# bobyqa + raised iterations: standard first-line remedy for glmer
# convergence warnings (Bolker et al. 2009; Harrison et al. 2018).
comp_glmm <- tryCatch(
  glmer(glmm_formula, data = glmm_data, family = binomial,
        control = glmerControl(optimizer = "bobyqa",
                               optCtrl = list(maxfun = 2e5))),
  error   = function(e) { cat("  GLMM error:", e$message, "\n"); NULL },
  warning = function(w) {
    cat("  GLMM warning:", conditionMessage(w), "\n")
    suppressWarnings(
      glmer(glmm_formula, data = glmm_data, family = binomial,
            control = glmerControl(optimizer = "bobyqa",
                                   optCtrl = list(maxfun = 2e5)))
    )
  }
)

# --- E4. DIAGNOSTICS + SUMMARY ----------------------------------------------
comp_glmm_file <- paste0("composition_glmm_", output_suffix, ".txt")
if (!is.null(comp_glmm)) {
  sink(comp_glmm_file)
  cat("COMPOSITION BINOMIAL GLMM —", output_suffix, "\n")
  cat("Formula:\n"); print(glmm_formula)
  cat("\nSingular fit:", isSingular(comp_glmm), "\n")
  cat("\n--- Random effects (variance partition) ---\n")
  print(VarCorr(comp_glmm))
  cat("\n--- Fixed effects (log-odds scale) ---\n")
  print(summary(comp_glmm)$coefficients)
  cat("\n--- Type III Wald test of fixed effects ---\n")
  print(car::Anova(comp_glmm, type = 3))
  sink()
  
  if (isSingular(comp_glmm))
    cat("  NOTE: singular fit — a higher-level variance collapsed to ~0.",
        "Odds ratios remain valid; interpret variance partition with caution",
        "(expected with few mice; Bolker et al. 2009).\n")
  
  # --- E5. ODDS RATIOS: Treatment within Location -------------------------
  # Mock vs RP odds of being motile, conditional on each Location. This is the
  # region-specific composition effect (the direct test of "RP -> more motile
  # in Core"). type="response" back-transforms logits to odds ratios.
  emm_trt_by_loc <- tryCatch(
    emmeans(comp_glmm, ~ Treatment | Location),
    error = function(e) { cat("  EMM failed:", e$message, "\n"); NULL }
  )
  if (!is.null(emm_trt_by_loc)) {
    or_tbl <- as.data.frame(
      pairs(emm_trt_by_loc, type = "response", infer = TRUE, adjust = "tukey")
    )
    write_csv(or_tbl,
              paste0("composition_glmm_OR_trt_by_location_", output_suffix, ".csv"))
    cat("\n--- Odds ratios: Treatment within Location ---\n")
    print(or_tbl)
    
    # Estimated motile probabilities per cell (for plotting / reporting)
    prob_tbl <- as.data.frame(
      emmeans(comp_glmm, ~ Treatment | Location, type = "response")
    )
    write_csv(prob_tbl,
              paste0("composition_glmm_probs_", output_suffix, ".csv"))
  }
  
  # --- E5b. ODDS RATIOS: Location within Treatment ------------------------
  # Location gradient (CoreEdge/Core, Edge/Core, Edge/CoreEdge) estimated
  # separately within Mock and RP. reverse = TRUE gives "relative to Core";
  # Tukey adjusts across the 3 Location contrasts within each arm.
  emm_loc_by_trt <- tryCatch(
    emmeans(comp_glmm, ~ Location | Treatment),
    error = function(e) { cat("  EMM (Location | Treatment) failed:", e$message, "\n"); NULL }
  )
  if (!is.null(emm_loc_by_trt)) {
    loc_or_tbl <- as.data.frame(
      pairs(emm_loc_by_trt, type = "response", infer = TRUE,
            adjust = "tukey", reverse = TRUE)
    )
    write_csv(loc_or_tbl,
              paste0("composition_glmm_OR_loc_by_treatment_", output_suffix, ".csv"))
    cat("\n--- Odds ratios: Location within Treatment ---\n")
    print(loc_or_tbl)
  }
  
  cat("\nSaved: composition_glmm_", output_suffix, ".txt\n", sep = "")
  cat("Saved: composition_glmm_OR_trt_by_location_", output_suffix, ".csv\n", sep = "")
  cat("Saved: composition_glmm_probs_", output_suffix, ".csv\n", sep = "")
  cat("Saved: composition_glmm_OR_loc_by_treatment_", output_suffix, ".csv\n", sep = "")
} else {
  cat("  GLMM did not fit — section skipped.\n")
}

# --- NOTE ON EXTENSIONS ------------------------------------------------------
# Location*Day: edit the RHS in E2. The three-way Treatment*Location*Day is
# the day-resolved composition test but will likely be singular at ~5 mice
# (as the speed Model 3 was). If it won't converge, fall back to fitting
# Day-stratified two-way models, or a per-ROI beta-binomial (Harrison 2018).




# =============================================================================
# SECTION F: DAY-STRATIFIED COMPOSITION GLMM (middle-ground Day diagnostic)
# =============================================================================
# APPEND after Section E (the composition GLMM) in 04_treatment_model_v2.R.
# Depends on `glmm_data` created in Section E (which carries the derived
# `motile` column) — Section E always runs first, so no re-derivation needed.
# Nothing above is modified.
#
# PURPOSE
#   The Section E model has Day as an ADDITIVE term, so it cannot say whether
#   RP's motile-shift differs between D1 and D2. A full Treatment*Location*Day
#   three-way GLMM is the formal test but is singularity-prone at ~5 mice on a
#   binary outcome. Instead we STRATIFY by Day: fit the SAME two-way model
#   (Treatment*Location) separately on the D1 and D2 subsets. If RP's Core
#   effect (OR ~2.75 in the pooled model) appears in D1 and washes out in D2,
#   that descriptively localises the effect to D1 — where the descriptive
#   composition data was stable — WITHOUT asking one model to estimate an
#   underpowered three-way interaction.
#
#   This is a DESCRIPTIVE localisation, not a formal test of the Day
#   difference (that would require the interaction term). Reported as such.
#
# EXPECTED: splitting by Day halves the clusters, so the full 4-level nesting
#   may go singular within a Day (esp. D2). Reported honestly, NOT auto-
#   simplified — a collapsed level is itself informative (Bolker et al. 2009).
#
# REFERENCES
#   Bolker BM et al. (2009) Trends Ecol Evol 24:127-135 (GLMM practical guide;
#     interaction terms need replication at the relevant level).
#   Gelman A & Stern H (2006) Am Stat 60:328 (a non-significant interaction
#     does not establish effects differ; stratification as description).
#   Agresti A (2013) Categorical Data Analysis 3rd ed., ch.13 (odds ratios).
# =============================================================================

cat("\n=== SECTION F: DAY-STRATIFIED COMPOSITION GLMM ===\n")

# --- F1. FIT TWO-WAY MODEL WITHIN ONE DAY -----------------------------------
# Same random structure as Section E. Returns list(model, singular, day) or
# NULL on failure. bobyqa + raised iterations as in Section E.
fit_day_glmm <- function(day_level) {
  d <- droplevels(subset(glmm_data, Day == day_level))
  cat("  Day", day_level, "- n tracks:", nrow(d),
      "| motile:", sum(d$motile),
      sprintf("(%.1f%%)\n", 100 * mean(d$motile)))
  
  m <- tryCatch(
    glmer(motile ~ Treatment * Location +
            (1 | MouseID/TumourID_global/SliceID_global/ROIID_global),
          data = d, family = binomial,
          control = glmerControl(optimizer = "bobyqa",
                                 optCtrl = list(maxfun = 2e5))),
    error   = function(e) { cat("    fit error:", e$message, "\n"); NULL },
    warning = function(w) {
      cat("    fit warning:", conditionMessage(w), "\n")
      suppressWarnings(
        glmer(motile ~ Treatment * Location +
                (1 | MouseID/TumourID_global/SliceID_global/ROIID_global),
              data = d, family = binomial,
              control = glmerControl(optimizer = "bobyqa",
                                     optCtrl = list(maxfun = 2e5)))
      )
    }
  )
  if (is.null(m)) return(NULL)
  sing <- isSingular(m)
  if (sing) cat("    NOTE: singular fit for Day", day_level,
                "- a variance level collapsed (expected with fewer clusters",
                "per Day); ORs valid, variance partition interpret with care.\n")
  list(model = m, singular = sing, day = day_level)
}

# --- F2. EXTRACT Mock/RP ODDS RATIOS WITHIN LOCATION, TAG WITH DAY -----------
or_within_day <- function(fit) {
  if (is.null(fit)) return(NULL)
  emm <- tryCatch(emmeans(fit$model, ~ Treatment | Location),
                  error = function(e) NULL)
  if (is.null(emm)) return(NULL)
  as.data.frame(pairs(emm, type = "response", infer = TRUE, adjust = "tukey")) %>%
    mutate(Day = fit$day, Singular = fit$singular, .before = 1)
}

# --- F3. RUN BOTH DAYS, STACK INTO ONE CSV ----------------------------------
day_levels <- levels(factor(glmm_data$Day))     # typically c("D1","D2")
day_fits   <- lapply(day_levels, fit_day_glmm)

or_stacked <- bind_rows(lapply(day_fits, or_within_day))

if (!is.null(or_stacked) && nrow(or_stacked) > 0) {
  write_csv(or_stacked,
            paste0("composition_glmm_OR_trt_by_location_BY_DAY_",
                   output_suffix, ".csv"))
  cat("\n--- Mock/RP odds ratios within Location, stratified by Day ---\n")
  print(or_stacked)
  cat("\nSaved: composition_glmm_OR_trt_by_location_BY_DAY_",
      output_suffix, ".csv\n", sep = "")
} else {
  cat("  No stratified ORs produced (both Day fits failed).\n")
}

# --- INTERPRETATION NOTE -----------------------------------------------------
# Compare the Core Mock/RP OR between D1 and D2 rows. If Core OR < 1 (RP more
# motile) and its CI excludes 1 in D1 but not D2, that descriptively places
# the RP composition shift in D1. Do NOT interpret a D1-vs-D2 difference as a
# tested effect — the interaction was not fitted here (Gelman & Stern 2006).

# =============================================================================
# --- SECTION E-DAY: COMPOSITION GLMM — Treatment*Day INTERACTION (Location covar)----
# =============================================================================
# APPEND after Section E (and Section F if present) in 04_treatment_model_v2.R.
# Reuses `glmm_data` from Section E (carries the derived `motile` column).
# Self-contained; nothing existing is modified.
#
# PURPOSE
#   Complements the Section E Treatment*Location model by asking the DAY
#   question instead: does the motile-proportion CHANGE from D1 to D2, and does
#   that change differ between Mock and RP? Location kept as an ADDITIVE
#   covariate (protects the Day estimate from Location-sampling confounds, same
#   logic as keeping Location in the speed EMMs). Does NOT address Core-
#   specificity (that is Section E's job).
#
# MODEL
#   motile ~ Treatment * Day + Location + (1 | MouseID/TumourID_global/
#            SliceID_global/ROIID_global),  family = binomial
#
# INTERPRETATION FRAMING (read BEFORE the numbers)
#   - Treatment main effect  : does RP shift motile odds overall (Day-averaged)?
#   - Day main effect        : does motile odds change D1->D2 overall?
#   - Treatment:Day interaction: does the D1->D2 CHANGE differ by treatment?
#       * ns interaction = NO EVIDENCE the arms differ over time (supports
#         "motility rose similarly in both arms"). It is an ABSENCE of evidence
#         for a difference, NOT evidence of a difference (Gelman & Stern 2006,
#         Am Stat 60:328).
#   - Conditional contrasts:
#       Mock/RP | Day  = treatment effect within each day
#       D1/D2  | Treatment = day change within each arm ("increase in both?")
#     A significant D1/D2 within BOTH Mock and RP = motility rose in both arms.
#
# REFERENCES: Agresti (2013) Categorical Data Analysis 3rd ed. ch.13;
#   Bolker et al. (2009) TREE 24:127; Gelman & Stern (2006) Am Stat 60:328;
#   emmeans conditional/marginal contrasts — Lenth, emmeans vignette.
# =============================================================================

cat("\n=== SECTION E-DAY: COMPOSITION GLMM (Treatment*Day + Location) ===\n")

if (!exists("glmm_data")) {
  cat("  Skipped — glmm_data not found (run Section E first).\n")
} else {
  
  library(lme4); library(emmeans)
  
  # --- ED1. FIT --------------------------------------------------------------
  comp_glmm_day <- tryCatch(
    glmer(motile ~ Treatment * Day + Location +
            (1 | MouseID/TumourID_global/SliceID_global/ROIID_global),
          data = glmm_data, family = binomial,
          control = glmerControl(optimizer = "bobyqa",
                                 optCtrl = list(maxfun = 2e5))),
    error   = function(e) { cat("  fit error:", e$message, "\n"); NULL },
    warning = function(w) {
      cat("  fit warning:", conditionMessage(w), "\n")
      suppressWarnings(
        glmer(motile ~ Treatment * Day + Location +
                (1 | MouseID/TumourID_global/SliceID_global/ROIID_global),
              data = glmm_data, family = binomial,
              control = glmerControl(optimizer = "bobyqa",
                                     optCtrl = list(maxfun = 2e5))))
    }
  )
  
  if (is.null(comp_glmm_day)) {
    cat("  GLMM did not fit — Section E-Day skipped.\n")
  } else {
    
    # --- ED2. DIAGNOSTICS + WALD TESTS --------------------------------------
    day_file <- paste0("composition_glmm_DAY_", output_suffix, ".txt")
    sink(day_file)
    cat("COMPOSITION GLMM Treatment*Day + Location —", output_suffix, "\n")
    cat("Singular fit:", isSingular(comp_glmm_day), "\n\n")
    cat("--- Random effects ---\n");         print(VarCorr(comp_glmm_day))
    cat("\n--- Fixed effects (log-odds) ---\n")
    print(summary(comp_glmm_day)$coefficients)
    cat("\n--- Type III Wald tests (main effects + interaction) ---\n")
    print(car::Anova(comp_glmm_day, type = 3))
    sink()
    if (isSingular(comp_glmm_day))
      cat("  NOTE: singular fit — ORs valid, variance partition cautious.\n")
    
    # --- ED3. MARGINAL contrasts (Day-averaged / Treatment-averaged) --------
    # Treatment marginal: RP vs Mock averaged over Day (and Location).
    emm_trt_marg <- emmeans(comp_glmm_day, ~ Treatment)
    or_trt_marg  <- as.data.frame(pairs(emm_trt_marg, type = "response",
                                        infer = TRUE)) %>%
      mutate(Contrast_type = "Treatment marginal (over Day, Location)")
    
    # Day marginal: D2 vs D1 averaged over Treatment (and Location).
    emm_day_marg <- emmeans(comp_glmm_day, ~ Day)
    or_day_marg  <- as.data.frame(pairs(emm_day_marg, type = "response",
                                        infer = TRUE)) %>%
      mutate(Contrast_type = "Day marginal (over Treatment, Location)")
    
    # --- ED4. CONDITIONAL contrasts -----------------------------------------
    # Mock/RP within each Day (treatment effect per day).
    emm_trt_by_day <- emmeans(comp_glmm_day, ~ Treatment | Day)
    or_trt_by_day  <- as.data.frame(pairs(emm_trt_by_day, type = "response",
                                          infer = TRUE, adjust = "tukey")) %>%
      mutate(Contrast_type = "Treatment within Day")
    
    # D1/D2 within each Treatment (day change per arm — "increase in both?").
    emm_day_by_trt <- emmeans(comp_glmm_day, ~ Day | Treatment)
    or_day_by_trt  <- as.data.frame(pairs(emm_day_by_trt, type = "response",
                                          infer = TRUE, adjust = "tukey")) %>%
      mutate(Contrast_type = "Day within Treatment")
    
    # --- ED5. STACK + SAVE ---------------------------------------------------
    or_all <- bind_rows(or_trt_marg, or_day_marg,
                        or_trt_by_day, or_day_by_trt)
    write_csv(or_all,
              paste0("composition_glmm_DAY_odds_ratios_", output_suffix, ".csv"))
    
    # Estimated motile probabilities per Treatment x Day cell (for plotting)
    prob_day <- as.data.frame(
      emmeans(comp_glmm_day, ~ Treatment | Day, type = "response"))
    write_csv(prob_day,
              paste0("composition_glmm_DAY_probs_", output_suffix, ".csv"))
    
    cat("\n--- Marginal + conditional odds ratios (see Contrast_type) ---\n")
    print(or_all)
    cat("\nSaved: composition_glmm_DAY_", output_suffix, ".txt\n", sep = "")
    cat("Saved: composition_glmm_DAY_odds_ratios_", output_suffix, ".csv\n", sep = "")
    cat("Saved: composition_glmm_DAY_probs_", output_suffix, ".csv\n", sep = "")
    
    cat("\nREAD: check the Treatment:Day interaction p in the .txt FIRST.\n")
    cat("If ns, the marginal Treatment/Day effects are the ones to report and\n")
    cat("the conditional contrasts describe (not test) the per-cell pattern.\n")
  }
}
# ==========================================================
#### MOTILE PROPORTION VS CYTOKINE CORRELATIONS #### 
##!!!! have to run section 1 to load tracks !!!
# ==========================================================
# =============================================================================
#  ---- SLICE-AVERAGED MOTILE PROPORTIONS  (for cytokine correlation) ------
#  Pools all tracks within each SliceID_global, computes motile fraction
#  (motile = GMM comp 1 or 3 with MaxDis > threshold = "1 motile"/"3 motile").
# =============================================================================
maxdis_motile_threshold <- 7                      # µm, matches composition analysis
motile_classes <- c("1 motile", "3 motile")

slice_props <- tracks %>%
  mutate(GMM_component_plot = case_when(
    GMM_component == 1 & MaxDis >  maxdis_motile_threshold ~ "1 motile",
    GMM_component == 1 & MaxDis <= maxdis_motile_threshold ~ "1",
    GMM_component == 3 & MaxDis >  maxdis_motile_threshold ~ "3 motile",
    GMM_component == 3 & MaxDis <= maxdis_motile_threshold ~ "3",
    TRUE ~ as.character(GMM_component)
  )) %>%
  group_by(SliceID_global,Day) %>%                    # add , Day here if cytokines are day-specific
  summarise(
    TumourID_global = first(TumourID_global),
    Treatment       = first(Treatment),
    n_total     = n(),
    n_motile    = sum(GMM_component_plot %in% motile_classes),
    motile_frac = n_motile / n_total,
    motile_pct  = round(100 * motile_frac, 1),
    .groups = "drop"
  ) %>%
  arrange(TumourID_global, SliceID_global)

print(slice_props, n = Inf)
write_csv(slice_props, paste0("slice_motile_proportions_", output_suffix, ".csv"))


##### Need to update the spreasheet below
slice_data <- read_csv("slice_data_with_cytokines_prop.csv")
stopifnot(nrow(slice_data) == nrow(slice_props))    

# =============================================================================
#  ---- SLICE-LEVEL CORRELATION: slice cytokines vs slice motile fraction ----
#  Exploratory. Spearman rho + p (cor.test).
# =============================================================================
library(dplyr); library(purrr); library(tibble)

cytokine_cols <- setdiff(names(slice_data), names(slice_props))
cat("Detected cytokine columns:\n"); print(cytokine_cols)      # <-- confirm these
# manual fallback: cytokine_cols <- c("CXCL10","CCL5", ...)

# slice_cor_data <- slice_data                       # READ ME: all slices !
slice_cor_data <- filter(slice_data, Treatment == "RP")   # READ ME: RP-only variant !

cor_results <- map_dfr(cytokine_cols, function(cy) {
  df <- slice_cor_data %>% filter(!is.na(.data[[cy]]), !is.na(motile_frac))
  if (nrow(df) < 3) return(tibble(cytokine = cy, n = nrow(df), rho = NA, p = NA))
  ct <- suppressWarnings(cor.test(df[[cy]], df$motile_frac, method = "spearman"))
  tibble(cytokine = cy, n = nrow(df), rho = unname(ct$estimate), p = ct$p.value)
}) %>% arrange(p)

print(cor_results, n = Inf)
write_csv(cor_results, paste0("slice_cytokine_motile_spearman_RPonly_", output_suffix, ".csv"))


# derive Mouse from TumourID_global ("M1_T1" -> "M1")
plot_df <- slice_cor_data %>%
  mutate(Mouse = sub("_.*$", "", TumourID_global)) %>%
  pivot_longer(all_of(cytokine_cols),
               names_to = "cytokine", values_to = "cyto_value") %>%
  filter(!is.na(cyto_value), !is.na(motile_frac))

# label table: one rho/p string per cytokine, from cor_results
lab_df <- cor_results %>%
  mutate(label = sprintf("rho = %.2f, p = %.3f (n = %d)", rho, p, n))

p_one <- list()                    
for (cy in cytokine_cols) {
  d   <- filter(plot_df, cytokine == cy)
  lab <- lab_df$label[lab_df$cytokine == cy]
  if (nrow(d) < 3) next
  
  p_one[[cy]] <- ggplot(d, aes(cyto_value, motile_pct, colour = Treatment, shape = Mouse)) +
    geom_smooth(aes(group = 1), method = "lm", se = FALSE,
                colour = "grey50", linewidth = 0.4, linetype = "dashed") +
    geom_point(size = 3, stroke = 0.8) +
    scale_colour_manual(values = c(Mock = "#2C6FBB", RP = "#C0392B")) +
    scale_shape_manual(values = c(16, 17, 15, 18, 8)) +
    annotate("text", x = Inf, y = Inf, label = lab,
             hjust = 1.05, vjust = 1.5, size = 3.5) +
    labs(x = paste(cy, "level"), y = "Motile fraction (%)",
         title = cy, colour = "Treatment", shape = "Mouse") +
    theme_bw(base_size = 11) +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold"))
  
  # ggsave(paste0("slice_", cy, "_vs_motile_", output_suffix, ".jpeg"),
  #        plot = p_one, width = 4.5, height = 4.5, dpi = 300)
}

# =============================================================================
#  ---- CORE-ONLY MOTILE PROPORTIONS per slice + slice cytokines ----
#  Pools only Core-region tracks within each slice (one value per slice).
#  Assumption: slice cytokine level applies equally across regions.
# =============================================================================
maxdis_motile_threshold <- 7
motile_classes <- c("1 motile", "3 motile")

core_props <- tracks %>%
  filter(Location == "Core") %>%                       # <-- EDIT 1: core only
  mutate(GMM_component_plot = case_when(
    GMM_component == 1 & MaxDis >  maxdis_motile_threshold ~ "1 motile",
    GMM_component == 1 & MaxDis <= maxdis_motile_threshold ~ "1",
    GMM_component == 3 & MaxDis >  maxdis_motile_threshold ~ "3 motile",
    GMM_component == 3 & MaxDis <= maxdis_motile_threshold ~ "3",
    TRUE ~ as.character(GMM_component)
  )) %>%
  group_by(SliceID_global, Day) %>%                         # <-- pool per slice, not per ROI
  summarise(
    TumourID_global = first(TumourID_global),
    Treatment       = first(Treatment),
    n_total_core  = n(),
    n_motile_core = sum(GMM_component_plot %in% motile_classes),
    core_motile_frac = n_motile_core / n_total_core,
    core_motile_pct  = round(100 * core_motile_frac, 1),
    .groups = "drop"
  )

# join slice cytokines (broadcast) by SliceID_global
cytokine_cols <- setdiff(names(slice_data), names(slice_props))   ### !!! NEEDS SLICE PROPS FROM ABOVE
core_data <- core_props %>%
  left_join(slice_data %>% select(SliceID_global, all_of(cytokine_cols)),
            by = "SliceID_global") %>%
  arrange(TumourID_global, SliceID_global)

stopifnot(nrow(core_data) == nrow(core_props))
cat("Core slices:", nrow(core_data),
    "| with cytokines:", sum(!is.na(core_data[[cytokine_cols[1]]])), "\n")
print(core_data, n = Inf)
write_csv(core_data, paste0("core_motile_proportions_cytokines_", output_suffix, ".csv"))

# =============================================================================
#  ---- CORE-ONLY CORRELATION: slice cytokines vs CORE motile fraction ----
#  Exploratory. Spearman rho + p (cor.test).
# =============================================================================
cytokine_cols <- setdiff(names(slice_data), names(slice_props))
cat("Detected cytokine columns:\n"); print(cytokine_cols)

#coreroi_cor_data <- core_data                                  # all slices
coreroi_cor_data <- filter(core_data, Treatment == "RP")     # RP-only variant

cor_results <- map_dfr(cytokine_cols, function(cy) {
  df <- coreroi_cor_data %>% filter(!is.na(.data[[cy]]), !is.na(core_motile_frac))
  if (nrow(df) < 3) return(tibble(cytokine = cy, n = nrow(df), rho = NA, p = NA))
  ct <- suppressWarnings(cor.test(df[[cy]], df$core_motile_frac, method = "spearman"))  # FIX 1
  tibble(cytokine = cy, n = nrow(df), rho = unname(ct$estimate), p = ct$p.value)
}) %>% arrange(p)

print(cor_results, n = Inf)
write_csv(cor_results, paste0("core_cytokine_motile_spearman_RPonly_", output_suffix, ".csv"))

# --- plot: use the SAME core data, not slice_cor_data --------------------
plot_df <- coreroi_cor_data %>%                                # FIX 2
  mutate(Mouse = sub("_.*$", "", TumourID_global)) %>%
  pivot_longer(all_of(cytokine_cols),
               names_to = "cytokine", values_to = "cyto_value") %>%
  filter(!is.na(cyto_value), !is.na(core_motile_frac))         # FIX 3a

lab_df <- cor_results %>%
  mutate(label = sprintf("rho = %.2f, p = %.3f (n = %d)", rho, p, n))

p_one_core <- list()     
for (cy in cytokine_cols) {
  d   <- filter(plot_df, cytokine == cy)
  lab <- lab_df$label[lab_df$cytokine == cy]
  if (nrow(d) < 3) next
  
  p_one_core[[cy]] <- ggplot(d, aes(cyto_value, core_motile_pct, colour = Treatment, shape = Mouse)) +  # FIX 3b
    geom_smooth(aes(group = 1), method = "lm", se = FALSE,
                colour = "grey50", linewidth = 0.4, linetype = "dashed") +
    geom_point(size = 3, stroke = 0.8) +
    scale_colour_manual(values = c(Mock = "#2C6FBB", RP = "#C0392B")) +
    scale_shape_manual(values = c(16, 17, 15, 18, 8)) +
    annotate("text", x = Inf, y = Inf, label = lab,
             hjust = 1.05, vjust = 1.5, size = 3.5) +
    labs(x = paste(cy, "level"), y = "Core motile fraction (%)",
         title = cy, colour = "Treatment", shape = "Mouse") +
    theme_bw(base_size = 11) +
    theme(legend.position = "right", plot.title = element_text(face = "bold"))
  
  # ggsave(paste0("core_", cy, "_vs_motile_", output_suffix, ".jpeg"),   # FIX 4: filename
  #        plot = p_one_core, width = 4.5, height = 4.5, dpi = 300)
}


# =============================================================================
#  TUMOUR VOLUME vs MOTILE PROPORTION — Option A (track-pooled per tumour)
#  One motile fraction per TumourID_global (all tracks pooled = track-weighted).
#  Exploratory: Spearman rho + p vs tumour volume. n = 5 tumours.
# =============================================================================
maxdis_motile_threshold <- 7
motile_classes <- c("1 motile", "3 motile")

tumour_props <- tracks %>%
  mutate(GMM_component_plot = case_when(
    GMM_component == 1 & MaxDis >  maxdis_motile_threshold ~ "1 motile",
    GMM_component == 1 & MaxDis <= maxdis_motile_threshold ~ "1",
    GMM_component == 3 & MaxDis >  maxdis_motile_threshold ~ "3 motile",
    GMM_component == 3 & MaxDis <= maxdis_motile_threshold ~ "3",
    TRUE ~ as.character(GMM_component)
  )) %>%
  group_by(TumourID_global) %>%                         # pool ALL tracks per tumour
  summarise(
    Mouse       = sub("_.*$", "", TumourID_global[1]),
    n_total     = n(),
    n_motile    = sum(GMM_component_plot %in% motile_classes),
    motile_frac = n_motile / n_total,
    motile_pct  = round(100 * motile_frac, 1),
    .groups = "drop"
  )

# --- attach tumour volume -------------------------------------------------
# EDIT: fill in volumes (mm^3) keyed to TumourID_global. Example:
volumes <- tibble::tribble(
  ~TumourID_global, ~volume_mm3,
  "M1_T1",  62,
  "M2_T1",  191,
  "M2_T2",  316,
  "M3_T1",  120,
  "M3_T2",  150,
  "M4_T1",  86,
  "M4_T2",  157,
  "M5_T1",  69
)
# volumes <- tibble::tribble(
#   ~TumourID_global, ~volume_mm3,
#   "M1_T1",  62,
#   "M2_T1",  316,
#   "M2_T2",  191,
#   "M3_T1",  150,
#   "M3_T2",  120,
#   "M4_T1",  157,
#   "M4_T2",  86,
#   "M5_T1",  69
# )

tumour_vol <- tumour_props %>%
  left_join(volumes, by = "TumourID_global")

stopifnot(nrow(tumour_vol) == nrow(tumour_props))       # no dup/drop from join
cat("Tumours:", nrow(tumour_vol),
    "| with volume:", sum(!is.na(tumour_vol$volume_mm3)), "\n")
print(tumour_vol, n = Inf)

# --- exploratory Spearman -------------------------------------------------
d_t  <- tumour_vol %>% filter(!is.na(volume_mm3), !is.na(motile_frac))
ct <- suppressWarnings(cor.test(d_t$volume_mm3, d_t$motile_frac, method = "spearman"))
cat(sprintf("\nSpearman: rho = %.2f, p = %.3f (n = %d) — exploratory\n",
            unname(ct$estimate), ct$p.value, nrow(d)))

write_csv(tumour_vol, paste0("tumour_volume_vs_motile_A_", output_suffix, ".csv"))

# --- build the rho/p label from the Spearman test already run -------------
lab <- sprintf("rho = %.2f, p = %.3f (n = %d)",
               unname(ct$estimate), ct$p.value, nrow(d_t))

p_vol <- ggplot(d_t, aes(volume_mm3, motile_pct, colour = Mouse, shape = Mouse)) +
  # geom_smooth(aes(group = 1), method = "lm", se = FALSE,
  #             colour = "grey50", linewidth = 0.4, linetype = "dashed") +
  geom_point(size = 3.5, stroke = 0.8) +
  scale_shape_manual(values = c(16, 17, 15, 18, 8)) +
  annotate("text", x = Inf, y = Inf, label = lab,
           hjust = 1.05, vjust = 1.5, size = 3.5) +
  labs(x = expression("Tumour volume (mm"^3*")"),
       y = "Motile fraction (%)",
       colour = "Mouse", shape = "Mouse") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")

print(p_vol)
ggsave(paste0("tumour_volume_vs_motile_A_", output_suffix, ".jpeg"),
       plot = p_vol, width = 5, height = 4.5, dpi = 300)