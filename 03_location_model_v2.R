# =============================================================================
# CELL TRACKING ANALYSIS - STEP 3: LOCATION MODEL
# =============================================================================
# Purpose: Test whether motility parameters differ by tumour location
#          (Core, CoreEdge, Edge) while controlling for Treatment and Day.
#
# Models fitted for each variable (raw + log):
#   ROI level:   variable ~ Location + Treatment + Day + (1|MouseID/TumourID_global/SliceID_global)
#   Track level: variable ~ Location + Treatment + Day + (1|MouseID/TumourID_global/SliceID_global/ROIID_global)
#
# Variables: Speed (raw + log), TrackLength (log only), Straightness (log only)
#
# Outputs (all suffixed with output_suffix):
#   location_model_summaries_<suffix>.txt     : full model summaries
#   location_emms_<suffix>.csv                : estimated marginal means per location
#   location_pairwise_<suffix>.csv            : pairwise comparisons (Tukey)
#   location_r2_<suffix>.csv                  : R² marginal/conditional
#   location_assumption_plots_<suffix>.pdf    : residual + QQ plots
#   location_emm_plots_<suffix>.pdf           : EMM dot plots + forest plots
#   location_effect_sizes_<suffix>.csv        : Cohen's d per effect
#   location_geometric_means_<suffix>.csv     : back-transformed group means
#   location_fold_changes_<suffix>.csv        : fold changes from log contrasts
#   location_model3_summaries_<suffix>.txt    : Model 3 (interaction) summaries
#   location_model3_lrt_<suffix>.csv          : Model 3 LRT results
#   location_model3_plots_<suffix>.pdf        : Model 3 EMM plots
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
gmm_keep <- c(1,3)

# Output suffix — auto-constructed
output_suffix <- if (length(gmm_keep) == 0) "all" else
                 paste0("GMM", paste(gmm_keep, collapse = "_"))

# ROI-level analysis gate
# Set TRUE only when gmm_keep is finalised AND roi_summary matches it.
roi_analysis_enabled <- TRUE

# --- Variable toggle ---
# Controls which variables enter the models. Mirror this in script 04.
# Screening (02) always runs on all variables regardless of this toggle.
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
maxdis_threshold <- 7  # µm

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

  # Alpha R² threshold — must match 02
  alpha_r2_threshold <- 0.5
  

# --- 0. PACKAGES --------------------------------------------------------------
# install.packages(c("lme4", "lmerTest", "emmeans", "effectsize",
#                    "performance", "tidyverse", "patchwork", "pbkrtest"))
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(effectsize)
library(performance)
library(patchwork)
library(pbkrtest)


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


# --- 1b. GMM and NA FILTER -----------------------------------------------------------

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
# Applied before model fitting; Alpha must be in (0, 2]
n_alpha_before <- sum(!is.na(tracks$alpha_30pct))
tracks <- tracks %>%
  mutate(alpha_30pct = ifelse(alpha_30pct <= 0 | alpha_30pct > 2.25,
                              NA_real_, alpha_30pct))
cat("Alpha bounds filter: ",
    n_alpha_before - sum(!is.na(tracks$alpha_30pct)),
    "values set to NA (outside 0-2 range)\n\n")



# --- 2. DEFINE VARIABLES ------------------------------------------------------

# Track level column names (from tracks_clean_v2.csv)
track_vars <- list(
  list(name = "Speed",        raw = "SpeedMean_ummin", log = "log_Speed",
       raw_label = "Speed (µm/min)",      log_label = "Log Speed (µm/min)",
       run_raw = TRUE),
  list(name = "TrackLength",  raw = NULL,              log = "log_TrackLength",
       raw_label = NULL,                 log_label = "Log Track Length (µm)",
       run_raw = FALSE),
  list(name = "Straightness", raw = "Track_Straightness", log = "log_Straightness",
       raw_label = "Straightness",        log_label = "Log Straightness",
       run_raw = TRUE),
  list(name = "MaxDis",       raw = "MaxDis",              log = "log_MaxDis",
       raw_label = "MaxDis (um)",          log_label = "Log MaxDis (um)",
       run_raw = TRUE),
  list(name = "TurningAngle", raw = "mean_turning_angle",  log = NULL,
       raw_label = "Mean Turning Angle (deg)", log_label = NULL,
       run_raw = TRUE),
  list(name = "Alpha",        raw = "alpha_30pct",         log = NULL,
       raw_label = "Alpha (MSD)",          log_label = NULL,
       run_raw = TRUE,
       filter_col = "msd_r2_30pct", filter_val = alpha_r2_threshold)
)

# ROI level column names (from roi_summary_v2.csv)
roi_vars <- list(
  list(name = "Speed",        raw = "Speed_mean",       log = "log_Speed_mean",
       raw_label = "Speed (µm/min)",      log_label = "Log Speed (µm/min)",
       run_raw = TRUE),
  list(name = "TrackLength",  raw = NULL,               log = "log_Length_mean",
       raw_label = NULL,                 log_label = "Log Track Length (µm)",
       run_raw = FALSE),
  list(name = "Straightness", raw = "Straightness_mean",  log = "log_Straight_mean",
       raw_label = "Straightness",        log_label = "Log Straightness",
       run_raw = TRUE),
  list(name = "MaxDis",       raw = "MaxDis_mean",         log = "log_MaxDis_mean",
       raw_label = "MaxDis (um)",          log_label = "Log MaxDis (um)",
       run_raw = TRUE),
  list(name = "TurningAngle", raw = "TurningAngle_mean",   log = NULL,
       raw_label = "Mean Turning Angle (deg)", log_label = NULL,
       run_raw = TRUE),
  list(name = "Alpha",        raw = "Alpha30_filt_mean",   log = NULL,
       raw_label = "Alpha (MSD, filtered)", log_label = NULL,
       run_raw = TRUE)
)

# Colours
location_colours <- c(
  "Core"     = "#E63946",
  "CoreEdge" = "#457B9D",
  "Edge"     = "#2A9D8F"
)


# --- 3. HELPER: FIT MODEL 2 (additive location model) ------------------------
# Fixed effects: Location + Treatment + Day 
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


fit_location_model <- function(data, col, random_formula) {
  formula <- as.formula(
    paste(col, "~ Location + Treatment + Day +", random_formula)
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

extract_emm <- function(model, var_name, level_name, transform) {

  # KR for ROI level (small n); Satterthwaite for Track level (large n)
  # Reference: Kenward & Roger (1997) Biometrics 53:983-997
  df_method <- if (level_name == "ROI") "kenward-roger" else "satterthwaite"

  emm <- emmeans(model, ~ Location, lmer.df = df_method)
  pw  <- pairs(emm, adjust = "tukey", infer = TRUE)

  emm_raw <- as.data.frame(emm)
  pw_raw  <- as.data.frame(pw)

  names(emm_raw)[names(emm_raw) == "emmean"] <- "emmean_val"
  names(emm_raw)[grepl("lower|LCL", names(emm_raw), ignore.case = TRUE)] <- "lower_CL"
  names(emm_raw)[grepl("upper|UCL", names(emm_raw), ignore.case = TRUE)] <- "upper_CL"
  names(pw_raw)[grepl("lower|LCL",  names(pw_raw),  ignore.case = TRUE)] <- "lower_CL"
  names(pw_raw)[grepl("upper|UCL",  names(pw_raw),  ignore.case = TRUE)] <- "upper_CL"

  emm_df <- emm_raw %>%
    mutate(Variable = var_name, Level = level_name, Transform = transform)

  pw_df <- pw_raw %>%
    mutate(
      Variable    = var_name,
      Level       = level_name,
      Transform   = transform,
      Significant = p.value < 0.05
    )

  list(emm = emm_df, pairwise = pw_df, emm_obj = emm, pw_obj = pw)
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


# --- 6. HELPER: EMM DOT PLOT --------------------------------------------------

make_emm_plot <- function(emm_obj, pw_obj, title, ylabel) {

  emm_raw <- as.data.frame(emm_obj)
  lower_col <- intersect(names(emm_raw), c("lower.CL", "asymp.LCL"))[1]
  upper_col <- intersect(names(emm_raw), c("upper.CL", "asymp.UCL"))[1]
  emm_df <- emm_raw %>%
    rename(mean = emmean, lower = !!lower_col, upper = !!upper_col)

  pw_raw <- as.data.frame(pw_obj)
  lower_col_pw <- intersect(names(pw_raw), c("lower.CL", "asymp.LCL"))[1]
  upper_col_pw <- intersect(names(pw_raw), c("upper.CL", "asymp.UCL"))[1]
  pw_df <- pw_raw %>%
    rename(lower = !!lower_col_pw, upper = !!upper_col_pw) %>%
    mutate(
      stars = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE            ~ "ns"
      )
    )

  # EMM dot plot with 95% CI
  p_emm <- ggplot(emm_df, aes(x = Location, y = mean, colour = Location)) +
    geom_point(size = 4) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15, linewidth = 0.8) +
    geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE),
                  width = 0.0, linewidth = 3, alpha = 0.5) +
    scale_colour_manual(values = location_colours) +
    labs(title = title, x = "Location", y = ylabel, 
         caption = "Point = EMM | Thick bar = ±SE | Thin bar = 95% CI") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 11))

  # Forest plot of pairwise contrasts
  p_forest <- ggplot(pw_df, aes(y = contrast, x = estimate)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(size = 3) +
    geom_errorbar(aes(xmin = lower, xmax = upper),
                  width = 0.2, linewidth = 0.8,
                  orientation = "y") +
    geom_text(aes(label = stars, x = upper),
              hjust = -0.3, size = 4) +
    labs(title = "Pairwise contrasts (Tukey)",
         x = "Estimate (difference)", y = NULL) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11))

  p_emm | p_forest
}


# --- 7. RUN MODEL 2 -----------------------------------------------------------

cat("\n=== FITTING LOCATION MODELS (MODEL 2 — ADDITIVE) ===\n")

all_emms         <- list()
all_pairwise     <- list()
all_r2           <- list()
all_models       <- list()
assumption_plots <- list()
emm_plots        <- list()

summary_file <- paste0("location_model_summaries_", output_suffix, ".txt")

sink(summary_file)
cat("LOCATION MODEL SUMMARIES — MODEL 2 (ADDITIVE)\n")
cat("===============================================\n")
cat("Fixed effects: Location + Treatment + Day\n")
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

      model <- fit_location_model(data_v, col, rf)
      if (is.null(model)) next

      all_models[[key]] <<- model

      singular <- isSingular(model)
      if (singular) cat("  NOTE: Singular fit for", title, "\n")

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

      all_r2[[key]] <<- tibble(
        Model          = key,
        Variable       = vp$name,
        Level          = level_name,
        Transform      = transform,
        Singular_fit   = singular,
        R2_marginal    = r2m,
        R2_conditional = r2c,
        R2_random      = round(r2c - r2m, 3)
      )

      sink(summary_file, append = TRUE)
      cat("\n", paste(rep("=", 60), collapse = ""), "\n")
      cat(title, "\n")
      cat(paste(rep("=", 60), collapse = ""), "\n")
      cat("Singular fit:", singular, "\n")
      cat("R2 marginal (fixed effects):  ", r2m, "\n")
      cat("R2 conditional (fixed+random):", r2c, "\n")
      cat("R2 random effects only:       ", round(r2c - r2m, 3), "\n\n")
      print(summary(model))
      cat("\nAnova (Type III):\n")
      print(anova(model))
      sink()

      emm_result <- tryCatch(
        extract_emm(model, vp$name, level_name, transform),
        error = function(e) { cat("  EMM failed:", e$message, "\n"); NULL }
      )

      if (!is.null(emm_result)) {
        all_emms[[key]]     <<- emm_result$emm
        all_pairwise[[key]] <<- emm_result$pairwise
        emm_plots[[key]]    <<- make_emm_plot(
          emm_obj = emm_result$emm_obj,
          pw_obj  = emm_result$pw_obj,
          title   = title,
          ylabel  = ylabel
        )
      }

      assumption_plots[[key]] <<- make_assumption_plots(model, title)
    }
  }
}

# ROI level — gated: only run when roi_analysis_enabled = TRUE
if (!roi_analysis_enabled) {
  cat("\nROI-level location models skipped (roi_analysis_enabled = FALSE).\n")
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


# --- 8. SAVE NUMERICAL RESULTS FOR MODEL 2 ------------------------------------

emm_table      <- bind_rows(all_emms)
pairwise_table <- bind_rows(all_pairwise) %>%
  mutate(
    p_display = case_when(
      is.na(p.value)      ~ "NA",
      p.value < 0.0001    ~ "< 0.0001",
      TRUE                ~ as.character(round(p.value, 4))
    ),
    stars = case_when(
      is.na(p.value)  ~ "NA",
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )

r2_table <- bind_rows(all_r2)

write_csv(emm_table,      paste0("location_emms_",      output_suffix, ".csv"))
write_csv(pairwise_table, paste0("location_pairwise_",  output_suffix, ".csv"))
write_csv(r2_table,       paste0("location_r2_",        output_suffix, ".csv"))

cat("\nSaved: location_model_summaries_", output_suffix, ".txt\n", sep = "")
cat("Saved: location_emms_",             output_suffix, ".csv\n", sep = "")
cat("Saved: location_pairwise_",         output_suffix, ".csv\n", sep = "")
cat("Saved: location_r2_",              output_suffix, ".csv\n", sep = "")

cat("\n--- R2 Summary ---\n")
print(r2_table, n = Inf)


# --- 9. SAVE ASSUMPTION PLOTS -------------------------------------------------

pdf(paste0("location_assumption_plots_", output_suffix, ".pdf"), width = 10, height = 4)
for (nm in names(assumption_plots)) {
  print(assumption_plots[[nm]])
}
dev.off()

cat("Saved: location_assumption_plots_", output_suffix, ".pdf\n", sep = "")


# --- 10. SAVE EMM PLOTS -------------------------------------------------------

pdf(paste0("location_emm_plots_", output_suffix, ".pdf"), width = 10, height = 5)

for (vp in roi_vars) {
  for (type in c("Raw", "Log")) {
    key_roi   <- paste("ROI",   vp$name, type, sep = "_")
    key_track <- paste("Track", vp$name, type, sep = "_")

    if (!is.null(emm_plots[[key_roi]]))   print(emm_plots[[key_roi]])
    if (!is.null(emm_plots[[key_track]])) print(emm_plots[[key_track]])
  }
}

dev.off()

cat("Saved: location_emm_plots_", output_suffix, ".pdf\n", sep = "")
cat("\nAll outputs saved. Review results before proceeding.\n")


# =============================================================================
# ---- SECTION A: EFFECT SIZES, GEOMETRIC MEANS, FOLD CHANGES ----
# =============================================================================
# Computed from stored Model 2 objects.
# Three outputs per model:
#   (1) Cohen's d using residual SD (sigma_resid) — primary estimate
#       and conservative d using total SD (sigma_total) — marginal effects only
#   (2) Geometric means — back-transformed EMMs (log variables only)
#   (3) Fold changes — exp(contrast estimate) between location pairs
#
# Reference: Cohen (1988) Statistical Power Analysis for the Behavioral Sciences
#            Nakagawa & Schielzeth (2013) Methods Ecol Evol 4:133-142 (R2)
# =============================================================================

cat("\n=== EFFECT SIZES (MODEL 2) ===\n")

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

  # --- Marginal EMMs: Location, Treatment, Day ---
  for (effect in c("Location", "Treatment", "Day")) {

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
    if (is.null(effect_df) || is.na(effect_df)) effect_df <- df.residual(model)

    eff_resid <- tryCatch(
      eff_size(emm_obj, sigma = sigma_resid, edf = effect_df),
      error = function(e) { cat("  eff_size failed:", e$message, "\n"); NULL }
    )
    eff_total <- tryCatch(
      eff_size(emm_obj, sigma = sigma_total, edf = effect_df),
      error = function(e) NULL
    )
    if (is.null(eff_resid)) next

    # Back-transformed group means
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

    # Fold changes from pairwise contrasts
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

  # --- Conditional: Location within each Treatment ---
  emm_loc_by_trt <- tryCatch(
    emmeans(model, ~ Location | Treatment,
            lmer.df = if (level_nm == "ROI") "kenward-roger" else "satterthwaite"),
    error = function(e) { cat("  EMM failed Location|Treatment:", e$message, "\n"); NULL }
  )

  if (!is.null(emm_loc_by_trt)) {

    eff_loc_by_trt <- tryCatch(
      eff_size(emm_loc_by_trt, sigma = sigma_resid,
               edf = aov_tbl["Location", "DenDF"]),
      error = function(e) NULL
    )
    if (!is.null(eff_loc_by_trt)) {
      effect_size_results[[paste(key, "Location_by_Treatment", sep = "_")]] <-
        as.data.frame(eff_loc_by_trt) %>%
        mutate(Model = key, Variable = var_nm, Level = level_nm,
               Transform = transform, Fixed_effect = "Location within Treatment",
               d_residual = round(effect.size, 3),
               sigma_residual = round(sigma_resid, 4),
               d_interpretation = case_when(
                 abs(effect.size) >= 0.8 ~ "Large",
                 abs(effect.size) >= 0.5 ~ "Medium",
                 abs(effect.size) >= 0.2 ~ "Small",
                 TRUE                    ~ "Negligible"
               ))
    }

    emm_loc_trt_raw <- as.data.frame(emm_loc_by_trt)
    names(emm_loc_trt_raw)[names(emm_loc_trt_raw) == "emmean"] <- "emmean_val"
    names(emm_loc_trt_raw)[grepl("lower|LCL", names(emm_loc_trt_raw), ignore.case = TRUE)] <- "lower_CL"
    names(emm_loc_trt_raw)[grepl("upper|UCL", names(emm_loc_trt_raw), ignore.case = TRUE)] <- "upper_CL"

    if (transform == "Log") {
      emm_loc_trt_raw <- emm_loc_trt_raw %>%
        mutate(
          geometric_mean     = round(exp(emmean_val), 3),
          geometric_mean_LCL = round(exp(lower_CL),   3),
          geometric_mean_UCL = round(exp(upper_CL),   3)
        )
      geomean_results[[paste(key, "Location_by_Treatment", "geomean", sep = "_")]] <-
        emm_loc_trt_raw %>%
        mutate(Model = key, Variable = var_nm, Level = level_nm,
               Transform = transform, Fixed_effect = "Location within Treatment")
    }

    pw_loc_by_trt <- tryCatch(
      as.data.frame(pairs(emm_loc_by_trt, infer = TRUE)),
      error = function(e) NULL
    )
    if (!is.null(pw_loc_by_trt) && transform == "Log") {
      names(pw_loc_by_trt)[grepl("lower|LCL", names(pw_loc_by_trt), ignore.case = TRUE)] <- "lower_CL"
      names(pw_loc_by_trt)[grepl("upper|UCL", names(pw_loc_by_trt), ignore.case = TRUE)] <- "upper_CL"
      pw_loc_by_trt <- pw_loc_by_trt %>%
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
      foldchange_results[[paste(key, "Location_by_Treatment", "foldchange", sep = "_")]] <-
        pw_loc_by_trt %>%
        mutate(Model = key, Variable = var_nm, Level = level_nm,
               Transform = transform, Fixed_effect = "Location within Treatment")
    }
  }

  # --- Conditional: Location within each Day ---
  emm_loc_by_day <- tryCatch(
    emmeans(model, ~ Location | Day,
            lmer.df = if (level_nm == "ROI") "kenward-roger" else "satterthwaite"),
    error = function(e) { cat("  EMM failed Location|Day:", e$message, "\n"); NULL }
  )

  if (!is.null(emm_loc_by_day)) {

    eff_loc_by_day <- tryCatch(
      eff_size(emm_loc_by_day, sigma = sigma_resid,
               edf = aov_tbl["Location", "DenDF"]),
      error = function(e) NULL
    )
    if (!is.null(eff_loc_by_day)) {
      effect_size_results[[paste(key, "Location_by_Day", sep = "_")]] <-
        as.data.frame(eff_loc_by_day) %>%
        mutate(Model = key, Variable = var_nm, Level = level_nm,
               Transform = transform, Fixed_effect = "Location within Day",
               d_residual = round(effect.size, 3),
               sigma_residual = round(sigma_resid, 4),
               d_interpretation = case_when(
                 abs(effect.size) >= 0.8 ~ "Large",
                 abs(effect.size) >= 0.5 ~ "Medium",
                 abs(effect.size) >= 0.2 ~ "Small",
                 TRUE                    ~ "Negligible"
               ))
    }

    emm_loc_day_raw <- as.data.frame(emm_loc_by_day)
    names(emm_loc_day_raw)[names(emm_loc_day_raw) == "emmean"] <- "emmean_val"
    names(emm_loc_day_raw)[grepl("lower|LCL", names(emm_loc_day_raw), ignore.case = TRUE)] <- "lower_CL"
    names(emm_loc_day_raw)[grepl("upper|UCL", names(emm_loc_day_raw), ignore.case = TRUE)] <- "upper_CL"

    if (transform == "Log") {
      emm_loc_day_raw <- emm_loc_day_raw %>%
        mutate(
          geometric_mean     = round(exp(emmean_val), 3),
          geometric_mean_LCL = round(exp(lower_CL),   3),
          geometric_mean_UCL = round(exp(upper_CL),   3)
        )
      geomean_results[[paste(key, "Location_by_Day", "geomean", sep = "_")]] <-
        emm_loc_day_raw %>%
        mutate(Model = key, Variable = var_nm, Level = level_nm,
               Transform = transform, Fixed_effect = "Location within Day")
    }

    pw_loc_by_day <- tryCatch(
      as.data.frame(pairs(emm_loc_by_day, infer = TRUE)),
      error = function(e) NULL
    )
    if (!is.null(pw_loc_by_day) && transform == "Log") {
      names(pw_loc_by_day)[grepl("lower|LCL", names(pw_loc_by_day), ignore.case = TRUE)] <- "lower_CL"
      names(pw_loc_by_day)[grepl("upper|UCL", names(pw_loc_by_day), ignore.case = TRUE)] <- "upper_CL"
      pw_loc_by_day <- pw_loc_by_day %>%
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
      foldchange_results[[paste(key, "Location_by_Day", "foldchange", sep = "_")]] <-
        pw_loc_by_day %>%
        mutate(Model = key, Variable = var_nm, Level = level_nm,
               Transform = transform, Fixed_effect = "Location within Day")
    }
  }
}

write_csv(bind_rows(effect_size_results), paste0("location_effect_sizes_",    output_suffix, ".csv"))
write_csv(bind_rows(geomean_results),     paste0("location_geometric_means_", output_suffix, ".csv"))
write_csv(bind_rows(foldchange_results),  paste0("location_fold_changes_",    output_suffix, ".csv"))

cat("Saved: location_effect_sizes_",    output_suffix, ".csv\n", sep = "")
cat("Saved: location_geometric_means_", output_suffix, ".csv\n", sep = "")
cat("Saved: location_fold_changes_",    output_suffix, ".csv\n", sep = "")


# =============================================================================
# ---- SECTION B: MODEL 3 — IS THE LOCATION EFFECT INVARIANT ACROSS TREATMENT? ----
# =============================================================================
# Formula: variable ~ Location + Treatment + Day + Treatment:Location + (1|random)
#
# The Treatment:Location interaction tests whether the spatial motility
# gradient (Core vs CoreEdge vs Edge) is consistent regardless of treatment,
# or whether RP modifies the spatial pattern.
#
# A significant interaction means: the location effect differs between
# Mock and RP — i.e. location is not an invariant property of the tumour
# microenvironment but is modified by the experimental intervention.
#
# Comparison: LRT and BIC between Model 2 (additive) and Model 3 (interaction)
# Uses REML=FALSE for LRT as required for nested fixed effect comparisons.
#
# Interpretation guide:
#   LRT p < 0.05 AND BIC_diff > 10  : Location effect modified by treatment
#   LRT p < 0.05 AND BIC_diff 0-10  : Weak evidence for modification
#   LRT p > 0.05 OR  BIC_diff < 0   : Location effect is invariant
# =============================================================================

cat("\n=== MODEL 3 — TREATMENT:LOCATION INTERACTION ===\n")
cat("Question: Is the location effect on motility invariant across treatment?\n\n")

fit_location_model3 <- function(data, col, random_formula) {
  formula <- as.formula(
    paste(col, "~ Location + Treatment + Day + Treatment:Location +",
          random_formula)
  )
  tryCatch(
    lmer(formula, data = data, REML = TRUE),
    error   = function(e) { cat("  Model 3 error:", e$message, "\n"); NULL },
    warning = function(w) {
      m <- lmer(formula, data = data, REML = TRUE)
      cat("  Model 3 warning:", w$message, "\n")
      m
    }
  )
}

fit_ml <- function(formula_str, data) {
  tryCatch(
    lmer(as.formula(formula_str), data = data, REML = FALSE),
    error   = function(e) NULL,
    warning = function(w) suppressWarnings(
      lmer(as.formula(formula_str), data = data, REML = FALSE)
    )
  )
}

model3_lrt_results <- list()
model3_plots       <- list()

model3_summary_file <- paste0("location_model3_summaries_", output_suffix, ".txt")

sink(model3_summary_file)
cat("MODEL 3 SUMMARIES — TREATMENT:LOCATION INTERACTION\n")
cat("====================================================\n")
cat("Formula: variable ~ Location + Treatment + Day + Treatment:Location + (1|random)\n")
cat("Output suffix:", output_suffix, "\n\n")
sink()

run_model3_for_level <- function(var_list, data, random_formula, level_name) {

  for (vp in var_list) {
    for (type in c("raw", "log")) {

      if (type == "raw" && !isTRUE(vp$run_raw)) next
      if (type == "log" && is.null(vp$log)) next

      col       <- if (type == "raw") vp$raw else vp$log
      ylabel    <- if (type == "raw") vp$raw_label else vp$log_label
      transform <- if (type == "raw") "Raw" else "Log"
      key       <- paste(level_name, vp$name, transform, sep = "_")
      title     <- paste(level_name, "—", vp$name, paste0("(", transform, ")"),
                         "— Model 3")

      cat("  Fitting Model 3:", key, "\n")

      m3 <- fit_location_model3(data, col, random_formula)
      if (is.null(m3)) next

      singular3 <- isSingular(m3)
      if (singular3) cat("  NOTE: Model 3 singular fit for", key, "\n")

      sink(model3_summary_file, append = TRUE)
      cat("\n", paste(rep("=", 60), collapse = ""), "\n")
      cat(title, "\n")
      cat(paste(rep("=", 60), collapse = ""), "\n")
      cat("Singular fit:", singular3, "\n\n")
      print(summary(m3))
      cat("\nAnova (Type III):\n")
      print(anova(m3))
      sink()

      f2 <- paste(col, "~ Location + Treatment + Day +", random_formula)
      f3 <- paste(col, "~ Location + Treatment + Day + Treatment:Location +",
                  random_formula)

      m2_ml <- fit_ml(f2, data)
      m3_ml <- fit_ml(f3, data)

      lrt_result <- if (!is.null(m2_ml) && !is.null(m3_ml)) {
        tryCatch(anova(m2_ml, m3_ml), error = function(e) NULL)
      } else NULL

      if (!is.null(lrt_result)) {
        bic_diff <- round(BIC(m2_ml) - BIC(m3_ml), 2)
        model3_lrt_results[[key]] <<- tibble(
          Model         = key,
          Variable      = vp$name,
          Level         = level_name,
          Transform     = transform,
          LRT_Chi2      = round(lrt_result[2, "Chisq"],      3),
          LRT_df        = lrt_result[2, "Df"],
          LRT_p         = round(lrt_result[2, "Pr(>Chisq)"], 4),
          LRT_p_display = ifelse(lrt_result[2, "Pr(>Chisq)"] < 0.0001,
                                 "< 0.0001",
                                 as.character(round(lrt_result[2, "Pr(>Chisq)"], 4))),
          BIC_model2    = round(BIC(m2_ml), 2),
          BIC_model3    = round(BIC(m3_ml), 2),
          BIC_diff      = bic_diff,
          BIC_decision  = case_when(
            bic_diff > 10 ~ "Model 3 preferred — location modified by treatment",
            bic_diff > 0  ~ "Weak evidence for Model 3",
            TRUE          ~ "Model 2 preferred — location effect invariant"
          ),
          Singular_model2 = isSingular(m2_ml),
          Singular_model3 = singular3
        )
      }

      # EMM plot for Model 3 — Location within each Treatment
      emm3 <- tryCatch(
        emmeans(m3, ~ Location | Treatment,
                lmer.df = if (level_name == "ROI") "kenward-roger" else "satterthwaite"),
        error = function(e) NULL
      )

      if (!is.null(emm3)) {
        emm3_df <- as.data.frame(emm3)
        lower_col <- intersect(names(emm3_df), c("lower.CL", "asymp.LCL"))[1]
        upper_col <- intersect(names(emm3_df), c("upper.CL", "asymp.UCL"))[1]
        emm3_df <- emm3_df %>%
          rename(mean = emmean, lower = !!lower_col, upper = !!upper_col)

        p_model3 <- ggplot(emm3_df,
                           aes(x = Location, y = mean, colour = Location)) +
          geom_point(size = 4) +
          geom_errorbar(aes(ymin = lower, ymax = upper),
                        width = 0.15, linewidth = 0.8) +
          facet_wrap(~ Treatment) +
          scale_colour_manual(values = location_colours) +
          labs(title = paste("Model 3 —", title),
               subtitle = "Location effect within each Treatment arm",
               x = "Location", y = ylabel) +
          theme_bw(base_size = 11) +
          theme(legend.position  = "none",
                strip.text       = element_text(face = "bold"),
                plot.title       = element_text(face = "bold", size = 11),
                plot.subtitle    = element_text(size = 10))

        model3_plots[[key]] <<- p_model3
      }
    }
  }
}

# ROI level — gated
if (!roi_analysis_enabled) {
  cat("ROI-level Model 3 skipped (roi_analysis_enabled = FALSE).\n\n")
} else {
  cat("\n--- ROI level ---\n")
  run_model3_for_level(
    var_list       = roi_vars[map_lgl(roi_vars, ~ isTRUE(var_include[.x$name]))],
    data           = roi,
    random_formula = "(1|MouseID/TumourID_global/SliceID_global)",
    level_name     = "ROI"
  )
}

# Track level — always runs
cat("\n--- Track level ---\n")
run_model3_for_level(
  var_list       = track_vars[map_lgl(track_vars, ~ isTRUE(var_include[.x$name]))],
  data           = tracks,
  random_formula = "(1|MouseID/TumourID_global/SliceID_global/ROIID_global)",
  level_name     = "Track"
)

# Save Model 3 results
lrt_table <- bind_rows(model3_lrt_results)
write_csv(lrt_table, paste0("location_model3_lrt_", output_suffix, ".csv"))

cat("\nSaved: location_model3_summaries_", output_suffix, ".txt\n", sep = "")
cat("Saved: location_model3_lrt_",        output_suffix, ".csv\n", sep = "")

cat("\n--- Model 3 LRT Summary ---\n")
print(lrt_table %>%
        select(Model, LRT_p_display, BIC_diff, BIC_decision,
               Singular_model2, Singular_model3),
      n = Inf)

pdf(paste0("location_model3_plots_", output_suffix, ".pdf"), width = 10, height = 5)
for (nm in names(model3_plots)) print(model3_plots[[nm]])
dev.off()
cat("Saved: location_model3_plots_", output_suffix, ".pdf\n", sep = "")

cat("\nAll outputs saved.\n")
cat("Next: review Model 3 LRT table — if Treatment:Location is significant,\n")
cat("the location effect is not invariant and should be reported as treatment-dependent.\n")
cat("Then review Script 04 as the interaction model.\n")

# =============================================================================
# SECTION C: MODEL 3 EFFECT SIZES (INTERACTION MODEL — Treatment:Location)
# =============================================================================
# APPEND THIS TO THE END of 03_location_model_v2.R (after the final
# cat("Then review Script 04 as the interaction model.\n") line, ~line 1007).
#
# PURPOSE
#   Section B (existing) fits Model 3 (Treatment:Location interaction) and
#   writes summaries + LRT + plots, but does NOT extract Cohen's d.
#   In the ADDITIVE Model 2, "Treatment within Location" contrasts are
#   statistically identical across locations (one Treatment coefficient,
#   no interaction) — so region-specific treatment effects can ONLY be read
#   from the interaction model. This section fills that gap.
#
# WHAT IT COMPUTES (Track-level, log Speed by default — see var filter below)
#   (1) Treatment within Location : Mock - RP | Core, | CoreEdge, | Edge
#       -> the region-specific treatment effect (your dilution story, Table B)
#   (2) Location within Treatment : Core-CoreEdge-Edge | Mock, | RP
#       -> whether the spatial gradient itself differs by arm
#   Both with d_residual (sigma_resid) and d_total (sigma_total), matching
#   the exact method used in Section A for Model 2.
#
# SELF-CONTAINED (Option 1): re-fits Model 3 internally via the existing
#   fit_location_model3() helper. Does not depend on objects from Section B
#   beyond helpers/params already defined earlier in the script. Nothing
#   above this section is modified.
#
# METHOD REFERENCES
#   Cohen J (1988) Statistical Power Analysis for the Behavioral Sciences,
#     2nd ed. — Cohen's d and 0.2/0.5/0.8 benchmarks.
#   eff_size() with sigma = residual vs total SD: Lenth R, emmeans package
#     vignette "Comparisons and contrasts"; total-SD scaling reflects
#     Nakagawa & Schielzeth (2013) Methods Ecol Evol 4:133-142 variance
#     partitioning (sigma_total = sqrt(sum of all VarCorr variance components)).
#   Kenward & Roger (1997) Biometrics 53:983-997 (ROI-level df);
#     Satterthwaite (1946) via lmerTest (Track-level df).
# =============================================================================

cat("\n=== SECTION C: MODEL 3 EFFECT SIZES (INTERACTION) ===\n")

# --- C0. Which variables to extract -----------------------------------------
# Speed only for now. To extend later, add variable names to this vector
# (must match vp$name in track_vars / roi_vars, e.g. "TrackLength",
# "Straightness"). Everything downstream loops over this filter.
model3_es_vars <- c("Speed")

effect_size_results_m3 <- list()
geomean_results_m3     <- list()
foldchange_results_m3  <- list()

# --- C1. Core worker: fit Model 3 for one variable/level and extract d -------
# Mirrors Section A structure exactly, but on the interaction model and with
# the two conditional slices that require the interaction to be meaningful.
extract_model3_effect_sizes <- function(var_list, data, random_formula,
                                        level_name) {
  
  df_method <- if (level_name == "ROI") "kenward-roger" else "satterthwaite"
  
  for (vp in var_list) {
    
    if (!vp$name %in% model3_es_vars) next   # variable gate (Speed only default)
    
    for (type in c("raw", "log")) {
      
      if (type == "raw" && !isTRUE(vp$run_raw)) next
      if (type == "log" && is.null(vp$log))     next
      
      col       <- if (type == "raw") vp$raw else vp$log
      transform <- if (type == "raw") "Raw" else "Log"
      key       <- paste(level_name, vp$name, transform, sep = "_")
      
      cat("  Model 3 effect sizes for:", key, "\n")
      
      # Re-fit Model 3 (self-contained; uses existing helper)
      m3 <- fit_location_model3(data, col, random_formula)
      if (is.null(m3)) { cat("    fit failed, skipping\n"); next }
      
      if (isSingular(m3))
        cat("    NOTE: singular fit — d_total valid but interpret variance",
            "partition with caution\n")
      
      # SDs for standardisation — identical definitions to Section A
      sigma_resid <- sigma(m3)
      sigma_total <- sqrt(sum(as.data.frame(VarCorr(m3))$vcov))
      
      aov_tbl <- tryCatch(anova(m3), error = function(e) NULL)
      
      # df to feed eff_size: use the interaction's DenDF where available,
      # else fall back to residual df. (eff_size edf affects the CI on d,
      # not d itself.)
      edf_use <- tryCatch(aov_tbl["Location:Treatment", "DenDF"],
                          error = function(e) df.residual(m3))
      if (is.null(edf_use) || is.na(edf_use)) edf_use <- df.residual(m3)
      
      # ---- (1) Treatment within Location : Mock - RP | each location -------
      emm_trt_by_loc <- tryCatch(
        emmeans(m3, ~ Treatment | Location, lmer.df = df_method),
        error = function(e) { cat("    EMM Treatment|Location failed:",
                                  e$message, "\n"); NULL }
      )
      
      if (!is.null(emm_trt_by_loc)) {
        eff_tbl_res <- tryCatch(
          eff_size(emm_trt_by_loc, sigma = sigma_resid, edf = edf_use),
          error = function(e) NULL
        )
        eff_tbl_tot <- tryCatch(
          eff_size(emm_trt_by_loc, sigma = sigma_total, edf = edf_use),
          error = function(e) NULL
        )
        if (!is.null(eff_tbl_res)) {
          es <- as.data.frame(eff_tbl_res) %>%
            mutate(
              Model          = key,
              Variable       = vp$name,
              Level          = level_name,
              Transform      = transform,
              Fixed_effect   = "Treatment within Location",
              sigma_residual = round(sigma_resid, 4),
              sigma_total    = round(sigma_total, 4),
              d_residual     = round(effect.size, 3),
              d_total        = if (!is.null(eff_tbl_tot))
                round(as.data.frame(eff_tbl_tot)$effect.size, 3) else NA,
              d_interpretation = case_when(
                abs(effect.size) >= 0.8 ~ "Large",
                abs(effect.size) >= 0.5 ~ "Medium",
                abs(effect.size) >= 0.2 ~ "Small",
                TRUE                    ~ "Negligible"
              )
            )
          effect_size_results_m3[[paste(key, "Treatment_by_Location",
                                        sep = "_")]] <<- es
        }
        
        # Fold change + geometric means (log only), same as Section A
        if (transform == "Log") {
          emm_raw <- as.data.frame(emm_trt_by_loc)
          names(emm_raw)[names(emm_raw) == "emmean"] <- "emmean_val"
          names(emm_raw)[grepl("lower|LCL", names(emm_raw), ignore.case = TRUE)] <- "lower_CL"
          names(emm_raw)[grepl("upper|UCL", names(emm_raw), ignore.case = TRUE)] <- "upper_CL"
          geomean_results_m3[[paste(key, "Treatment_by_Location", "geomean",
                                    sep = "_")]] <<- emm_raw %>%
            mutate(geometric_mean = round(exp(emmean_val), 3),
                   geometric_mean_LCL = round(exp(lower_CL), 3),
                   geometric_mean_UCL = round(exp(upper_CL), 3),
                   Model = key, Variable = vp$name, Level = level_name,
                   Transform = transform,
                   Fixed_effect = "Treatment within Location")
          
          pw <- tryCatch(as.data.frame(pairs(emm_trt_by_loc, infer = TRUE)),
                         error = function(e) NULL)
          if (!is.null(pw)) {
            names(pw)[grepl("lower|LCL", names(pw), ignore.case = TRUE)] <- "lower_CL"
            names(pw)[grepl("upper|UCL", names(pw), ignore.case = TRUE)] <- "upper_CL"
            foldchange_results_m3[[paste(key, "Treatment_by_Location",
                                         "foldchange", sep = "_")]] <<- pw %>%
              mutate(fold_change = round(exp(estimate), 3),
                     fold_change_LCL = round(exp(lower_CL), 3),
                     fold_change_UCL = round(exp(upper_CL), 3),
                     Model = key, Variable = vp$name, Level = level_name,
                     Transform = transform,
                     Fixed_effect = "Treatment within Location")
          }
        }
      }
      
      # ---- (2) Location within Treatment : gradient | each arm -------------
      emm_loc_by_trt <- tryCatch(
        emmeans(m3, ~ Location | Treatment, lmer.df = df_method),
        error = function(e) { cat("    EMM Location|Treatment failed:",
                                  e$message, "\n"); NULL }
      )
      
      if (!is.null(emm_loc_by_trt)) {
        eff_res <- tryCatch(
          eff_size(emm_loc_by_trt, sigma = sigma_resid, edf = edf_use),
          error = function(e) NULL
        )
        eff_tot <- tryCatch(
          eff_size(emm_loc_by_trt, sigma = sigma_total, edf = edf_use),
          error = function(e) NULL
        )
        if (!is.null(eff_res)) {
          es2 <- as.data.frame(eff_res) %>%
            mutate(
              Model          = key,
              Variable       = vp$name,
              Level          = level_name,
              Transform      = transform,
              Fixed_effect   = "Location within Treatment",
              sigma_residual = round(sigma_resid, 4),
              sigma_total    = round(sigma_total, 4),
              d_residual     = round(effect.size, 3),
              d_total        = if (!is.null(eff_tot))
                round(as.data.frame(eff_tot)$effect.size, 3) else NA,
              d_interpretation = case_when(
                abs(effect.size) >= 0.8 ~ "Large",
                abs(effect.size) >= 0.5 ~ "Medium",
                abs(effect.size) >= 0.2 ~ "Small",
                TRUE                    ~ "Negligible"
              )
            )
          effect_size_results_m3[[paste(key, "Location_by_Treatment",
                                        sep = "_")]] <<- es2
        }
      }
    }
  }
}

# --- C2. Run: ROI level (gated) then Track level ----------------------------
if (roi_analysis_enabled) {
  extract_model3_effect_sizes(
    var_list       = roi_vars[map_lgl(roi_vars, ~ isTRUE(var_include[.x$name]))],
    data           = roi,
    random_formula = "(1|MouseID/TumourID_global/SliceID_global)",
    level_name     = "ROI"
  )
} else {
  cat("  ROI-level Model 3 effect sizes skipped (roi_analysis_enabled = FALSE).\n")
}

extract_model3_effect_sizes(
  var_list       = track_vars[map_lgl(track_vars, ~ isTRUE(var_include[.x$name]))],
  data           = tracks,
  random_formula = "(1|MouseID/TumourID_global/SliceID_global/ROIID_global)",
  level_name     = "Track"
)

# --- C3. Write outputs (schema-identical to Section A files) ----------------
if (length(effect_size_results_m3) > 0) {
  write_csv(bind_rows(effect_size_results_m3),
            paste0("location_model3_effect_sizes_", output_suffix, ".csv"))
  cat("Saved: location_model3_effect_sizes_", output_suffix, ".csv\n", sep = "")
}
if (length(geomean_results_m3) > 0) {
  write_csv(bind_rows(geomean_results_m3),
            paste0("location_model3_geometric_means_", output_suffix, ".csv"))
  cat("Saved: location_model3_geometric_means_", output_suffix, ".csv\n", sep = "")
}
if (length(foldchange_results_m3) > 0) {
  write_csv(bind_rows(foldchange_results_m3),
            paste0("location_model3_fold_changes_", output_suffix, ".csv"))
  cat("Saved: location_model3_fold_changes_", output_suffix, ".csv\n", sep = "")
}

cat("\nSection C complete. Re-run once per subpopulation (change gmm_keep)\n")
cat("to generate location_model3_effect_sizes_<suffix>.csv for each tier.\n")
