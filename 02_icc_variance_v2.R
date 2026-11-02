# =============================================================================
# CELL TRACKING ANALYSIS - STEP 2: ICC & VARIANCE DECOMPOSITION
# =============================================================================
# Purpose: Quantify non-independence at each hierarchical level at both
#          track level and ROI level to guide random effects structure
#          for the main location model (03_location_model.R).
#          Includes variable screening (Section 2b) to assess all candidate
#          dependent variables before including them in scripts 03 and 04.
#
# Track level hierarchy: Track -> ROI -> Slice -> Tumour -> Mouse
# ROI level hierarchy:   ROI   -> Slice -> Tumour -> Mouse
#
# Variables: Speed, TrackLength, Straightness, MaxDis, TurningAngle, Alpha
#
# Outputs (all suffixed with output_suffix):
#   variable_screening_summary_<suffix>.csv        : distribution stats per variable
#   variable_screening_distributions_<suffix>.pdf  : histograms + density per variable
#   variable_correlation_heatmap_<suffix>.pdf      : Spearman correlation heatmap
#   icc_track_results_<suffix>.csv                 : variance components + ICC (track)
#   icc_roi_results_<suffix>.csv                   : variance components + ICC (ROI)
#   lrt_track_results_<suffix>.csv                 : LRT per level (track)
#   lrt_roi_results_<suffix>.csv                   : LRT per level (ROI)
#   icc_barcharts_<suffix>.pdf                     : variance decomposition bar charts
#
# Inputs: tracks_clean_v2.csv, roi_summary_v2.csv (from 02e_batch_revised.R)
#         GMM_component column required if gmm_keep != c() (from 02d)
# =============================================================================


# --- PARAMETERS ---------------------------------------------------------------

data_dir          <- "Z:\\lightmicro-inst\\scharalambous\\Resident_cell_motility_tr_v3\\ROI files"  # <-- UPDATE PATH
tracks_clean_path <- "Z:\\lightmicro-inst\\scharalambous\\Resident_cell_motility_tr_v3\\ROI files\\Complete data - retrospective fitting G=5 VEE (from v2)\\tracks_clean_v3_alldata_GMM_retrospective.csv"
roi_summary_path  <- "roi_summary_v3_GMM_1_MaxDis7.csv"

# GMM filter
# Requires GMM_component column (run 02d first). Set c() for full dataset.
gmm_keep <- c(1)

# Output suffix — auto-constructed from gmm_keep
output_suffix <- if (length(gmm_keep) == 0) "all" else
                 paste0("GMM", paste(gmm_keep, collapse = "_"))

# ROI-level analysis gate
# Set TRUE only when gmm_keep is finalised AND roi_summary matches it.
roi_analysis_enabled <- TRUE

# --- Variable toggle ---
# Controls which variables are passed to ICC (Section 5/6) and should be
# mirrored in track_vars / roi_vars in scripts 03 and 04.
# Set to FALSE to exclude a variable from modelling; screening always runs on all.
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
maxdis_threshold <- 7 # µm


# --- 0. PACKAGES --------------------------------------------------------------
# install.packages(c("lme4", "lmerTest", "performance", "tidyverse"))
library(tidyverse)
library(lme4)
library(lmerTest)
library(performance)


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
# run_log   : FALSE = raw only (no log ICC). Omit to default TRUE.
# filter_col / filter_val : applied only for that variable before ICC fitting.
#   Used for alpha_30pct, which requires msd_r2_30pct >= 0.5 for valid estimates.

var_pairs <- list(
  list(
    name       = "Speed",
    raw_col    = "SpeedMean_ummin",    raw_label = "Speed (raw, µm/min)",
    log_col    = "log_Speed",          log_label = "Speed (log)"
  ),
  list(
    name       = "TrackLength",
    raw_col    = "Track_Length",       raw_label = "Track Length (raw, µm)",
    log_col    = "log_TrackLength",    log_label = "Track Length (log)"
  ),
  list(
    name       = "Straightness",
    raw_col    = "Track_Straightness", raw_label = "Straightness (raw)",
    log_col    = "log_Straightness",   log_label = "Straightness (log)"
  ),
  list(
    name       = "MaxDis",
    raw_col    = "MaxDis",             raw_label = "MaxDis (raw, µm)",
    log_col    = "log_MaxDis",         log_label = "MaxDis (log)"
  ),
  list(
    name       = "TurningAngle",
    raw_col    = "mean_turning_angle", raw_label = "Mean Turning Angle (°)",
    log_col    = NULL,                 log_label = NULL,
    run_log    = FALSE
  ),
  list(
    name       = "Alpha",
    raw_col    = "alpha_30pct",        raw_label = "Alpha (MSD, R2≥0.5)",
    log_col    = NULL,                 log_label = NULL,
    run_log    = FALSE,
    filter_col = "msd_r2_30pct",       filter_val = 0.5
  )
)

# ROI-level column names (from roi_summary_v2.csv)
# Alpha uses Alpha30_filt_mean — already filtered to R2>=0.5 by 02e.
var_pairs_roi <- list(
  list(
    name    = "Speed",
    raw_col = "Speed_mean",         raw_label = "Speed (raw, µm/min)",
    log_col = "log_Speed_mean",     log_label = "Speed (log)"
  ),
  list(
    name    = "TrackLength",
    raw_col = "Length_mean",        raw_label = "Track Length (raw, µm)",
    log_col = "log_Length_mean",    log_label = "Track Length (log)"
  ),
  list(
    name    = "Straightness",
    raw_col = "Straightness_mean",  raw_label = "Straightness (raw)",
    log_col = "log_Straight_mean",  log_label = "Straightness (log)"
  ),
  list(
    name    = "MaxDis",
    raw_col = "MaxDis_mean",        raw_label = "MaxDis (raw, µm)",
    log_col = "log_MaxDis_mean",    log_label = "MaxDis (log)"
  ),
  list(
    name    = "TurningAngle",
    raw_col = "TurningAngle_mean",  raw_label = "Mean Turning Angle (°)",
    log_col = NULL,                 log_label = NULL,
    run_log = FALSE
  ),
  list(
    name    = "Alpha",
    raw_col = "Alpha30_filt_mean",  raw_label = "Alpha (MSD, R2≥0.5)",
    log_col = NULL,                 log_label = NULL,
    run_log = FALSE
  )
)


# =============================================================================
# --- 2b. VARIABLE SCREENING ---------------------------------------------------
# =============================================================================
# Runs before ICC to assess all candidate dependent variables.
# Three diagnostics:
#
#   1. Distribution summary: mean, median, SD, skewness, kurtosis, % NA
#      Skewness > |1| = substantially skewed; consider transformation.
#      Kurtosis > 3   = heavy tails relative to normal.
#
#   2. Distribution plots: histogram + density overlay per variable.
#
#   3. Spearman correlation matrix (collinearity check):
#      Variables are modelled separately as outcomes — collinearity here means
#      two variables capture the same biological signal. Pairs with |r| > 0.7
#      flagged as potentially redundant.
#      Reference: Dormann et al. (2013) Ecography 36:27-46
#
# alpha_30pct: internally filtered to msd_r2_30pct >= 0.5 for all checks.
#              n used vs n total is reported.
# =============================================================================

cat("\n=== SECTION 2b: VARIABLE SCREENING ===\n")

# Candidate variables for screening (includes existing + new)
screen_vars <- list(
  list(name = "Speed",          col = "SpeedMean_ummin",    label = "Speed (µm/min)",
       filter_col = NULL,        filter_val = NULL),
  list(name = "log_Speed",      col = "log_Speed",          label = "Log Speed",
       filter_col = NULL,        filter_val = NULL),
  list(name = "TrackLength",    col = "log_TrackLength",    label = "Log Track Length (µm)",
       filter_col = NULL,        filter_val = NULL),
  list(name = "Straightness",   col = "Track_Straightness", label = "Straightness",
       filter_col = NULL,        filter_val = NULL),
  list(name = "log_Straightness", col = "log_Straightness", label = "Log Straightness",
       filter_col = NULL,         filter_val = NULL),
  list(name = "MaxDis",         col = "MaxDis",             label = "MaxDis (µm)",
       filter_col = NULL,        filter_val = NULL),
  list(name = "log_MaxDis",     col = "log_MaxDis",         label = "Log MaxDis (µm)",
       filter_col = NULL,        filter_val = NULL),
  list(name = "TurningAngle",   col = "mean_turning_angle", label = "Mean Turning Angle (°)",
       filter_col = NULL,        filter_val = NULL),
  list(name = "Alpha",          col = "alpha_30pct",        label = "Alpha (MSD)",
       filter_col = "msd_r2_30pct", filter_val = 0.5)
)

# Helper: return clean values with optional row filter applied
get_screen_vals <- function(data, col, filter_col = NULL, filter_val = NULL) {
  if (!col %in% names(data)) return(NULL)
  if (!is.null(filter_col) && filter_col %in% names(data)) {
    data <- data %>% filter(.data[[filter_col]] >= filter_val)
  }
  data[[col]][!is.na(data[[col]])]
}

# Manual skewness and excess kurtosis (no additional package required)
skewness_fn <- function(x) {
  n <- length(x); mu <- mean(x); s <- sd(x)
  (n / ((n - 1) * (n - 2))) * sum(((x - mu) / s)^3)
}
kurtosis_fn <- function(x) {
  n <- length(x); mu <- mean(x); s <- sd(x)
  ((n * (n + 1)) / ((n - 1) * (n - 2) * (n - 3))) *
    sum(((x - mu) / s)^4) -
    (3 * (n - 1)^2) / ((n - 2) * (n - 3))
}

# --- 1. Distribution summary --------------------------------------------------

cat("\n--- Distribution summary (track level) ---\n")

screening_summary <- map_dfr(screen_vars, function(vp) {

  if (!vp$col %in% names(tracks)) {
    cat("  WARNING:", vp$col, "not found in tracks — skipping\n")
    return(NULL)
  }

  n_total <- nrow(tracks)
  x       <- get_screen_vals(tracks, vp$col, vp$filter_col, vp$filter_val)
  n_used  <- length(x)

  filter_note <- if (!is.null(vp$filter_col)) {
    paste0("filtered (", vp$filter_col, " >= ", vp$filter_val,
           "); n = ", n_used, " / ", n_total)
  } else {
    paste0("all; n = ", n_used)
  }

  tibble(
    Variable  = vp$name,
    Column    = vp$col,
    n_used    = n_used,
    n_total   = n_total,
    pct_used  = round(100 * n_used / n_total, 1),
    mean      = round(mean(x), 4),
    median    = round(median(x), 4),
    sd        = round(sd(x), 4),
    min       = round(min(x), 4),
    max       = round(max(x), 4),
    skewness  = round(skewness_fn(x), 3),
    kurtosis  = round(kurtosis_fn(x), 3),
    Filter    = filter_note
  )
})

print(screening_summary, n = Inf)
write_csv(screening_summary,
          paste0("variable_screening_summary_", output_suffix, ".csv"))
cat("Saved: variable_screening_summary_", output_suffix, ".csv\n", sep = "")


# --- 2. Distribution plots ----------------------------------------------------

cat("\n--- Distribution plots ---\n")

pdf(paste0("variable_screening_distributions_", output_suffix, ".pdf"),
    width = 8, height = 4)

for (vp in screen_vars) {

  if (!vp$col %in% names(tracks)) next

  n_total <- nrow(tracks)
  x       <- get_screen_vals(tracks, vp$col, vp$filter_col, vp$filter_val)
  n_used  <- length(x)

  subtitle_text <- if (!is.null(vp$filter_col)) {
    paste0("Filtered: ", vp$filter_col, " >= ", vp$filter_val,
           "  |  n = ", n_used, " / ", n_total)
  } else {
    paste0("n = ", n_used)
  }

  p <- ggplot(tibble(y = x), aes(x = y)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 60, fill = "#457B9D", colour = "white",
                   linewidth = 0.2, alpha = 0.8) +
    geom_density(colour = "#E63946", linewidth = 0.9) +
    labs(
      title    = paste("Distribution —", vp$label),
      subtitle = subtitle_text,
      x        = vp$label,
      y        = "Density"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10, colour = "grey40")
    )

  print(p)
}

dev.off()
cat("Saved: variable_screening_distributions_", output_suffix, ".pdf\n", sep = "")


# --- 3. Spearman correlation heatmap ------------------------------------------
# alpha_30pct set to NA where R2 < 0.5 so pairwise correlations use
# only valid rows for that variable, consistent with the ICC approach.

cat("\n--- Spearman correlations (|r| > 0.7 = potentially redundant) ---\n")

cor_data <- tracks
if ("msd_r2_30pct" %in% names(cor_data)) {
  cor_data <- cor_data %>%
    mutate(alpha_30pct = ifelse(msd_r2_30pct >= 0.5, alpha_30pct, NA_real_))
}

valid_vars  <- screen_vars[map_lgl(screen_vars, ~ .x$col %in% names(cor_data))]
cor_cols    <- map_chr(valid_vars, "col")
cor_labels  <- map_chr(valid_vars, "label")

cor_matrix <- cor_data %>%
  select(all_of(cor_cols)) %>%
  rename_with(~ cor_labels, all_of(cor_cols)) %>%
  cor(method = "spearman", use = "pairwise.complete.obs")

cor_long <- as.data.frame(as.table(cor_matrix)) %>%
  rename(Var1 = Var1, Var2 = Var2, r = Freq) %>%
  filter(Var1 != Var2) %>%
  mutate(
    r       = round(r, 3),
    abs_r   = abs(r),
    flagged = abs_r > 0.7
  ) %>%
  arrange(desc(abs_r))

cat("High correlations (|r| > 0.7):\n")
if (any(cor_long$flagged)) {
  cor_long %>%
    filter(flagged) %>%
    select(Var1, Var2, r) %>%
    distinct() %>%
    as.data.frame() %>%
    print()
} else {
  cat("  None found — no variables exceed the |r| > 0.7 threshold.\n")
}

pdf(paste0("variable_correlation_heatmap_", output_suffix, ".pdf"),
    width = 8, height = 7)

p_cor <- ggplot(cor_long, aes(x = Var1, y = Var2, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(
    aes(label  = sprintf("%.2f", r),
        colour = flagged),
    size = 3
  ) +
  scale_fill_gradient2(
    low      = "#2C7BB6",
    mid      = "white",
    high     = "#D7191C",
    midpoint = 0,
    limits   = c(-1, 1),
    name     = "Spearman r"
  ) +
  scale_colour_manual(
    values = c("FALSE" = "black", "TRUE" = "white"),
    guide  = "none"
  ) +
  labs(
    title    = "Spearman correlation — candidate outcome variables (track level)",
    subtitle = paste0("Bold white = |r| > 0.7 (potentially redundant)  |  suffix: ",
                      output_suffix),
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y      = element_text(size = 9),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 10, colour = "grey40"),
    legend.position  = "right"
  )

print(p_cor)
dev.off()
cat("Saved: variable_correlation_heatmap_", output_suffix, ".pdf\n\n", sep = "")


# --- 3. HELPER: FIT NULL MODEL AND EXTRACT VARIANCE COMPONENTS ----------------

fit_icc <- function(data, col, label, random_formula) {

  full_formula <- as.formula(paste(col, "~", "1 +", random_formula))

  model <- tryCatch(
    lmer(full_formula, data = data, REML = TRUE),
    error = function(e) {
      cat("  Model failed for", label, ":", e$message, "\n")
      return(NULL)
    },
    warning = function(w) {
      m <- lmer(full_formula, data = data, REML = TRUE)
      cat("  Warning for", label, ":", w$message, "\n")
      return(m)
    }
  )

  if (is.null(model)) return(NULL)

  vc        <- as.data.frame(VarCorr(model))
  total_var <- sum(vc$vcov)

  var_table <- vc %>%
    mutate(
      Variable     = label,
      Level = case_when(
        str_detect(grp, "ROIID")   ~ "ROI",
        str_detect(grp, "SliceID") ~ "Slice",
        str_detect(grp, "Tumour")  ~ "Tumour",
        grp == "MouseID"           ~ "Mouse",
        grp == "Residual"          ~ "Residual",
        TRUE                       ~ grp
      ),
      Variance     = round(vcov, 6),
      Pct_Variance = round(100 * vcov / total_var, 1),
      ICC          = round(vcov / total_var, 3)
    ) %>%
    select(Variable, Level, Variance, Pct_Variance, ICC)

  list(model = model, var_table = var_table)
}


# --- 4. HELPER: LIKELIHOOD RATIO TESTS ----------------------------------------
# Tests whether removing each random effect level significantly worsens fit.
# Uses REML=FALSE as required for LRT comparisons between nested random structures.
# Reference: Bates et al. (2015) J Stat Softw 67(1):1-48

run_lrt <- function(data, col, label, level) {

  if (level == "track") {
    m_full      <- paste(col, "~ 1 + (1|MouseID/TumourID_global/SliceID_global/ROIID_global)")
    m_no_roi    <- paste(col, "~ 1 + (1|MouseID/TumourID_global/SliceID_global)")
    m_no_slice  <- paste(col, "~ 1 + (1|MouseID/TumourID_global)")
    m_no_tumour <- paste(col, "~ 1 + (1|MouseID)")
    m_no_mouse  <- paste(col, "~ 1")
    levels_to_test   <- c("ROI", "Slice", "Tumour", "Mouse")
    formulas_reduced <- list(m_no_roi, m_no_slice, m_no_tumour, m_no_mouse)
  } else {
    m_full      <- paste(col, "~ 1 + (1|MouseID/TumourID_global/SliceID_global)")
    m_no_slice  <- paste(col, "~ 1 + (1|MouseID/TumourID_global)")
    m_no_tumour <- paste(col, "~ 1 + (1|MouseID)")
    m_no_mouse  <- paste(col, "~ 1")
    levels_to_test   <- c("Slice", "Tumour", "Mouse")
    formulas_reduced <- list(m_no_slice, m_no_tumour, m_no_mouse)
  }

  fit_model <- function(f) {
    tryCatch(
      lmer(as.formula(f), data = data, REML = FALSE),
      error = function(e) NULL
    )
  }

  full_model <- fit_model(m_full)
  if (is.null(full_model)) {
    full_model <- lm(as.formula(paste(col, "~ 1")), data = data)
  }

  results <- map2_dfr(levels_to_test, formulas_reduced, function(lvl, f) {
    reduced <- fit_model(f)
    if (is.null(reduced)) {
      reduced <- lm(as.formula(paste(col, "~ 1")), data = data)
    }

    lrt <- tryCatch(anova(reduced, full_model), error = function(e) NULL)

    if (is.null(lrt)) {
      return(tibble(Variable = label, Level = lvl,
                    LRT_Chi2 = NA, LRT_df = NA, LRT_p = NA, Significant = NA))
    }

    tibble(
      Variable    = label,
      Level       = lvl,
      LRT_Chi2    = round(lrt[2, "Chisq"], 3),
      LRT_df      = lrt[2, "Df"],
      LRT_p       = round(lrt[2, "Pr(>Chisq)"], 4),
      LRT_p_display = ifelse(lrt[2, "Pr(>Chisq)"] < 0.0001, "< 0.0001",
                             as.character(round(lrt[2, "Pr(>Chisq)"], 4))),
      BIC_full    = round(BIC(full_model), 1),
      BIC_reduced = round(BIC(reduced), 1),
      BIC_diff    = round(BIC(reduced) - BIC(full_model), 1),
      BIC_decision = case_when(
        is.na(BIC(reduced) - BIC(full_model)) ~ "Cannot assess",
        BIC(reduced) - BIC(full_model) > 10   ~ "Strong: keep level",
        BIC(reduced) - BIC(full_model) > 0    ~ "Weak: consider keeping",
        TRUE                                   ~ "Prefer simpler model"
      ),
      Significant = lrt[2, "Pr(>Chisq)"] < 0.05
    )
  })

  results
}


# --- 5. RUN AT TRACK LEVEL ----------------------------------------------------
# run_log defaults to TRUE when the field is absent.
# filter_col / filter_val applied per-variable before ICC fitting.
# For alpha_30pct: ICC runs on the msd_r2_30pct >= 0.5 subset; n is reported.

cat("=== TRACK LEVEL ICC ===\n")

track_icc_list <- list()
track_lrt_list <- list()

for (vp in var_pairs[map_lgl(var_pairs, ~ isTRUE(var_include[.x$name]))]) {
  for (type in c("raw", "log")) {

    # Skip log when run_log = FALSE or log_col is NULL
    if (type == "log" && (isFALSE(vp$run_log) || is.null(vp$log_col))) next

    col   <- if (type == "raw") vp$raw_col   else vp$log_col
    label <- if (type == "raw") vp$raw_label else vp$log_label

    # Apply variable-specific filter (e.g. alpha_30pct requires R2 >= 0.5)
    if (!is.null(vp$filter_col)) {
      n_all  <- nrow(tracks)
      data_v <- tracks %>% filter(.data[[vp$filter_col]] >= vp$filter_val)
      cat("  Processing:", label,
          "| filtered:", vp$filter_col, ">=", vp$filter_val,
          "| n (post-R2 filter):", nrow(data_v), "/", n_all,
          "| n (model-ready, excl. NA):", sum(!is.na(data_v[[col]])), "\n")
    } else {
      data_v <- tracks
      cat("  Processing:", label, "\n")
    }

    result <- fit_icc(
      data           = data_v,
      col            = col,
      label          = label,
      random_formula = "(1|MouseID/TumourID_global/SliceID_global/ROIID_global)"
    )
    if (!is.null(result)) track_icc_list[[label]] <- result$var_table

    lrt <- tryCatch(
      run_lrt(data_v, col, label, level = "track"),
      error = function(e) { cat("  LRT failed:", e$message, "\n"); NULL }
    )
    if (!is.null(lrt)) track_lrt_list[[label]] <- lrt
  }
}

icc_track <- bind_rows(track_icc_list)
lrt_track <- bind_rows(track_lrt_list)


# --- 6. RUN AT ROI LEVEL (gated) ----------------------------------------------

if (!roi_analysis_enabled) {
  cat("\nROI-level ICC skipped (roi_analysis_enabled = FALSE).\n")
  cat("Set roi_analysis_enabled <- TRUE once final roi_summary is confirmed.\n\n")
  icc_roi <- tibble()
  lrt_roi <- tibble()
} else {

  cat("\n=== ROI LEVEL ICC ===\n")

  roi_icc_list <- list()
  roi_lrt_list <- list()

  for (vp in var_pairs_roi[map_lgl(var_pairs_roi, ~ isTRUE(var_include[.x$name]))]) {
    for (type in c("raw", "log")) {

      if (type == "log" && (isFALSE(vp$run_log) || is.null(vp$log_col))) next

      col   <- if (type == "raw") vp$raw_col   else vp$log_col
      label <- if (type == "raw") vp$raw_label else vp$log_label

      cat("  Processing:", label, "\n")

      result <- fit_icc(
        data           = roi,
        col            = col,
        label          = label,
        random_formula = "(1|MouseID/TumourID_global/SliceID_global)"
      )
      if (!is.null(result)) roi_icc_list[[label]] <- result$var_table

      lrt <- tryCatch(
        run_lrt(roi, col, label, level = "roi"),
        error = function(e) { cat("  LRT failed:", e$message, "\n"); NULL }
      )
      if (!is.null(lrt)) roi_lrt_list[[label]] <- lrt
    }
  }

  icc_roi <- bind_rows(roi_icc_list)
  lrt_roi <- bind_rows(roi_lrt_list)

}


# --- 7. PRINT RESULTS ---------------------------------------------------------

cat("\n--- Track level ICC ---\n")
print(icc_track, n = Inf)
cat("\n--- Track level LRT ---\n")
print(lrt_track, n = Inf)

if (roi_analysis_enabled) {
  cat("\n--- ROI level ICC ---\n")
  print(icc_roi, n = Inf)
  cat("\n--- ROI level LRT ---\n")
  print(lrt_roi, n = Inf)
}


# --- 8. DECISION GUIDE --------------------------------------------------------

print_decision <- function(icc_data, level_name) {
  cat("\n--- Random effects decision guide (", level_name, ") ---\n", sep = "")
  cat("ICC thresholds: >0.3 must include | 0.1-0.3 should include | <0.1 consider dropping\n\n")

  icc_data %>%
    filter(Level != "Residual") %>%
    group_by(Level) %>%
    summarise(Mean_ICC = round(mean(ICC, na.rm = TRUE), 3), .groups = "drop") %>%
    mutate(
      Recommendation = case_when(
        Mean_ICC >= 0.3 ~ "MUST include — high clustering",
        Mean_ICC >= 0.1 ~ "SHOULD include — moderate clustering",
        TRUE            ~ "Consider dropping — low clustering"
      )
    ) %>%
    print()
}

print_decision(icc_track, "Track level")
if (roi_analysis_enabled && nrow(icc_roi) > 0) print_decision(icc_roi, "ROI level")


# --- 9. SAVE RESULTS ----------------------------------------------------------

write_csv(icc_track, paste0("icc_track_results_", output_suffix, ".csv"))
write_csv(lrt_track, paste0("lrt_track_results_", output_suffix, ".csv"))
cat("\nSaved: icc_track_results_", output_suffix, ".csv\n", sep = "")
cat("Saved: lrt_track_results_", output_suffix, ".csv\n", sep = "")

if (roi_analysis_enabled && nrow(icc_roi) > 0) {
  write_csv(icc_roi, paste0("icc_roi_results_", output_suffix, ".csv"))
  write_csv(lrt_roi, paste0("lrt_roi_results_", output_suffix, ".csv"))
  cat("Saved: icc_roi_results_", output_suffix, ".csv\n", sep = "")
  cat("Saved: lrt_roi_results_", output_suffix, ".csv\n", sep = "")
}


# --- 10. BAR CHARTS -----------------------------------------------------------
# One chart per variable in var_pairs.
# Panels: Track Raw | Track Log (when run_log = TRUE) | ROI Raw | ROI Log
# ROI panels only appear when roi_analysis_enabled = TRUE.
# Variables with run_log = FALSE (TurningAngle, Alpha) show Track Raw only
# (or Track Raw + ROI Raw when gate is open).

level_colours <- c(
  "Mouse"    = "#2C7BB6",
  "Tumour"   = "#ABD9E9",
  "Slice"    = "#FDAE61",
  "ROI"      = "#A1D99B",
  "Residual" = "#D7191C"
)

make_icc_chart <- function(vp) {

  # Find the matching roi entry by name
  vp_roi <- Filter(function(x) x$name == vp$name, var_pairs_roi)
  vp_roi <- if (length(vp_roi) > 0) vp_roi[[1]] else NULL
  
  rl <- vp$raw_label
  ll <- vp$log_label   # NULL for variables with run_log = FALSE

  cat("\nBuilding chart for:", vp$name, "\n")

  # --- Track data ---
  track_raw_data <- icc_track[icc_track$Variable == rl, ]
  track_raw_data$Analysis  <- "Track level"
  track_raw_data$Transform <- "Raw"

  track_log_data <- if (!is.null(ll)) {
    d <- icc_track[icc_track$Variable == ll, ]
    d$Analysis  <- "Track level"
    d$Transform <- "Log"
    d
  } else NULL

  # --- ROI data (only if gate open) ---
  roi_raw_data <- if (roi_analysis_enabled && nrow(icc_roi) > 0) {
    d <- icc_roi[icc_roi$Variable == vp_roi$raw_label, ]
    if (nrow(d) == 0) return(NULL)
    d$Analysis  <- "ROI level"
    d$Transform <- "Raw"
    d
  } else NULL

  roi_log_data <- if (roi_analysis_enabled && nrow(icc_roi) > 0 && !is.null(ll)) {
    d <- icc_roi[icc_roi$Variable == vp_roi$log_label, ]
    if (nrow(d) == 0) return(NULL)
    d$Analysis  <- "ROI level"
    d$Transform <- "Log"
    d
  } else NULL
  
  # Build panel factor levels dynamically
  panel_levels <- "Track level\nRaw"
  if (!is.null(track_log_data) && nrow(track_log_data) > 0)
    panel_levels <- c(panel_levels, "Track level\nLog")
  if (!is.null(roi_raw_data)   && nrow(roi_raw_data) > 0)
    panel_levels <- c(panel_levels, "ROI level\nRaw")
  if (!is.null(roi_log_data)   && nrow(roi_log_data) > 0)
    panel_levels <- c(panel_levels, "ROI level\nLog")

  plot_data <- bind_rows(
    track_raw_data, track_log_data, roi_raw_data, roi_log_data
  ) %>%
    mutate(
      Panel = paste(Analysis, Transform, sep = "\n"),
      Panel = factor(Panel, levels = panel_levels),
      Level = factor(Level, levels = c("Mouse", "Tumour", "Slice", "ROI", "Residual"))
    )

  subtitle_text <- if (!is.null(vp$filter_col)) {
    paste0("alpha filtered to ", vp$filter_col, " >= ", vp$filter_val)
  } else NULL

  ggplot(plot_data, aes(x = Panel, y = Pct_Variance, fill = Level)) +
    geom_bar(stat = "identity", position = "stack", width = 0.6, na.rm = TRUE) +
    geom_text(
      aes(label = ifelse(!is.na(Pct_Variance) & Pct_Variance >= 5,
                         paste0(Pct_Variance, "%"), "")),
      position = position_stack(vjust = 0.5),
      size = 3, colour = "white", fontface = "bold", na.rm = TRUE
    ) +
    scale_fill_manual(values = level_colours, drop = FALSE) +
    scale_y_continuous(limits = c(0, 105),
                       labels = function(x) paste0(x, "%")) +
    labs(
      title    = paste("Variance decomposition —", vp$name),
      subtitle = subtitle_text,
      x        = NULL,
      y        = "% Variance",
      fill     = "Level"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title       = element_text(face = "bold", size = 13),
      plot.subtitle    = element_text(size = 10, colour = "grey40"),
      axis.text.x      = element_text(size = 9),
      legend.position  = "right",
      panel.grid.minor = element_blank()
    )
}

pdf(paste0("icc_barcharts_", output_suffix, ".pdf"), width = 8, height = 5)
for (vp in var_pairs[map_lgl(var_pairs, ~ isTRUE(var_include[.x$name]))]) {
  p <- make_icc_chart(vp)
  print(p)
}
dev.off()

cat("Saved: icc_barcharts_", output_suffix, ".pdf\n", sep = "")
cat("\nNext step: review screening outputs and ICC results,\n")
cat("then update track_vars / roi_vars in 03 and 04 as needed.\n")
