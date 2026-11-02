# =============================================================================
# CELL TRACKING ANALYSIS - STEP 2b: EXPLORATORY VISUALISATION
# =============================================================================
# Purpose: Visualise motility parameters across Location, Treatment and Day
#          before formal modelling. Supports all candidate variables from 02
#          with the same toggle and GMM filter used in scripts 02, 03, 04.
#
# Layout per parameter:
#   - 2x2 facet grid: Treatment (Mock/RP) columns, Day (D1/D2) rows
#   - Within each panel: 3 location groups on x-axis
#   - Violin + jittered dots + mean±SD overlay
#   - Two colour versions per variable: coloured by Location, coloured by Mouse
#
# Output (suffixed with output_suffix):
#   exploratory_plots_<suffix>.pdf
#     Track level — Location coloured  (all toggled variables, raw then log)
#     Track level — Mouse coloured
#     ROI level   — Location coloured  (gated by roi_analysis_enabled)
#     ROI level   — Mouse coloured     (gated by roi_analysis_enabled)
#
# Inputs: tracks_clean_v2.csv, roi_summary_v2.csv (from 02e_batch_revised.R)
#         GMM_component column required if gmm_keep != c() (from 02d)
# =============================================================================


# --- PARAMETERS ---------------------------------------------------------------

data_dir          <- "Z:\\lightmicro-inst\\scharalambous\\Resident_cell_motility_tr_v3\\ROI files"  # <-- UPDATE PATH
tracks_clean_path <- "Z:\\lightmicro-inst\\scharalambous\\Resident_cell_motility_tr_v3\\ROI files\\Complete data - retrospective fitting G=5 VEE (from v2)\\tracks_clean_v3_alldata_GMM_retrospective.csv"
roi_summary_path  <- "roi_summary_v3_GMM_1_MaxDis7.csv"


# GMM filter — must match what was used in 02 and 03/04
gmm_keep <- c(1)

output_suffix <- if (length(gmm_keep) == 0) "all" else
                 paste0("GMM", paste(gmm_keep, collapse = "_"))

# ROI-level gate — set TRUE only when roi_summary matches the final dataset
roi_analysis_enabled <- TRUE

# Alpha R² threshold — must match 02
alpha_r2_threshold <- 0.5

# MaxDis displacement filter
# Applied after GMM filter within retained components.
# Removes near-stationary tracks that never meaningfully displaced
# from their starting position regardless of GMM component assignment.
# Threshold: MaxDis > 1.1 µm (equivalent to log_MaxDis > 0.095)
maxdis_threshold <- 7  # µm

# --- Variable toggle ----------------------------------------------------------
# Mirrors var_include from scripts 02, 03, 04.
# Each TRUE entry produces plots for both raw and log versions where available.
var_include <- c(
  "Speed"        = TRUE,
  "TrackLength"  = TRUE,
  "Straightness" = TRUE,
  "MaxDis"       = FALSE,
  "TurningAngle" = TRUE,
  "Alpha"        = TRUE
)


# --- 0. PACKAGES --------------------------------------------------------------
# install.packages(c("tidyverse", "patchwork"))
library(tidyverse)
library(patchwork)


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


# --- 2. COLOUR PALETTES -------------------------------------------------------

location_colours <- c(
  "Core"     = "#E63946",
  "CoreEdge" = "#457B9D",
  "Edge"     = "#2A9D8F"
)

# 5-mouse palette — extend if more mice are added
mouse_colours <- c(
  "M1" = "#457B9D",
  "M2" = "#E63946",
  "M3" = "#F4A261",
  "M4" = "#2A9D8F",
  "M5" = "#6A4C93"  ##"6A0572"
)


# --- 3. DEFINE VARIABLES ------------------------------------------------------
# Each entry has raw and log versions where applicable.
# run_raw / run_log control which transforms are plotted.
# filter_col / filter_val applied before plotting (alpha_30pct only).
# Column names follow tracks_clean_v2.csv and roi_summary_v2.csv conventions.

track_vars <- list(
  list(name = "Speed",
       raw_col = "SpeedMean_ummin",    raw_label = "Speed (µm/min)",
       log_col = "log_Speed",          log_label = "Log Speed (µm/min)",
       run_raw = TRUE,  run_log = TRUE,
       filter_col = NULL, filter_val = NULL),
  list(name = "TrackLength",
       raw_col = "Track_Length",       raw_label = "Track Length (µm)",
       log_col = "log_TrackLength",    log_label = "Log Track Length (µm)",
       run_raw = FALSE, run_log = TRUE,
       filter_col = NULL, filter_val = NULL),
  list(name = "Straightness",
       raw_col = "Track_Straightness", raw_label = "Straightness",
       log_col = "log_Straightness",   log_label = "Log Straightness",
       run_raw = TRUE,  run_log = TRUE,
       filter_col = NULL, filter_val = NULL),
  list(name = "MaxDis",
       raw_col = "MaxDis",             raw_label = "MaxDis (µm)",
       log_col = "log_MaxDis",         log_label = "Log MaxDis (µm)",
       run_raw = TRUE,  run_log = TRUE,
       filter_col = NULL, filter_val = NULL),
  list(name = "TurningAngle",
       raw_col = "mean_turning_angle", raw_label = "Mean Turning Angle (°)",
       log_col = NULL,                 log_label = NULL,
       run_raw = TRUE,  run_log = FALSE,
       filter_col = NULL, filter_val = NULL),
  list(name = "Alpha",
       raw_col = "alpha_30pct",        raw_label = "Alpha (MSD)",
       log_col = NULL,                 log_label = NULL,
       run_raw = TRUE,  run_log = FALSE,
       filter_col = "msd_r2_30pct",    filter_val = alpha_r2_threshold)
)

roi_vars <- list(
  list(name = "Speed",
       raw_col = "Speed_mean",         raw_label = "Speed (µm/min)",
       log_col = "log_Speed_mean",     log_label = "Log Speed (µm/min)",
       run_raw = TRUE,  run_log = TRUE,
       filter_col = NULL, filter_val = NULL),
  list(name = "TrackLength",
       raw_col = "Length_mean",        raw_label = "Track Length (µm)",
       log_col = "log_Length_mean",    log_label = "Log Track Length (µm)",
       run_raw = FALSE, run_log = TRUE,
       filter_col = NULL, filter_val = NULL),
  list(name = "Straightness",
       raw_col = "Straightness_mean",  raw_label = "Straightness",
       log_col = "log_Straight_mean",  log_label = "Log Straightness",
       run_raw = TRUE,  run_log = TRUE,
       filter_col = NULL, filter_val = NULL),
  list(name = "MaxDis",
       raw_col = "MaxDis_mean",        raw_label = "MaxDis (µm)",
       log_col = "log_MaxDis_mean",    log_label = "Log MaxDis (µm)",
       run_raw = TRUE,  run_log = TRUE,
       filter_col = NULL, filter_val = NULL),
  list(name = "TurningAngle",
       raw_col = "TurningAngle_mean",  raw_label = "Mean Turning Angle (°)",
       log_col = NULL,                 log_label = NULL,
       run_raw = TRUE,  run_log = FALSE,
       filter_col = NULL, filter_val = NULL),
  list(name = "Alpha",
       raw_col = "Alpha30_filt_mean",  raw_label = "Alpha (MSD, filtered)",
       log_col = NULL,                 log_label = NULL,
       run_raw = TRUE,  run_log = FALSE,
       filter_col = NULL, filter_val = NULL)   # already filtered in roi_summary
)


# --- 4. PLOT FUNCTION ---------------------------------------------------------
# Makes one faceted figure for a single parameter.
#
# data         : tracks or roi
# col          : column name to plot on y-axis
# ylabel       : y-axis label
# title        : plot title
# point_size   : jitter dot size (smaller for dense track-level data)
# alpha        : dot transparency
# colour_var   : "Location" or "MouseID"
# colour_palette: named colour vector matching levels of colour_var
# subtitle     : optional subtitle (used to note alpha filter)

make_facet_plot <- function(data, col, ylabel, title,
                            point_size = 1, alpha = 0.3,
                            colour_var = "Location",
                            colour_palette = location_colours,
                            subtitle = NULL) {

  plot_data <- data %>%
    rename(y = all_of(col)) %>%
    filter(!is.na(y), !is.na(Location), !is.na(Treatment), !is.na(Day))

  summary_data <- plot_data %>%
    group_by(Location, Treatment, Day) %>%
    summarise(
      mean_y = mean(y, na.rm = TRUE),
      sd_y   = sd(y,   na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      ymin = mean_y - sd_y,
      ymax = mean_y + sd_y
    )

  ggplot() +
    geom_violin(
      data      = plot_data,
      aes(x = Location, y = y),
      fill      = "grey85",
      colour    = "grey50",
      linewidth = 0.4,
      alpha     = 0.6,
      trim      = FALSE,
      scale     = "width"
    ) +
    geom_jitter(
      data   = plot_data,
      aes(x = Location, y = y, colour = .data[[colour_var]]),
      width  = 0.2,
      size   = point_size,
      alpha  = alpha
    ) +
    geom_point(
      data   = summary_data,
      aes(x = Location, y = mean_y),
      shape  = 21,
      size   = 3,
      fill   = "white",
      colour = "black",
      stroke = 0.8
    ) +
    geom_errorbar(
      data      = summary_data,
      aes(x = Location, ymin = ymin, ymax = ymax),
      width     = 0.15,
      colour    = "black",
      linewidth = 0.6
    ) +
    facet_grid(Day ~ Treatment) +
    scale_colour_manual(values = colour_palette, name = colour_var) +
    labs(
      title    = title,
      subtitle = subtitle,
      x        = "Location",
      y        = ylabel,
      colour   = colour_var
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 12),
      plot.subtitle    = element_text(size = 10, colour = "grey40"),
      strip.background = element_rect(fill = "grey90", colour = "grey60"),
      strip.text       = element_text(face = "bold", size = 10),
      axis.text.x      = element_text(angle = 30, hjust = 1, size = 9),
      legend.position  = "bottom",
      panel.grid.minor = element_blank()
    )
}


# --- 5. BUILD AND SAVE PLOTS --------------------------------------------------
# Loop over all toggled variables and both transforms.
# Alpha at track level: filters to msd_r2_30pct >= alpha_r2_threshold rows.
# ROI level: gated by roi_analysis_enabled.

pdf(paste0("exploratory_plots_", output_suffix, ".pdf"), width = 8, height = 7)

# Helper: build both Location and Mouse versions for one variable/transform
print_pair <- function(data, col, ylabel, level_name,
                       point_size, alpha, subtitle = NULL) {

  p_loc <- make_facet_plot(
    data           = data,
    col            = col,
    ylabel         = ylabel,
    title          = paste(level_name, "—", ylabel, "| by Location"),
    point_size     = point_size,
    alpha          = alpha,
    colour_var     = "Location",
    colour_palette = location_colours,
    subtitle       = subtitle
  )

  p_mouse <- make_facet_plot(
    data           = data,
    col            = col,
    ylabel         = ylabel,
    title          = paste(level_name, "—", ylabel, "| by Mouse"),
    point_size     = point_size,
    alpha          = alpha,
    colour_var     = "MouseID",
    colour_palette = mouse_colours,
    subtitle       = subtitle
  )

  print(p_loc)
  print(p_mouse)
}

# --- Track level --------------------------------------------------------------
cat("\nBuilding track-level plots...\n")

active_track <- track_vars[map_lgl(track_vars, ~ isTRUE(var_include[.x$name]))]

for (vp in active_track) {

  # Prepare data — apply alpha filter if needed
  if (!is.null(vp$filter_col)) {
    data_v   <- tracks %>% filter(.data[[vp$filter_col]] >= vp$filter_val)
    n_used   <- nrow(data_v)
    n_total  <- nrow(tracks)
    subtitle <- paste0("Filtered: ", vp$filter_col, " >= ", vp$filter_val,
                       "  |  n = ", n_used, " / ", n_total)
    cat("  Alpha filter applied:", n_used, "/", n_total, "tracks\n")
  } else {
    data_v   <- tracks
    subtitle <- NULL
  }

  # Raw version
  if (isTRUE(vp$run_raw) && !is.null(vp$raw_col)) {
    cat("  Plotting:", vp$raw_label, "(track raw)\n")
    print_pair(data_v, vp$raw_col, vp$raw_label, "Track",
               point_size = 0.8, alpha = 0.2, subtitle = subtitle)
  }

  # Log version
  if (isTRUE(vp$run_log) && !is.null(vp$log_col)) {
    cat("  Plotting:", vp$log_label, "(track log)\n")
    print_pair(data_v, vp$log_col, vp$log_label, "Track",
               point_size = 0.8, alpha = 0.2, subtitle = subtitle)
  }
}

# --- ROI level ----------------------------------------------------------------
if (!roi_analysis_enabled) {
  cat("\nROI-level plots skipped (roi_analysis_enabled = FALSE).\n")
  cat("Set roi_analysis_enabled <- TRUE once final roi_summary is confirmed.\n")
} else {
  cat("\nBuilding ROI-level plots...\n")

  active_roi <- roi_vars[map_lgl(roi_vars, ~ isTRUE(var_include[.x$name]))]

  for (vp in active_roi) {

    if (isTRUE(vp$run_raw) && !is.null(vp$raw_col)) {
      cat("  Plotting:", vp$raw_label, "(ROI raw)\n")
      print_pair(roi, vp$raw_col, vp$raw_label, "ROI",
                 point_size = 1.8, alpha = 0.6)
    }

    if (isTRUE(vp$run_log) && !is.null(vp$log_col)) {
      cat("  Plotting:", vp$log_label, "(ROI log)\n")
      print_pair(roi, vp$log_col, vp$log_label, "ROI",
                 point_size = 1.8, alpha = 0.6)
    }
  }
}

dev.off()

cat("\nSaved: exploratory_plots_", output_suffix, ".pdf\n", sep = "")
cat("Each variable produces 2 pages (Location coloured, Mouse coloured)\n")
cat("per transform (raw and/or log) at each level (Track; ROI if gated open)\n")
