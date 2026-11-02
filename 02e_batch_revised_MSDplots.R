# =============================================================================
# CELL TRACKING ANALYSIS - STEP 02e BATCH (REVISED): TURNING ANGLE & MSD
# =============================================================================
# Purpose: Batch process all XML track files, computing per-track turning
#          angle and MSD metrics, joining to merged_tracks CSVs, combining
#          into one large CSV, applying duration + QC filters, and producing
#          a clean track-level CSV and ROI-level summary.
#
# Processing pipeline:
#   STAGE 1 — Per ROI (batch loop):
#     a. Locate _merged_tracks.csv and _Tracks.xml in same subfolder
#     b. Load merged_tracks CSV (raw, unfiltered) — 4 header rows
#     c. Parse XML — all particles, no filter
#     d. Flag count mismatch if XML particles != CSV rows
#     e. Apply duration filter: keep tracks > 10% of max TrackDuration in ROI
#     f. Compute turning angle per track from XML coordinates
#     g. Compute MSD per track — alpha_30pct + msd_r2_30pct (weighted,
#        lag <= 30% of max track duration among surviving tracks)
#     h. Join metrics to merged_tracks CSV positionally (track i = row i)
#     i. Save per-ROI joined CSV to same subfolder as source files
#     j. Accumulate into combined list
#
#   STAGE 2 — Combine:
#     Bind all per-ROI results into one large CSV
#
#   STAGE 3 — Filter step 1:
#     Retain columns of interest only — save as filtered CSV
#
#   STAGE 4 — Quality control checks:
#     Flag tracks with implausible speed, length, or straightness values.
#     Save flagged_tracks.csv for inspection.
#     Removal controlled by remove_flagged toggle in Section 1.
#
#   SECTION 8b — Derived columns:
#     Compute inv_outreach_ratio, Translocating flag, log_TrackDuration.
#     Save tracks_clean_v2.csv.
#
#   SECTION 9 — ROI-level summary:
#     One row per ROI. Raw + transformed + new metric summaries.
#     Alpha summaries restricted to tracks with msd_r2_30pct >= 0.5.
#     Save roi_summary_v2.csv.
#
# Outputs (working directory):
#   turning_msd_combined_raw.csv  : all ROIs, all columns (Stage 2)
#   turning_msd_filtered_2.csv    : columns of interest only (Stage 3)
#   flagged_tracks.csv            : tracks failing QC flags (Stage 4)
#   tracks_clean_v2.csv           : final clean track-level CSV (Section 8b)
#   roi_summary_v2.csv            : ROI-level summary (Section 9)
#   missing_xml_report.csv        : ROIs with CSV but no XML
#   batch_processing_log.txt      : per-ROI processing summary
#   batch_processing_log.csv      : machine-readable log
#
# Per-ROI outputs (in experiment subfolder alongside source files):
#   <prefix>_turning_msd.csv      : merged_tracks + metrics for that ROI
#
# Optional MSD curve PDFs (working directory):
#   msd_<ROIID_global>.pdf        : one per requested ROI
#
# Inputs:
#   experimental_file_list.csv    : metadata master list
#   *_merged_tracks.csv           : unfiltered track CSVs (4 header rows)
#   *_Tracks.xml                  : renamed XML files (prefix matches CSV)
# =============================================================================



# --- 0. PACKAGES --------------------------------------------------------------
# install.packages(c("xml2", "tidyverse"))

library(xml2)
library(tidyverse)


# --- 1. CONFIGURATION ---------------------------------------------------------

setwd("Z:\\lightmicro-inst\\scharalambous\\Resident_cell_motility_tr_v3\\ROI files")

metadata_file <- "experimental_file_list_newly_analysed_JUN26.csv"
data_folder   <- "Z:\\lightmicro-inst\\scharalambous\\Resident_cell_motility_tr_v3\\ROI files"
#tracks_clean_file <- "tracks_clean.csv"

maxdis_threshold <- 1.1          # µm — minimum MaxDis to flag as Translocating  # <----- this was used previously before GMM was used for classification
output_clean     <- "tracks_clean_v3_newdata.csv"   # CHANGE: output filename

# --- QC thresholds (used in Stage 4 / Section 8) ---
speed_max_ummin <- 3000   # µm/min — upper plausibility bound   ###  a bit ridiculous but ensures all are included
remove_flagged  <- TRUE  # Set TRUE after reviewing flagged_tracks.csv

# MSD fitting thresholds
min_lags_30pct <- 3     # minimum unique lags for 30% truncated regression

# Required columns from merged_tracks CSV (row 1 = column names, skip rows 2-4)
required_cols <- c(
  "TRACK_INDEX",
  "TRACK_MEAN_SPEED",
  "TOTAL_DISTANCE_TRAVELED",
  "MAX_DISTANCE_TRAVELED",
  "CONFINEMENT_RATIO",
  "TRACK_START",
  "TRACK_DURATION"
)

friendly_names <- c(
  "ID",
  "Track_Speed_Mean",
  "Track_Length",
  "MaxDis",
  "Track_Straightness",
  "t0",
  "TrackDuration"
)

# Columns to retain in Stage 3 (filter step 1)
# Includes: track features, metadata, computed metrics
keep_cols_stage3 <- c(
  # Metadata
  "MouseID", "TumourID_global", "SliceID_global", "ROIID_global",
  "Location", "Treatment", "Day",
  "Filepath", "Filename",
  # Track features (raw)
  "ID", "Track_Speed_Mean", "Track_Length", "MaxDis",
  "Track_Straightness", "TrackDuration", "t0",
  # Computed metrics
  "mean_turning_angle", "alpha_30pct", "msd_r2_30pct"
)

# Match key for Stage 4 join to tracks_clean.csv
# Uses raw pre-conversion values present in both files
match_key <- c(
  "ROIID_global",
  "Track_Speed_Mean",   # = SpeedMean_ums in tracks_clean
  "Track_Length",       # = TrackLength   in tracks_clean
  "MaxDis",
  "TrackDuration"       # = TrackDuration  in tracks_clean
)

# Optional MSD curve PDFs — specify ROIID_global values
# Leave empty c() to skip all PDFs
msd_pdf_rois       <- c("M1_T1_RP_S4_ROI3")
msd_pdf_max_tracks <- 20

# # tracks_clean.csv column name mapping for Stage 4 match key
# # Left = name in combined CSV | Right = name in tracks_clean.csv
# tracks_clean_key_map <- c(
#   "Track_Speed_Mean" = "SpeedMean_ums",
#   "Track_Length"     = "TrackLength",
#   "MaxDis"           = "MaxDis",
#   "TrackDuration"    = "TrackDuration"
# )


# --- 2. LOAD METADATA ---------------------------------------------------------

cat("=== 02e BATCH (REVISED): TURNING ANGLE & MSD ===\n\n")
cat("Loading metadata:", metadata_file, "\n")

metadata <- read_csv(metadata_file, show_col_types = FALSE) %>%
  select(Filename, MouseID, TumourID, Day,
         Treatment, SliceID, ROIID, Location) %>%
  mutate(
    Filename  = paste0(tools::file_path_sans_ext(Filename), ".csv"),
    Location  = factor(Location,  levels = c("Core", "CoreEdge", "Edge")),
    Treatment = factor(Treatment, levels = c("Mock", "RP")),
    Day       = factor(Day,       levels = c("D1", "D2")),
    MouseID   = factor(MouseID),
    TumourID  = factor(TumourID),
    TumourID_global = factor(paste(MouseID, TumourID, sep = "_")),
    SliceID_global  = paste(MouseID, TumourID, Treatment, SliceID, sep = "_"),
    ROIID_global    = paste(MouseID, TumourID, Treatment, SliceID, ROIID,
                            sep = "_")
  )

cat("Metadata loaded:", nrow(metadata), "ROIs\n\n")


# --- 3. LOCATE FILES ----------------------------------------------------------
# For each ROI in metadata:
#   merged_tracks CSV: matched by Filename column in metadata
#   XML:              derived by stripping _merged_tracks.csv -> _Tracks.xml

cat("Locating merged_tracks CSVs and XML files...\n")

all_csvs <- list.files(data_folder, pattern = "_tracks\\.csv$",
                       full.names = TRUE, recursive = TRUE) %>%
  .[!grepl("tracks_clean|flagged_tracks|turning_msd", .)]

all_xmls <- list.files(data_folder, pattern = "_Tracks\\.xml$",
                       full.names = TRUE, recursive = TRUE)

cat("Merged-tracks CSVs found:", length(all_csvs), "\n")
cat("XML files found:         ", length(all_xmls), "\n\n")

# XML lookup: prefix (everything before _Tracks.xml) -> full path
xml_lookup <- setNames(
  all_xmls,
  basename(all_xmls) %>% str_remove("_Tracks\\.xml$")
)

# Attach file paths to metadata
metadata <- metadata %>%
  rowwise() %>%
  mutate(
    # Full path to merged_tracks CSV
    csv_path = {
      match <- all_csvs[basename(all_csvs) == Filename]
      if (length(match) == 1) match else NA_character_
    },
    # XML prefix = CSV filename stripped of _merged_tracks.csv
    xml_prefix = tools::file_path_sans_ext(Filename) %>%
      str_remove("\\s*-\\s*Position.*$") %>%   # strip " - Position 1 [16] 3DMontage..." 
      str_remove("\\s*[-_]\\s*merged.*$"),      # existing strip kept for safety
    # Full path to XML
    xml_path = {
      match <- xml_lookup[xml_prefix]
      if (!is.na(match) && length(match) == 1) unname(match)
      else NA_character_
    },
    csv_exists = !is.na(csv_path) & file.exists(csv_path),
    xml_exists = !is.na(xml_path) & file.exists(xml_path)
  ) %>%
  ungroup()

# Report missing XMLs
missing_xml <- metadata %>%
  filter(csv_exists & !xml_exists) %>%
  select(MouseID, TumourID_global, Day, Treatment,
         SliceID_global, ROIID_global, Filename, xml_prefix)

cat("ROIs with CSV but missing XML:", nrow(missing_xml), "\n")
if (nrow(missing_xml) > 0) {
  cat("Missing XMLs:\n")
  walk(missing_xml$xml_prefix, ~ cat("  ", .x, "_Tracks.xml\n"))
  write_csv(missing_xml, "missing_xml_report.csv")
  cat("Saved: missing_xml_report.csv\n")
}

to_process <- metadata %>% filter(csv_exists & xml_exists)
cat("\nROIs to process:", nrow(to_process), "\n\n")


# --- 4. CORE FUNCTIONS --------------------------------------------------------

# 4a. Load merged_tracks CSV
load_merged_tracks <- function(csv_path, roi_meta) {

  # Row 1 = column names; rows 2-4 = sub-headers/units; data from row 5
  col_names <- read.csv(csv_path, nrows = 1, header = FALSE,
                        check.names = TRUE) %>%
    unlist() %>% as.character()

  raw <- read.csv(csv_path, skip = 4, header = FALSE,
                  fill = TRUE, check.names = FALSE)
  colnames(raw) <- col_names

  # Standardise to UPPER_SNAKE_CASE
  colnames(raw) <- colnames(raw) %>%
    str_trim() %>%
    str_replace_all("[. ]+", "_") %>%
    str_to_upper()

  # Check required columns
  missing_cols <- setdiff(required_cols, colnames(raw))
  if (length(missing_cols) > 0) {
    warning(paste("Missing columns in", basename(csv_path), ":",
                  paste(missing_cols, collapse = ", ")))
    return(NULL)
  }

  raw %>%
    select(all_of(required_cols)) %>%
    setNames(friendly_names) %>%
    mutate(
      Filepath        = csv_path,
      Filename        = basename(csv_path),
      MouseID         = as.character(roi_meta$MouseID),
      TumourID_global = as.character(roi_meta$TumourID_global),
      SliceID_global  = roi_meta$SliceID_global,
      ROIID_global    = roi_meta$ROIID_global,
      Location        = as.character(roi_meta$Location),
      Treatment       = as.character(roi_meta$Treatment),
      Day             = as.character(roi_meta$Day),
      track_index     = row_number()
    )
}

# 4b. Parse XML — no duration filter, all particles
parse_xml_tracks <- function(xml_path) {
  doc              <- read_xml(xml_path)
  frame_interval_s <- as.numeric(xml_attr(doc, "frameInterval"))
  particle_nodes   <- xml_find_all(doc, ".//particle")

  tracks_raw <- lapply(seq_along(particle_nodes), function(i) {
    node      <- particle_nodes[[i]]
    det_nodes <- xml_find_all(node, "detection")
    n         <- length(det_nodes)
    if (n == 0) return(NULL)
    tibble(
      particle_index = i,
      n_spots        = n,
      t              = as.integer(xml_attr(det_nodes, "t")),
      x              = as.numeric(xml_attr(det_nodes, "x")),
      y              = as.numeric(xml_attr(det_nodes, "y")),
      z              = as.numeric(xml_attr(det_nodes, "z"))
    ) %>%
      arrange(t) %>%
      mutate(time_s = t * frame_interval_s)
  })

  list(
    tracks           = Filter(Negate(is.null), tracks_raw),
    frame_interval_s = frame_interval_s,
    n_raw            = length(particle_nodes)
  )
}

# 4c. Mean turning angle (unsigned degrees)
compute_turning_angle <- function(track_df) {
  n <- nrow(track_df)
  if (n < 3) return(NA_real_)
  dx      <- diff(track_df$x)
  dy      <- diff(track_df$y)
  n_steps <- length(dx)
  if (n_steps < 2) return(NA_real_)

  angles_deg <- sapply(seq_len(n_steps - 1), function(i) {
    v1 <- c(dx[i],   dy[i])
    v2 <- c(dx[i+1], dy[i+1])
    mag1 <- sqrt(sum(v1^2))
    mag2 <- sqrt(sum(v2^2))
    if (mag1 == 0 || mag2 == 0) return(NA_real_)
    cos_theta <- max(-1, min(1, sum(v1 * v2) / (mag1 * mag2)))
    acos(cos_theta) * 180 / pi
  })
  mean(angles_deg, na.rm = TRUE)
}

# 4d. MSD — 30% truncated weighted regression only
compute_msd_30pct <- function(track_df, max_dur_file) {

  n <- nrow(track_df)
  empty <- list(alpha_30pct = NA_real_, D_30pct = NA_real_,
                n_lags_30pct = 0L, msd_r2_30pct = NA_real_,
                fit_quality_30pct = "insufficient spots")
  if (n < 2) return(empty)

  pairs <- combn(n, 2, simplify = FALSE)
  msd_pairs <- do.call(rbind, lapply(pairs, function(idx) {
    i <- idx[1]; j <- idx[2]
    data.frame(
      lag_s   = track_df$time_s[j] - track_df$time_s[i],
      sq_disp = (track_df$x[j] - track_df$x[i])^2 +
                (track_df$y[j] - track_df$y[i])^2
    )
  }))

  msd_curve <- msd_pairs %>%
    group_by(lag_s) %>%
    summarise(MSD = mean(sq_disp, na.rm = TRUE),
              n_pairs = n(), .groups = "drop") %>%
    arrange(lag_s) %>%
    filter(lag_s > 0, MSD > 0)

  # 30% truncation based on max track duration in this file
  cutoff_s     <- 0.30 * max_dur_file
  msd_curve_30 <- msd_curve %>% filter(lag_s <= cutoff_s)
  n_lags       <- nrow(msd_curve_30)

  if (n_lags < min_lags_30pct) {
    return(list(alpha_30pct = NA_real_, D_30pct = NA_real_,
                n_lags_30pct = n_lags,
                msd_r2_30pct = NA_real_,
                fit_quality_30pct = paste0("insufficient lags (", n_lags, ")")))
  }

  fit <- tryCatch(
    lm(log(MSD) ~ log(lag_s), data = msd_curve_30, weights = n_pairs),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(list(alpha_30pct = NA_real_, D_30pct = NA_real_,
                n_lags_30pct = n_lags, msd_r2_30pct = NA_real_,
                fit_quality_30pct = "fit failed"))
  }

  list(
    alpha_30pct       = coef(fit)[["log(lag_s)"]],
    D_30pct           = exp(coef(fit)[["(Intercept)"]]) / 4,
    n_lags_30pct      = n_lags,
    msd_r2_30pct      = summary(fit)$r.squared,
    fit_quality_30pct = case_when(
      summary(fit)$r.squared >= 0.7 ~ "good",
      summary(fit)$r.squared >= 0.5 ~ "moderate",
      TRUE                          ~ "poor"
    )
  )
}


# --- 5. BATCH PROCESSING LOOP -------------------------------------------------

cat("=== STAGE 1: BATCH PROCESSING ===\n\n")

all_results <- list()
log_entries <- list()


# To reprocess specific ROIs only, list their ROIID_global values here.
# Leave empty c() to process all available ROIs.

reprocess_rois <- c("M1_T1_D2_MOCK_S3_ROI_1","M1_T1_D2_MOCK_S3_ROI_3","M1_T1_D2_MOCK_S4_ROI_2","M1_T1_D2_RP_S4_ROI_2","M1_T1_D2_RP_S4_ROI_3")

# Example: reprocess_rois <- c("M3_T2_RP_S2_ROI3")

sink("batch_processing_log_2.txt")
cat("BATCH PROCESSING LOG\n====================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("ROIs to process:", nrow(to_process), "\n")
cat("ROIs missing XML:", nrow(missing_xml), "\n\n")
sink()

to_process_run <- if (length(reprocess_rois) > 0) {
  to_process %>% filter(ROIID_global %in% reprocess_rois)
} else {
  to_process
}

for (row_idx in seq_len(nrow(to_process_run)))  {

  roi_meta <- to_process_run[row_idx, ]
  roi_id   <- roi_meta$ROIID_global

  cat(sprintf("[%d/%d] %s\n", row_idx, nrow(to_process), roi_id))

  log_entry <- list(
    ROIID_global         = roi_id,
    status               = "processing",
    n_csv_rows           = NA, n_xml_particles = NA,
    count_mismatch       = FALSE, n_angle_computed = NA,
    n_alpha_computed     = NA, n_removed_duration = NA,
    error_message        = NA_character_
  )

  tryCatch({

    # --- Load merged_tracks CSV ---
    csv_data <- load_merged_tracks(roi_meta$csv_path, roi_meta)
    if (is.null(csv_data)) stop("Failed to load merged_tracks CSV")
    log_entry$n_csv_rows <- nrow(csv_data)
    cat(sprintf("  CSV rows: %d\n", nrow(csv_data)))

    # --- Parse XML ---
    parsed <- parse_xml_tracks(roi_meta$xml_path)
    log_entry$n_xml_particles <- length(parsed$tracks)
    cat(sprintf("  XML particles: %d\n", length(parsed$tracks)))

    # --- Flag mismatch ---
    if (length(parsed$tracks) != nrow(csv_data)) {
      cat(sprintf("  WARNING: count mismatch — XML %d vs CSV %d\n",
                  length(parsed$tracks), nrow(csv_data)))
      log_entry$count_mismatch <- TRUE
    }

    # --- Duration filter: keep tracks > 10% of max TrackDuration in this ROI ---
    # Threshold computed on all raw tracks before any subsetting.
    # Mirrors the filter applied upstream to _Allstats.csv files in script 01.
    dur_threshold    <- 0.10 * max(csv_data$TrackDuration, na.rm = TRUE)
    keep_idx         <- which(csv_data$TrackDuration >= dur_threshold)
    n_removed_dur    <- nrow(csv_data) - length(keep_idx)
    log_entry$n_removed_duration <- n_removed_dur
    
    cat(sprintf("  Duration threshold: %.1f s | Removed: %d | Kept: %d\n",
                dur_threshold, n_removed_dur, length(keep_idx)))
    
    csv_data_filtered <- csv_data[keep_idx, ]
    
    # Positional XML alignment after duration filter.
    # keep_idx holds original integer row positions from the unfiltered CSV.
    # These are used to select the corresponding XML particles, preserving the
    # positional CSV<->XML correspondence that the batch loop relies on.
    # Capped at the number of XML particles to handle any pre-existing mismatch.
    valid_xml_idx <- keep_idx[keep_idx <= length(parsed$tracks)]
    n_matched     <- min(length(valid_xml_idx), nrow(csv_data_filtered))
    tracks_use    <- parsed$tracks[valid_xml_idx[seq_len(n_matched)]]
    csv_use       <- csv_data_filtered[seq_len(n_matched), ]
    
    # Max duration for 30% MSD cutoff — computed on duration-filtered tracks.
    # Since the filter removes only short tracks, max(pre-filter) == max(post-filter)
    # numerically, but using the filtered set is correct by construction.
    max_dur_file  <- max(csv_use$TrackDuration, na.rm = TRUE)
    
    # --- Compute turning angles ---
    turning_angles <- sapply(tracks_use, compute_turning_angle)
    log_entry$n_angle_computed <- sum(!is.na(turning_angles))

    # --- Compute MSD ---
    msd_results <- lapply(tracks_use, function(tr) {
      compute_msd_30pct(tr, max_dur_file = max_dur_file)
    })
    log_entry$n_alpha_computed <- sum(
      !is.na(sapply(msd_results, function(x) x$alpha_30pct))
    )

    cat(sprintf("  Turning angle: %d | Alpha 30pct: %d\n",
                log_entry$n_angle_computed, log_entry$n_alpha_computed))

    # --- Join metrics to CSV positionally ---
    metrics_df <- tibble(
      track_index        = seq_len(n_matched),
      mean_turning_angle = turning_angles,
      alpha_30pct        = sapply(msd_results, function(x) x$alpha_30pct),
      D_30pct_um2_s      = sapply(msd_results, function(x) x$D_30pct),
      n_lags_30pct       = sapply(msd_results, function(x) x$n_lags_30pct),
      msd_r2_30pct       = sapply(msd_results, function(x) x$msd_r2_30pct),
      fit_quality_30pct  = sapply(msd_results,
                                  function(x) x$fit_quality_30pct)
    )

    roi_joined <- csv_use %>%
      bind_cols(metrics_df %>% select(-track_index))

    # --- Save per-ROI CSV ---
    roi_csv_out <- file.path(
      dirname(roi_meta$csv_path),
      paste0(roi_meta$xml_prefix, "_turning_msd.csv")
    )
    write_csv(roi_joined, roi_csv_out)
    cat(sprintf("  Saved: %s\n", basename(roi_csv_out)))

    all_results[[roi_id]] <- roi_joined
    log_entry$status      <- "success"

  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
    log_entry$status        <<- "error"
    log_entry$error_message <<- e$message
  })

  log_entries[[row_idx]] <- log_entry

  sink("batch_processing_log.txt", append = TRUE)
  cat(sprintf("\n[%d/%d] %s — %s\n",
              row_idx, nrow(to_process), roi_id, log_entry$status))
  if (!is.na(log_entry$error_message)) {
    cat("  ERROR:", log_entry$error_message, "\n")
  } else {
    cat(sprintf("  CSV: %d | XML: %d | Mismatch: %s\n",
                log_entry$n_csv_rows, log_entry$n_xml_particles,
                log_entry$count_mismatch))
    cat(sprintf("  Turning angle: %d | Alpha 30pct: %d\n",
                log_entry$n_angle_computed, log_entry$n_alpha_computed))
  }
  sink()

  cat("\n")
}


# --- 6. STAGE 2: COMBINE ALL ROIs ---------------------------------------------
# Rebuild all_results from saved per-ROI turning_msd CSVs
# Use this if all_results was overwritten by a single-ROI reprocess run

merge_metadata_only <- FALSE  # TRUE  = only ROIs present in metadata file
# FALSE = all _turning_msd.csv files found in data_folder

cat("Rebuilding all_results from saved per-ROI CSVs...\n")

all_roi_csvs <- list.files(data_folder,
                           pattern   = "_turning_msd\\.csv$",
                           full.names = TRUE,
                           recursive  = TRUE)

# --- Deduplicate: same ROI may have two files (short XML-style vs long
#     Position-style name). Keep the shorter filename for each ROI key.
all_roi_csvs <- tibble(path = all_roi_csvs) %>%
  mutate(
    roi_key = basename(path) %>%
      str_remove("\\s*-\\s*Position.*$") %>%   # strip long Position block
      str_remove("_turning_msd\\.csv$")         # strip suffix → bare ROI key
  ) %>%
  group_by(roi_key) %>%
  slice_min(nchar(basename(path)), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  pull(path)

cat("Turning MSD CSVs after deduplication:", length(all_roi_csvs), "\n")

cat("Per-ROI CSVs found:", length(all_roi_csvs), "\n")

all_results <- lapply(all_roi_csvs, function(f) {
  tryCatch(
    read_csv(f, show_col_types = FALSE),
    error = function(e) {
      cat("  Failed to load:", basename(f), "\n")
      NULL
    }
  )
})

all_results <- Filter(Negate(is.null), all_results)
names(all_results) <- sapply(all_results, function(x) x$ROIID_global[1])
cat("ROIs successfully loaded:", length(all_results), "\n\n")

if (merge_metadata_only) {
  before_filter  <- length(all_results)
  all_results    <- all_results[names(all_results) %in% to_process$ROIID_global]
  cat(sprintf("Metadata-only mode: %d ROIs matched to metadata | %d excluded\n\n",
              length(all_results), before_filter - length(all_results)))
} else {
  cat("Full mode: merging all found ROIs regardless of metadata.\n\n")
}

cat("=== STAGE 2: COMBINING ALL ROIs ===\n")

if (length(all_results) == 0) {
  stop("No results to combine — check batch_processing_log.txt")
}

combined_raw <- bind_rows(all_results)

# Speed conversion and log transforms applied once on full combined dataset
combined_raw <- combined_raw %>%
  mutate(
    # Speed conversion: µm/s -> µm/min
    SpeedMean_ummin = Track_Speed_Mean * 60,
    
    # Log transforms (1e-6 offset guards against exact zeros)
    log_Speed          = log(SpeedMean_ummin    + 1e-6),
    log_TrackLength    = log(Track_Length        + 1e-6),
    log_MaxDis         = log(MaxDis             + 1e-6),
    log_Straightness   = log(Track_Straightness + 1e-6),
    logit_Straightness = log((Track_Straightness + 1e-6) /
                               (1 - Track_Straightness + 1e-6))
  )

write_csv(combined_raw, "turning_msd_combined_raw_v3_alldata.csv")
cat("Saved: turning_msd_combined_raw_v3_alldata.csv\n")
cat("Total tracks:", nrow(combined_raw), "\n")
cat("ROIs:        ", n_distinct(combined_raw$ROIID_global), "\n\n")


# --- 7. STAGE 3: FILTER STEP 1 — RETAIN COLUMNS OF INTEREST ------------------

# Selects columns of interest from the combined raw CSV.
# Log transforms and speed conversion were already applied in Stage 2,
# so all derived columns are available here for selection.
#
# Columns retained:
#   Metadata        : MouseID, TumourID_global, SliceID_global, ROIID_global,
#                     Location, Treatment, Day, Filepath, Filename
#   Track features  : ID, Track_Speed_Mean (µm/s), SpeedMean_ummin (µm/min),
#                     Track_Length, MaxDis, Track_Straightness,
#                     TrackDuration, t0
#   Log transforms  : log_Speed, log_TrackLength, log_MaxDis,
#                     log_Straightness, logit_Straightness
#   Computed metrics: mean_turning_angle, alpha_30pct, msd_r2_30pct

cat("=== STAGE 3: FILTER STEP 1 — COLUMN SELECTION ===\n")

keep_cols_stage3 <- c(
  # Metadata
  "MouseID", "TumourID_global", "SliceID_global", "ROIID_global",
  "Location", "Treatment", "Day",
  "Filepath", "Filename",
  
  # Track features (raw)
  "ID",
  "Track_Speed_Mean",     # raw speed µm/s (original from CSV)
  "SpeedMean_ummin",      # converted speed µm/min
  "Track_Length",
  "MaxDis",
  "Track_Straightness",
  "TrackDuration",
  "t0",
  
  # Log/logit transforms
  "log_Speed",
  "log_TrackLength",
  "log_MaxDis",
  "log_Straightness",
  "logit_Straightness",
  
  # Computed metrics from XML
  "mean_turning_angle",
  "alpha_30pct",
  "msd_r2_30pct",
  "fit_quality_30pct"
)

# Check which requested columns are actually present
available_cols <- intersect(keep_cols_stage3, names(combined_raw))
missing_keep   <- setdiff(keep_cols_stage3, names(combined_raw))

if (length(missing_keep) > 0) {
  cat("WARNING: following requested columns not found and will be skipped:\n")
  walk(missing_keep, ~ cat("  ", .x, "\n"))
}

filtered_stage3 <- combined_raw %>%
  select(all_of(available_cols))

write_csv(filtered_stage3, "turning_msd_filtered_v3_alldata.csv")
cat("Saved: turning_msd_filtered_v3_alldata.csv\n")
cat("Columns retained:", ncol(filtered_stage3), "\n")
cat("Tracks:          ", nrow(filtered_stage3), "\n\n")

# --- 8. STAGE 4: FILTER STEP 2 — MATCH TO tracks_clean.csv ------------------

# OLD 
#cat("=== STAGE 4: FILTER STEP 2 — MATCH TO tracks_clean.csv ===\n")
# 
# cat("Loading:", tracks_clean_file, "\n")
# tracks_clean <- read_csv(tracks_clean_file, show_col_types = FALSE)
# cat("tracks_clean rows:", nrow(tracks_clean), "\n")
# 
# # Rename tracks_clean key columns to match combined CSV naming
# # so the join uses consistent column names on both sides
# tracks_clean_keys <- tracks_clean %>%
#   select(ROIID_global, ID) %>%
#   mutate(in_tracks_clean = TRUE)
# 
# final <- filtered_stage3 %>%
#   left_join(tracks_clean_keys,
#             by = c("ROIID_global", "ID")) %>%
#   filter(in_tracks_clean == TRUE) %>%
#   select(-in_tracks_clean)
# 
# 
# cat("Tracks before match:", nrow(filtered_stage3), "\n")
# cat("Tracks after match: ", nrow(final), "\n")
# cat("Tracks removed:     ", nrow(filtered_stage3) - nrow(final),
#     "(not in tracks_clean — duration filter, flags etc.)\n\n")
# 
# # Warn if substantial drop
# pct_removed <- round(
#   (nrow(filtered_stage3) - nrow(final)) / nrow(filtered_stage3) * 100, 1
# )
# if (pct_removed > 20) {
#   cat("WARNING:", pct_removed,
#       "% of tracks removed in Stage 4 — check match key alignment\n")
#   cat("Verify that Track_Speed_Mean, Track_Length, MaxDis, TrackDuration\n")
#   cat("have identical values in both files (no rounding differences)\n\n")
# }
# 
# write_csv(final, "turning_msd_final_2.csv")
# cat("Saved: turning_msd_final.csv\n")
# cat("Columns:", ncol(final), "\n")
# cat("Tracks: ", nrow(final), "\n\n")

# --- 8. STAGE 4: QUALITY CONTROL CHECKS --------------------------------------
# Flags tracks with implausible values for speed, length, or straightness.
# Duration filtering was already applied per-ROI in the batch loop (Section 5).
# Mirrors Section 7 of script 01.

cat("=== STAGE 4: QUALITY CONTROL CHECKS ===\n\n")

filtered_stage3 <- filtered_stage3 %>%
  mutate(
    flag_speed        = SpeedMean_ummin <= 0 | SpeedMean_ummin > speed_max_ummin,
    flag_length       = Track_Length <= 0,
    flag_straightness = Track_Straightness < 0 | Track_Straightness > 1
  )

# --- Distribution summaries ---
cat("Tracks per Location:\n");  print(table(filtered_stage3$Location,        useNA = "ifany"))
cat("Tracks per Treatment:\n"); print(table(filtered_stage3$Treatment,       useNA = "ifany"))
cat("Tracks per Day:\n");       print(table(filtered_stage3$Day,             useNA = "ifany"))
cat("Tracks per TumourID:\n");  print(table(filtered_stage3$TumourID_global, useNA = "ifany"))

# --- NA audit on key columns ---
key_cols_qc <- c("SpeedMean_ummin", "Track_Length", "Track_Straightness",
                 "MaxDis", "TrackDuration")
cat("\nNA counts in key columns:\n")
print(filtered_stage3 %>% summarise(across(all_of(key_cols_qc), ~ sum(is.na(.)))))

# --- Flag counts ---
cat("\nFlagged tracks:\n")
cat("  Speed out of range (0 –", speed_max_ummin, "µm/min):",
    sum(filtered_stage3$flag_speed,        na.rm = TRUE), "\n")
cat("  Length <= 0:                         ",
    sum(filtered_stage3$flag_length,       na.rm = TRUE), "\n")
cat("  Straightness out of [0, 1]:          ",
    sum(filtered_stage3$flag_straightness, na.rm = TRUE), "\n")

# --- Save flagged rows for inspection ---
flagged_rows <- filtered_stage3 %>%
  filter(flag_speed | flag_length | flag_straightness) %>%
  select(ROIID_global, ID, SpeedMean_ummin, Track_Length, Track_Straightness,
         flag_speed, flag_length, flag_straightness)

if (nrow(flagged_rows) > 0) {
  cat("\nFlagged rows for inspection:\n")
  print(flagged_rows, n = Inf)
  write_csv(flagged_rows, "flagged_tracks.csv")
  cat("\nFlagged rows saved to: flagged_tracks.csv\n")
  cat("Review, then set remove_flagged <- TRUE in Section 1 and re-run.\n")
} else {
  cat("\nNo flagged rows found.\n")
}

# --- Removal step ---
# Controlled by remove_flagged in Section 1.
# Leave FALSE on first run; set TRUE after reviewing flagged_tracks.csv.
if (remove_flagged) {
  n_before        <- nrow(filtered_stage3)
  filtered_stage3 <- filtered_stage3 %>%
    filter(!flag_speed, !flag_length, !flag_straightness)
  cat("\nRemoved", n_before - nrow(filtered_stage3), "flagged tracks.",
      nrow(filtered_stage3), "tracks remaining.\n")
} else {
  cat("\nFlagged tracks retained (remove_flagged = FALSE in Section 1).\n")
}

# Drop flag columns before passing downstream
filtered_stage3 <- filtered_stage3 %>% select(-starts_with("flag_"))

# Pass to Section 8b (variable name kept consistent with rest of script)
final <- filtered_stage3

# =============================================================================
# ---- SECTION 8b — COMPUTE inv_outreach_ratio AND SAVE ----
# =============================================================================
# Adds three columns to final and saves to output_clean:
#   inv_outreach_ratio : Track_Length / MaxDis — tortuosity relative to
#                        furthest excursion. 1 = straight line to max point,
#                        >1 = meandering. Differs from Straightness because
#                        MaxDis >= net displacement.
#   Translocating      : TRUE if MaxDis >= maxdis_threshold (Section 1)
#   log_TrackDuration  : natural log of TrackDuration (for PCA use)

tracks_clean <- final %>%
  mutate(
    inv_outreach_ratio = if_else(MaxDis > 0, Track_Length / MaxDis, NA_real_),
    Translocating      = MaxDis >= maxdis_threshold,
    log_TrackDuration  = log(TrackDuration),
    Location           = factor(Location,       levels = c("Core", "CoreEdge", "Edge")),
    Treatment          = factor(Treatment,       levels = c("Mock", "RP")),
    Day                = factor(Day,             levels = c("D1", "D2")),
    MouseID            = factor(MouseID),
    TumourID_global    = factor(TumourID_global),
    SliceID_global     = factor(SliceID_global),
    ROIID_global       = factor(ROIID_global)
  )

# --- Summary ---
cat("--- inv_outreach_ratio summary ---\n")
cat("  MaxDis = 0 (ratio NA):                        ",
    sum(final$MaxDis == 0, na.rm = TRUE), "\n")
cat("  Non-translocating (MaxDis <", maxdis_threshold, "µm):",
    sum(!tracks_clean$Translocating, na.rm = TRUE), "\n")
cat("  Ratio range, all tracks:          [",
    round(min(tracks_clean$inv_outreach_ratio, na.rm = TRUE), 3), ",",
    round(max(tracks_clean$inv_outreach_ratio, na.rm = TRUE), 3), "]\n")
cat("  Ratio range, translocating only:  [",
    round(min(tracks_clean$inv_outreach_ratio[tracks_clean$Translocating],
              na.rm = TRUE), 3), ",",
    round(max(tracks_clean$inv_outreach_ratio[tracks_clean$Translocating],
              na.rm = TRUE), 3), "]\n\n")

write_csv(tracks_clean, output_clean)                   #<------ Make sure you change the name of the output_clean in section 1
cat("Saved:", output_clean, "\n")
cat("Rows:", nrow(tracks_clean), "| Columns:", ncol(tracks_clean), "\n")

# ---- 9a. ROI-LEVEL SUMMARY -----------------------------------------------------
# One row per ROI.
# Follows the structure of script 01 Section 9 and extends it with the
# additional metrics computed in this pipeline.
#
# Alpha handling:
#   Raw mean/median/SD are not computed because negative alpha values can arise
#   from (a) poor log-log fit on short tracks or (b) genuine confined/oscillatory
#   motion. Mixing these without distinguishing source is misleading.
#   Instead:
#     - Filtered summaries use only tracks with msd_r2_30pct >= 0.5
#       (moderate or better fit; Qian et al. 1991, Biophys J 60:910-921)
#     - Fit quality is reported as both counts and proportions per ROI
#     - prop_negative_alpha uses only tracks with a computed alpha as denominator
#       (tracks where alpha is NA due to insufficient lags are excluded)

roi_summary <- tracks_clean %>%
  group_by(MouseID, TumourID_global, SliceID_global, ROIID_global,
           Location, Treatment, Day) %>%
  summarise(
    n_tracks = n(),
    
    # --- Speed (µm/min) ---
    Speed_mean   = mean(SpeedMean_ummin,   na.rm = TRUE),
    Speed_median = median(SpeedMean_ummin, na.rm = TRUE),
    Speed_sd     = sd(SpeedMean_ummin,     na.rm = TRUE),
    
    # --- Track Length (µm) ---
    Length_mean   = mean(Track_Length,   na.rm = TRUE),
    Length_median = median(Track_Length, na.rm = TRUE),
    Length_sd     = sd(Track_Length,     na.rm = TRUE),
    
    # --- Straightness (0–1) ---
    Straightness_mean   = mean(Track_Straightness,   na.rm = TRUE),
    Straightness_median = median(Track_Straightness, na.rm = TRUE),
    Straightness_sd     = sd(Track_Straightness,     na.rm = TRUE),
    
    # --- MaxDis (µm) ---
    MaxDis_mean   = mean(MaxDis,   na.rm = TRUE),
    MaxDis_median = median(MaxDis, na.rm = TRUE),
    MaxDis_sd     = sd(MaxDis,     na.rm = TRUE),
    
    # --- Inv outreach ratio (Track_Length / MaxDis) ---
    # 1 = straight line to furthest point; >1 = meandering path
    InvOutreach_mean   = mean(inv_outreach_ratio,   na.rm = TRUE),
    InvOutreach_median = median(inv_outreach_ratio, na.rm = TRUE),
    InvOutreach_sd     = sd(inv_outreach_ratio,     na.rm = TRUE),
    
    # --- Turning angle (degrees, unsigned mean per track) ---
    TurningAngle_mean   = mean(mean_turning_angle,   na.rm = TRUE),
    TurningAngle_median = median(mean_turning_angle, na.rm = TRUE),
    TurningAngle_sd     = sd(mean_turning_angle,     na.rm = TRUE),
    
    # --- MSD anomalous diffusion exponent α — filtered summaries only ---
    # Restricted to tracks with msd_r2_30pct >= 0.5 (moderate or better).
    # Reference: Qian et al. (1991) Biophys J 60:910-921;
    #            Saxton & Jacobson (1997) Annu Rev Biophys Biomol Struct 26:373-399
    Alpha30_filt_n      = sum(!is.na(alpha_30pct) & msd_r2_30pct >= 0.5),
    Alpha30_filt_mean   = mean(  alpha_30pct[!is.na(alpha_30pct) & msd_r2_30pct >= 0.5], na.rm = TRUE),
    Alpha30_filt_median = median(alpha_30pct[!is.na(alpha_30pct) & msd_r2_30pct >= 0.5], na.rm = TRUE),
    Alpha30_filt_sd     = sd(    alpha_30pct[!is.na(alpha_30pct) & msd_r2_30pct >= 0.5], na.rm = TRUE),
    
    # --- Proportion of negative alpha (among tracks with a computed alpha) ---
    # Denominator excludes NA alpha (insufficient lags) — not the same as a
    # negative fit. High prop_negative_alpha may indicate confined/oscillatory ROI.
    n_alpha_computed    = sum(!is.na(alpha_30pct)),
    n_negative_alpha    = sum(alpha_30pct < 0, na.rm = TRUE),
    prop_negative_alpha = sum(alpha_30pct < 0, na.rm = TRUE) /
      sum(!is.na(alpha_30pct)),
    
    # --- Fit quality breakdown (counts and proportions) ---
    # Proportions are out of n_tracks (all tracks in ROI, not just fitted ones),
    # so they reflect overall data quality per ROI.
    n_fit_good        = sum(fit_quality_30pct == "good",     na.rm = TRUE),
    n_fit_moderate    = sum(fit_quality_30pct == "moderate", na.rm = TRUE),
    n_fit_poor        = sum(fit_quality_30pct == "poor",     na.rm = TRUE),
    n_fit_NA          = sum(is.na(alpha_30pct)),
    prop_fit_good     = sum(fit_quality_30pct == "good",     na.rm = TRUE) / n(),
    prop_fit_moderate = sum(fit_quality_30pct == "moderate", na.rm = TRUE) / n(),
    prop_fit_poor     = sum(fit_quality_30pct == "poor",     na.rm = TRUE) / n(),
    prop_fit_NA       = sum(is.na(alpha_30pct)) / n(),
    
    # --- Proportion of translocating tracks ---
    prop_Translocating = mean(Translocating, na.rm = TRUE),
    
    # --- Log/logit transform summaries ---
    # Mean and SD of track-level transformed values (not transform of averages).
    log_Speed_mean        = mean(log_Speed,          na.rm = TRUE),
    log_Speed_sd          = sd(log_Speed,            na.rm = TRUE),
    log_Length_mean       = mean(log_TrackLength,    na.rm = TRUE),
    log_Length_sd         = sd(log_TrackLength,      na.rm = TRUE),
    log_MaxDis_mean       = mean(log_MaxDis,         na.rm = TRUE),
    log_MaxDis_sd         = sd(log_MaxDis,           na.rm = TRUE),
    log_Straight_mean     = mean(log_Straightness,   na.rm = TRUE),
    log_Straight_sd       = sd(log_Straightness,     na.rm = TRUE),
    logit_Straight_mean   = mean(logit_Straightness, na.rm = TRUE),
    logit_Straight_sd     = sd(logit_Straightness,   na.rm = TRUE),
    
    .groups = "drop"
  )

cat("\n--- ROI summary ---\n")
cat("Total ROIs:", nrow(roi_summary), "\n")
cat("ROIs per Location:\n");  print(table(roi_summary$Location))
cat("ROIs per Treatment:\n"); print(table(roi_summary$Treatment))

write_csv(roi_summary, "roi_summary_v3_alldata.csv")
cat("Saved: roi_summary_v3_alldata.csv\n\n")

# ---- 9b. BATCH SUMMARY ---------------------------------------------------------

log_df <- bind_rows(lapply(log_entries, as_tibble))
write_csv(log_df, "batch_processing_log.csv")

n_success  <- sum(log_df$status == "success",      na.rm = TRUE)
n_error    <- sum(log_df$status == "error",         na.rm = TRUE)
n_mismatch <- sum(log_df$count_mismatch == TRUE,    na.rm = TRUE)

cat("=== BATCH SUMMARY ===\n")
cat("ROIs attempted:     ", nrow(to_process),  "\n")
cat("ROIs succeeded:     ", n_success,          "\n")
cat("ROIs errored:       ", n_error,            "\n")
cat("ROIs with mismatch: ", n_mismatch,         "\n")
cat("ROIs missing XML:   ", nrow(missing_xml),
    "(will process when XMLs added)\n\n")

if (n_error > 0) {
  cat("Failed ROIs:\n")
  log_df %>%
    filter(status == "error") %>%
    select(ROIID_global, error_message) %>%
    print()
}

cat("Saved: batch_processing_log.csv\n\n")

# ---- 9c. ROI-LEVEL SUMMARY — STANDALONE RE-RUN (commented out by default) ----
# Use this section to re-run the ROI summary on a filtered subset of
# tracks_clean_v2.csv after GMM components have been assigned by script 02d.
#
# Steps:
#   1. Set the parameters below
#   2. Uncomment the entire section (Ctrl+Shift+C in RStudio selects block)
#   3. Run
#   4. Re-comment when done

 # --- Parameters ---
 tracks_clean_path    <- "Z:\\lightmicro-inst\\scharalambous\\Resident_cell_motility_tr_v3\\ROI files\\Complete data - retrospective fitting G=5 VEE (from v2)\\tracks_clean_v3_alldata_GMM_retrospective.csv" # path to GMM-annotated file
 gmm_keep             <- c(1)                # components to keep; c() = no filter (full dataset)
 roi_summary_filtered_output <- "roi_summary_v3_GMM_1_MaxDis7.csv"  # output filename — edit to match gmm_keep

 # --- Load ---
 cat("=== SECTION 9b: ROI SUMMARY RE-RUN ===\n\n")
 tracks_rerun <- read_csv(tracks_clean_path, show_col_types = FALSE) %>%
   mutate(
     Location        = factor(Location,        levels = c("Core", "CoreEdge", "Edge")),
     Treatment       = factor(Treatment,        levels = c("Mock", "RP")),
     Day             = factor(Day,              levels = c("D1", "D2")),
     MouseID         = factor(MouseID),
     TumourID_global = factor(TumourID_global),
     SliceID_global  = factor(SliceID_global),
     ROIID_global    = factor(ROIID_global)
   )
 cat("Loaded:", nrow(tracks_rerun), "tracks from", tracks_clean_path, "\n")

 # --- GMM filter ---
 if (length(gmm_keep) > 0) {
   if (!"GMM_component" %in% names(tracks_rerun)) {
     stop("GMM_component column not found in ", tracks_clean_path,
          " — ensure script 02d has been run first.")
   }
   n_before     <- nrow(tracks_rerun)
   tracks_rerun <- tracks_rerun %>% filter(GMM_component %in% gmm_keep)
   cat(sprintf("GMM filter applied: kept components %s | %d -> %d tracks | %d removed\n\n",
               paste(gmm_keep, collapse = ", "),
               n_before, nrow(tracks_rerun), n_before - nrow(tracks_rerun)))
 } else {
   cat("No GMM filter applied — using full dataset.\n\n")
 }

 # --- MaxDis filter ---
 # Removes near-stationary tracks within retained GMM components.
 # Must match maxdis_threshold used in scripts 02, 03, 04, 02b.
 maxdis_threshold  <- 7  # µm — update if threshold changes
 n_before_maxdis   <- nrow(tracks_rerun)
 tracks_rerun_excl <- tracks_rerun %>% filter(MaxDis <= maxdis_threshold)
 tracks_rerun      <- tracks_rerun %>% filter(MaxDis >  maxdis_threshold)
 cat(sprintf("MaxDis filter: MaxDis > %.1f µm | %d -> %d tracks | %d removed\n\n",
             maxdis_threshold,
             n_before_maxdis, nrow(tracks_rerun),
             n_before_maxdis - nrow(tracks_rerun)))
 write_csv(
   tracks_rerun_excl %>%
     select(any_of(c(
       "ROIID_global", "SliceID_global", "TumourID_global", "MouseID",
       "Location", "Treatment", "Day", "GMM_component",
       "MaxDis", "log_MaxDis", "SpeedMean_ummin",
       "Track_Length", "Track_Straightness", "TrackDuration"
     ))),
   paste0("excluded_tracks_maxdis_",
          if (length(gmm_keep) == 0) "all" else
            paste0("GMM", paste(gmm_keep, collapse = "_")), ".csv")
 )

 # Complete-case filter
 # Removes tracks with NA in any core variable column.
 # Alpha is excluded from this filter — it runs on its own subset via filter_col.
 # The 20 excluded tracks have no turning angle due to complete immobility.
 core_cols <- c("SpeedMean_ummin", "log_Speed", "log_TrackLength",
                "Track_Straightness", "log_Straightness",
                "MaxDis", "log_MaxDis", "mean_turning_angle")

 # --- 1c. COMPLETE-CASE FILTER -------------------------------------------------
 n_before <- nrow(tracks_rerun)
 tracks   <- tracks_rerun %>% filter(if_all(any_of(core_cols), ~ !is.na(.)))
 cat("Complete-case filter: ", n_before - nrow(tracks_rerun),
     "tracks removed (NA in core columns)\n")
 cat("Tracks retained:", nrow(tracks_rerun), "\n\n")


 # Alpha bounds filter — remove physically impossible values
 # Applied before model fitting; Alpha must be in (0, 2]
 n_alpha_before <- sum(!is.na(tracks_rerun$alpha_30pct))
 tracks_rerun <- tracks_rerun %>%
   mutate(alpha_30pct = ifelse(alpha_30pct <= 0 | alpha_30pct > 2.25,
                               NA_real_, alpha_30pct))
 cat("Alpha bounds filter: ",
     n_alpha_before - sum(!is.na(tracks_rerun$alpha_30pct)),
     "values set to NA (outside 0-2 range)\n\n")



 # --- ROI summary (identical logic to Section 9) ---
 roi_summary_filtered <- tracks_rerun %>%
   group_by(MouseID, TumourID_global, SliceID_global, ROIID_global,
            Location, Treatment, Day) %>%
   summarise(
     n_tracks = n(),
     Speed_mean   = mean(SpeedMean_ummin,   na.rm = TRUE),
     Speed_median = median(SpeedMean_ummin, na.rm = TRUE),
     Speed_sd     = sd(SpeedMean_ummin,     na.rm = TRUE),
     Length_mean   = mean(Track_Length,   na.rm = TRUE),
     Length_median = median(Track_Length, na.rm = TRUE),
     Length_sd     = sd(Track_Length,     na.rm = TRUE),
     Straightness_mean   = mean(Track_Straightness,   na.rm = TRUE),
     Straightness_median = median(Track_Straightness, na.rm = TRUE),
     Straightness_sd     = sd(Track_Straightness,     na.rm = TRUE),
     MaxDis_mean   = mean(MaxDis,   na.rm = TRUE),
     MaxDis_median = median(MaxDis, na.rm = TRUE),
     MaxDis_sd     = sd(MaxDis,     na.rm = TRUE),
     InvOutreach_mean   = mean(inv_outreach_ratio,   na.rm = TRUE),
     InvOutreach_median = median(inv_outreach_ratio, na.rm = TRUE),
     InvOutreach_sd     = sd(inv_outreach_ratio,     na.rm = TRUE),
     TurningAngle_mean   = mean(mean_turning_angle,   na.rm = TRUE),
     TurningAngle_median = median(mean_turning_angle, na.rm = TRUE),
     TurningAngle_sd     = sd(mean_turning_angle,     na.rm = TRUE),
     Alpha30_filt_n      = sum(!is.na(alpha_30pct) & msd_r2_30pct >= 0.5),
     Alpha30_filt_mean   = mean(  alpha_30pct[!is.na(alpha_30pct) & msd_r2_30pct >= 0.5], na.rm = TRUE),
     Alpha30_filt_median = median(alpha_30pct[!is.na(alpha_30pct) & msd_r2_30pct >= 0.5], na.rm = TRUE),
     Alpha30_filt_sd     = sd(    alpha_30pct[!is.na(alpha_30pct) & msd_r2_30pct >= 0.5], na.rm = TRUE),
     n_alpha_computed    = sum(!is.na(alpha_30pct)),
     n_negative_alpha    = sum(alpha_30pct < 0, na.rm = TRUE),
     prop_negative_alpha = sum(alpha_30pct < 0, na.rm = TRUE) / sum(!is.na(alpha_30pct)),
     n_fit_good        = sum(fit_quality_30pct == "good",     na.rm = TRUE),
     n_fit_moderate    = sum(fit_quality_30pct == "moderate", na.rm = TRUE),
     n_fit_poor        = sum(fit_quality_30pct == "poor",     na.rm = TRUE),
     n_fit_NA          = sum(is.na(alpha_30pct)),
     prop_fit_good     = sum(fit_quality_30pct == "good",     na.rm = TRUE) / n(),
     prop_fit_moderate = sum(fit_quality_30pct == "moderate", na.rm = TRUE) / n(),
     prop_fit_poor     = sum(fit_quality_30pct == "poor",     na.rm = TRUE) / n(),
     prop_fit_NA       = sum(is.na(alpha_30pct)) / n(),
     prop_Translocating    = mean(Translocating,      na.rm = TRUE),
     log_Speed_mean        = mean(log_Speed,          na.rm = TRUE),
     log_Speed_sd          = sd(log_Speed,            na.rm = TRUE),
     log_Length_mean       = mean(log_TrackLength,    na.rm = TRUE),
     log_Length_sd         = sd(log_TrackLength,      na.rm = TRUE),
     log_MaxDis_mean       = mean(log_MaxDis,         na.rm = TRUE),
     log_MaxDis_sd         = sd(log_MaxDis,           na.rm = TRUE),
     log_Straight_mean     = mean(log_Straightness,   na.rm = TRUE),
     log_Straight_sd       = sd(log_Straightness,     na.rm = TRUE),
     logit_Straight_mean   = mean(logit_Straightness, na.rm = TRUE),
     logit_Straight_sd     = sd(logit_Straightness,   na.rm = TRUE),
     .groups = "drop"
   )

 cat("--- Filtered ROI summary ---\n")
 cat("Total ROIs:", nrow(roi_summary_filtered), "\n")
 cat("ROIs per Location:\n");  print(table(roi_summary_filtered$Location))
 cat("ROIs per Treatment:\n"); print(table(roi_summary_filtered$Treatment))

 write_csv(roi_summary_filtered, roi_summary_filtered_output)
 cat("Saved:", roi_summary_filtered_output, "\n\n")

# --- 10. OPTIONAL MSD CURVE PDFs ----------------------------------------------
# One PDF per ROI listed in msd_pdf_rois.
# Re-processes that ROI to recover msd_curve objects.

if (length(msd_pdf_rois) == 0) {
  cat("No MSD PDFs requested.\n")
  cat("Add ROI IDs to msd_pdf_rois at top of script to generate PDFs.\n")

} else {

  cat("=== GENERATING MSD CURVE PDFs ===\n")

  for (roi_id in msd_pdf_rois) {

    if (!roi_id %in% names(all_results)) {
      cat("Skipping", roi_id, "— not in processed results\n")
      next
    }

    cat("Generating PDF for:", roi_id, "\n")

    roi_meta <- to_process %>% filter(ROIID_global == roi_id)
    if (nrow(roi_meta) == 0) next
    roi_meta <- roi_meta[1, ]

    csv_data    <- load_merged_tracks(roi_meta$csv_path, roi_meta)
    parsed      <- parse_xml_tracks(roi_meta$xml_path)
    n_matched   <- min(length(parsed$tracks), nrow(csv_data))
    tracks_use  <- parsed$tracks[seq_len(n_matched)]
    csv_use     <- csv_data[seq_len(n_matched), ]
    max_dur_file <- max(csv_use$TrackDuration, na.rm = TRUE)
    cutoff_30   <- 0.30 * max_dur_file
    n_plot      <- min(msd_pdf_max_tracks, n_matched)

    # Recompute full MSD curves for plotting
    msd_full <- lapply(tracks_use, function(tr) {
      n <- nrow(tr)
      if (n < 2) return(NULL)
      pairs <- combn(n, 2, simplify = FALSE)
      msd_pairs <- do.call(rbind, lapply(pairs, function(idx) {
        i <- idx[1]; j <- idx[2]
        data.frame(
          lag_s   = tr$time_s[j] - tr$time_s[i],
          sq_disp = (tr$x[j] - tr$x[i])^2 + (tr$y[j] - tr$y[i])^2
        )
      }))
      msd_pairs %>%
        group_by(lag_s) %>%
        summarise(MSD = mean(sq_disp, na.rm = TRUE),
                  n_pairs = n(), .groups = "drop") %>%
        arrange(lag_s) %>%
        filter(lag_s > 0, MSD > 0)
    })

    msd_metrics <- lapply(tracks_use, function(tr) {
      compute_msd_30pct(tr, max_dur_file)
    })

    pdf_path <- paste0("msd_", roi_id, ".pdf")
    pdf(pdf_path, width = 10, height = 6)

    for (i in seq_len(n_plot)) {

      curve_df <- msd_full[[i]]
      if (is.null(curve_df) || nrow(curve_df) == 0) next

      alpha_30 <- msd_metrics[[i]]$alpha_30pct
      D_30     <- msd_metrics[[i]]$D_30pct
      r2_30    <- msd_metrics[[i]]$msd_r2_30pct
      quality  <- msd_metrics[[i]]$fit_quality_30pct
      n_spots  <- nrow(tracks_use[[i]])
      track_id <- if (!is.null(csv_data$ID[i])) csv_data$ID[i] else i

      quality_colour <- case_when(
        quality == "good"     ~ "#2A9D8F",
        quality == "moderate" ~ "#F4A261",
        TRUE                  ~ "#E63946"
      )

      subtitle_str <- paste0(
        "30%: α = ", ifelse(is.na(alpha_30), "NA", round(alpha_30, 3)),
        "  D = ",    ifelse(is.na(D_30),     "NA", round(D_30,     4)),
        " µm²/s  R² = ", ifelse(is.na(r2_30), "NA", round(r2_30, 3)),
        "  [", quality, "]  | Spots: ", n_spots
      )

      if (!is.na(alpha_30)) {
        tau_max_30   <- min(max(curve_df$lag_s), cutoff_30)
        tau_seq_30   <- seq(min(curve_df$lag_s), tau_max_30, length.out = 100)
        fitted_30_df <- tibble(lag_s = tau_seq_30,
                               MSD   = 4 * D_30 * tau_seq_30^alpha_30)
      }

      p <- ggplot(curve_df, aes(x = lag_s, y = MSD)) +
        geom_vline(xintercept = cutoff_30,
                   linetype = "dotted", colour = "grey70", linewidth = 0.6) +
        geom_point(size = 2.5, colour = quality_colour) +
        geom_line(colour = "grey60", linewidth = 0.5) +
        { if (!is.na(alpha_30))
            geom_line(data = fitted_30_df, aes(x = lag_s, y = MSD),
                      colour = "#E63946", linewidth = 0.9, linetype = "dotted")
        } +
        annotate("text", x = cutoff_30, y = Inf,
                 label = "30% cutoff", hjust = -0.1, vjust = 1.5,
                 size = 3, colour = "grey50") +
        scale_x_log10() +
        scale_y_log10() +
        annotation_logticks(sides = "bl", linewidth = 0.3) +
        labs(
          title    = paste0("Track ID: ", track_id,
                            "  (index ", i, ")  —  ", roi_id),
          subtitle = subtitle_str,
          x = "Time lag τ (seconds)", y = "MSD (µm²)"
        ) +
        theme_bw(base_size = 11) +
        theme(
          plot.title       = element_text(face = "bold", size = 11),
          plot.subtitle    = element_text(size = 8.5, colour = "grey30"),
          panel.grid.minor = element_blank()
        )

      print(p)
    }

    dev.off()
    cat("  Saved:", pdf_path, "\n")
  }
}

cat("\n=== COMPLETE ===\n")
cat("Key outputs:\n")
cat("  turning_msd_combined_raw.csv : all ROIs, all columns\n")
cat("  turning_msd_filtered.csv     : columns of interest\n")
cat("  turning_msd_final.csv        : matched to tracks_clean.csv\n")
cat("  batch_processing_log.csv     : per-ROI processing summary\n")
cat("\nNext step: join turning_msd_final.csv to tracks_clean.csv\n")
cat("Join keys: ROIID_global + Track_Speed_Mean + Track_Length +",
    "MaxDis + TrackDuration\n")
