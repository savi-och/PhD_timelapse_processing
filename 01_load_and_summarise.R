# =============================================================================
# CELL TRACKING ANALYSIS - STEP 1: DATA LOADING & SUMMARISATION
# =============================================================================
# Hierarchy: Track -> ROI -> Slice -> Tumour -> Mouse
# Fixed effects (later): (Size), ROI Location, Treatment, Day
# Random effects (later): (1|Mouse/Tumour/Slice/ROI)
#
# Filename format (on disk):
#   M2_T1_D1_Mock_S1_ROI_1 - Position 1 [25] 3DMontage...._Allstats.xlsx
#
# Metadata CSV columns:
#   Filename | MouseID | TumourID | Day | Treatment | SliceID | ROIID | Location | Comments
#   (Comments column is ignored on load)
# =============================================================================
# ----- clear environment --------
rm(list=ls())
gc()


# --- 0. PACKAGES --------------------------------------------------------------
# install.packages(c("tidyverse", "readxl", "lme4", "lmerTest",
#                    "emmeans", "performance", "ggrain", "patchwork"))

install.packages(c("lme4","lmerTest","readxl","emmeans"))
library(tidyverse)
library(readxl)


# --- 1. CONFIGURATION ---------------------------------------------------------

# Folder containing all xlsx files (searches subfolders automatically)
setwd("Z:\\lightmicro-inst\\scharalambous\\Resident_Cell_Cum_Tracking_analysis\\ROI files")
data_folder <- "Z:\\lightmicro-inst\\scharalambous\\Resident_Cell_Cum_Tracking_analysis\\ROI files"   # <-- CHANGE THIS to your actual path

# Your metadata CSV file
metadata_file <- "experimental_file_list.csv" # <-- CHANGE THIS if needed


# --- 2. LOAD METADATA CSV -----------------------------------------------------
# Reads your CSV and drops the Comments column.
# Filename column must exactly match the filename on disk (including the long suffix).

metadata <- read_csv(metadata_file, show_col_types = FALSE) %>%
  mutate(Filename = paste0(tools::file_path_sans_ext(Filename), ".csv")) %>%
  select(Filename, MouseID, TumourID, Day,
         Treatment, SliceID, ROIID, Location) %>%
  mutate(
    # Ordered factors - reference categories: Core, Mock, D1
    Location  = factor(Location,  levels = c("Core", "CoreEdge", "Edge")),
    Treatment = factor(Treatment, levels = c("Mock", "RP")),
    Day       = factor(Day,       levels = c("D1", "D2")),
    MouseID   = factor(MouseID),
    TumourID  = factor(TumourID),
    # Globally unique IDs to prevent collisions across mice/tumours
    TumourID_global = factor(paste(MouseID, TumourID, sep = "_")),
    SliceID_global  = paste(MouseID, TumourID, Treatment, SliceID, sep = "_"),
    ROIID_global    = paste(MouseID, TumourID, Treatment, SliceID, ROIID, sep = "_")
  )

cat("Metadata loaded:", nrow(metadata), "ROIs\n")
cat("Mice:   ", levels(metadata$MouseID),        "\n")
cat("Tumours:", levels(metadata$TumourID_global), "\n")
cat("Locations:\n"); print(table(metadata$Location))
cat("Treatment:\n");print(table(metadata$Treatment))


# --- 3. LOAD ALL CSV FILES ---------------------------------------------------

all_files <- list.files(data_folder, pattern = "\\.csv$",
                        full.names = TRUE, recursive = TRUE)

cat("\ncsv files found on disk:", length(all_files), "\n")

if (length(all_files) == 0) {
  stop(paste("No xlsx files found in:", data_folder,
             "\nCheck your data_folder path."))
}

load_roi_file <- function(filepath) {
  filepath <- normalizePath(filepath, winslash = "/", mustWork = FALSE)
   tryCatch({
    read_csv(filepath) %>%
      rename_with(~ str_replace_all(., "[. ]", "_")) %>%
      mutate(
        Filepath        = filepath,
        Filename        = basename(filepath)  # used to join with metadata
      )
  }, error = function(e) {
    warning(paste("Failed to load:", filepath, "\n  Error:", e$message))
    return(NULL)
  })
}

raw_list <- map(all_files, load_roi_file)

# Report any files that failed to load
failed <- all_files[map_lgl(raw_list, is.null)]
if (length(failed) > 0) {
  cat("\nWARNING - Failed to load", length(failed), "file(s):\n")
  walk(failed, ~ cat("  ", .x, "\n"))
}

raw_data <- bind_rows(compact(raw_list))
cat("Raw tracks loaded:", nrow(raw_data), "tracks from",
    n_distinct(raw_data$Filename), "files\n")


# --- 4. CHECK FILES MATCH METADATA --------------------------------------------
# Every file on disk should have a matching row in metadata, and vice versa.

files_on_disk    <- basename(all_files)
files_in_meta    <- metadata$Filename

# Files on disk but missing from metadata CSV
not_in_meta <- setdiff(files_on_disk, files_in_meta)
if (length(not_in_meta) > 0) {
  cat("\nWARNING -", length(not_in_meta),
      "file(s) on disk have no metadata entry:\n")
  walk(not_in_meta, ~ cat("  ", .x, "\n"))
}

# Rows in metadata CSV but no matching file on disk
not_on_disk <- setdiff(files_in_meta, files_on_disk)
if (length(not_on_disk) > 0) {
  cat("\nWARNING -", length(not_on_disk),
      "metadata row(s) have no matching file on disk:\n")
  walk(not_on_disk, ~ cat("  ", .x, "\n"))
}

if (length(not_in_meta) == 0 & length(not_on_disk) == 0) {
  cat("File check passed: all files match metadata\n")
}


# --- 5. JOIN METADATA TO TRACKS -----------------------------------------------
# Drop any metadata columns that already exist in raw_data (from csv headers)
# to prevent .x / .y duplicate column suffixes

meta_cols <- setdiff(names(metadata), "Filename")  # all metadata cols except join key

tracks <- raw_data %>%
  select(-any_of(meta_cols)) %>%          # drop duplicates from raw data first
  left_join(metadata, by = "Filename")

# Check for any tracks with missing metadata after join
missing_meta <- tracks %>% filter(is.na(Location)) %>% distinct(Filename)
if (nrow(missing_meta) > 0) {
  cat("\nWARNING -", nrow(missing_meta),
      "file(s) did not join to metadata:\n")
  print(missing_meta$Filename)
}

cat("Join complete:", nrow(tracks), "tracks with full metadata\n")

# --- 6. RENAME & CONVERT KEY COLUMNS ------------------------------------------
# If unsure of your column names after loading, run: names(tracks)
# Adjust left-hand side names below to match your xlsx headers exactly.

tracks <- tracks %>%
  rename(
    TrackLength   = Track_Length,        # µm
    SpeedMean_ums = Track_Speed_Mean,    # µm/s (original)
    Straightness  = Track_Straightness,  # 0-1
    MaxDis        = MaxDis,              # µm
    TrackDuration = Track_Duration       # seconds
  ) %>%
  mutate(
    SpeedMean_ummin = SpeedMean_ums * 60  # convert to µm/min
  )

# Sanity check on conversion
cat("\nSpeed conversion check (first 5 rows):\n")
print(tracks %>% select(Filename, SpeedMean_ums, SpeedMean_ummin) %>% head(5))


# --- 7. QUALITY CHECKS --------------------------------------------------------

cat("\n--- Quality checks ---\n")
cat("Tracks per Location:\n");  print(table(tracks$Location,        useNA = "ifany"))
cat("Tracks per Treatment:\n"); print(table(tracks$Treatment,       useNA = "ifany"))
cat("Tracks per Day:\n");       print(table(tracks$Day,             useNA = "ifany"))
cat("Tracks per Tumour:\n");    print(table(tracks$TumourID_global, useNA = "ifany"))

# NA check on key columns
key_cols <- c("SpeedMean_ummin", "TrackLength", "Straightness",
              "MaxDis", "TrackDuration")
cat("\nNA counts in key columns:\n")
print(tracks %>% summarise(across(all_of(key_cols), ~ sum(is.na(.)))))

# Flag implausible values - adjust thresholds to your cell type
tracks <- tracks %>%
  mutate(
    flag_speed        = SpeedMean_ummin <= 0 | SpeedMean_ummin > 3000,
    flag_length       = TrackLength <= 0,
    flag_straightness = Straightness < 0 | Straightness > 1
  )

cat("\nFlagged tracks:\n")
cat("  Speed out of range (0-3000 µm/min):", sum(tracks$flag_speed,        na.rm = TRUE), "\n")
cat("  Length <= 0:                        ", sum(tracks$flag_length,       na.rm = TRUE), "\n")
cat("  Straightness out of 0-1:            ", sum(tracks$flag_straightness, na.rm = TRUE), "\n")

# Print flagged rows so you can inspect them before removal
flagged_rows <- tracks %>%
  filter(flag_speed | flag_length | flag_straightness) %>%
  select(Filename, ID, SpeedMean_ummin, TrackLength, Straightness,
         flag_speed, flag_length, flag_straightness)

if (nrow(flagged_rows) > 0) {
  cat("\nFlagged rows for inspection:\n")
  print(flagged_rows, n = Inf)
  write_csv(flagged_rows, "flagged_tracks.csv")
  cat("\nFlagged rows saved to: flagged_tracks.csv\n")
}

# --- REMOVAL STEP -------------------------------------------------------------
# After inspecting flagged_tracks.csv, set remove_flagged <- TRUE to remove them.
# Leave as FALSE until you are satisfied with the flags.

remove_flagged <- TRUE   # <-- CHANGE TO TRUE once you have reviewed

if (remove_flagged) {
  n_before <- nrow(tracks)
  tracks <- tracks %>%
    filter(!flag_speed, !flag_length, !flag_straightness)
  n_removed <- n_before - nrow(tracks)
  cat("\nRemoved", n_removed, "flagged tracks.",
      nrow(tracks), "tracks remaining.\n")
} else {
  cat("\nFlagged tracks retained (remove_flagged = FALSE).\n")
  cat("Review flagged_tracks.csv then set remove_flagged <- TRUE to remove.\n")
}


# --- 8. TRACK-LEVEL TRANSFORMATIONS -------------------------------------------
# Transform FIRST at track level, then summarise to ROI level.
# This is the correct order - average the transformed values, not transform the averages.
#
# Speed, TrackLength, MaxDis -> log transform (right-skewed continuous)
# Straightness (0-1 bounded) -> both log and logit provided; remove the one you don't need
# Duration                   -> no transformation (not biologically meaningful to transform)
#
# Zero check: log(0) = -Inf, so we check for zeros before deciding on offset

cat("\n--- Zero checks (log transform safety) ---\n")
cat("Zeros in SpeedMean_ummin:", sum(tracks$SpeedMean_ummin == 0, na.rm = TRUE), "\n")
cat("Zeros in TrackLength:    ", sum(tracks$TrackLength     == 0, na.rm = TRUE), "\n")
cat("Zeros in MaxDis:         ", sum(tracks$MaxDis          == 0, na.rm = TRUE), "\n")
cat("Zeros in Straightness:   ", sum(tracks$Straightness    == 0, na.rm = TRUE), "\n")
cat("Ones  in Straightness:   ", sum(tracks$Straightness    == 1, na.rm = TRUE), "\n")

tracks <- tracks %>%
  mutate(
    # Log transforms (1e-6 offset guards against exact zeros)
    log_Speed       = log(SpeedMean_ummin + 1e-6),
    log_TrackLength = log(TrackLength     + 1e-6),
    log_MaxDis      = log(MaxDis          + 1e-6),
    
    # Straightness: log transform (handles 0s via offset)
    log_Straightness   = log(Straightness + 1e-6),
    
    # Straightness: logit transform (standard for 0-1 bounded variables)
    # 1e-6 offset handles exact 0s and 1s (avoids -Inf / +Inf)
    logit_Straightness = log((Straightness + 1e-6) / (1 - Straightness + 1e-6))
  )


# --- 9. ROI-LEVEL SUMMARY -----------------------------------------------------
# One row per ROI.
# Raw variables: mean, median, SD
# Transformed variables: mean and SD of the already-transformed track-level values
# Duration: kept in raw data only, not summarised

roi_summary <- tracks %>%
  group_by(
    MouseID, TumourID, TumourID_global,
    SliceID_global, ROIID_global,
    Location, Treatment, Day
  ) %>%
  summarise(
    n_tracks = n(),
    
    # --- Raw summaries ---
    # Speed (µm/min)
    Speed_mean   = mean(SpeedMean_ummin, na.rm = TRUE),
    Speed_median = median(SpeedMean_ummin, na.rm = TRUE),
    Speed_sd     = sd(SpeedMean_ummin,   na.rm = TRUE),
    
    # Track Length (µm)
    Length_mean   = mean(TrackLength, na.rm = TRUE),
    Length_median = median(TrackLength, na.rm = TRUE),
    Length_sd     = sd(TrackLength,   na.rm = TRUE),
    
    # Straightness (0-1)
    Straightness_mean   = mean(Straightness, na.rm = TRUE),
    Straightness_median = median(Straightness, na.rm = TRUE),
    Straightness_sd     = sd(Straightness,   na.rm = TRUE),
    
    # MaxDis (µm)
    MaxDis_mean   = mean(MaxDis, na.rm = TRUE),
    MaxDis_median = median(MaxDis, na.rm = TRUE),
    MaxDis_sd     = sd(MaxDis,   na.rm = TRUE),
    
    # --- Transformed summaries (mean and SD of track-level transformed values) ---
    log_Speed_mean      = mean(log_Speed,          na.rm = TRUE),
    log_Speed_sd        = sd(log_Speed,            na.rm = TRUE),
    
    log_Length_mean     = mean(log_TrackLength,    na.rm = TRUE),
    log_Length_sd       = sd(log_TrackLength,      na.rm = TRUE),
    
    log_MaxDis_mean     = mean(log_MaxDis,         na.rm = TRUE),
    log_MaxDis_sd       = sd(log_MaxDis,           na.rm = TRUE),
    
    log_Straight_mean   = mean(log_Straightness,   na.rm = TRUE),
    log_Straight_sd     = sd(log_Straightness,     na.rm = TRUE),
    
    logit_Straight_mean = mean(logit_Straightness, na.rm = TRUE),
    logit_Straight_sd   = sd(logit_Straightness,   na.rm = TRUE),
    
    .groups = "drop"
  )

cat("\n--- ROI summary ---\n")
cat("Total ROIs:", nrow(roi_summary), "\n")
cat("ROIs per Location:\n");  print(table(roi_summary$Location))
cat("ROIs per Treatment:\n"); print(table(roi_summary$Treatment))


# --- 10. SAVE -----------------------------------------------------------------

write_csv(tracks,      "tracks_clean.csv")
write_csv(roi_summary, "roi_summary.csv")

cat("\nSaved: tracks_clean.csv  (", nrow(tracks),     "rows — one row per track)\n")
cat("Saved: roi_summary.csv   (", nrow(roi_summary), "rows — one row per ROI)\n")
cat("\nNext step: run 02_icc_variance.R\n")

