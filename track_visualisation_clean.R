# =============================================================================
# CELL TRACKING ANALYSIS - TRACK VISUALISATION (STANDALONE)
# =============================================================================
# Purpose: Interactively select one or more renamed XML files and visualise
#          their tracks as 2D plots (raw coordinates and origin-normalised).
#          Track endpoints are labelled with the track ID from the matching
#          Allstats CSV file.
#
# Two plots produced per XML file:
#   Raw        : tracks in original x,y coordinate space (y-axis inverted
#                to match microscopy convention — origin top-left)
#   Normalised : all tracks shifted so first detection = (0,0), symmetric
#                axes, useful for comparing displacement patterns
#
# Track labelling:
#   Each track is labelled at its endpoint with the ID from the Allstats CSV.
#   Matching is positional after applying the same duration filter (>10% of
#   max track duration) to both XML and CSV — consistent with 02e pipeline.
#
# Usage:
#   1. Run the script
#   2. A file picker dialog opens — select one or more renamed XML files
#   3. The script automatically finds the matching Allstats CSV in the same
#      folder (same prefix, ends in Allstats.csv)
#   4. Plots are saved to a PDF named after the XML file(s)
#      and displayed in the RStudio viewer
#
# Naming convention assumed (after renaming):
#   XML : M3_T1_D1_Mock_S1_ROI_1_Tracks.xml
#   CSV : M3_T1_D1_Mock_S1_ROI_1_...Allstats.csv
#   The script matches by stripping "_Tracks.xml" and finding a CSV in the
#   same folder whose name starts with the same prefix.
#
# Inputs:
#   XML files selected interactively
#   Matching Allstats CSV files (auto-located in same folder)
#
# Outputs:
#   track_plots_<timestamp>.pdf : all plots for selected files
# =============================================================================


# --- 0. PACKAGES --------------------------------------------------------------
# install.packages(c("xml2", "tidyverse", "patchwork"))

library(xml2)
library(tidyverse)
library(patchwork)

#install.packages("colorspace")
library(colorspace)


# --- 1. CONFIGURATION ---------------------------------------------------------

# Duration filter — must match 02e pipeline exactly
duration_filter_pct <- 0.01  # tracks with duration <= 10% of max are excluded - Comment: here I edited to make sure there is more alignment between the csvs and xmls

# Plot appearance
track_alpha    <- 0.6          # track line transparency
track_linewidth <- 0.6         # track line width
label_size     <- 2          # endpoint label text size
point_size     <- 0.8          # endpoint dot size

# PDF output
pdf_width  <- 10
pdf_height <- 8
ncol_wrap  <- 2                # plots per row in patchwork layout

subsample_seed <- 123   # fixed seed for reproducible track subsampling- relevant for section 9


# --- Location metadata (used by Section 11) ---
location_metadata_csv <- "Z:\\lightmicro-inst\\scharalambous\\Resident_cell_motility_tr_v3\\ROI files\\experimental_file_list_new_COMPLETE.csv"   # CHANGE: full path to location metadata CSV
location_col          <- "Location"  # CHANGE: column name for location label
#         (e.g. "Edge", "CoreEdge", "Core")
meta_mouse_col        <- "MouseID"     # CHANGE: column in metadata matching MouseID
#         values must match pooled_all (e.g. "M1")
meta_slice_col        <- "SliceID"     # CHANGE: column in metadata matching Slice
#         values must match ROIID_global (e.g. "S3")
meta_roi_col          <- "ROIID"       # CHANGE: column in metadata matching ROI
#         values must match ROIID_global (e.g. "ROI1")
meta_tumour_col <- "TumourID"   # column in metadata matching TumourID (e.g. "T1")
meta_treatment_col <- "Treatment"  # column in metadata matching Treatment (e.g. "Mock", "RP")


# --- GMM component plotting ---------------------------------------------------
# Set use_gmm_components <- TRUE to colour tracks by GMM component
# Requires a gmm_component_assignments_<model>_Space<X>.csv file from 02d/6e

use_gmm_components <- TRUE        # gate: set TRUE to enable component colouring

# CHANGE THE FOLLOWING CSV BELOW!
gmm_assignment_csv <- "Z:\\lightmicro-inst\\scharalambous\\Resident_cell_motility_tr_v3\\ROI files\\GMM modelling - complete dataset - new fitting G=5 VEE\\gmm_component_assignments_VEE_SpaceB_NEWFITVEE.csv" # full path to assignment CSV from 6e
min_prob_threshold <- 0.80         # threshold separating high/low confidence pages

# Per-ROI PDF gate — set TRUE to build and save individual ROI plot pages
# (raw + high-confidence component + low-confidence component per XML file)
save_per_roi_pdf <- TRUE
label_component_tracks <- TRUE   # show TrackID at track endpoint in component plots

# Summary plot appearance — used by Section 9 (per-mouse) and Section 10 (all-mice)
summary_track_alpha     <- 0.40   # transparency for pooled summary plots
summary_track_linewidth <- 0.25   # line width for pooled summary plots

# --- READ ME !!! - RELEVANT TO SECTION 09B ----
maxdis_filter_threshold <- 0        # CHANGE: µm — keeps tracks with MaxDis >= this value in filter_components
filter_components       <- c()    # components to apply MaxDis filter to
tracks_clean_csv        <- "Z:\\lightmicro-inst\\scharalambous\\Resident_cell_motility_tr_v3\\ROI files\\Complete data - retrospective fitting G=5 VEE (from v2)\\tracks_clean_v3_alldata_GMM_retrospective.csv"         # CHANGE: full path to tracks_clean.csv
clean_trackid_col       <- "ID"  # TrackID column name in tracks_clean.csv
clean_roi_col           <- "ROIID_global" # ROI column name in tracks_clean.csv
rds_folder              <- "Z:\\lightmicro-inst\\scharalambous\\Resident_cell_motility_tr_v3\\ROI files\\Track Visuals - RETROSPECTIVE - GMM VEE - 5 COMPONENTS"         # folder containing accumulated_tracks_M*.rds

# --- 2. HELPER: FIND MATCHING CSV ---------------------------------------------
# Given an XML filepath, finds the matching Allstats CSV in the same folder.
# Strips "_Tracks.xml" suffix and looks for a CSV starting with the same prefix.

find_matching_csv <- function(xml_path, prefer_merged_tracks = FALSE) {
  # prefer_merged_tracks = TRUE  : used in the GMM path (Section 8) to select
  #                                the merged_tracks CSV for TRACK_INDEX bridging
  # prefer_merged_tracks = FALSE : default; selects Allstats CSV for raw ID labelling

  folder   <- dirname(xml_path)
  xml_name <- basename(xml_path)

  # Strip known XML suffix variants
  prefix <- xml_name %>%
    str_remove("_Tracks\\.xml$") %>%
    str_remove("\\.xml$")

  # Find all CSVs in same folder whose name starts with this prefix
  csv_files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)
  matching  <- csv_files[str_starts(str_to_lower(basename(csv_files)),
                                    str_to_lower(prefix))]

  if (length(matching) == 0) {
    cat("  WARNING: no matching CSV found for", xml_name, "\n")
    return(NULL)
  }

  if (length(matching) == 1) return(matching)

  # Multiple matches — resolve by context
  if (prefer_merged_tracks) {
    preferred <- matching[str_detect(basename(matching), "merged_tracks")]
    if (length(preferred) == 1) return(preferred)
  } else {
    preferred <- matching[str_detect(basename(matching), "Allstats")]
    if (length(preferred) == 1) return(preferred)
  }

  # Still ambiguous — fall back to first with warning
  cat("  WARNING: multiple CSV matches for", xml_name,
      "\u2014 using first:\n  ", basename(matching[1]), "\n")
  return(matching[1])
}
# --- 2b. HELPER: LOAD GMM ASSIGNMENTS -----------------------------------------
# Loads component assignment CSV exported from Script 02d Section 6e.
# Builds a match key from ROIID_global + Day that mirrors the XML filename
# structure for unambiguous joining.
#
# Match key construction from XML filename:
#   M1_T1_D2_Mock_S3_ROI_1_Tracks.xml
#   -> Mouse=M1, Tumour=T1, Day=D2, Treatment=Mock, Slice=S3, ROI=ROI1
#   -> ROIID_global = M1_T1_Mock_S3_ROI1  (Treatment title-cased, ROI_ removed)
#   -> Day = D2
#
# Case normalisation: Treatment is title-cased in both the CSV and the
# parsed filename to handle Mock vs MOCK inconsistency across files.

load_gmm_assignments <- function(csv_path) {
  df <- read_csv(csv_path, show_col_types = FALSE) %>%
    mutate(
      # Normalise Treatment embedded in ROIID_global to title case
      
      # ROIID format: Mouse_Tumour_Treatment_Slice_ROI
      # ROIID_norm = str_replace(
      #   ROIID_global,
      #   "(M\\d+_T\\d+_)(\\w+)(_S\\d+_ROI\\d+)",
      #   function(m) {
      #     parts <- str_match(m, "(M\\d+_T\\d+_)(\\w+)(_S\\d+_ROI\\d+)")
      #     paste0(parts[2], str_to_title(parts[3]), parts[4])
      #   }
      # ),
      
      ROIID_norm = paste0(
        str_extract(ROIID_global, "^M\\d+_T\\d+_"),
        {
          trt <- str_extract(ROIID_global, "(?<=T\\d_)\\w+(?=_S\\d+)")
          if_else(str_to_upper(trt) == "RP", "RP", str_to_title(trt))
        },
        str_extract(ROIID_global, "_S\\d+_ROI\\d+$")
      ),
      Day = as.character(Day)
    )
  cat("GMM assignments loaded:", nrow(df), "tracks\n")
  cat("Components:", paste(sort(unique(df$GMM_component)), collapse = ", "), "\n")
  cat("ROIs covered:", n_distinct(df$ROIID_global), "\n\n")
  df
}

# Parse ROIID_global from XML filename
# Returns list(roiid = "M1_T1_Mock_S3_ROI1", day = "D2")
parse_roiid_from_xml <- function(xml_filename) {
  # Strip path and suffix
  name <- basename(xml_filename) %>%
    str_remove("_Tracks\\.xml$") %>%
    str_remove("\\.xml$")
  
  # Expected pattern: Mouse_Tumour_Day_Treatment_Slice_ROI_N
  parts <- str_match(name,
                     "^(M\\d+)_(T\\d+)_(D\\d+)_([A-Za-z]+)_(S\\d+)_ROI_(\\d+)$")
  
  if (is.na(parts[1])) {
    cat("  WARNING: could not parse ROIID from filename:", xml_filename, "\n")
    return(NULL)
  }
  
  mouse     <- parts[2]
  tumour    <- parts[3]
  day       <- parts[4]
  treatment <- if (str_to_upper(parts[5]) == "RP") "RP" else str_to_title(parts[5])
  slice     <- parts[6]
  roi_num   <- parts[7]
  
  roiid <- paste0(mouse, "_", tumour, "_", treatment, "_",
                  slice, "_ROI", roi_num)
  
  list(roiid = roiid, day = day, mouse = mouse)
}


# --- 3. HELPER: PARSE XML -----------------------------------------------------

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

  tracks_raw        <- Filter(Negate(is.null), tracks_raw)
  frame_interval_s  <- frame_interval_s

  list(tracks = tracks_raw,
       frame_interval_s = frame_interval_s,
       n_raw = length(tracks_raw))
}


# --- 4. HELPER: APPLY DURATION FILTER -----------------------------------------

apply_duration_filter <- function(tracks_raw, frame_interval_s) {

  durations <- sapply(tracks_raw, function(tr) {
    (max(tr$t) - min(tr$t)) * frame_interval_s
  })

  max_dur   <- max(durations)
  threshold <- duration_filter_pct * max_dur
  keep      <- which(durations > threshold)

  list(
    tracks    = tracks_raw[keep],
    durations = durations[keep],
    max_dur   = max_dur,
    threshold = threshold,
    n_removed = length(tracks_raw) - length(keep)
  )
}


# --- 5. HELPER: LOAD IDs FROM CSV ---------------------------------------------

load_track_ids <- function(csv_path, duration_filter_pct) {
  
  # Read column names from row 1, data from row 5 onwards
  # (rows 2-4 are verbose labels, units, and empty — skipped)
  col_names <- read.csv(csv_path, nrows = 1, header = FALSE,
                        check.names = TRUE) %>%
    unlist() %>% as.character()
  
  raw <- read.csv(csv_path, skip = 4, header = FALSE,
                  fill = TRUE, check.names = FALSE)
  colnames(raw) <- col_names
  
  # Standardise to UPPER_SNAKE_CASE for robust matching
  colnames(raw) <- colnames(raw) %>%
    str_trim() %>%
    str_replace_all("[. ]+", "_") %>%
    str_to_upper()
  
  # TRACK_INDEX is the TrackMate ID that matches GMM assignment CSV ID column
  if (!"TRACK_INDEX" %in% names(raw)) {
    cat("  WARNING: TRACK_INDEX not found in CSV — using sequential IDs\n")
    cat("  Available columns:", paste(names(raw)[1:6], collapse = ", "), "\n")
    return(tibble(row_index = seq_len(nrow(raw)),
                  TrackID   = seq_len(nrow(raw))))
  }
  
  tibble(
    row_index = seq_len(nrow(raw)),
    TrackID   = as.integer(raw$TRACK_INDEX)
  )
}
# --- 6. HELPER: BUILD TRACK DATAFRAME -----------------------------------------
# Combines filtered XML tracks with IDs from CSV into one flat dataframe
# with raw and normalised coordinates.

build_track_df <- function(filtered_tracks, track_ids) {

  n_xml <- length(filtered_tracks)
  n_csv <- nrow(track_ids)

  if (n_xml != n_csv) {
    cat("  WARNING: track count mismatch after filtering —",
        "XML:", n_xml, "| CSV:", n_csv, "\n")
    cat("  Using sequential index as ID\n")
    track_ids <- tibble(row_index = seq_len(n_xml),
                        TrackID   = seq_len(n_xml))
  }

  # Flatten all tracks into one dataframe with ID attached
  track_df <- lapply(seq_along(filtered_tracks), function(i) {
    tr <- filtered_tracks[[i]]
    tr$TrackID    <- track_ids$TrackID[i]
    tr$track_seq  <- i
    tr
  }) %>%
    bind_rows() %>%
    # Normalise coordinates: shift each track so first detection = (0,0)
    group_by(TrackID) %>%
    mutate(
      x_norm = x - first(x),
      y_norm = y - first(y)
    ) %>%
    ungroup()

  track_df
}


# --- 7. HELPER: MAKE PLOTS ----------------------------------------------------

# Colour for tracks — single neutral colour per file
track_colour <- "#457B9D"

make_raw_plot <- function(track_df, file_label) {

  # Endpoint of each track for labelling
  endpoints <- track_df %>%
    group_by(TrackID) %>%
    slice_tail(n = 1) %>%
    ungroup()

  ggplot(track_df, aes(x = x, y = y, group = TrackID)) +
    geom_path(colour = track_colour, alpha = track_alpha,
              linewidth = track_linewidth) +
    geom_point(data = endpoints, aes(x = x, y = y),
               colour = track_colour, size = point_size) +
    geom_text(data = endpoints, aes(x = x, y = y, label = TrackID),
              size = label_size, colour = "grey20",
              hjust = -0.3, vjust = 0.5) +
    scale_y_reverse() +
    coord_fixed() +
    labs(
      title    = paste0(file_label, " — Raw coordinates"),
      subtitle = paste0(n_distinct(track_df$TrackID), " tracks"),
      x        = "X (µm)",
      y        = "Y (µm)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title    = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 9, colour = "grey40"),
      panel.grid.minor = element_blank(),
      legend.position  = "none"
    )
}

make_norm_plot <- function(track_df, file_label) {

  # Endpoint of each normalised track for labelling
  endpoints_norm <- track_df %>%
    group_by(TrackID) %>%
    slice_tail(n = 1) %>%
    ungroup()

  # Symmetric axis limits based on max excursion
  max_axis <- max(abs(c(track_df$x_norm, track_df$y_norm)),
                  na.rm = TRUE) %>%
    ceiling()
  if (!is.finite(max_axis) || max_axis == 0) max_axis <- 1

  ggplot(track_df, aes(x = x_norm, y = y_norm, group = TrackID)) +
    geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.4) +
    geom_path(colour = track_colour, alpha = track_alpha,
              linewidth = track_linewidth) +
    geom_point(data = endpoints_norm, aes(x = x_norm, y = y_norm),
               colour = track_colour, size = point_size) +
    geom_text(data = endpoints_norm,
              aes(x = x_norm, y = y_norm, label = TrackID),
              size = label_size, colour = "grey20",
              hjust = -0.3, vjust = 0.5) +
    # Mark origin
    
    # geom_point(aes(x = 0, y = 0), colour = "black",
    #            shape = 3, size = 2, stroke = 0.8,
    #            inherit.aes = FALSE) +
    annotate("point", x = 0, y = 0, colour = "black",
             shape = 3, size = 2, stroke = 0.8)+
    coord_equal(xlim = c(-max_axis, max_axis),
                ylim = c(-max_axis, max_axis),
                expand = FALSE) +
    labs(
      title    = paste0(file_label, " — Normalised (origin = start)"),
      subtitle = paste0(n_distinct(track_df$TrackID), " tracks",
                        "  |  + marks origin"),
      x        = "X displacement (µm)",
      y        = "Y displacement (µm)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title    = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 9, colour = "grey40"),
      panel.grid.minor = element_blank(),
      legend.position  = "none"
    )
}

make_component_plot <- function(track_df, file_label, subtitle_extra) {
  
  if (nrow(track_df) == 0) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5,
                      label = paste("No tracks —", subtitle_extra),
                      size = 4, colour = "grey40") +
             theme_void())
  }
  
  comp_pal <- c("1" = "#457B9D", "2" = "#E63946",
                "3" = "#2A9D8F", "4" = "#F4A261", "5" = "#6A0572")
  
  max_axis <- max(abs(c(track_df$x_norm, track_df$y_norm)),
                  na.rm = TRUE) %>% ceiling()
  if (!is.finite(max_axis) || max_axis == 0) max_axis <- 1
  
  endpoints_norm <- track_df %>%
    group_by(TrackID) %>%
    slice_tail(n = 1) %>%
    ungroup()
  
  ggplot(track_df,
         aes(x = x_norm, y = y_norm, group = TrackID,
             colour = factor(GMM_component))) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.4) +
    geom_path(alpha = track_alpha, linewidth = track_linewidth) +
    geom_point(data = endpoints_norm,
               aes(x = x_norm, y = y_norm,
                   colour = factor(GMM_component)),
               size = point_size) +
    # geom_point(aes(x = 0, y = 0), colour = "black",
    #            shape = 3, size = 2, stroke = 0.8,
    #            inherit.aes = FALSE) +
    {if (label_component_tracks)
      geom_text(data        = endpoints_norm,
                aes(x       = x_norm,
                    y       = y_norm,
                    label   = TrackID),
                colour      = "black",
                size        = label_size,
                hjust       = -0.2,
                vjust       = 0,
                show.legend = FALSE)} +
  
    annotate("point", x = 0, y = 0, colour = "black",
             shape = 3, size = 2, stroke = 0.8)+
    scale_colour_manual(values = comp_pal, name = "Component",
                        na.value = "grey70") +
    facet_wrap(~ paste0("Comp ", GMM_component),
               nrow = 2, ncol = 3) +
    coord_equal(xlim = c(-max_axis, max_axis),
                ylim = c(-max_axis, max_axis),
                expand = FALSE) +
    labs(
      title    = paste0(file_label, " — Normalised tracks by GMM component"),
      subtitle = subtitle_extra,
      x        = "X displacement (µm)",
      y        = "Y displacement (µm)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title       = element_text(face = "bold", size = 11),
      plot.subtitle    = element_text(size = 9, colour = "grey40"),
      panel.grid.minor = element_blank(),
      legend.position  = "right",
      strip.text       = element_text(face = "bold", size = 9),
      panel.spacing    = unit(0.5, "lines")   
    )
}

# --- 8. MAIN: FILE SELECTION AND PROCESSING -----------------------------------

cat("=== TRACK VISUALISATION ===\n")
cat("Select one or more XML track files in the dialog...\n\n")

# Load GMM assignments if enabled
if (use_gmm_components) {
  if (gmm_assignment_csv == "" || !file.exists(gmm_assignment_csv))
    stop("gmm_assignment_csv path missing or file not found. Set in Section 1.")
  gmm_assignments <- load_gmm_assignments(gmm_assignment_csv)
}

xml_paths <- choose.files(
  caption = "Select XML track file(s)",
  filters = matrix(c("XML files", "*.xml", "All files", "*.*"),
                   nrow = 2, byrow = TRUE)
)

if (length(xml_paths) == 0) stop("No files selected.")

cat("Files selected:", length(xml_paths), "\n\n")

# Output PDF filename
timestamp  <- format(Sys.time(), "%Y%m%d_%H%M%S")
pdf_output <- paste0("track_plots_", timestamp, ".pdf")

if (save_per_roi_pdf) {
  all_plots  <- list()
  plot_index <- 1
  plot_steps <- c()
}
all_tracks <- list()   # accumulates track_df across all XMLs for summary plots

for (xml_path in xml_paths) {

  xml_name   <- basename(xml_path)
  file_label <- str_remove(xml_name, "_Tracks\\.xml$") %>%
                str_remove("\\.xml$")

  cat("Processing:", xml_name, "\n")

  # --- Find matching CSV (only needed when use_gmm_components = FALSE) ---
  csv_path <- if (!use_gmm_components) find_matching_csv(xml_path) else NULL

  # --- Parse XML ---
  parsed <- parse_xml_tracks(xml_path)
  cat("  XML tracks (raw):", parsed$n_raw, "\n")

  # --- Apply duration filter ---
  filtered <- apply_duration_filter(parsed$tracks, parsed$frame_interval_s)
  cat("  Duration filter: removed", filtered$n_removed, "tracks",
      "(threshold:", round(filtered$threshold, 1), "s)\n")
  cat("  XML tracks (filtered):", length(filtered$tracks), "\n")

  # --- Load IDs ---
  if (use_gmm_components) {
    parsed_roi <- parse_roiid_from_xml(xml_path)
    
    if (!is.null(parsed_roi)) {
      roi_assignments <- gmm_assignments %>%
        filter(ROIID_norm == parsed_roi$roiid,
               Day        == parsed_roi$day)
      
      cat("  GMM assignments found:", nrow(roi_assignments), "tracks",
          "| ROIID:", parsed_roi$roiid, "| Day:", parsed_roi$day, "\n")
      
      # Load TRACK_INDEX from merged_tracks CSV to bridge XML position → TrackMate ID
      csv_path <- find_matching_csv(xml_path, prefer_merged_tracks = TRUE)
      
      if (!is.null(csv_path)) {
        track_ids <- load_track_ids(csv_path, duration_filter_pct)
        cat("  Merged tracks CSV loaded:", nrow(track_ids), "rows\n")
        
        if (nrow(track_ids) != length(filtered$tracks)) {
          cat("  WARNING: XML tracks (", length(filtered$tracks),
              ") != CSV rows (", nrow(track_ids),
              ") — positional match unreliable\n")
        }
      } else {
        cat("  WARNING: no merged_tracks CSV found — cannot map IDs\n")
        track_ids <- tibble(
          row_index = seq_along(filtered$tracks),
          TrackID   = seq_along(filtered$tracks)
        )
      }
      
    } else {
      cat("  Filename parse failed — using sequential IDs\n")
      track_ids <- tibble(
        row_index = seq_along(filtered$tracks),
        TrackID   = seq_along(filtered$tracks)
      )
      roi_assignments <- NULL
    }
    
  } else {
    if (!is.null(csv_path)) {
      cat("  CSV found:", basename(csv_path), "\n")
      track_ids <- load_track_ids(csv_path, duration_filter_pct)
      cat("  CSV tracks (filtered):", nrow(track_ids), "\n")
    } else {
      cat("  No CSV found — using sequential IDs\n")
      track_ids <- tibble(
        row_index = seq_along(filtered$tracks),
        TrackID   = seq_along(filtered$tracks)
      )
    }
    roi_assignments <- NULL
  }

  # --- Build combined dataframe ---
  track_df <- build_track_df(filtered$tracks, track_ids)
  
  
  # --- Attach GMM component if enabled ---
  if (use_gmm_components && !is.null(roi_assignments)) {
    
    comp_cols <- roi_assignments %>%
      select(ID, GMM_component, GMM_max_prob, starts_with("P_comp"))
    
    n_before <- n_distinct(track_df$TrackID)
    
    track_df <- track_df %>%
      left_join(comp_cols, by = c("TrackID" = "ID"))
    
    n_matched <- sum(!is.na(track_df$GMM_component[
      !duplicated(track_df$TrackID)]))
    
    cat("  GMM components matched:", n_matched, "/", n_before, "tracks\n")
    
    if (n_before - n_matched > 0)
      cat("  WARNING:", n_before - n_matched,
          "tracks unmatched — check ID alignment\n")
    
  } else if (use_gmm_components && is.null(roi_assignments)) {
    cat("  GMM component colouring skipped — filename parse failed\n")
    track_df$GMM_component <- NA_integer_
    track_df$GMM_max_prob  <- NA_real_
  }
  
  # --- Accumulate tracks for summary plots ---
  if (use_gmm_components && !is.null(parsed_roi)) {
    all_tracks[[length(all_tracks) + 1]] <- track_df %>%
      mutate(
        MouseID     = parsed_roi$mouse,
        ROIID_global = parsed_roi$roiid
      )
  }

  # --- Make plots (gated by save_per_roi_pdf) ---
  if (save_per_roi_pdf) {

    if (use_gmm_components && "GMM_component" %in% names(track_df) &&
        any(!is.na(track_df$GMM_component))) {

      # Three pages: raw coordinates, high-confidence components, low-confidence
      p_raw       <- make_raw_plot(track_df, file_label)
      p_comp_high <- make_component_plot(
        track_df %>% filter(GMM_max_prob >= min_prob_threshold),
        file_label,
        paste0("GMM components (max_prob ≥ ", min_prob_threshold, ")")
      )
      p_comp_low <- make_component_plot(
        track_df %>% filter(GMM_max_prob < min_prob_threshold),
        file_label,
        paste0("GMM components (max_prob < ", min_prob_threshold, ")")
      )

      all_plots[[plot_index]]     <- p_raw
      all_plots[[plot_index + 1]] <- p_comp_high
      all_plots[[plot_index + 2]] <- p_comp_low
      plot_index <- plot_index + 3
      plot_steps <- c(plot_steps, 3)

    } else {
      # Fallback path: GMM colouring unavailable for this file
      tryCatch({
        p_raw  <- make_raw_plot(track_df,  file_label)
        p_norm <- make_norm_plot(track_df, file_label)
        all_plots[[plot_index]]     <- p_raw
        all_plots[[plot_index + 1]] <- p_norm
        plot_index <<- plot_index + 2
        plot_steps <<- c(plot_steps, 2)
      }, error = function(e) {
        cat("  WARNING: fallback plot failed for", file_label,
            "—", conditionMessage(e), "\n")
        cat("  File skipped in PDF output.\n")
      })
    }

    cat("  Plots created\n\n")
  }
}


# --- 9. SAVE AND DISPLAY ------------------------------------------------------

if (save_per_roi_pdf && length(all_plots) > 0) {

  cat("Saving per-ROI PDF to:", pdf_output, "\n")

  pdf(pdf_output, width = pdf_width, height = pdf_height)

  if (use_gmm_components) {
    idx <- 1
    for (step in plot_steps) {
      for (j in 0:(step - 1)) print(all_plots[[idx + j]])
      idx <- idx + step
    }
  } else {
    for (i in seq(1, length(all_plots), by = 2)) {
      print(all_plots[[i]] | all_plots[[i + 1]])
    }
  }

  dev.off()

  # Display first plot in RStudio viewer
  if (use_gmm_components && length(all_plots) >= 2) {
    print(all_plots[[2]])   # high-confidence component page
  } else if (length(all_plots) >= 2) {
    print(all_plots[[1]] | all_plots[[2]])
  }

  cat("Saved:", pdf_output, "\n")

} else if (save_per_roi_pdf) {
  cat("No plots accumulated — per-ROI PDF skipped.\n")
}

cat("\n=== Done ===\n")
cat("Per-ROI PDF: gated by save_per_roi_pdf in Section 1\n")
cat("Summary PDFs: built in the block below (Section 9 summary)\n")

# --- Summary component plots (pooled across all XMLs) -------------------------

if (use_gmm_components && length(all_tracks) > 0) {
  
  cat("\nBuilding summary component plots...\n")
  
  # Combine all accumulated tracks
  pooled <- bind_rows(all_tracks) %>%
    filter(!is.na(GMM_component)) %>%
    mutate(
      GMM_component = factor(GMM_component),
      MouseID       = factor(MouseID)
    )
  
  cat("Total tracks pooled:", n_distinct(pooled %>% select(TrackID, ROIID_global)),
      "\n")
  cat("Mice present:", paste(levels(pooled$MouseID), collapse = ", "), "\n")
  
  # Shared symmetric axis limit across all components and mice
  max_axis <- max(abs(c(pooled$x_norm, pooled$y_norm)), na.rm = TRUE) %>%
    ceiling()
  
  # Track counts per component × mouse (used for subsample reporting)
  mouse_counts <- pooled %>%
    group_by(GMM_component, MouseID) %>%
    summarise(n_tracks = n_distinct(TrackID), .groups = "drop")

  # Subsample helper — max N tracks per component per mouse
  subsample_tracks <- function(df, n_max = 200) {
    df %>%
      group_by(GMM_component, MouseID, TrackID) %>%
      mutate(keep = cur_group_id()) %>%
      ungroup() %>%
      group_by(GMM_component, MouseID) %>%
      mutate(track_sample = {
        set.seed(subsample_seed)
        TrackID %in% sample(unique(TrackID),
                            min(n_max, n_distinct(TrackID)))
      }) %>%
      ungroup() %>%
      filter(track_sample) %>%
      select(-track_sample, -keep)
  }
  
  # Per-component normalised track plot function.
  # Tracks coloured by MouseID using fixed_colour palette (mouse_pal).
  # summary_track_linewidth and summary_track_alpha from Section 1.
  make_summary_plot <- function(df, title_str, subtitle_str, fixed_colour) {

    if (nrow(df) == 0) {
      return(ggplot() +
               annotate("text", x = 0.5, y = 0.5,
                        label = "No tracks", size = 4, colour = "grey40") +
               theme_void())
    }

    ggplot(df,
           aes(x = x_norm, y = y_norm,
               group = interaction(TrackID, ROIID_global))) +
      geom_path(aes(colour = MouseID),
                linewidth = summary_track_linewidth,
                alpha     = summary_track_alpha) +
      scale_colour_manual(values = fixed_colour, name = "Mouse") +
      geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
      geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
      annotate("point", x = 0, y = 0, colour = "black",
               shape = 3, size = 1.5, stroke = 0.6) +
      facet_wrap(~ paste0("Comp ", GMM_component), nrow = 2, ncol = 3) +
      coord_equal(xlim = c(-max_axis, max_axis),
                  ylim = c(-max_axis, max_axis),
                  expand = FALSE) +
      labs(title    = title_str,
           subtitle = subtitle_str,
           x        = "X displacement (µm)",
           y        = "Y displacement (µm)") +
      theme_bw(base_size = 10) +
      theme(plot.title       = element_text(face = "bold", size = 12),
            plot.subtitle    = element_text(size = 9, colour = "grey40"),
            panel.grid.minor = element_blank(),
            strip.text       = element_text(face = "bold", size = 9),
            legend.position  = "right")
  }
  
  # --- Save accumulated tracks to disk -----------------------------------------
  # Preserves all track data across runs so all-mice summary can be built
  # by loading and combining per-mouse saved files rather than re-running.
  # File named by mouse(es) present in this run and timestamp.
  

    pooled_save <- bind_rows(all_tracks)
    mice_label  <- paste(sort(unique(pooled_save$MouseID)), collapse = "_")
    tracks_file <- paste0("accumulated_tracks_", mice_label, "_",
                          timestamp, ".rds")
    
    saveRDS(pooled_save, tracks_file)
    cat("Saved accumulated tracks:", tracks_file, "\n")
    cat("  Rows:", nrow(pooled_save), "\n")
    cat("  Tracks:", n_distinct(
      pooled_save %>% select(TrackID, ROIID_global)), "\n")
    cat("  Mice:", paste(unique(pooled_save$MouseID), collapse = ", "), "\n\n")
  
  
  
  # Mouse colour palette for per-mouse plots
  mouse_pal <- c("M1" = "#457B9D", "M2" = "#E63946", "M3" = "#F4A261", "M4" = "#2A9D8F","M5" = "#6A4C93")
  
  # --- Per-mouse PDFs ---
  for (mouse in levels(pooled$MouseID)) {
    
    mouse_data <- pooled %>% filter(MouseID == mouse)
    mouse_file <- paste0("component_summary_", mouse, "_",
                         timestamp, ".pdf")
    
    high_sub <- subsample_tracks(
      mouse_data %>% filter(GMM_max_prob >= min_prob_threshold))
    low_sub  <- subsample_tracks(
      mouse_data %>% filter(GMM_max_prob <  min_prob_threshold))

    p_high <- make_summary_plot(
      high_sub,
      title_str    = paste0(mouse, " — Normalised tracks by component"),
      subtitle_str = paste0("max_prob ≥ ", min_prob_threshold,
                            " | subsampled ≤ 200 tracks/component"),
      fixed_colour = mouse_pal
    )

    p_low <- make_summary_plot(
      low_sub,
      title_str    = paste0(mouse, " — Normalised tracks by component"),
      subtitle_str = paste0("max_prob < ", min_prob_threshold,
                            " | subsampled ≤ 200 tracks/component"),
      fixed_colour = mouse_pal
    )
    
    pdf(mouse_file, width = 12, height = 10)
    print(p_high)
    print(p_low)
    dev.off()
    cat("Saved:", mouse_file, "\n")
  }

} else {
  cat("\nSummary plots skipped — use_gmm_components = FALSE or no tracks accumulated.\n")
}


# =============================================================================
# SECTION 9b — STANDALONE BLOCK: MaxDis FILTER ON SELECTED COMPONENTS
# =============================================================================
# Loads all mouse RDS files, joins MaxDis and GMM_component from
# tracks_clean.csv using TrackID + ROIID_global, then removes tracks from
# filter_components where MaxDis < maxdis_filter_threshold.
#
# IMPORTANT: filtering uses GMM_component_csv (from tracks_clean.csv) not the
# RDS GMM_component — these may differ if built from different GMM runs.
# Plots and saved RDS retain the original RDS GMM_component for coordinates.
#
# Section 1 additions required:
#   maxdis_filter_threshold <- 7
#   filter_components       <- c(1, 3)   # in tracks_clean.csv numbering
#   tracks_clean_csv        <- ""
#   clean_trackid_col       <- "ID"
#   clean_roi_col           <- "ROIID_global"
#
# Outputs:
#   component_summary_filtered_<mouse>_<timestamp>.pdf
#   accumulated_tracks_<mouse>_filtered.rds
# =============================================================================

# --- Clear stale objects from previous runs ----------------------------------
rm(list = intersect(ls(), c(
  "pooled_9b", "pooled_filt", "pooled_sub",
  "maxdis_df", "max_axis",
  "before_counts", "after_counts", "summary_counts",
  "rds_list", "rds_files", "mice_filt",
  "subsample_filt", "timestamp_filt",
  "n_total", "n_matched"
)))

# --- Standalone preamble: fallback values if Section 1 was not run -----------
if (!exists("summary_track_alpha"))     summary_track_alpha     <- 0.40
if (!exists("summary_track_linewidth")) summary_track_linewidth <- 0.25
if (!exists("min_prob_threshold"))      min_prob_threshold      <- 0.80
if (!exists("subsample_seed"))          subsample_seed          <- 123
if (!exists("rds_folder"))              rds_folder              <- ""  # CHANGE
if (!exists("mouse_pal")) {
  mouse_pal <- c("M1" = "#457B9D", "M2" = "#E63946", "M3" = "#F4A261",
                 "M4" = "#2A9D8F", "M5" = "#6A0572")
}

# --- Validate Section 1 parameters -------------------------------------------
required_9b <- c("maxdis_filter_threshold", "filter_components",
                 "tracks_clean_csv", "clean_trackid_col",
                 "clean_roi_col", "rds_folder")
missing_9b <- required_9b[!sapply(required_9b, exists)]
if (length(missing_9b) > 0)
  stop("Missing from environment \u2014 set in Section 1:\n  ",
       paste(missing_9b, collapse = ", "))
if (nchar(rds_folder) == 0)
  stop("rds_folder is empty. Set path in Section 1 or standalone preamble.")
if (!file.exists(tracks_clean_csv))
  stop("tracks_clean_csv not found. Set path in Section 1.")

cat("=== Section 9b: MaxDis filter ===\n")
cat("    Components (CSV numbering):", paste(filter_components, collapse = ", "), "\n")
cat("    Threshold:", maxdis_filter_threshold, "\u00b5m\n\n")

# --- Local plot function (self-contained, no dependency on Section 9) --------
# Normalised track plot faceted by GMM_component, coloured by MouseID.
# Uses max_axis (computed in 9b-iv) and summary_track_* from preamble.
make_summary_plot_9b <- function(df, title_str, subtitle_str, colour_pal) {
  
  if (nrow(df) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = "No tracks", size = 4, colour = "grey40") +
        theme_void()
    )
  }
  
  ggplot(df,
         aes(x = x_norm, y = y_norm,
             group = interaction(TrackID, ROIID_global))) +
    geom_path(aes(colour = MouseID),
              linewidth = summary_track_linewidth,
              alpha     = summary_track_alpha) +
    scale_colour_manual(values = colour_pal, name = "Mouse") +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    annotate("point", x = 0, y = 0, colour = "black",
             shape = 3, size = 1.5, stroke = 0.6) +
    facet_wrap(~ paste0("Comp ", GMM_component), nrow = 2, ncol = 3) +
    coord_equal(xlim = c(-max_axis, max_axis),
                ylim = c(-max_axis, max_axis),
                expand = FALSE) +
    labs(title    = title_str,
         subtitle = subtitle_str,
         x        = "X displacement (\u00b5m)",
         y        = "Y displacement (\u00b5m)") +
    theme_bw(base_size = 10) +
    theme(plot.title       = element_text(face = "bold", size = 12),
          plot.subtitle    = element_text(size = 9, colour = "grey40"),
          panel.grid.minor = element_blank(),
          strip.text       = element_text(face = "bold", size = 9),
          legend.position  = "right")
}

# --- 9b-i. Load RDS files ----------------------------------------------------
rds_files <- list.files(rds_folder,
                        pattern    = "^accumulated_tracks_M\\d+.*\\.rds$",
                        full.names = TRUE) %>%
  .[!grepl("_filtered", basename(.))]

if (length(rds_files) == 0)
  stop("No unfiltered RDS files found in rds_folder:\n  ", rds_folder)

cat("RDS files found:", length(rds_files), "\n")
invisible(lapply(rds_files, function(f) cat(" ", basename(f), "\n")))

pooled_9b <- bind_rows(lapply(rds_files, readRDS)) %>%
  select(-any_of("MaxDis"))  # drop pre-existing MaxDis to avoid column conflict

cat("Unique tracks loaded:",
    nrow(distinct(pooled_9b, TrackID, ROIID_global)), "\n\n")

# --- 9b-ii. Load MaxDis + GMM_component_csv and join -------------------------
# GMM_component_csv kept separate from RDS GMM_component because the two may
# have been built from different GMM runs with shuffled component numbering.
# MaxDis filtering and plots both use the CSV component numbering.
maxdis_df <- read_csv(tracks_clean_csv, show_col_types = FALSE) %>%
  select(all_of(c(clean_trackid_col, clean_roi_col,
                  "MaxDis", "GMM_component"))) %>%
  rename(TrackID           = all_of(clean_trackid_col),
         ROIID_global      = all_of(clean_roi_col),
         GMM_component_csv = GMM_component) %>%
  mutate(TrackID      = as.numeric(TrackID),
         ROIID_global = as.character(ROIID_global)) %>%
  distinct(TrackID, ROIID_global, .keep_all = TRUE)

pooled_9b <- pooled_9b %>%
  left_join(maxdis_df, by = c("TrackID", "ROIID_global"))

n_total   <- nrow(distinct(pooled_9b, TrackID, ROIID_global))
n_matched <- pooled_9b %>%
  filter(!is.na(MaxDis)) %>%
  distinct(TrackID, ROIID_global) %>%
  nrow()
cat("MaxDis join:", n_matched, "/", n_total, "tracks matched\n")
if (n_total - n_matched > 0)
  cat("WARNING:", n_total - n_matched,
      "tracks have no MaxDis \u2014 kept (not penalised for join miss)\n\n")

# --- 9b-iii. Apply filter and replace GMM_component with CSV version ---------
# filter_components uses CSV component numbering.
# NA MaxDis = unknown → kept regardless.
# GMM_component is overwritten with GMM_component_csv so plots and saved RDS
# use the same component numbering as tracks_clean.csv throughout.

before_counts <- pooled_9b %>%
  distinct(TrackID, ROIID_global, GMM_component_csv, MouseID) %>%
  count(MouseID, GMM_component_csv, name = "n_before")

pooled_filt <- pooled_9b %>%
  filter(
    !(GMM_component_csv %in% filter_components &
        !is.na(MaxDis) & MaxDis < maxdis_filter_threshold)
  ) %>%
  mutate(GMM_component = factor(GMM_component_csv))

pooled_sub <- pooled_9b %>%
  filter(
    GMM_component_csv %in% filter_components &
      !is.na(MaxDis) & MaxDis < maxdis_filter_threshold
  ) %>%
  mutate(GMM_component = factor(GMM_component_csv))

after_counts <- pooled_filt %>%
  distinct(TrackID, ROIID_global, GMM_component_csv, MouseID) %>%
  count(MouseID, GMM_component_csv, name = "n_after")

summary_counts <- before_counts %>%
  left_join(after_counts, by = c("MouseID", "GMM_component_csv")) %>%
  mutate(n_after   = replace_na(n_after, 0L),
         n_removed = n_before - n_after) %>%
  arrange(MouseID, GMM_component_csv)

cat("--- Filter summary (CSV component numbering) ---\n")
cat("Threshold:", maxdis_filter_threshold, "\u00b5m",
    "| Components:", paste(filter_components, collapse = ", "), "\n")
print(summary_counts, n = Inf)
cat("\n")

# --- 9b-iv. Global axis limit ------------------------------------------------
max_axis <- max(abs(c(pooled_filt$x_norm, pooled_filt$y_norm)),
                na.rm = TRUE) %>% ceiling()

# --- 9b-v. Subsample helper --------------------------------------------------
subsample_filt <- function(df, n_max = 200) {
  df %>%
    group_by(GMM_component) %>%
    mutate(track_sample = {
      set.seed(subsample_seed)
      TrackID %in% sample(unique(TrackID), min(n_max, n_distinct(TrackID)))
    }) %>%
    ungroup() %>%
    filter(track_sample) %>%
    select(-track_sample)
}

# --- 9b-vi. Per-mouse PDFs ---------------------------------------------------
# Page 1: translocating tracks (MaxDis >= threshold, all components)
# Page 2: non-translocating tracks (MaxDis < threshold, filter_components only)
timestamp_filt <- format(Sys.time(), "%Y%m%d_%H%M%S")
mice_filt      <- sort(unique(pooled_filt$MouseID))

invisible(lapply(mice_filt, function(mouse) {
  
  mouse_above <- pooled_filt %>%
    filter(MouseID == mouse, !is.na(GMM_component))
  mouse_below <- pooled_sub %>%
    filter(MouseID == mouse, !is.na(GMM_component))
  
  if (nrow(mouse_above) == 0) {
    cat("  Skipping", mouse, "\u2014 no tracks after filtering\n")
    return(invisible(NULL))
  }
  
  above_sub <- subsample_filt(mouse_above)
  below_sub <- if (nrow(mouse_below) > 0) subsample_filt(mouse_below) else mouse_below
  
  subtitle_base <- paste0(
    "CSV comps ", paste(filter_components, collapse = "/"),
    " | MaxDis \u2265 ", maxdis_filter_threshold, " \u00b5m",
    " | \u2264 200 tracks/component"
  )
  
  p_above <- make_summary_plot_9b(
    above_sub,
    title_str    = paste0(mouse, " \u2014 Translocating (MaxDis \u2265 ",
                          maxdis_filter_threshold, " \u00b5m)"),
    subtitle_str = subtitle_base,
    colour_pal   = mouse_pal
  )
  
  p_below <- make_summary_plot_9b(
    below_sub,
    title_str    = paste0(mouse, " \u2014 Non-translocating (MaxDis < ",
                          maxdis_filter_threshold, " \u00b5m)"),
    subtitle_str = paste0("CSV comps ", paste(filter_components, collapse = "/"),
                          " only | \u2264 200 tracks/component"),
    colour_pal   = mouse_pal
  )
  
  filt_pdf <- paste0("component_summary_filtered_", mouse,
                     "_", timestamp_filt, ".pdf")
  pdf(filt_pdf, width = 12, height = 10)
  print(p_above)
  print(p_below)
  dev.off()
  cat("Saved:", filt_pdf, "\n")
}))

# --- 9b-vii. Save filtered RDS per mouse -------------------------------------
# MaxDis and GMM_component_csv dropped — RDS structure matches originals.
# GMM_component now reflects CSV numbering (overwritten in 9b-iii).
cat("\n")
invisible(lapply(mice_filt, function(mouse) {
  
  mouse_filt <- pooled_filt %>%
    filter(MouseID == mouse) %>%
    select(-any_of(c("MaxDis", "GMM_component_csv")))
  
  rds_out <- file.path(rds_folder,
                       paste0("accumulated_tracks_", mouse, "_filtered.rds"))
  saveRDS(mouse_filt, rds_out)
  cat("Saved:", basename(rds_out),
      "| Unique tracks:",
      nrow(distinct(mouse_filt, TrackID, ROIID_global)), "\n")
}))

cat("\n=== Section 9b complete ===\n")
cat("Filtered RDS files saved to:", rds_folder, "\n")
cat("To use in Sections 10/11: update rds_folder pattern to include '_filtered'\n")

  # =============================================================================
  # ---- SECTION 10 — STANDALONE BLOCK: ALL-MICE SUMMARY PDF----
  # =============================================================================
  # Run ONLY after all three mouse runs are complete.
  # Highlight entire block and press Ctrl+Enter.
  # Do NOT source the whole script — it will re-open the file dialog.
  #
  # Requires on disk:
  #   accumulated_tracks_M1_<timestamp>.rds
  #   accumulated_tracks_M2_<timestamp>.rds
  #   accumulated_tracks_M3_<timestamp>.rds
  # =============================================================================
  
  # --- Load RDS files ---
  # summary_track_alpha, summary_track_linewidth, min_prob_threshold, and
  # subsample_seed are read from Section 1. If running this block in isolation
  # without having run the main script, define them above this block manually.

if (!exists("summary_track_alpha"))     summary_track_alpha     <- 0.40
if (!exists("summary_track_linewidth")) summary_track_linewidth <- 0.25
if (!exists("min_prob_threshold"))      min_prob_threshold      <- 0.80
if (!exists("subsample_seed"))          subsample_seed          <- 123

  rds_folder <- "Z:\\lightmicro-inst\\scharalambous\\Resident_cell_motility_tr_v3\\ROI files\\Track Visuals - NEW FIT - GMM VEE - 5 COMPONENTS"
  
  rds_files <- list.files(
    path       = rds_folder,
    pattern    = "^accumulated_tracks_M\\d+_.*\\.rds$",
    full.names = TRUE
  )
  
  cat("RDS files found:", length(rds_files), "\n")
  print(basename(rds_files))
  
  if (length(rds_files) == 0) stop("No RDS files found — run all three mouse runs first.")
  
  pooled_all <- bind_rows(lapply(rds_files, readRDS)) %>%
    filter(!is.na(GMM_component)) %>%
    mutate(
      GMM_component = factor(GMM_component),
      MouseID       = factor(MouseID)
    )
  
  cat("Total rows:", nrow(pooled_all), "\n")
  cat("Mice:", paste(levels(pooled_all$MouseID), collapse = ", "), "\n")
  cat("Components:", paste(levels(pooled_all$GMM_component), collapse = ", "), "\n")
  cat("Unique tracks:",
      n_distinct(pooled_all %>% select(TrackID, ROIID_global)), "\n\n")
  
  # --- Helper objects ---
  
  max_axis_all <- max(abs(c(pooled_all$x_norm, pooled_all$y_norm)),
                      na.rm = TRUE) %>% ceiling()
  
  mouse_pal <- c("M1" = "#457B9D", "M2" = "#E63946", "M3" = "#F4A261", "M4" = "#2A9D8F","M5" = "#6A4C93")
  
  subsample_all <- function(df, n_max = 200) {
    df %>%
      group_by(GMM_component, MouseID) %>%
      mutate(track_sample = {
        set.seed(subsample_seed)
        TrackID %in% sample(unique(TrackID), min(n_max, n_distinct(TrackID)))
      }) %>%
      ungroup() %>%
      filter(track_sample) %>%
      select(-track_sample)
  }
  
  # --- Plot function ---
  
  make_summary_plot_all <- function(df, title_str, subtitle_str, free_axis = FALSE) {
    
    if (nrow(df) == 0) {
      return(ggplot() +
               annotate("text", x = 0.5, y = 0.5,
                        label = "No tracks", size = 4, colour = "grey40") +
               theme_void())
    }
    
    ggplot(df,
           aes(x = x_norm, y = y_norm,
               group = interaction(TrackID, ROIID_global),
               colour = MouseID)) +
      geom_path(linewidth = summary_track_linewidth,
                alpha     = summary_track_alpha) +
      geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
      geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
      annotate("point", x = 0, y = 0, colour = "black",
               shape = 3, size = 1.5, stroke = 0.6) +
      scale_colour_manual(values = mouse_pal, name = "Mouse") +
      facet_wrap(~ paste0("Comp ", GMM_component),
                 nrow   = 2, ncol = 3,
                 scales = if (free_axis) "free" else "fixed") +
      {if (free_axis)
        coord_cartesian(expand = FALSE)
        else
          coord_equal(xlim   = c(-max_axis_all, max_axis_all),
                      ylim   = c(-max_axis_all, max_axis_all),
                      expand = FALSE)} +
      labs(title    = title_str,
           subtitle = subtitle_str,
           x        = "X displacement (µm)",
           y        = "Y displacement (µm)") +
      theme_bw(base_size = 10) +
      theme(plot.title       = element_text(face = "bold", size = 12),
            plot.subtitle    = element_text(size = 9, colour = "grey40"),
            panel.grid.minor = element_blank(),
            strip.text       = element_text(face = "bold", size = 9),
            legend.position  = "right")
  }
  
  # --- Subsample and build plots ---
  
  high_sub_all <- subsample_all(
    pooled_all %>% filter(GMM_max_prob >= min_prob_threshold))
  low_sub_all  <- subsample_all(
    pooled_all %>% filter(GMM_max_prob <  min_prob_threshold))
  
  p_high_all <- make_summary_plot_all(
    high_sub_all,
    title_str    = "All mice — Normalised tracks by component",
    subtitle_str = paste0("max_prob ≥ ", min_prob_threshold,
                          " | ≤ 200 tracks/component/mouse | shared axis")
  )
  
  p_high_all_free <- make_summary_plot_all(
    high_sub_all,
    title_str    = "All mice — Normalised tracks by component (free axis)",
    subtitle_str = paste0("max_prob ≥ ", min_prob_threshold,
                          " | ≤ 200 tracks/component/mouse | free axis per component"),
    free_axis    = TRUE
  )
  
  p_low_all <- make_summary_plot_all(
    low_sub_all,
    title_str    = "All mice — Normalised tracks by component",
    subtitle_str = paste0("max_prob < ", min_prob_threshold,
                          " | ≤ 200 tracks/component/mouse | shared axis")
  )
  
  p_low_all_free <- make_summary_plot_all(
    low_sub_all,
    title_str    = "All mice — Normalised tracks by component (free axis)",
    subtitle_str = paste0("max_prob < ", min_prob_threshold,
                          " | ≤ 200 tracks/component/mouse | free axis per component"),
    free_axis    = TRUE
  )
  
  # --- Save PDF ---
  
  timestamp_summary <- format(Sys.time(), "%Y%m%d_%H%M%S")
  all_file <- paste0("component_summary_all_mice_", timestamp_summary, ".pdf")
  
  pdf(all_file, width = 12, height = 10)
  print(p_high_all)
  print(p_high_all_free)
  print(p_low_all)
  print(p_low_all_free)
  dev.off()
  
  cat("Saved:", all_file, "\n")
  cat("=== All-mice summary complete ===\n")
  
  
  
  
  # =============================================================================
  # ---- SECTION 11 — STANDALONE BLOCK: LOCATION-BASED TRACK VISUALISATION -----
  # =============================================================================
  # Run ONLY after Section 10 has completed. pooled_all, mouse_pal,
  # summary_track_alpha, summary_track_linewidth, min_prob_threshold, and
  # subsample_seed must all be present in the environment.
  #
  # Loads a separate metadata CSV (one row per ROI) containing a categorical
  # location label (e.g. "cortex", "edge", "centre"). Joins to pooled_all via
  # a composite key: Mouse + Slice + ROI, parsed from ROIID_global.
  #
  # make_location_plot() accepts a facet_var argument, so the same function
  # can later be reused to facet by "GMM_component", "Day", Treatment, or any
  # other column present in pooled_loc — change facet_var in the call.
  #
  # Outputs per run:
  #   location_tracks_all_mice_<timestamp>.pdf    — all mice, faceted by location,
  #                                                 coloured by MouseID (4 pages:
  #                                                 high/low conf × shared/free axis)
  #   location_tracks_<mouse>_<timestamp>.pdf     — one file per mouse, faceted by
  #                                                 location (4 pages each)
  # =============================================================================
  
  # --- 11a. Validate environment dependencies -----------------------------------
  
  required_objects <- c("pooled_all", "mouse_pal", "summary_track_alpha",
                        "summary_track_linewidth", "min_prob_threshold",
                        "subsample_seed", "location_metadata_csv",
                        "location_col", "meta_mouse_col",
                        "meta_slice_col", "meta_roi_col","meta_tumour_col", "meta_treatment_col")
  
  missing_obj <- required_objects[!sapply(required_objects, exists)]
  if (length(missing_obj) > 0)
    stop("Missing from environment — run Section 10 first and check Section 1 config:\n  ",
         paste(missing_obj, collapse = ", "))
  
  if (!file.exists(location_metadata_csv))
    stop("location_metadata_csv not found. Set path in Section 1.")
  
  # --- 11b. Load and validate metadata CSV -------------------------------------
  # One row per ROI. Must contain location_col, meta_mouse_col,
  # meta_slice_col, meta_roi_col as specified in Section 1.
  loc_meta <- read_csv(location_metadata_csv,
                       locale         = locale(encoding = "windows-1252"),
                       show_col_types = FALSE) %>%
    select(-starts_with("...")) %>%
    mutate(across(where(is.character),
                  ~iconv(., from = "windows-1252", to = "UTF-8", sub = "?")))
  
  required_cols <- c(location_col, meta_mouse_col, meta_slice_col, meta_roi_col)
  missing_cols  <- setdiff(required_cols, colnames(loc_meta))
  if (length(missing_cols) > 0)
    stop("Missing columns in metadata CSV: ", paste(missing_cols, collapse = ", "),
         "\n  Available columns: ", paste(colnames(loc_meta), collapse = ", "))
  
  # Coerce key columns to character for safe paste-based joining,
  # then build a composite join key: Mouse_Slice_ROI (e.g. "M1_S3_ROI1")
  loc_meta <- loc_meta %>%
    mutate(
      across(all_of(c(meta_mouse_col, meta_slice_col, meta_roi_col)),
             as.character),
      # In loc_meta (Section 11b) — build key to match ROIID_global exactly
      .loc_key = paste0(
        .data[[meta_mouse_col]],     "_",
        .data[[meta_tumour_col]],    "_",
        .data[[meta_treatment_col]], "_",
        .data[[meta_slice_col]],     "_",
        .data[[meta_roi_col]]
      )
      
    )
  
  cat("=== Section 11: Location-based track visualisation ===\n\n")
  cat("Location metadata loaded:", nrow(loc_meta), "ROIs\n")
  cat("Locations present:",
      paste(sort(unique(loc_meta[[location_col]])), collapse = ", "), "\n\n")
  
  # --- 11c. Parse join key from ROIID_global and left-join ---------------------
  # ROIID_global format: M1_T1_Mock_S3_ROI1
  # Extracts Mouse (^M\d+), Slice (S\d+), ROI (ROI\d+) via regex.
  # These three fields are pasted into the same composite key used above.
  # NOTE: values in meta_mouse_col / meta_slice_col / meta_roi_col must match
  # this extracted format exactly (e.g. "M1", "S3", "ROI1" — not "1" or "roi1").
  
  pooled_loc <- pooled_all %>%
    mutate(.loc_key = ROIID_global) %>%
    left_join(
      loc_meta %>% select(.loc_key, all_of(location_col)),
      by           = ".loc_key",
      relationship = "many-to-many"
    ) %>%
    select(-.loc_key) %>%
    mutate(across(all_of(location_col), ~ as.factor(as.character(.x))))
  
  # --- Match coverage report ---------------------------------------------------
  
  n_rois_total   <- n_distinct(pooled_loc$ROIID_global)
  n_rois_matched <- pooled_loc %>%
    filter(!is.na(.data[[location_col]])) %>%
    distinct(ROIID_global) %>%
    nrow()
  
  cat("ROI location match:", n_rois_matched, "/", n_rois_total, "ROIs matched\n")
  
  if (n_rois_total - n_rois_matched > 0) {
    unmatched_rois <- pooled_loc %>%
      filter(is.na(.data[[location_col]])) %>%
      distinct(ROIID_global) %>%
      pull(ROIID_global)
    cat("Unmatched ROIs (check meta_slice_col / meta_roi_col format):\n")
    lapply(unmatched_rois, function(r) cat("  ", r, "\n"))
  }
  
  pooled_loc <- pooled_loc %>% filter(!is.na(.data[[location_col]]))
  
  cat("Unique tracks after join:",
      nrow(distinct(pooled_loc, TrackID, ROIID_global)), "\n\n")
  
  # --- 11d. Shared axis limit --------------------------------------------------
  # Computed from the location-joined subset so the axis reflects the actual
  # data range, consistent with Section 10's max_axis_all approach.
  
  max_axis_loc <- max(abs(c(pooled_loc$x_norm, pooled_loc$y_norm)),
                      na.rm = TRUE) %>% ceiling()
  
  # --- 11e. Subsample helper ---------------------------------------------------
  # Same logic as Section 10's subsample_all(): samples up to n_max unique
  # TrackIDs within each combination of group_vars, using subsample_seed for
  # reproducibility. group_vars is a character vector of column name(s).
  
  subsample_loc <- function(df, group_vars, n_max = 200) {
    df %>%
      group_by(across(all_of(group_vars))) %>%
      mutate(track_sample = {
        set.seed(subsample_seed)
        TrackID %in% sample(unique(TrackID), min(n_max, n_distinct(TrackID)))
      }) %>%
      ungroup() %>%
      filter(track_sample) %>%
      select(-track_sample)
  }
  
  # --- 11f. Plot function -------------------------------------------------------
  # Produces a normalised-track plot with flexible faceting and colouring.
  #
  # Arguments:
  #   df           — track data (rows = detections), must contain x_norm, y_norm,
  #                  TrackID, ROIID_global, and the columns named by facet_var
  #                  and colour_var
  #   title_str    — plot title
  #   subtitle_str — plot subtitle
  #   facet_var    — column to facet panels by. Defaults to location_col, but
  #                  can be swapped to "GMM_component", "Day", etc. at the call
  #                  site without changing this function
  #   colour_var   — column mapped to the colour aesthetic (default "MouseID")
  #   colour_pal   — named character vector: names = levels of colour_var,
  #                  values = hex colour codes
  #   free_axis    — if TRUE each panel scales independently (coord_cartesian);
  #                  if FALSE all panels share a symmetric square axis
  #                  (coord_equal using max_axis_loc)
  
  make_location_plot <- function(df,
                                 title_str,
                                 subtitle_str,
                                 facet_var   = location_col,
                                 colour_var  = "MouseID",
                                 colour_pal  = mouse_pal,
                                 free_axis   = FALSE) {
    
    if (nrow(df) == 0) {
      return(
        ggplot() +
          annotate("text", x = 0.5, y = 0.5,
                   label = "No tracks", size = 4, colour = "grey40") +
          theme_void()
      )
    }
    
    # Auto-size grid: max 3 columns, rows fill as needed
    n_facets <- n_distinct(df[[facet_var]])
    ncol_loc <- min(n_facets, 3)
    nrow_loc <- ceiling(n_facets / ncol_loc)
    
    ggplot(df,
           aes(x      = x_norm,
               y      = y_norm,
               group  = interaction(TrackID, ROIID_global),
               colour = .data[[colour_var]])) +
      geom_path(linewidth = summary_track_linewidth,
                alpha     = summary_track_alpha) +
      geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
      geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
      annotate("point", x = 0, y = 0, colour = "black",
               shape = 3, size = 1.5, stroke = 0.6) +
      scale_colour_manual(values = colour_pal, name = colour_var) +
      facet_wrap(vars(.data[[facet_var]]),
                 nrow   = nrow_loc,
                 ncol   = ncol_loc,
                 scales = if (free_axis) "free" else "fixed") +
      {if (free_axis)
        coord_cartesian(expand = FALSE)
        else
          coord_equal(xlim   = c(-max_axis_loc, max_axis_loc),
                      ylim   = c(-max_axis_loc, max_axis_loc),
                      expand = FALSE)} +
      labs(title    = title_str,
           subtitle = subtitle_str,
           x        = "X displacement (\u00b5m)",
           y        = "Y displacement (\u00b5m)") +
      theme_bw(base_size = 10) +
      theme(plot.title       = element_text(face = "bold", size = 12),
            plot.subtitle    = element_text(size = 9,  colour = "grey40"),
            panel.grid.minor = element_blank(),
            strip.text       = element_text(face = "bold", size = 9),
            legend.position  = "right")
  }
  
  # --- 11g. All-mice PDF: faceted by location, coloured by MouseID -------------
  # Four pages: high confidence × {shared, free} axis
  #             low  confidence × {shared, free} axis
  # Subsampled to 200 tracks per location × mouse combination.
  
  timestamp_loc <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  high_loc_all <- subsample_loc(
    pooled_loc %>% filter(GMM_max_prob >= min_prob_threshold),
    group_vars = c(location_col, "MouseID")
  )
  low_loc_all <- subsample_loc(
    pooled_loc %>% filter(GMM_max_prob < min_prob_threshold),
    group_vars = c(location_col, "MouseID")
  )
  
  p_loc_high_all <- make_location_plot(
    high_loc_all,
    title_str    = paste0("All mice \u2014 Normalised tracks by ", location_col),
    subtitle_str = paste0("max_prob \u2265 ", min_prob_threshold,
                          " | \u2264 200 tracks/location/mouse | shared axis")
  )
  
  p_loc_high_all_free <- make_location_plot(
    high_loc_all,
    title_str    = paste0("All mice \u2014 Normalised tracks by ", location_col,
                          " (free axis)"),
    subtitle_str = paste0("max_prob \u2265 ", min_prob_threshold,
                          " | \u2264 200 tracks/location/mouse | free axis per panel"),
    free_axis    = TRUE
  )
  
  p_loc_low_all <- make_location_plot(
    low_loc_all,
    title_str    = paste0("All mice \u2014 Normalised tracks by ", location_col),
    subtitle_str = paste0("max_prob < ", min_prob_threshold,
                          " | \u2264 200 tracks/location/mouse | shared axis")
  )
  
  p_loc_low_all_free <- make_location_plot(
    low_loc_all,
    title_str    = paste0("All mice \u2014 Normalised tracks by ", location_col,
                          " (free axis)"),
    subtitle_str = paste0("max_prob < ", min_prob_threshold,
                          " | \u2264 200 tracks/location/mouse | free axis per panel"),
    free_axis    = TRUE
  )
  
  all_loc_file <- paste0("location_tracks_all_mice_", timestamp_loc, ".pdf")
  pdf(all_loc_file, width = 12, height = 10)
  print(p_loc_high_all)
  print(p_loc_high_all_free)
  print(p_loc_low_all)
  print(p_loc_low_all_free)
  dev.off()
  cat("Saved:", all_loc_file, "\n")
  
  # --- 11h. Per-mouse PDFs: faceted by location, single mouse colour -----------
  # One PDF per mouse (four pages: high/low conf × shared/free axis).
  # Colour is a single flat colour taken from mouse_pal for that mouse,
  # consistent with per-mouse plots in Section 9.
  # Subsampled to 200 tracks per location (no MouseID grouping since single mouse).
  
  mice_in_loc <- levels(droplevels(pooled_loc$MouseID))
  
  lapply(mice_in_loc, function(mouse) {
    
    mouse_loc_data <- pooled_loc %>% filter(MouseID == mouse)
    
    if (nrow(mouse_loc_data) == 0) {
      cat("  Skipping", mouse, "— no tracks after location join\n")
      return(invisible(NULL))
    }
    
    high_loc_m <- subsample_loc(
      mouse_loc_data %>% filter(GMM_max_prob >= min_prob_threshold),
      group_vars = location_col
    )
    low_loc_m <- subsample_loc(
      mouse_loc_data %>% filter(GMM_max_prob < min_prob_threshold),
      group_vars = location_col
    )
    
    # Single entry named palette so scale_colour_manual maps correctly
    mouse_single_pal <- setNames(mouse_pal[[mouse]], mouse)
    
    p_loc_high_m <- make_location_plot(
      high_loc_m,
      title_str    = paste0(mouse, " \u2014 Normalised tracks by ", location_col),
      subtitle_str = paste0("max_prob \u2265 ", min_prob_threshold,
                            " | \u2264 200 tracks/location | shared axis"),
      colour_pal   = mouse_single_pal
    )
    
    p_loc_high_m_free <- make_location_plot(
      high_loc_m,
      title_str    = paste0(mouse, " \u2014 Normalised tracks by ", location_col,
                            " (free axis)"),
      subtitle_str = paste0("max_prob \u2265 ", min_prob_threshold,
                            " | \u2264 200 tracks/location | free axis per panel"),
      colour_pal   = mouse_single_pal,
      free_axis    = TRUE
    )
    
    p_loc_low_m <- make_location_plot(
      low_loc_m,
      title_str    = paste0(mouse, " \u2014 Normalised tracks by ", location_col),
      subtitle_str = paste0("max_prob < ", min_prob_threshold,
                            " | \u2264 200 tracks/location | shared axis"),
      colour_pal   = mouse_single_pal
    )
    
    p_loc_low_m_free <- make_location_plot(
      low_loc_m,
      title_str    = paste0(mouse, " \u2014 Normalised tracks by ", location_col,
                            " (free axis)"),
      subtitle_str = paste0("max_prob < ", min_prob_threshold,
                            " | \u2264 200 tracks/location | free axis per panel"),
      colour_pal   = mouse_single_pal,
      free_axis    = TRUE
    )
    
    mouse_loc_file <- paste0("location_tracks_", mouse, "_",
                             timestamp_loc, ".pdf")
    pdf(mouse_loc_file, width = 12, height = 10)
    print(p_loc_high_m)
    print(p_loc_high_m_free)
    print(p_loc_low_m)
    print(p_loc_low_m_free)
    dev.off()
    cat("Saved:", mouse_loc_file, "\n")
  })
  
  cat("\n=== Section 11: Location-based visualisation complete ===\n")
  
  
  
  # =============================================================================
  # ---- SUB-ANALYSIS: EXTRACT INITIAL TRACK POSITIONS (t0, x, y) ----
  # =============================================================================
  # Interactively select renamed XML files. For each XML:
  #   1. Parse all particles (no duration filter — every raw track kept)
  #   2. Match positionally to the merged_tracks CSV in the same folder to get
  #      TRACK_INDEX (the TrackMate ID that bridges to tracks_clean.csv)
  #   3. Extract the FIRST detection (t0, x, y, z) of each track
  #   4. Join GMM_component and MaxDis from tracks_clean.csv via ID + ROIID_global
  #   5. Save one CSV per ROI alongside the XML: <prefix>_initial_positions.csv
  #
  # File matching: merged_tracks CSV located by prefix BEFORE "_merged_[Tt]racks".
  # t0 is the first detection time per particle (NOT necessarily 0).
  #
  # Reuses parse_xml_tracks(), load_track_ids(), parse_roiid_from_xml() and the
  # Section 1 config from track_visualisation_clean.R — run that script's
  # Sections 0-6 first so those functions and tracks_clean_csv are in scope.
  # =============================================================================
  
  # --- Config specific to this sub-analysis ------------------------------------
  # tracks_clean_csv, clean_trackid_col, clean_roi_col inherited from Section 1.
  init_pos_suffix <- "_initial_positions.csv"   # output filename suffix
  
  # --- Helper: find merged_tracks CSV by prefix before _merged_[Tt]racks -------
  find_merged_tracks_csv <- function(xml_path) {
    folder   <- dirname(xml_path)
    xml_name <- basename(xml_path)
    
    # Prefix = XML name stripped of _Tracks.xml
    prefix <- xml_name %>%
      str_remove("_Tracks\\.xml$") %>%
      str_remove("\\.xml$")
    
    csv_files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)
    
    # Match CSVs whose name (before _merged_tracks / _merged_Tracks) equals prefix
    csv_prefix <- basename(csv_files) %>%
      str_remove("\\s*-\\s*Position.*$") %>%       # strip long Position block
      str_remove("_merged_[Tt]racks.*$") %>%       # strip merged suffix
      str_remove("\\.csv$")
    
    matching <- csv_files[str_to_lower(csv_prefix) == str_to_lower(prefix) &
                            str_detect(basename(csv_files), "merged_registered_[Tt]racks")]
    
    if (length(matching) == 0) {
      cat("  WARNING: no merged_tracks CSV for", xml_name, "\n")
      return(NULL)
    }
    if (length(matching) > 1) {
      # Prefer shortest filename (XML-style over long Position-style)
      matching <- matching[which.min(nchar(basename(matching)))]
    }
    matching
  }
  
  # --- Load tracks_clean.csv once (GMM_component + MaxDis source) ---------------
  clean_lookup <- read_csv(tracks_clean_csv, show_col_types = FALSE) %>%
    select(all_of(c(clean_trackid_col, clean_roi_col,
                    "GMM_component", "MaxDis"))) %>%
    rename(TrackID      = all_of(clean_trackid_col),
           ROIID_global = all_of(clean_roi_col)) %>%
    mutate(TrackID      = as.numeric(TrackID),
           ROIID_global = as.character(ROIID_global)) %>%
    distinct(TrackID, ROIID_global, .keep_all = TRUE)
  
  cat("tracks_clean lookup loaded:", nrow(clean_lookup), "tracks\n\n")
  
  # --- Interactive file selection ----------------------------------------------
  # xml_paths <- choose.files(
  #   caption = "Select renamed XML files for initial position extraction",
  #   filters = matrix(c("XML files", "*.xml"), 1, 2, byrow = TRUE)
  # )
  # if (length(xml_paths) == 0) stop("No files selected.")
  # cat("Selected", length(xml_paths), "XML file(s)\n\n")
  
  xml_paths <- list.files(
    "Z:\\lightmicro-inst\\scharalambous\\Resident_cell_motility_tr_v3\\ROI files\\EXP112",
    pattern = "_Tracks\\.xml$", full.names = TRUE, recursive = TRUE
  )
  cat("Found", length(xml_paths), "XMLs\n")
  
  # --- Process each XML --------------------------------------------------------
  invisible(lapply(xml_paths, function(xml_path) {
    
    file_label <- basename(xml_path)
    cat("Processing:", file_label, "\n")
    
    # 1. Parse XML — keep ALL particles (no duration filter)
    parsed <- parse_xml_tracks(xml_path)
    n_xml  <- length(parsed$tracks)
    cat("  XML tracks (raw):", n_xml, "\n")
    if (n_xml == 0) { cat("  Skipped — no tracks\n\n"); return(invisible(NULL)) }
    
    # 2. Match to merged_tracks CSV for TRACK_INDEX (positional)
    csv_path <- find_merged_tracks_csv(xml_path)
    if (is.null(csv_path)) { cat("  Skipped — no CSV\n\n"); return(invisible(NULL)) }
    
    track_ids <- load_track_ids(csv_path, duration_filter_pct)
    cat("  CSV rows:", nrow(track_ids), "\n")
    
    if (nrow(track_ids) != n_xml) {
      cat("  WARNING: XML (", n_xml, ") != CSV (", nrow(track_ids),
          ") — positional match may be unreliable\n")
    }
    
    # 3. Parse ROIID_global from filename
    roi_info <- parse_roiid_from_xml(xml_path)
    if (is.null(roi_info)) { cat("  Skipped — ROIID parse failed\n\n"); return(invisible(NULL)) }
    roiid <- roi_info$roiid
    
    # 4. Extract FIRST detection per track, attach TrackID positionally
    n_match <- min(n_xml, nrow(track_ids))
    init_pos <- lapply(seq_len(n_match), function(i) {
      tr <- parsed$tracks[[i]]
      first_row <- tr[1, ]   # tracks already arranged by t in parse_xml_tracks
      tibble(
        ROIID_global   = roiid,
        TrackID        = track_ids$TrackID[i],
        particle_index = i,
        t0             = first_row$t,
        time_s0        = first_row$time_s,
        x              = first_row$x,
        y              = first_row$y,
        z              = first_row$z
      )
    }) %>% bind_rows()
    
    # 5. Join GMM_component + MaxDis from tracks_clean.csv
    init_pos <- init_pos %>%
      left_join(clean_lookup, by = c("TrackID", "ROIID_global"))
    
    n_gmm <- sum(!is.na(init_pos$GMM_component))
    cat("  GMM/MaxDis matched:", n_gmm, "/", nrow(init_pos), "tracks\n")
    
    # 6. Save one CSV per ROI alongside the XML
    prefix   <- file_label %>% str_remove("_Tracks\\.xml$") %>% str_remove("\\.xml$")
    out_path <- file.path(dirname(xml_path), paste0(prefix, init_pos_suffix))
    write_csv(init_pos, out_path)
    cat("  Saved:", basename(out_path), "\n\n")
  }))
  
  cat("=== Initial position extraction complete ===\n")