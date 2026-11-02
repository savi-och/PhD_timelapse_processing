# =============================================================================
# CELL TRACKING ANALYSIS - STEP 05: AREA-MOTILITY CORRELATION
# =============================================================================
# Purpose: Correlate source cell area (from ImageJ macro D_AREA) with
#          ROI-level motility parameters for RP-treated ROIs.
#
# Key change from original:
#   roi_summary_path points to any pre-computed ROI summary — unfiltered
#   (roi_summary_v2.csv) or GMM/MaxDis filtered (roi_summary_v2_GMM1_2.csv).
#   All outputs are suffixed with output_suffix to distinguish populations.
#
# Join strategy (two-stage fallback):
#   Primary:  match on filename prefix (everything before first " - ")
#   Fallback: match on middle block (between first " - " and "_merged_*")
#
# Motility variables: defined via motility_vars list, toggled by var_include.
#
# Correlation method: toggled by correlation_method ("pearson" or "spearman")
#
# Outputs (all suffixed with output_suffix):
#   correlation_results_<suffix>.csv
#   correlation_log_log_<suffix>.pdf
#   correlation_log_linear_<suffix>.pdf
#   correlation_results_per_mouse_<suffix>.csv
#   correlation_per_mouse_log_log_<suffix>.pdf
#   correlation_per_mouse_log_linear_<suffix>.pdf
#
# Inputs:
#   tracks_clean_v2.csv           (for Filename → ROIID_global mapping only)
#   roi_summary_v2.csv            (or filtered equivalent — set roi_summary_path)
#   Grand_source_cell_areas.csv   (from ImageJ Macro D_AREA)
# =============================================================================


# --- PARAMETERS ---------------------------------------------------------------

data_dir          <- "Z:\\lightmicro-inst\\scharalambous\\Resident_cell_motility_tr_v3\\ROI files"  # <-- UPDATE

# ROI summary to use for correlations.
# Point to any pre-computed filtered ROI summary to run on a specific population:
#   "roi_summary_v2.csv"            — full unfiltered population
#   "roi_summary_v2_GMM1_2.csv"     — GMM components 1 and 2 only
#   "roi_summary_v2_GMM1_2_MaxDis.csv" — GMM + MaxDis filtered
# The output_suffix should match the population being analysed.
roi_summary_path  <- "roi_summary_v3_GMM_1_3_MaxDis7.csv"
tracks_clean_path <- "Z:\\lightmicro-inst\\scharalambous\\Resident_cell_motility_tr_v3\\ROI files\\Complete data - retrospective fitting G=5 VEE (from v2)\\tracks_clean_v3_alldata_GMM_retrospective.csv"  # used for Filename mapping only
area_path         <- "Grand_source_cell_areas_v2_minus4581outlier.csv"

# Output suffix — update to match the population in roi_summary_path
# e.g. "all", "GMM1_2", "GMM1_2_MaxDis"
output_suffix     <- "1_MaxDis7"

# Correlation method: "pearson" or "spearman"
# Pearson:  parametric, assumes normality, appropriate for log-transformed data
# Spearman: rank-based, robust to outliers, appropriate for median values
correlation_method <- "pearson"

# Total ROI area in µm² used to standardise source cell area
ROI_AREA_TOTAL <- 435600.392

# --- Variable toggle ----------------------------------------------------------
# Mirrors var_include in scripts 02, 03, 04.
# Only variables set to TRUE and present in roi_summary are plotted.
var_include <- c(
  "Speed"        = TRUE,
  "TrackLength"  = TRUE,
  "Straightness" = TRUE,
  "MaxDis"       = TRUE,
  "TurningAngle" = TRUE,
  "Alpha"        = TRUE
)


# --- 0. PACKAGES --------------------------------------------------------------
library(tidyverse)
library(patchwork)


# --- 1. LOAD DATA -------------------------------------------------------------

setwd(data_dir)

# tracks_clean used only for Filename → ROIID_global mapping (Section 3)
tracks <- read_csv(tracks_clean_path, show_col_types = FALSE) %>%
  select(-any_of("___9")) %>%
  mutate(ROIID_global = as.character(ROIID_global))

roi_summary <- read_csv(roi_summary_path, show_col_types = FALSE) %>%
  mutate(
    ROIID_global    = factor(ROIID_global),
    TumourID_global = factor(TumourID_global),
    SliceID_global  = factor(SliceID_global),
    MouseID         = factor(MouseID),
    Treatment       = factor(Treatment, levels = c("Mock", "RP")),
    Day             = factor(Day,       levels = c("D1", "D2")),
    Location        = factor(Location,  levels = c("Core", "CoreEdge", "Edge"))
  )

area_raw <- read_csv(area_path, show_col_types = FALSE)

cat("tracks_clean rows:    ", nrow(tracks),      "\n")
cat("roi_summary rows:     ", nrow(roi_summary),  "\n")
cat("area CSV rows:        ", nrow(area_raw),     "\n")
cat("roi_summary source:   ", roi_summary_path,   "\n")
cat("output_suffix:        ", output_suffix,      "\n\n")


# --- 2. GET FILENAME → ROIID_global MAPPING -----------------------------------

filename_map <- read_csv(tracks_clean_path, show_col_types = FALSE) %>%
  select(ROIID_global, Filename) %>%
  mutate(
    ROIID_global = as.character(ROIID_global),
    Filename     = as.character(Filename)
  ) %>%
  filter(!is.na(ROIID_global), !is.na(Filename)) %>%
  distinct()

cat("Unique ROIID → Filename pairs:", nrow(filename_map), "\n")

tbl  <- table(filename_map$ROIID_global)
dups <- tbl[tbl > 1]
if (length(dups) > 0) {
  cat("WARNING: ROIIDs mapping to more than one Filename:\n")
  print(dups)
} else {
  cat("ROIID to Filename mapping is one-to-one\n")
}


# --- 3. JOIN FILENAME ONTO ROI_SUMMARY ----------------------------------------

roi <- roi_summary %>%
  left_join(filename_map, by = "ROIID_global")

n_missing <- sum(is.na(roi$Filename))
if (n_missing > 0) {
  cat("WARNING:", n_missing, "ROI(s) have no matching Filename\n")
  print(roi %>% filter(is.na(Filename)) %>% select(ROIID_global))
}


# --- 4. JOIN KEY HELPERS ------------------------------------------------------

make_keys <- function(filename, suffix_pattern) {
  stem   <- str_remove(filename, suffix_pattern)
  stem   <- trimws(stem)
  prefix <- if_else(str_detect(stem, " - "),
                    str_extract(stem, "^.*?(?= - )"), stem)
  middle <- if_else(str_detect(stem, " - "),
                    str_extract(stem, "(?<= - ).*$"), NA_character_)
  tibble(prefix_key = trimws(prefix), middle_key = trimws(middle))
}

csv_keys <- make_keys(
  roi$Filename,
  suffix_pattern = "(_| - )merged(_registered)?_tracks\\.csv$"   #_Allstats
)
roi <- roi %>%
  mutate(prefix_key_csv = csv_keys$prefix_key,
         middle_key_csv = csv_keys$middle_key)

tif_keys <- make_keys(
  area_raw$Filename,
  suffix_pattern = "(_| - )merged_registered\\.tif$"
)
area_keyed <- area_raw %>%
  rename(Filename_tif = Filename) %>%
  mutate(prefix_key_tif = tif_keys$prefix_key,
         middle_key_tif = tif_keys$middle_key)


# --- 5. JOIN AREA TO ROI_SUMMARY ----------------------------------------------

joined_primary <- roi %>%
  left_join(
    area_keyed %>% select(prefix_key_tif, middle_key_tif,
                          Filename_tif, N_Particles, Total_Area_um2),
    by = c("prefix_key_csv" = "prefix_key_tif")
  )

matched   <- joined_primary %>% filter(!is.na(Total_Area_um2))
unmatched <- joined_primary %>% filter(is.na(Total_Area_um2))

cat("\nPrimary join (prefix key):  ", nrow(matched),
    "matched |", nrow(unmatched), "unmatched\n")

if (nrow(unmatched) > 0) {
  fallback_area   <- area_keyed %>%
    filter(!is.na(middle_key_tif)) %>%
    select(middle_key_tif, Filename_tif, N_Particles, Total_Area_um2)

  joined_fallback <- unmatched %>%
    select(-Filename_tif, -N_Particles, -Total_Area_um2) %>%
    left_join(fallback_area, by = c("middle_key_csv" = "middle_key_tif"))

  newly_matched   <- joined_fallback %>% filter(!is.na(Total_Area_um2))
  still_unmatched <- joined_fallback %>% filter(is.na(Total_Area_um2))

  cat("Fallback join (middle key):", nrow(newly_matched),
      "additional matched |", nrow(still_unmatched), "still unmatched\n")

  if (nrow(still_unmatched) > 0) {
    cat("\nWARNING — unmatched ROIs:\n")
    print(still_unmatched %>%
            select(ROIID_global, Filename, prefix_key_csv, middle_key_csv))
  }

  roi_area <- bind_rows(matched, joined_fallback)
} else {
  roi_area <- matched
}

cat("Total matched:", sum(!is.na(roi_area$Total_Area_um2)),
    "/", nrow(roi_area), "ROIs\n")


# --- 6. STANDARDISE AREA ------------------------------------------------------

roi_area <- roi_area %>%
  mutate(
    Area_std     = Total_Area_um2 / ROI_AREA_TOTAL,
    log_Area_std = log(Area_std + 1e-9)
  )


# --- 7. FILTER TO RP ONLY -----------------------------------------------------

roi_rp <- roi_area %>%
  filter(
    Treatment      == "RP",
    !is.na(Total_Area_um2),
    Total_Area_um2  > 0
  )

cat("\nRP ROIs available for correlation:", nrow(roi_rp), "\n")
cat("By Day:\n");      print(table(roi_rp$Day))
cat("By Location:\n"); print(table(roi_rp$Location))


# --- 8. DEFINE MOTILITY VARIABLES ---------------------------------------------
# Column names follow roi_summary_v2.csv conventions.
# Variables toggled by var_include — only included if TRUE AND column present.
# Add or remove entries here; toggle via var_include in parameters block.
# !!! READ ME !!! - CHANGE BETWEEN MEDIAN/MEAN AND LOG/NON-LOG DEPENDING ON DISTRIBUTIONS


all_motility_vars <- list(
  list(name   = "Speed",
       col    = "Speed_mean",     sd_col = "Speed_sd",   #log_Speed_mean, log_Speed_sd (based on distribution), Speed_median, Speed_sd
       label  = "Speed (µm/min)"),
  list(name   = "TrackLength",
       col    = "log_Length_mean",    sd_col = "log_Length_sd",   #log_Length_mean, log_Length_sd (based on distribution), Length_median, Length_sd
       label  = "Log Track Length (µm)"),
  list(name   = "Straightness",
       col    = "Straightness_mean", sd_col = "Straightness_sd",
       label  = "Straightness"),
  list(name   = "MaxDis",
       col    = "log_MaxDis_mean",    sd_col = "log_MaxDis_sd",    #log_MaxDis_mean, log_MaxDis_sd
       label  = "Log Maximum Displacement (µm)"),
  list(name   = "TurningAngle",
       col    = "TurningAngle_mean", sd_col = "TurningAngle_sd",
       label  = "Mean Turning Angle (°)"),
  list(name   = "Alpha",
       col    = "Alpha30_filt_mean", sd_col = "Alpha30_filt_sd",
       label  = "Alpha (MSD, filtered)")
)

# Apply toggle and check column presence
motility_vars <- all_motility_vars[map_lgl(all_motility_vars, function(vp) {
  if (!isTRUE(var_include[vp$name])) return(FALSE)
  if (!vp$col %in% names(roi_rp)) {
    cat("NOTE: column", vp$col, "not found in roi_summary — skipping", vp$name, "\n")
    return(FALSE)
  }
  TRUE
})]

cat("\nVariables included in correlation:\n")
cat(paste(map_chr(motility_vars, "name"), collapse = ", "), "\n\n")


# --- 9. CORRELATION FUNCTION --------------------------------------------------
# Single function handling both Pearson and Spearman via correlation_method.
# Pearson  : r, 95% CI available
# Spearman : rho, CI = NA (not analytically available for rank correlation)

run_correlations <- function(data, group_label) {

  results <- list()

  for (vp in motility_vars) {
    for (area_type in c("log", "linear")) {

      area_col <- if (area_type == "log") "log_Area_std" else "Area_std"

      d <- data %>%
        filter(!is.na(.data[[vp$col]]),
               !is.na(.data[[area_col]]),
               is.finite(.data[[area_col]]))

      if (nrow(d) < 3) {
        cat("  Skipping:", vp$col, "(", group_label, area_type,
            ") — insufficient data (n =", nrow(d), ")\n")
        next
      }

      ct <- cor.test(d[[area_col]], d[[vp$col]],
                     method = correlation_method,
                     exact  = FALSE)

      stat_val <- round(ct$estimate, 3)
      names(stat_val) <- NULL

      results[[paste(vp$col, area_type, group_label, sep = "__")]] <- tibble(
        Group          = group_label,
        Variable       = vp$col,
        Variable_label = vp$label,
        Area_scale     = area_type,
        Method         = correlation_method,
        n              = nrow(d),
        stat           = stat_val,
        stat_name      = if (correlation_method == "pearson") "r" else "rho",
        p_value        = round(ct$p.value, 4),
        p_display      = if_else(ct$p.value < 0.0001, "< 0.0001",
                                 as.character(round(ct$p.value, 4))),
        CI_lower       = if (!is.null(ct$conf.int)) round(ct$conf.int[1], 3)
                         else NA_real_,
        CI_upper       = if (!is.null(ct$conf.int)) round(ct$conf.int[2], 3)
                         else NA_real_,
        Significant    = ct$p.value < 0.05
      )
    }
  }

  bind_rows(results)
}

corr_pooled <- run_correlations(roi_rp,                         "Pooled")
corr_D1     <- run_correlations(roi_rp %>% filter(Day == "D1"), "D1")
corr_D2     <- run_correlations(roi_rp %>% filter(Day == "D2"), "D2")

corr_all <- bind_rows(corr_pooled, corr_D1, corr_D2)

write_csv(corr_all, paste0("correlation_results_", output_suffix, ".csv"))
cat("\nSaved: correlation_results_", output_suffix, ".csv\n", sep = "")

cat("\n--- Correlation Summary ---\n")
print(
  corr_all %>% select(Group, Variable, Area_scale, Method,
                      n, stat, stat_name, p_display, Significant),
  n = Inf
)


# --- 10. COLOUR CONFIGURATION -------------------------------------------------

day_colours <- c("D1" = "#457B9D", "D2" = "#E63946")

# 5-mouse shape palette — extend if more mice added
mouse_shapes <- c("M1" = 16, "M2" = 17, "M3" = 15, "M4" = 18, "M5" = 8)


# --- 11. SCATTER PLOT FUNCTION ------------------------------------------------

make_scatter <- function(data, vp, area_col, area_label,
                         corr_table, group_label, day_col = NULL) {

  area_type <- if_else(str_detect(area_col, "log"), "log", "linear")
  corr_row  <- corr_table %>%
    filter(Group == group_label,
           Variable == vp$col,
           Area_scale == area_type)

  stat_name <- if (correlation_method == "pearson") "r" else "rho"
  r_label   <- if (nrow(corr_row) == 1) {
    paste0(stat_name, " = ", corr_row$stat,
           ",  p ",
           if_else(corr_row$p_display == "< 0.0001",
                   "< 0.0001",
                   paste0("= ", corr_row$p_display)),
           "  (n = ", corr_row$n, ")")
  } else { "" }

  plot_data <- data %>%
    filter(!is.na(.data[[vp$col]]),
           !is.na(.data[[area_col]]),
           is.finite(.data[[area_col]]))

  p <- ggplot(plot_data, aes(x = .data[[area_col]], y = .data[[vp$col]]))

  if (vp$sd_col %in% names(plot_data)) {
    p <- p +
      geom_errorbar(
        aes(ymin = .data[[vp$col]] - .data[[vp$sd_col]],
            ymax = .data[[vp$col]] + .data[[vp$sd_col]]),
        # aes(ymin = pmax(0, .data[[vp$col]] - .data[[vp$sd_col]]),
        #     ymax = .data[[vp$col]] + .data[[vp$sd_col]]),
         width = 0, colour = "grey70", linewidth = 0.35, alpha = 0.7
      )
  }

  show_line <- nrow(corr_row) == 1 && isTRUE(corr_row$Significant)
  if (show_line) {
    p <- p +
      geom_smooth(method = "lm", formula = y ~ x,
                  colour = "black", linewidth = 0.8,
                  fill = "grey80", alpha = 0.3)
  }

  if (!is.null(day_col)) {
    p <- p +
      geom_point(aes(colour = .data[[day_col]], shape = MouseID),
                 size = 2.5, alpha = 0.85) +
      scale_colour_manual(values = day_colours, name = "Day") +
      scale_shape_manual(values  = mouse_shapes, name = "Mouse")
  } else {
    p <- p +
      geom_point(aes(shape = MouseID),
                 size = 2.5, alpha = 0.85, colour = "#2C7BB6") +
      scale_shape_manual(values = mouse_shapes, name = "Mouse")
  }

  p <- p +
    annotate("text",
             x = Inf, y = -Inf,
             label    = r_label,
             hjust    = 1.05, vjust = -0.5,
             size     = 3.2, fontface = "italic",
             colour   = if (nrow(corr_row) == 1 &&
                            isTRUE(corr_row$Significant)) "black" else "grey50") +
    labs(
      title    = paste0(group_label, "  —  ", vp$label),
      subtitle = if_else(show_line,
                         "Line of best fit shown (p < 0.05)",
                         "No line fitted (p ≥ 0.05)"),
      x = area_label,
      y = vp$label
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold",    size = 11),
      plot.subtitle    = element_text(colour = "grey50", size = 9),
      panel.grid.minor = element_blank(),
      legend.position  = if (!is.null(day_col)) "right" else "none"
    )

  p
}


# --- 12. BUILD AND SAVE PLOTS -------------------------------------------------

groups <- list(
  list(label   = "Pooled",
       data    = roi_rp,
       day_col = "Day"),
  list(label   = "D1",
       data    = roi_rp %>% filter(Day == "D1"),
       day_col = NULL),
  list(label   = "D2",
       data    = roi_rp %>% filter(Day == "D2"),
       day_col = NULL)
)

pdf(paste0("correlation_log_log_", output_suffix, ".pdf"), width = 7, height = 5)
for (grp in groups) {
  for (vp in motility_vars) {
    print(make_scatter(
      data        = grp$data, vp = vp,
      area_col    = "log_Area_std",
      area_label  = "Log Standardised Source Cell Area",
      corr_table  = corr_all,
      group_label = grp$label,
      day_col     = grp$day_col
    ))
  }
}
dev.off()
cat("Saved: correlation_log_log_", output_suffix, ".pdf\n", sep = "")

pdf(paste0("correlation_log_linear_", output_suffix, ".pdf"), width = 7, height = 5)
for (grp in groups) {
  for (vp in motility_vars) {
    print(make_scatter(
      data        = grp$data, vp = vp,
      area_col    = "Area_std",
      area_label  = "Standardised Source Cell Area (fraction of ROI)",
      corr_table  = corr_all,
      group_label = grp$label,
      day_col     = grp$day_col
    ))
  }
}
dev.off()
cat("Saved: correlation_log_linear_", output_suffix, ".pdf\n", sep = "")

cat("\nAll pooled outputs saved.\n")


# --- 13. PER-MOUSE CORRELATIONS -----------------------------------------------

mouse_ids <- levels(factor(roi_rp$MouseID))

corr_mouse <- map_dfr(mouse_ids, function(m) {
  run_correlations(
    data        = roi_rp %>% filter(MouseID == m),
    group_label = m
  )
})

write_csv(corr_mouse,
          paste0("correlation_results_per_mouse_", output_suffix, ".csv"))
cat("\nSaved: correlation_results_per_mouse_", output_suffix, ".csv\n", sep = "")

cat("\n--- Per-Mouse Correlation Summary ---\n")
print(
  corr_mouse %>% select(Group, Variable, Area_scale, n,
                        stat, stat_name, p_display, Significant),
  n = Inf
)


# --- 14. PER-MOUSE SCATTER PLOT FUNCTION --------------------------------------

make_scatter_mouse <- function(vp, area_col, area_label, corr_table) {

  plots <- map(mouse_ids, function(m) {

    d <- roi_rp %>%
      filter(MouseID == m,
             !is.na(.data[[vp$col]]),
             !is.na(.data[[area_col]]),
             is.finite(.data[[area_col]]))

    corr_row <- corr_table %>%
      filter(Group == m, Variable == vp$col,
             Area_scale == if_else(str_detect(area_col, "log"), "log", "linear"))

    stat_name <- if (correlation_method == "pearson") "r" else "rho"
    r_label   <- if (nrow(corr_row) == 1) {
      paste0(stat_name, " = ", corr_row$stat,
             ",  p ",
             if_else(corr_row$p_display == "< 0.0001", "< 0.0001",
                     paste0("= ", corr_row$p_display)),
             "  (n = ", corr_row$n, ")")
    } else { "" }

    p <- ggplot(d, aes(x = .data[[area_col]], y = .data[[vp$col]]))

    if (vp$sd_col %in% names(d)) {
      p <- p +
        geom_errorbar(
          aes(ymin = .data[[vp$col]] - .data[[vp$sd_col]],
              ymax = .data[[vp$col]] + .data[[vp$sd_col]]),
          # aes(ymin = pmax(0, .data[[vp$col]] - .data[[vp$sd_col]]),
          #     ymax = .data[[vp$col]] + .data[[vp$sd_col]]),
          width = 0, colour = "grey70", linewidth = 0.35, alpha = 0.7
        )
    }

    n_finite <- sum(is.finite(d[[area_col]]) & is.finite(d[[vp$col]]))
    if (n_finite >= 3) {
      p <- p +
        tryCatch(
          geom_smooth(method = "lm", formula = y ~ x,
                      colour = "black", linewidth = 0.8, se = FALSE),
          error = function(e) {
            message("  geom_smooth failed for ", m, ": ", e$message); NULL
          }
        )
    }

    p <- p +
      geom_point(aes(colour = Day), size = 2.5, alpha = 0.85) +
      scale_colour_manual(values = day_colours, name = "Day") +
      annotate("text",
               x = Inf, y = -Inf, label = r_label,
               hjust = 1.05, vjust = -0.5,
               size = 3, fontface = "italic", colour = "grey40") +
      labs(
        title    = m,
        subtitle = paste0("n = ", nrow(d), " ROIs  (exploratory only)"),
        x        = area_label,
        y        = vp$label
      ) +
      theme_bw(base_size = 11) +
      theme(
        plot.title       = element_text(face = "bold", size = 11),
        plot.subtitle    = element_text(colour = "grey50", size = 8),
        panel.grid.minor = element_blank(),
        legend.position  = "right"
      )

    p
  })

  wrap_plots(plots, nrow = 2) +
    plot_annotation(
      title   = paste0(vp$label,
                       "  —  Per-mouse (RP only)  |  suffix: ", output_suffix),
      caption = "Line is directional only — interpret p-values cautiously given small n per mouse",
      theme   = theme(
        plot.title   = element_text(face = "bold", size = 12),
        plot.caption = element_text(colour = "grey50", size = 8)
      )
    )
}


# --- 15. SAVE PER-MOUSE PLOTS -------------------------------------------------

pdf(paste0("correlation_per_mouse_log_log_", output_suffix, ".pdf"),
    width = 14, height = 9)
for (vp in motility_vars) {
  print(make_scatter_mouse(
    vp         = vp,
    area_col   = "log_Area_std",
    area_label = "Log Standardised Source Cell Area",
    corr_table = corr_mouse
  ))
}
dev.off()
cat("Saved: correlation_per_mouse_log_log_", output_suffix, ".pdf\n", sep = "")

pdf(paste0("correlation_per_mouse_log_linear_", output_suffix, ".pdf"),
    width = 14, height = 9)
for (vp in motility_vars) {
  print(make_scatter_mouse(
    vp         = vp,
    area_col   = "Area_std",
    area_label = "Standardised Source Cell Area (fraction of ROI)",
    corr_table = corr_mouse
  ))
}
dev.off()
cat("Saved: correlation_per_mouse_log_linear_", output_suffix, ".pdf\n", sep = "")

cat("\nAll outputs saved.\n")
cat("To run on a filtered population: update roi_summary_path and output_suffix\n")
cat("in the parameters block and re-run.\n")
