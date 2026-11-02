# Extracting Motility metric statistics from TrackMate (CellposeSAM)
## !!!! REMEMBER TO RENAME YOUR FILES WITH THE ..._ROI_... NAME !!!!
## !!!!CHANGE DURATION FOR 6 FRAMES (FOR 61 FRAMES)- TO CAPTURE ONLY TRACKST GREATER THAN 10% OF THE MOVIE DURATION

# ----- clear environment --------
rm(list=ls())
gc()

# Library upload
library(tidyverse)

######## 1) COLLECT ALL THE DATA NEEDED TO COMPARE MOTILITY METRICS #######

# Set working directory as the folder containing all the input CSV files
setwd("Z:\\lightmicro-inst\\scharalambous\\EXP25N107\\MIPs\\Registered\\Track Analysis (New_and_oldset_jun)\\D2\\MOCK") 
wdir <- "Z:\\lightmicro-inst\\scharalambous\\EXP25N107\\MIPs\\Registered\\Track Analysis (New_and_oldset_jun)\\D2\\MOCK" 

# Get all Statistics.csv files
stat_files <- list.files(path = wdir, pattern = "_merged_tracks\\.csv$", full.names = TRUE)  ### change name according to file naming convention
cat("Found", length(stat_files), "Statistics files\n")

####### ---- FOR PER ROI based analysis ---- 
# Get all Statistics.csv files
stat_files <- list.files(path = wdir, pattern = "_merged_tracks_Allstats\\.csv$", full.names = TRUE)
cat("Found", length(stat_files), "Statistics files\n")

#######

# Check if any files were found
if (length(stat_files) == 0) {
  stop("No *_Statistics.csv files found in the working directory")
}

# Print the files that will be processed
print(stat_files)

######################################

# Iteration for each Statistics file
for (filepath in stat_files) {
  
  cat("\nProcessing file:", basename(filepath), "\n")
  
  # Read the CSV file (skip first 4 rows, use row 1 as headers)
  raw <- read.csv(filepath, skip = 4, header = FALSE, fill = TRUE, check.names = TRUE)
  
  # Read the header separately from row 1
  headers <- read.csv(filepath, nrows = 1, header = FALSE, stringsAsFactors = FALSE)
  colnames(raw) <- as.character(headers[1, ])
  
  # Define required columns by name
  required_cols <- c("TRACK_DURATION","TRACK_INDEX", "TRACK_MEAN_SPEED", "TOTAL_DISTANCE_TRAVELED", 
                     "MAX_DISTANCE_TRAVELED", "CONFINEMENT_RATIO")
  
  # Check if all required columns exist
  if (!all(required_cols %in% colnames(raw))) {
    cat("Warning: Not all required columns found in", basename(filepath), "\n")
    cat("Available columns:", colnames(raw), "\n")
    next
  }
  
  # Extract only the required columns
  Allstats <- raw[, required_cols, drop = FALSE]
  
  # Rename columns to match your desired output format
  colnames(Allstats) <- c("Track.Duration","ID", "Track.Speed.Mean", "Track.Length", "MaxDis", "Track.Straightness")
  
  # Reorder columns to: ID, Track.Length, Track.Speed.Mean, Track.Straightness, MaxDis
  Allstats <- Allstats[, c("ID", "Track.Length", "Track.Speed.Mean", "Track.Straightness", "MaxDis", "Track.Duration")]
  
  # Filter out rows with missing or blank TrackIDs and Tracks with duration lower than 6 time frames (Time frame = 30-35 sec)
  Allstats <- Allstats[!is.na(Allstats$ID) & Allstats$ID != "", , drop = FALSE]
  # Dynamically filter: keep tracks > 10% of the max Track.Duration in this file
  duration_threshold <- 0.10 * max(Allstats$Track.Duration, na.rm = TRUE)
  Allstats <- Allstats[Allstats$Track.Duration > duration_threshold, , drop = FALSE]
  # Arrange by TrackID in order
  Allstats <- Allstats[order(Allstats$ID), , drop = FALSE]
  
  # Check the output
  cat("Processed", nrow(Allstats), "tracks\n")
  print(head(Allstats))
  
  # Create output filename
  base_name <- tools::file_path_sans_ext(basename(filepath))
  output_filename <- paste0(base_name, "_Allstats.csv")
  output_path <- file.path(wdir, output_filename)
  
  # Write the output CSV file
  write.csv(Allstats, output_path, row.names = FALSE)
  cat("Saved:", output_filename, "\n")
  
  cat("Processed", nrow(Allstats), "tracks (duration threshold:", round(duration_threshold, 1), "s)\n")
}

cat("\nAll files processed successfully!\n")

######### 2) MERGE ALLSTATS CSVS WITH COMMON PREFIX ##########

# 1) Collect only Allstats CSVs (skip already-merged outputs)
files <- list.files(wdir, pattern = "Allstats\\.csv$", full.names = TRUE)
files <- files[!grepl("\\.merged\\.csv$", files, ignore.case = TRUE)]

# 2) Build a grouping key: everything before `_ROI_<number>`
base_names <- basename(files)
keys <- sub("_ROI_\\d+.*$", "", base_names)   # e.g. "slice1", "slice2", RP, MOCK

# 3) Split file paths by key
groups <- split(files, keys)

# 4) For each group, read all files and row-bind them
for (k in names(groups)) {
  dfs <- lapply(groups[[k]], function(f) read.csv(f, stringsAsFactors = FALSE))
  
  # rbind assumes same columns; if columns might differ, use the dplyr variant below
  merged <- do.call(rbind, dfs)
  
  # 5) Write one merged file per group
  out <- file.path(wdir, paste0(k, "_Statistics_Allstats.merged.csv"))
  write.csv(merged, out, row.names = FALSE)
  
  # keep it in the workspace as k.merged
  assign(paste0(k, ".merged"), merged, envir = .GlobalEnv)
}


##############
# Separate four populations: 
# - Static cells (max displacement<7um, track length<13um),
# - Wobbling cells (max displacement<7um, track length>=13um),
# - Migrating cells (7um<=max displacement<25um, track length<40um)
# - Long migrating cells (max displacement>=25um, track length>=40um)

## Get the merged Allstats files and read them
files_Allstats <- list.files(wdir, pattern = "Allstats.merged.csv", full.names = FALSE) # Get all merged Allstats.merged files
files_Allstats
ldf_Allstats <- sapply(files_Allstats, read.csv, simplify = FALSE, USE.NAMES = TRUE) # Read all merged Allstats files and keep their names in the list
ldf_Allstats # Check the output
length(ldf_Allstats)

## Subset static, wobbling, migrating and long migrating populations, and iterate this for each sample (each merged Allstats.merged file)
for (i in 1:length(ldf_Allstats)){
  st <- ldf_Allstats[[i]][(ldf_Allstats[[i]]$MaxDis<7) & (ldf_Allstats[[i]]$Track.Length<13), ] # subset static cells
  colnames(st) <- c("ID", "Track.Length", "Track.Speed.Mean", "Track.Straightness", "MaxDis","Track.Duration")
  wob <- ldf_Allstats[[i]][(ldf_Allstats[[i]]$MaxDis<7) & (ldf_Allstats[[i]]$Track.Length>=13),] # subset wobbling cells
  colnames(wob) <- c("ID", "Track.Length", "Track.Speed.Mean", "Track.Straightness", "MaxDis","Track.Duration")
  mig <- ldf_Allstats[[i]][(ldf_Allstats[[i]]$MaxDis>=7) & ((ldf_Allstats[[i]]$MaxDis<25) | (ldf_Allstats[[i]]$Track.Length<40)),] # Subset migrating cells
  colnames(mig) <- c("ID", "Track.Length", "Track.Speed.Mean", "Track.Straightness", "MaxDis","Track.Duration")
  mig2 <- ldf_Allstats[[i]][(ldf_Allstats[[i]]$MaxDis>=7), ] # Subset migrating + long migrating cells together
  colnames(mig2) <- c("ID", "Track.Length", "Track.Speed.Mean", "Track.Straightness", "MaxDis","Track.Duration")
  Lmig <- ldf_Allstats[[i]][(ldf_Allstats[[i]]$MaxDis>=25) & (ldf_Allstats[[i]]$Track.Length>=40),] # Subset long migrating cells
  colnames(Lmig) <- c("ID", "Track.Length", "Track.Speed.Mean", "Track.Straightness", "MaxDis","Track.Duration")
  Prop.types <- data.frame((nrow(st)*100/nrow(ldf_Allstats[[i]])), (nrow(wob)*100/nrow(ldf_Allstats[[i]])), (nrow(mig)*100/nrow(ldf_Allstats[[i]])), (nrow(Lmig)*100/nrow(ldf_Allstats[[i]]))) # Calculate the porportion of each behaviour category
  colnames(Prop.types) <- c("St", "Wob", "Mig", "Lmig")
  
  # Set proper file names
  name.st <- gsub("_Allstats.merged.csv", "_st", names(ldf_Allstats[i]))
  name.wob <- gsub("_Allstats.merged.csv", "_wob", names(ldf_Allstats[i]))
  name.mig <- gsub("_Allstats.merged.csv", "_mig", names(ldf_Allstats[i]))
  name.mig2 <- gsub("_Allstats.merged.csv", "_mig2", names(ldf_Allstats[i]))
  name.Lmig <- gsub("_Allstats.merged.csv", "_Lmig", names(ldf_Allstats[i]))
  name.Prop.types <- gsub("_Allstats.merged.csv", "_Proportions", names(ldf_Allstats[i]))
  
  # Write .csv output for each behaviour category
  write.csv(st, paste(name.st, ".csv", sep = ""), row.names = FALSE)
  write.csv(wob, paste(name.wob, ".csv", sep = ""), row.names = FALSE)
  write.csv(mig, paste(name.mig, ".csv", sep = ""), row.names = FALSE)
  write.csv(mig2, paste(name.mig2, ".csv", sep = ""), row.names = FALSE)
  write.csv(Lmig, paste(name.Lmig, ".csv", sep = ""), row.names = FALSE)
  write.csv(Prop.types, paste(name.Prop.types, ".csv", sep = ""), row.names = FALSE)
  assign(name.st, st)
  assign(name.wob, wob)
  assign(name.mig, mig)
  assign(name.mig2, mig2)
  assign(name.Lmig, Lmig)
  assign(name.Prop.types, Prop.types)
}


# Group all samples together into one file for each statistic of interest    ######
## Create a list of data frames containing all statistics of interest for all samples
files.all <- list.files(wdir, pattern = "Allstats\\.merged\\.csv$", full.names = FALSE) # Get all Allstats.merged files  #### CHANGE here for ROI-based analysis FROM Allstats\\.merged\\.csv$ to Allstats\\.csv$
files.all

ldf.all <- sapply(files.all, read.csv, simplify = FALSE, USE.NAMES = TRUE) # Read the files and store output in a list
ldf.all <- lapply(ldf.all, as.data.frame, col.names = c("ID", "Track.Length", "Track.Speed.Mean", "Track.Straightness", "MaxDis","Track.Duration")) # Ensure each element of the list is a data frame and assign clear column names

# Step 1: Clean data - remove NAs and zeros
ldf.all <- lapply(ldf.all, function(df) {
  df <- df[complete.cases(df), , drop = FALSE]
  df <- df[df$Track.Length != 0, , drop = FALSE]
  df <- df[df$Track.Speed.Mean != 0, , drop = FALSE]
  df <- df[df$Track.Straightness != 0, , drop = FALSE]
  df
})

# Step 2: Calculate maximum number of rows
max_len <- max(vapply(ldf.all, nrow, integer(1))) # nrow is a function that returns the number of rows

# Step 3: Create helper function
create_padded_df <- function(stat_column_name) {
  padded_list <- lapply(seq_along(ldf.all), function(i) {
    v <- ldf.all[[i]][[stat_column_name]]
    length(v) <- max_len  # Pads with NA
    return(v)
  })
  df <- as.data.frame(padded_list)
  colnames(df) <- sub("_Allstats\\.merged\\.csv$", "", names(ldf.all))
  return(df)
}

# Step 4: Create all padded data frames
all.Length <- create_padded_df("Track.Length")
all.Speed.Mean <- create_padded_df("Track.Speed.Mean")
all.Straightness <- create_padded_df("Track.Straightness")

# Check outputs
head(all.Length)
head(all.Speed.Mean)
head(all.Straightness)

#######################################################################
####### ---- ROI based analysis  ---- #-> ONLY AFTER SLICE ANALYSIS
#######################################################################
#### Calculate means and median of indivual ROIS from Allstats files #### !!! DO THIS ONCE YOU COMPLETE THE PER SLICE ANALYSIS
# Create a list of data frames containing all statistics of interest for all samples
files.all <- list.files(wdir, pattern = "Allstats\\.csv$", full.names = FALSE) # Get all Allstats.merged files
files.all

ldf.all <- sapply(files.all, read.csv, simplify = FALSE, USE.NAMES = TRUE) # Read the files and store output in a list
ldf.all <- lapply(ldf.all, as.data.frame, col.names = c("ID", "Track.Length", "Track.Speed.Mean", "Track.Straightness", "MaxDis","Track.Duration")) # Ensure each element of the list is a data frame and assign clear column names

# Step 1: Clean data - remove NAs and zeros
ldf.all <- lapply(ldf.all, function(df) {
  df <- df[complete.cases(df), , drop = FALSE]
  df <- df[df$Track.Length != 0, , drop = FALSE]
  df <- df[df$Track.Speed.Mean != 0, , drop = FALSE]
  df <- df[df$Track.Straightness != 0, , drop = FALSE]
  df
})

# Step 2: Calculate maximum number of rows
max_len <- max(vapply(ldf.all, nrow, integer(1))) # nrow is a function that returns the number of rows

# Step 3: Create helper function
create_padded_df <- function(stat_column_name) {
  padded_list <- lapply(seq_along(ldf.all), function(i) {
    v <- ldf.all[[i]][[stat_column_name]]
    length(v) <- max_len  # Pads with NA
    return(v)
  })
  df <- as.data.frame(padded_list)
  colnames(df) <- sub("_Allstats\\.merged\\.csv$", "", names(ldf.all))
  return(df)
}

# Step 4: Create all padded data frames
all.Length <- create_padded_df("Track.Length")
all.Speed.Mean <- create_padded_df("Track.Speed.Mean")
all.Straightness <- create_padded_df("Track.Straightness")

# Check outputs
head(all.Length)
head(all.Speed.Mean)
head(all.Straightness)
###################################################################################################

##############  Individual populations
############################### Static ###########################

files.st <- list.files(wdir, pattern = "*st.csv", full.names = FALSE)
ldf.st <- sapply(files.st, read.csv, simplify = FALSE, USE.NAMES = TRUE)
ldf.st <- lapply(ldf.st, as.data.frame, col.names = c("ID", "Track.Length", "Track.Speed.Mean", "Track.Straightness", "MaxDis","Track.Duration"))

# Step 1: Clean data - remove NAs and zeros
ldf.st <- lapply(ldf.st, function(df) {
  df <- df[complete.cases(df), , drop = FALSE]
  df <- df[df$Track.Length != 0, , drop = FALSE]
  df <- df[df$Track.Speed.Mean != 0, , drop = FALSE]
  df <- df[df$Track.Straightness != 0, , drop = FALSE]
  df
})

# Step 2: Calculate maximum number of rows
max_len <- max(vapply(ldf.st, nrow, integer(1))) # nrow is a function that returns the number of rows

# Step 3: Create helper function
create_padded_df <- function(stat_column_name) {
  padded_list <- lapply(seq_along(ldf.st), function(i) {
    v <- ldf.st[[i]][[stat_column_name]]
    length(v) <- max_len  # Pads with NA
    return(v)
  })
  df <- as.data.frame(padded_list)
  colnames(df) <- sub("_st.csv$", "", names(ldf.st))
  return(df)
}

# Step 4: Create all padded data frames
st.Length <- create_padded_df("Track.Length")
st.Speed.Mean <- create_padded_df("Track.Speed.Mean")
st.Straightness <- create_padded_df("Track.Straightness")

# Check outputs
head(st.Length)
head(st.Speed.Mean)
head(st.Straightness)


############################### Wobbling ###########################

files.wob <- list.files(wdir, pattern = "*wob.csv", full.names = FALSE)
ldf.wob <- sapply(files.wob, read.csv, simplify = FALSE, USE.NAMES = TRUE)
ldf.wob <- lapply(ldf.wob, as.data.frame, col.names = c("ID", "Track.Length", "Track.Speed.Mean", "Track.Straightness", "MaxDis","Track.Duration"))

# Step 1: Clean data - remove NAs and zeros
ldf.wob <- lapply(ldf.st, function(df) {
  df <- df[complete.cases(df), , drop = FALSE]
  df <- df[df$Track.Length != 0, , drop = FALSE]
  df <- df[df$Track.Speed.Mean != 0, , drop = FALSE]
  df <- df[df$Track.Straightness != 0, , drop = FALSE]
  df
})

# Step 2: Calculate maximum number of rows
max_len <- max(vapply(ldf.wob, nrow, integer(1))) # nrow is a function that returns the number of rows

# Step 3: Create helper function
create_padded_df <- function(stat_column_name) {
  padded_list <- lapply(seq_along(ldf.wob), function(i) {
    v <- ldf.wob[[i]][[stat_column_name]]
    length(v) <- max_len  # Pads with NA
    return(v)
  })
  df <- as.data.frame(padded_list)
  colnames(df) <- sub("_wob.csv$", "", names(ldf.wob))
  return(df)
}

# Step 4: Create all padded data frames
wob.Length <- create_padded_df("Track.Length")
wob.Speed.Mean <- create_padded_df("Track.Speed.Mean")
wob.Straightness <- create_padded_df("Track.Straightness")

# Check outputs
head(wob.Length)
head(wob.Speed.Mean)
head(wob.Straightness)

######################################################## Migrating ########################


## Create a data frame containing the Length statistic for only migrating cells for all samples.
## ---------------- LOAD & CLEAN BASE ---------------- ##
files.mig <- list.files(wdir, pattern = "_mig\\.csv$", full.names = FALSE)
ldf.mig <- sapply(files.mig, read.csv, simplify = FALSE, USE.NAMES = TRUE)
ldf.mig <- lapply(ldf.mig, as.data.frame, col.names = c("ID", "Track.Length", "Track.Speed.Mean", "Track.Straightness", "MaxDis","Track.Duration"))

# Step 1: Clean data - remove NAs and zeros
ldf.mig <- lapply(ldf.mig, function(df) {
  df <- df[complete.cases(df), , drop = FALSE]
  df <- df[df$Track.Length != 0, , drop = FALSE]
  df <- df[df$Track.Speed.Mean != 0, , drop = FALSE]
  df <- df[df$Track.Straightness != 0, , drop = FALSE]
  df
})

# Step 2: Calculate maximum number of rows
max_len <- max(vapply(ldf.mig, nrow, integer(1))) # nrow is a function that returns the number of rows

# Step 3: Create helper function
create_padded_df <- function(stat_column_name) {
  padded_list <- lapply(seq_along(ldf.mig), function(i) {
    v <- ldf.mig[[i]][[stat_column_name]]
    length(v) <- max_len  # Pads with NA
    return(v)
  })
  df <- as.data.frame(padded_list)
  colnames(df) <- sub("_mig.csv$", "", names(ldf.mig))
  return(df)
}

# Step 4: Create all padded data frames
mig.Length <- create_padded_df("Track.Length")
mig.Speed.Mean <- create_padded_df("Track.Speed.Mean")
mig.Straightness <- create_padded_df("Track.Straightness")

# Check outputs
head(mig.Length)
head(mig.Speed.Mean)
head(mig.Straightness)


####################################################  Migrating + long Migrating #######################

## Create a data frame containing the Length statistic for only migrating + Long migrating cells for all samples.
files.mig2 <- list.files(wdir, pattern = "_mig2\\.csv$", full.names = FALSE)
ldf.mig2 <- sapply(files.mig2, read.csv, simplify = FALSE, USE.NAMES = TRUE)
ldf.mig2 <- lapply(ldf.mig2, as.data.frame, col.names = c("ID", "Track.Length", "Track.Speed.Mean", "Track.Straightness", "MaxDis","Track.Duration"))

# Step 1: Clean data - remove NAs and zeros
ldf.mig2 <- lapply(ldf.mig2, function(df) {
  df <- df[complete.cases(df), , drop = FALSE]
  df <- df[df$Track.Length != 0, , drop = FALSE]
  df <- df[df$Track.Speed.Mean != 0, , drop = FALSE]
  df <- df[df$Track.Straightness != 0, , drop = FALSE]
  df
})

# Step 2: Calculate maximum number of rows
max_len <- max(vapply(ldf.mig2, nrow, integer(1))) # nrow is a function that returns the number of rows

# Step 3: Create helper function
create_padded_df <- function(stat_column_name) {
  padded_list <- lapply(seq_along(ldf.mig2), function(i) {
    v <- ldf.mig2[[i]][[stat_column_name]]
    length(v) <- max_len  # Pads with NA
    return(v)
  })
  df <- as.data.frame(padded_list)
  colnames(df) <- sub("_mig2.csv$", "", names(ldf.mig2))
  return(df)
}

# Step 4: Create all padded data frames
mig2.Length <- create_padded_df("Track.Length")
mig2.Speed.Mean <- create_padded_df("Track.Speed.Mean")
mig2.Straightness <- create_padded_df("Track.Straightness")

# Check outputs
head(mig2.Length)
head(mig2.Speed.Mean)
head(mig2.Straightness)

###############################################  Long Migrating ########################################

## Create a data frame containing the Length statistic for only Long migrating cells for all samples.

files.Lmig <- list.files(wdir, pattern = "_Lmig\\.csv$", full.names = FALSE)
ldf.Lmig <- sapply(files.Lmig, read.csv, simplify = FALSE, USE.NAMES = TRUE)
ldf.Lmig <- lapply(ldf.Lmig, as.data.frame, col.names = c("ID", "Track.Length", "Track.Speed.Mean", "Track.Straightness", "MaxDis","Track.Duration"))


# Step 1: Clean data - remove NAs and zeros
ldf.Lmig <- lapply(ldf.Lmig, function(df) {
  df <- df[complete.cases(df), , drop = FALSE]
  df <- df[df$Track.Length != 0, , drop = FALSE]
  df <- df[df$Track.Speed.Mean != 0, , drop = FALSE]
  df <- df[df$Track.Straightness != 0, , drop = FALSE]
  df
})

# Step 2: Calculate maximum number of rows
max_len <- max(vapply(ldf.Lmig, nrow, integer(1))) # nrow is a function that returns the number of rows

# Step 3: Create helper function
create_padded_df <- function(stat_column_name) {
  padded_list <- lapply(seq_along(ldf.Lmig), function(i) {
    v <- ldf.Lmig[[i]][[stat_column_name]]
    length(v) <- max_len  # Pads with NA
    return(v)
  })
  df <- as.data.frame(padded_list)
  colnames(df) <- sub("_Lmig.csv$", "", names(ldf.Lmig))
  return(df)
}

# Step 4: Create all padded data frames
Lmig.Length <- create_padded_df("Track.Length")
Lmig.Speed.Mean <- create_padded_df("Track.Speed.Mean")
Lmig.Straightness <- create_padded_df("Track.Straightness")

# Check outputs
head(Lmig.Length)
head(Lmig.Speed.Mean)
head(Lmig.Straightness)

###################### Save CSVs for individual categories ##############

###################### Helper Functions ##############

# Function to calculate summary statistics (mean and median)
calculate_summary <- function(df) {
  data.frame(
    Mean = colMeans(df, na.rm = TRUE),
    Median = apply(df, 2, median, na.rm = TRUE)
  )
}

# Function to save a data frame as CSV
save_csv <- function(df, name) {
  write.csv(df, paste0(name, ".csv"), row.names = FALSE)
}

###################### Save Raw Data CSVs ##############

# List of all raw data frames
raw_data_list <- list(
  all.Length = all.Length,
  all.Speed.Mean = all.Speed.Mean,
  all.Straightness = all.Straightness,
  st.Length = st.Length,
  st.Speed.Mean = st.Speed.Mean,
  st.Straightness = st.Straightness,
  wob.Length = wob.Length,
  wob.Speed.Mean = wob.Speed.Mean,
  wob.Straightness = wob.Straightness,
  mig.Length = mig.Length,
  mig.Speed.Mean = mig.Speed.Mean,
  mig.Straightness = mig.Straightness,
  mig2.Length = mig2.Length,
  mig2.Speed.Mean = mig2.Speed.Mean,
  mig2.Straightness = mig2.Straightness,
  Lmig.Length = Lmig.Length,
  Lmig.Speed.Mean = Lmig.Speed.Mean,
  Lmig.Straightness = Lmig.Straightness
)

# Save all raw data CSVs
for (name in names(raw_data_list)) {
  save_csv(raw_data_list[[name]], name)
}

###################### Calculate and Save Summary Statistics ##############

# Calculate summaries for all data frames
all.Length.mean <- calculate_summary(all.Length)
all.Speed.Mean.mean <- calculate_summary(all.Speed.Mean)
all.Straightness.mean <- calculate_summary(all.Straightness)

st.Length.mean <- calculate_summary(st.Length)
st.Speed.Mean.mean <- calculate_summary(st.Speed.Mean)
st.Straightness.mean <- calculate_summary(st.Straightness)

wob.Length.mean <- calculate_summary(wob.Length)
wob.Speed.Mean.mean <- calculate_summary(wob.Speed.Mean)
wob.Straightness.mean <- calculate_summary(wob.Straightness)

mig.Length.mean <- calculate_summary(mig.Length)
mig.Speed.Mean.mean <- calculate_summary(mig.Speed.Mean)
mig.Straightness.mean <- calculate_summary(mig.Straightness)

mig2.Length.mean <- calculate_summary(mig2.Length)
mig2.Speed.Mean.mean <- calculate_summary(mig2.Speed.Mean)
mig2.Straightness.mean <- calculate_summary(mig2.Straightness)

Lmig.Length.mean <- calculate_summary(Lmig.Length)
Lmig.Speed.Mean.mean <- calculate_summary(Lmig.Speed.Mean)
Lmig.Straightness.mean <- calculate_summary(Lmig.Straightness)

# List of all summary data frames
summary_data_list <- list(
  all.Length.mean = all.Length.mean,
  all.Speed.Mean.mean = all.Speed.Mean.mean,
  all.Straightness.mean = all.Straightness.mean,
  st.Length.mean = st.Length.mean,
  st.Speed.Mean.mean = st.Speed.Mean.mean,
  st.Straightness.mean = st.Straightness.mean,
  wob.Length.mean = wob.Length.mean,
  wob.Speed.Mean.mean = wob.Speed.Mean.mean,
  wob.Straightness.mean = wob.Straightness.mean,
  mig.Length.mean = mig.Length.mean,
  mig.Speed.Mean.mean = mig.Speed.Mean.mean,
  mig.Straightness.mean = mig.Straightness.mean,
  mig2.Length.mean = mig2.Length.mean,
  mig2.Speed.Mean.mean = mig2.Speed.Mean.mean,
  mig2.Straightness.mean = mig2.Straightness.mean,
  Lmig.Length.mean = Lmig.Length.mean,
  Lmig.Speed.Mean.mean = Lmig.Speed.Mean.mean,
  Lmig.Straightness.mean = Lmig.Straightness.mean
)

# Save all summary CSVs
for (name in names(summary_data_list)) {
  save_csv(summary_data_list[[name]], name)
}

# Check some outputs
all.Length.mean
all.Speed.Mean.mean
all.Straightness.mean

############################## USE ONLY FOR PER ROI BASIS ANALYSIS #############################
# Calculate summaries for all data frames
all.Length.mean <- calculate_summary(all.Length)
all.Speed.Mean.mean <- calculate_summary(all.Speed.Mean)
all.Straightness.mean <- calculate_summary(all.Straightness)

# List of all summary data frames FOR ROIS
summary_data_list <- list(
  ROI.all.Length.mean = all.Length.mean,
  ROI.all.Speed.Mean.mean = all.Speed.Mean.mean,
  ROI.all.Straightness.mean = all.Straightness.mean)

# Save all summary CSVs
for (name in names(summary_data_list)) {
  save_csv(summary_data_list[[name]], name)
}

# Check some outputs
all.Length.mean
all.Speed.Mean.mean
all.Straightness.mean

####################################

###### Combine proportions files
files.prop <- list.files(wdir, pattern = "*Proportions.csv", full.names = FALSE) # Get proportion files
files.prop
ldf.prop <- sapply(files.prop, read.csv, simplify = FALSE, USE.NAMES = TRUE) # Read proportion files
ldf.prop
Proportions <- do.call("rbind", ldf.prop) # Merge proportion files together into one data frame
Proportions # Check output

write.csv(Proportions, "Proportions.csv", row.names = FALSE) # Write output as a .csv file
