# The aim of this script is to facilitate and automate the aggregation of movement metrics produced by Imaris.
# The script will extract the statistics of interest from each film and enable the splitting of behaviour categories
# described in Lafor?ts et al. iScience. 2023 Mar 29;26(4):106514.
# The script also produces the following outputs:
#   - files for each film containing all statistics for individual cells
#   - files for each sample (pooling two film of the same sample together) containing all statistics for individual cells
#   - A file aggregating all samples for each statistic of interest, both with individual cell value and sample mean. This is done for all cells but also for each behaviour categories
#   - A files reporting the proportions of each behaviour categories in all samples

rm(list=ls())
gc()

# Install and call the required R packages.
#install.packages("tidyverse")
library(tidyverse)

# Set working directory as the folder containing all the input folders. Replace with your working directory's file path
setwd("Z:\\lightmicro-inst\\scharalambous\\EXP25N106_MIPs_t\\Outputs_celltracking\\D2\\RP") 
wdir <- "Z:\\lightmicro-inst\\scharalambous\\EXP25N106_MIPs_t\\Outputs_celltracking\\D2\\RP" 
# It is useful to have the working directory path saved so it can be called


# Get all film folders
folders <- list.dirs(path = wdir)
folders <- folders[-1]
folders # Check the all input folders are loaded


# Iteration for each film folder: read the four required statistic files, merge them according to Track (cell) ID and write the output as .csv
for (folder in folders) { # iterates the script in all input folders
  
  ## Load .csv files containing the statistics of interest for one film, with only the track number and stat columns
  dir.name <- basename(folder) # Gets the input folder name
  filepath <- paste(wdir, dir.name, sep = "/") # Creates a generic file path for input files
  cat("Processing folder:", folder, "\n")
  
  ## Calculate the maximum displacement and create a .csv file for it
  # Find the file you want (exclude Track_ and Position_[XYZ])
  filepath.Pos <- list.files(filepath, pattern = "_Position\\.csv$", full.names = TRUE)
  filepath.Pos <- filepath.Pos[!grepl("Track_|Position_[XYZ]", filepath.Pos)]
  
  cat("Found position files:", length(filepath.Pos), "\n")
  
  if (length(filepath.Pos) == 0) {
    cat("No valid Position file in", filepath, "\n")
    next
  }
  
  filepath.Pos <- filepath.Pos[1]  # choose one if multiple
  print(filepath.Pos)
  

  
  # Read everything, then keep columns by NAME
  raw <- read.csv(filepath.Pos, skip = 3, header = TRUE, fill = TRUE, check.names = TRUE)
  # Names now like: Position.X, Position.Y, Position.Z, ..., Time, TrackID, ID
  View (raw)
  need <- c("Position.X", "Position.Y", "Time", "TrackID")
  stopifnot(all(need %in% names(raw)))
  
  Positions <- raw[need]
  
  # Drop missing/blank TrackID rows
  Positions <- Positions[!is.na(Positions$TrackID) & Positions$TrackID != "", , drop = FALSE]
  
  # Use the actual IDs
  IDs <- unique(Positions$TrackID)
  
  # Arrange (optional)
  Positions.arranged <- Positions[order(Positions$TrackID, Positions$Time), , drop = FALSE]
  
  # Easiest: split by TrackID (no loop needed)
  listx <- split(Positions.arranged, Positions.arranged$TrackID)


  # Create a list identical to listx but with only t0. This is used as the initial position for calculating the displacement at each time point
  list.t0 <- list()
    for (i in seq_along(IDs)) { 
      list.t0[[i]] <- subset(Positions.arranged, TrackID == IDs[i])
      list.t0[[i]] <- list.t0[[i]][1,] # first time point row
      list.t0[[i]] <- list.t0[[i]][rep(seq_len(nrow(list.t0[[i]])), nrow(listx[[i]])), ] 
  }
  
  # Calculate the displacement of each cell at each time point
    for (i in seq_along(IDs)) { 
    listx[[i]]$Dis <- sqrt((listx[[i]][1]-list.t0[[i]][1])^2 + (listx[[i]][2]-list.t0[[i]][2])^2)
  }
  
  # Creates a list with only the maximum displacement observed for each track
  list.DisMax <- list()
  for (i in seq_along(IDs)) { 
    list.DisMax[[i]] <- data.frame(ID = IDs[i],MaxDis = max(listx[[i]]$Dis, na.rm = TRUE))
  }

  
  TMaxDis <- do.call("rbind", list.DisMax) # Binds all the maximum displacement data frame in the list into one.
  
  write.csv(TMaxDis, paste(filepath, "__Max_Displacement.csv", sep = "/"), row.names = FALSE) # write the file in the corresponding input folder

  ## Create filepaths for each required statistic file
  filepath.Length <- list.files(filepath, pattern = "_Track_Length\\.csv$", full.names = TRUE)
  filepath.Speed <- list.files(filepath, pattern ="_Track_Speed_Mean\\.csv$", full.names = TRUE)
  filepath.Speed <- filepath.Speed[!grepl("_Track_Speed_Mean_[XYZ]", filepath.Speed)]
  filepath.Straightness <- list.files(filepath, pattern="_Track_Straightness\\.csv$", full.names = TRUE)
  
  ## Read all required statistic files
  TLength <- read.csv(filepath.Length, skip = 3, fill = TRUE, header = TRUE)[,c(1,4)]
  TSpeed <- read.csv(filepath.Speed, skip = 3, fill = TRUE, header = TRUE)[,c(1,4)]
  TStraightness <- read.csv(filepath.Straightness, skip = 3,  fill = TRUE, header = TRUE)[,c(1,4)]
  
  ## Merge the required statistic files according to the ID columns
  Allstats <- Reduce(function(x, y) merge(x, y, by = "ID", all=TRUE), list(TLength, TSpeed, TStraightness, TMaxDis)) # The name "Allstats" refers all to the statistics of interest used to split the different behavioural categories, but it does not actually contain all the statistics produced by Imaris
  head(Allstats) # Check the Allstats output
  dim(Allstats)
  
  ## Write a .csv file containing all statistics of interest for the corresponding film
  filename1 <- paste(dir.name, "Allstats.csv", sep = "_")
  write.csv(Allstats, paste(wdir, filename1, sep = "/"), row.names = FALSE)

}
########################## Remove non paired input folders from working directory #############################

# 1) Collect only Allstats CSVs (skip already-merged outputs)
files <- list.files(wdir, pattern = "Allstats\\.csv$", full.names = TRUE)
files <- files[!grepl("\\.merged\\.csv$", files, ignore.case = TRUE)]

# 2) Build a grouping key: everything before `_ROI_<number>`
base_names <- basename(files)
keys <- sub("_ROI_\\d+.*$", "", base_names)   # e.g. "slice1", "slice2", ...

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
  
  # (optional) keep it in the workspace as k.merged
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
  colnames(st) <- c("ID", "Track.Length", "Track.Speed.Mean", "Track.Straightness", "Max.Displacement.Length")
  wob <- ldf_Allstats[[i]][(ldf_Allstats[[i]]$MaxDis<7) & (ldf_Allstats[[i]]$Track.Length>=13),] # subset wobbling cells
  colnames(wob) <- c("ID", "Track.Length", "Track.Speed.Mean", "Track.Straightness", "Max.Displacement.Length")
  mig <- ldf_Allstats[[i]][(ldf_Allstats[[i]]$MaxDis>=7) & ((ldf_Allstats[[i]]$MaxDis<25) | (ldf_Allstats[[i]]$Track.Length<40)),] # Subset migrating cells
  colnames(mig) <- c("ID", "Track.Length", "Track.Speed.Mean", "Track.Straightness", "Max.Displacement.Length")
  mig2 <- ldf_Allstats[[i]][(ldf_Allstats[[i]]$MaxDis>=7), ] # Subset migrating + long migrating cells together
  colnames(mig2) <- c("ID", "Track.Length", "Track.Speed.Mean", "Track.Straightness", "Max.Displacement.Length")
  Lmig <- ldf_Allstats[[i]][(ldf_Allstats[[i]]$MaxDis>=25) & (ldf_Allstats[[i]]$Track.Length>=40),] # Subset long migrating cells
  colnames(Lmig) <- c("ID", "Track.Length", "Track.Speed.Mean", "Track.Straightness", "Max.Displacement.Length")
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
############################################################### Grouping all populations into the three categories

# Group all samples together into one file for each statistic of interest
## Create a list of data frames containing all statistics of interest for all samples
files.all <- list.files(wdir, pattern = "Allstats\\.merged\\.csv$", full.names = FALSE) # Get all Allstats.merged files
files.all

ldf.all <- sapply(files.all, read.csv, simplify = FALSE, USE.NAMES = TRUE) # Read the files and store output in a list
ldf.all <- lapply(ldf.all, as.data.frame, col.names = c("ID", "Track.Length", "Track.Speed.Mean", "Track.Straightness", "MaxDis")) # Ensure each element of the list is a data frame and assign clear column names

 ldf.all <- lapply(ldf.all, function(df) {
  # drop any rows with NA (in any column)
  df <- df[complete.cases(df), , drop = FALSE]
 # drop rows where Track.Length is 0
  df <- df[df$Track.Length != 0, , drop = FALSE]
  # drop rows where Track.Length is 0
  df <- df[df$Track.Speed.Mean != 0, , drop = FALSE]
 # drop rows where Track.Length is 0
  df <- df[df$Track.Straightness != 0, , drop = FALSE]
  df
 })
ldf.all # Check output

## Create a data frame containing the Length statistic for all samples and all cells.
### As all films do not contain the same number of cells (tracks) we need to create a base data frame where the first column is very long
### We use the Length statistic of the first sample as the base
all.Length <- data.frame(ldf.all[[1]]$Track.Length)
head(all.Length)
colnames(all.Length)[1] <- gsub("_Allstats.merged.csv", "", names(ldf.all[1])) # Rename the column with just the sample's name
dim(all.Length)
head(all.Length)

### Expand the data frame (just one column for now) with NA values long enough that it can be merged with very long columns

#NA.vec <- rep(NA, 10000) # Creates a vector with 10 000 "NA" values (10000 is long enough unless one or more of your analysed samples include more than 10000 track)
#NA.df <- as.data.frame(NA.vec) # Transform the vector to a data frame (one column)
#colnames(NA.df)[1] <- colnames(all.Length)[1] # Set the NA data frame column name
#dim(NA.df)
#all.Length <- rbind(all.Length, NA.df) # Bind the two data frames so that the NA values are now below the stat values
#all.Length
#dim(all.Length) # Check output

# ## Compute how long the final table should be (max number of rows across samples)
# max_len <- max(vapply(ldf.all, nrow, integer(1)))
# 
# ## Pad the base column up to max_len with NAs
# nr <- nrow(all.Length)
# if (nr < max_len) {
#   all.Length[[1]] <- c(all.Length[[1]], rep(NA, max_len - nr))
# }
# 
# for (i in 2:length(ldf.all)) {
#   v <- ldf.all[[i]]$Track.Length
#   length(v) <- max_len  # pads with NA if shorter, truncates if somehow longer
#   sample_name <- sub("_Allstats\\.merged\\.csv$", "", names(ldf.all)[i])
#   all.Length[[sample_name]] <- v  # add with the final name directly
# }

max_len <- max(vapply(ldf.all, nrow, integer(1)))

# Step 2: Pad all columns to the same length
padded_list <- lapply(seq_along(ldf.all), function(i) {
  v <- ldf.all[[i]]$Track.Length
  length(v) <- max_len  # safely pad with NA
  return(v)
})

# Step 3: Combine into a single data frame
all.Length <- as.data.frame(padded_list)
colnames(all.Length) <- sub("_Allstats\\.merged\\.csv$", "", names(ldf.all))
head(all.Length) # Check output


## Create a data frame containing the Mean Speed statistic for all samples and all cells. 
all.Speed.Mean <- data.frame(ldf.all[[1]]$Track.Speed.Mean)
head(all.Speed.Mean)
colnames(all.Speed.Mean)[1] <- gsub("_Allstats.merged.csv", "", names(ldf.all[1]))
dim(all.Speed.Mean)
head(all.Speed.Mean)


# ## Pad the base column up to max_len with NAs
# nr <- nrow(all.Speed.Mean)
# if (nr < max_len) {
#   all.Speed.Mean[[1]] <- c(all.Speed.Mean[[1]], rep(NA, max_len - nr))
# }
# 
# for (i in 2:length(ldf.all)) {
#   v <- ldf.all[[i]]$Track.Speed.Mean
#   length(v) <- max_len  # pads with NA if shorter, truncates if somehow longer
#   sample_name <- sub("_Allstats\\.merged\\.csv$", "", names(ldf.all)[i])
#   all.Speed.Mean[[sample_name]] <- v  # add with the final name directly
# }

# Step 2: Pad all columns to the same length
padded_list <- lapply(seq_along(ldf.all), function(i) {
  v <- ldf.all[[i]]$Track.Speed.Mean
  length(v) <- max_len  # safely pad with NA
  return(v)
})

# Step 3: Combine into a single data frame
all.Speed.Mean <- as.data.frame(padded_list)
colnames(all.Speed.Mean) <- sub("_Allstats\\.merged\\.csv$", "", names(ldf.all))
head(all.Speed.Mean) # Check output

### Create a data frame containing the Straightness statistic for all samples and all cells.
all.Straightness <- data.frame(ldf.all[[1]]$Track.Straightness)
head(all.Straightness)
colnames(all.Straightness)[1] <- gsub("_Allstats.merged.csv", "", names(ldf.all[1]))
dim(all.Straightness)
head(all.Straightness)


# ## Pad the base column up to max_len with NAs
# nr <- nrow(all.Straightness)
# if (nr < max_len) {
#   all.Straightness[[1]] <- c(all.Straightness[[1]], rep(NA, max_len - nr))
# }
# 
# for (i in 2:length(ldf.all)) {
#   v <- ldf.all[[i]]$Track.Straightness
#   length(v) <- max_len  # pads with NA if shorter, truncates if somehow longer
#   sample_name <- sub("_Allstats\\.merged\\.csv$", "", names(ldf.all)[i])
#   all.Straightness[[sample_name]] <- v  # add with the final name directly
# }

# Step 2: Pad all columns to the same length
padded_list <- lapply(seq_along(ldf.all), function(i) {
  v <- ldf.all[[i]]$Track.Straightness
  length(v) <- max_len  # safely pad with NA
  return(v)
})

# Step 3: Combine into a single data frame
all.Straightness <- as.data.frame(padded_list)
colnames(all.Straightness) <- sub("_Allstats\\.merged\\.csv$", "", names(ldf.all))
head(all.Straightness) # Check output

##############  Individual populations
############################### Static ###########################

## ---------------- LOAD & CLEAN BASE ---------------- ##
files.st <- list.files(wdir, pattern = "*st.csv", full.names = FALSE)
ldf.st <- sapply(files.st, read.csv, simplify = FALSE, USE.NAMES = TRUE)
ldf.st <- lapply(ldf.st, as.data.frame, col.names = c("ID", "Track.Length", "Track.Speed.Mean", "Track.Straightness", "Track.Displacement.Length"))

# Global cleaning:
#  - remove NAs
#  - drop rows where ANY of Length, Speed.Mean, or Straightness = 0
ldf.st.cleaned <- lapply(ldf.st, function(df) {
  df <- df[complete.cases(df), , drop = FALSE]
  df <- df[
    df$Track.Length != 0 &
      df$Track.Speed.Mean != 0 &
      df$Track.Straightness != 0,
    , drop = FALSE
  ]
  df
})

## ---------------- LENGTH ---------------- ##
ldf.st.Length <- ldf.st.cleaned

max_len <- max(vapply(ldf.st.Length, nrow, integer(1)))
padded_list <- lapply(seq_along(ldf.st.Length), function(i) {
  v <- ldf.st.Length[[i]]$Track.Length
  length(v) <- max_len
  return(v)
})
st.Length <- as.data.frame(padded_list)
colnames(st.Length) <- sub("_st.csv$", "", names(ldf.st.Length))
head(st.Length)

## ---------------- SPEED ---------------- ##
ldf.st.Speed <- ldf.st.cleaned

max_len <- max(vapply(ldf.st.Speed, nrow, integer(1)))
padded_list <- lapply(seq_along(ldf.st.Speed), function(i) {
  v <- ldf.st.Speed[[i]]$Track.Speed.Mean
  length(v) <- max_len
  return(v)
})
st.Speed.Mean <- as.data.frame(padded_list)
colnames(st.Speed.Mean) <- sub("_st.csv$", "", names(ldf.st.Speed))
head(st.Speed.Mean)

## ---------------- STRAIGHTNESS ---------------- ##
ldf.st.Straightness <- ldf.st.cleaned

max_len <- max(vapply(ldf.st.Straightness, nrow, integer(1)))
padded_list <- lapply(seq_along(ldf.st.Straightness), function(i) {
  v <- ldf.st.Straightness[[i]]$Track.Straightness
  length(v) <- max_len
  return(v)
})
st.Straightness <- as.data.frame(padded_list)
colnames(st.Straightness) <- sub("_st.csv$", "", names(ldf.st.Straightness))
head(st.Straightness)


############################################################# Wobbling #######################################

## Create a data frame containing the Length statistic for only wobbling cells for all samples.
## ---------------- LOAD & CLEAN BASE ---------------- ##
files.wob <- list.files(wdir, pattern = "*wob.csv", full.names = FALSE)
ldf.wob <- sapply(files.wob, read.csv, simplify = FALSE, USE.NAMES = TRUE)
ldf.wob <- lapply(ldf.wob, as.data.frame, 
                  col.names = c("ID", "Track.Length", "Track.Speed.Mean", 
                                "Track.Straightness", "Track.Displacement.Length"))

# Global cleaning:
#  - remove NAs
#  - drop rows where ANY of Length, Speed.Mean, or Straightness = 0
ldf.wob.cleaned <- lapply(ldf.wob, function(df) {
  df <- df[complete.cases(df), , drop = FALSE]
  df <- df[
    df$Track.Length != 0 &
      df$Track.Speed.Mean != 0 &
      df$Track.Straightness != 0,
    , drop = FALSE
  ]
  df
})

## ---------------- LENGTH ---------------- ##
ldf.wob.Length <- ldf.wob.cleaned

max_len <- max(vapply(ldf.wob.Length, nrow, integer(1)))
padded_list <- lapply(seq_along(ldf.wob.Length), function(i) {
  v <- ldf.wob.Length[[i]]$Track.Length
  length(v) <- max_len
  return(v)
})
wob.Length <- as.data.frame(padded_list)
colnames(wob.Length) <- sub("_wob.csv$", "", names(ldf.wob.Length))
head(wob.Length)

## ---------------- SPEED ---------------- ##
ldf.wob.Speed <- ldf.wob.cleaned

max_len <- max(vapply(ldf.wob.Speed, nrow, integer(1)))
padded_list <- lapply(seq_along(ldf.wob.Speed), function(i) {
  v <- ldf.wob.Speed[[i]]$Track.Speed.Mean
  length(v) <- max_len
  return(v)
})
wob.Speed.Mean <- as.data.frame(padded_list)
colnames(wob.Speed.Mean) <- sub("_wob.csv$", "", names(ldf.wob.Speed))
head(wob.Speed.Mean)

## ---------------- STRAIGHTNESS ---------------- ##
ldf.wob.Straightness <- ldf.wob.cleaned

max_len <- max(vapply(ldf.wob.Straightness, nrow, integer(1)))
padded_list <- lapply(seq_along(ldf.wob.Straightness), function(i) {
  v <- ldf.wob.Straightness[[i]]$Track.Straightness
  length(v) <- max_len
  return(v)
})
wob.Straightness <- as.data.frame(padded_list)
colnames(wob.Straightness) <- sub("_wob.csv$", "", names(ldf.wob.Straightness))
head(wob.Straightness)


######################################################## Migrating ########################


## Create a data frame containing the Length statistic for only migrating cells for all samples.
## ---------------- LOAD & CLEAN BASE ---------------- ##
files.mig <- list.files(wdir, pattern = "_mig\\.csv$", full.names = FALSE)
ldf.mig <- sapply(files.mig, read.csv, simplify = FALSE, USE.NAMES = TRUE)
ldf.mig <- lapply(ldf.mig, as.data.frame, 
                  col.names = c("ID", "Track.Length", "Track.Speed.Mean", 
                                "Track.Straightness", "Track.Displacement.Length"))

# Global cleaning:
#  - remove NAs
#  - drop rows where ANY of Length, Speed.Mean, or Straightness = 0
ldf.mig.cleaned <- lapply(ldf.mig, function(df) {
  df <- df[complete.cases(df), , drop = FALSE]
  df <- df[
    df$Track.Length != 0 &
      df$Track.Speed.Mean != 0 &
      df$Track.Straightness != 0,
    , drop = FALSE
  ]
  df
})

## ---------------- LENGTH ---------------- ##
ldf.mig.Length <- ldf.mig.cleaned

max_len <- max(vapply(ldf.mig.Length, nrow, integer(1)))
padded_list <- lapply(seq_along(ldf.mig.Length), function(i) {
  v <- ldf.mig.Length[[i]]$Track.Length
  length(v) <- max_len
  return(v)
})
mig.Length <- as.data.frame(padded_list)
colnames(mig.Length) <- sub("_mig.csv$", "", names(ldf.mig.Length))
head(mig.Length)

## ---------------- SPEED ---------------- ##
ldf.mig.Speed <- ldf.mig.cleaned

max_len <- max(vapply(ldf.mig.Speed, nrow, integer(1)))
padded_list <- lapply(seq_along(ldf.mig.Speed), function(i) {
  v <- ldf.mig.Speed[[i]]$Track.Speed.Mean
  length(v) <- max_len
  return(v)
})
mig.Speed.Mean <- as.data.frame(padded_list)
colnames(mig.Speed.Mean) <- sub("_mig.csv$", "", names(ldf.mig.Speed))
head(mig.Speed.Mean)

## ---------------- STRAIGHTNESS ---------------- ##
ldf.mig.Straightness <- ldf.mig.cleaned

max_len <- max(vapply(ldf.mig.Straightness, nrow, integer(1)))
padded_list <- lapply(seq_along(ldf.mig.Straightness), function(i) {
  v <- ldf.mig.Straightness[[i]]$Track.Straightness
  length(v) <- max_len
  return(v)
})
mig.Straightness <- as.data.frame(padded_list)
colnames(mig.Straightness) <- sub("_mig.csv$", "", names(ldf.mig.Straightness))
head(mig.Straightness)


####################################################  Migrating + long Migrating #######################

## Create a data frame containing the Length statistic for only migrating + Long migrating cells for all samples.
## ---------------- LOAD & CLEAN BASE ---------------- ##
files.mig2 <- list.files(wdir, pattern = "_mig2\\.csv$", full.names = FALSE)
ldf.mig2 <- sapply(files.mig2, read.csv, simplify = FALSE, USE.NAMES = TRUE)
ldf.mig2 <- lapply(ldf.mig2, as.data.frame, 
                   col.names = c("ID", "Track.Length", "Track.Speed.Mean", 
                                 "Track.Straightness", "Track.Displacement.Length"))

# Global cleaning:
#  - remove NAs
#  - drop rows where ANY of Length, Speed.Mean, or Straightness = 0
ldf.mig2.cleaned <- lapply(ldf.mig2, function(df) {
  df <- df[complete.cases(df), , drop = FALSE]
  df <- df[
    df$Track.Length != 0 &
      df$Track.Speed.Mean != 0 &
      df$Track.Straightness != 0,
    , drop = FALSE
  ]
  df
})

## ---------------- LENGTH ---------------- ##
ldf.mig2.Length <- ldf.mig2.cleaned

max_len <- max(vapply(ldf.mig2.Length, nrow, integer(1)))
padded_list <- lapply(seq_along(ldf.mig2.Length), function(i) {
  v <- ldf.mig2.Length[[i]]$Track.Length
  length(v) <- max_len
  return(v)
})
mig2.Length <- as.data.frame(padded_list)
colnames(mig2.Length) <- sub("_mig2.csv$", "", names(ldf.mig2.Length))
head(mig2.Length)

## ---------------- SPEED ---------------- ##
ldf.mig2.Speed <- ldf.mig2.cleaned

max_len <- max(vapply(ldf.mig2.Speed, nrow, integer(1)))
padded_list <- lapply(seq_along(ldf.mig2.Speed), function(i) {
  v <- ldf.mig2.Speed[[i]]$Track.Speed.Mean
  length(v) <- max_len
  return(v)
})
mig2.Speed.Mean <- as.data.frame(padded_list)
colnames(mig2.Speed.Mean) <- sub("_mig2.csv$", "", names(ldf.mig2.Speed))
head(mig2.Speed.Mean)

## ---------------- STRAIGHTNESS ---------------- ##
ldf.mig2.Straightness <- ldf.mig2.cleaned

max_len <- max(vapply(ldf.mig2.Straightness, nrow, integer(1)))
padded_list <- lapply(seq_along(ldf.mig2.Straightness), function(i) {
  v <- ldf.mig2.Straightness[[i]]$Track.Straightness
  length(v) <- max_len
  return(v)
})
mig2.Straightness <- as.data.frame(padded_list)
colnames(mig2.Straightness) <- sub("_mig2.csv$", "", names(ldf.mig2.Straightness))
head(mig2.Straightness)


###############################################  Long Migrating ########################################

## Create a data frame containing the Length statistic for only Long migrating cells for all samples.
## ---------------- LOAD & CLEAN BASE ---------------- ##
files.Lmig <- list.files(wdir, pattern = "_Lmig\\.csv$", full.names = FALSE)
ldf.Lmig <- sapply(files.Lmig, read.csv, simplify = FALSE, USE.NAMES = TRUE)
ldf.Lmig <- lapply(ldf.Lmig, as.data.frame, 
                   col.names = c("ID", "Track.Length", "Track.Speed.Mean", 
                                 "Track.Straightness", "Track.Displacement.Length"))

# Global cleaning:
#  - remove NAs
#  - drop rows where ANY of Length, Speed.Mean, or Straightness = 0
ldf.Lmig.cleaned <- lapply(ldf.Lmig, function(df) {
  df <- df[complete.cases(df), , drop = FALSE]
  df <- df[
    df$Track.Length != 0 &
      df$Track.Speed.Mean != 0 &
      df$Track.Straightness != 0,
    , drop = FALSE
  ]
  df
})

## ---------------- LENGTH ---------------- ##
ldf.Lmig.Length <- ldf.Lmig.cleaned

max_len <- max(vapply(ldf.Lmig.Length, nrow, integer(1)))
padded_list <- lapply(seq_along(ldf.Lmig.Length), function(i) {
  v <- ldf.Lmig.Length[[i]]$Track.Length
  length(v) <- max_len
  return(v)
})
Lmig.Length <- as.data.frame(padded_list)
colnames(Lmig.Length) <- sub("_Lmig.csv$", "", names(ldf.Lmig.Length))
head(Lmig.Length)

## ---------------- SPEED ---------------- ##
ldf.Lmig.Speed <- ldf.Lmig.cleaned

max_len <- max(vapply(ldf.Lmig.Speed, nrow, integer(1)))
padded_list <- lapply(seq_along(ldf.Lmig.Speed), function(i) {
  v <- ldf.Lmig.Speed[[i]]$Track.Speed.Mean
  length(v) <- max_len
  return(v)
})
Lmig.Speed.Mean <- as.data.frame(padded_list)
colnames(Lmig.Speed.Mean) <- sub("_Lmig.csv$", "", names(ldf.Lmig.Speed))
head(Lmig.Speed.Mean)

## ---------------- STRAIGHTNESS ---------------- ##
ldf.Lmig.Straightness <- ldf.Lmig.cleaned

max_len <- max(vapply(ldf.Lmig.Straightness, nrow, integer(1)))
padded_list <- lapply(seq_along(ldf.Lmig.Straightness), function(i) {
  v <- ldf.Lmig.Straightness[[i]]$Track.Straightness
  length(v) <- max_len
  return(v)
})
Lmig.Straightness <- as.data.frame(padded_list)
colnames(Lmig.Straightness) <- sub("_Lmig.csv$", "", names(ldf.Lmig.Straightness))
head(Lmig.Straightness)


############

# Write all outputs as .csv files
write.csv(all.Length, "all.Length.csv", row.names = FALSE)
write.csv(all.Speed.Mean, "all.Speed.Mean.csv", row.names = FALSE)
write.csv(all.Straightness, "all.Straightness.csv", row.names = FALSE)
write.csv(st.Length, "st.Length.csv", row.names = FALSE)
write.csv(st.Speed.Mean, "st.Speed.Mean.csv", row.names = FALSE)
write.csv(st.Straightness, "st.Straightness.csv", row.names = FALSE)
write.csv(wob.Length, "wob.Length.csv", row.names = FALSE)
write.csv(wob.Speed.Mean, "wob.Speed.Mean.csv", row.names = FALSE)
write.csv(wob.Straightness, "wob.Straightness.csv", row.names = FALSE)
write.csv(mig.Length, "mig.Length.csv", row.names = FALSE)
write.csv(mig.Speed.Mean, "mig.Speed.Mean.csv", row.names = FALSE)
write.csv(mig.Straightness, "mig.Straightness.csv", row.names = FALSE)
write.csv(mig2.Length, "mig2.Length.csv", row.names = FALSE)
write.csv(mig2.Speed.Mean, "mig2.Speed.Mean.csv", row.names = FALSE)
write.csv(mig2.Straightness, "mig2.Straightness.csv", row.names = FALSE)
write.csv(Lmig.Length, "Lmig.Length.csv", row.names = FALSE)
write.csv(Lmig.Speed.Mean, "Lmig.Speed.Mean.csv", row.names = FALSE)
write.csv(Lmig.Straightness, "Lmig.Straightness.csv", row.names = FALSE)

############

# Organize data to have one data point per sample (the mean) for all cells
all.Length.mean <- data.frame(Mean=colMeans(all.Length, na.rm = TRUE)) # Creates a data frame with the mean value for each samples
all.Length.mean

all.Speed.Mean.mean <- data.frame(Mean=colMeans(all.Speed.Mean, na.rm = TRUE)) # Creates a data frame with the mean value for each samples
all.Speed.Mean.mean

all.Straightness.mean <- data.frame(Mean=colMeans(all.Straightness, na.rm = TRUE)) # Creates a data frame with the mean value for each samples
all.Straightness.mean

# Organize data to have one data point per sample (the mean) for static cells
st.Length.mean <- data.frame(Mean=colMeans(st.Length, na.rm = TRUE)) # Creates a data frame with the mean value for each samples
st.Length.mean

st.Speed.Mean.mean <- data.frame(Mean=colMeans(st.Speed.Mean, na.rm = TRUE)) # Creates a data frame with the mean value for each samples
st.Speed.Mean.mean

st.Straightness.mean <- data.frame(Mean=colMeans(st.Straightness, na.rm = TRUE)) # Creates a data frame with the mean value for each samples
st.Straightness.mean

# Organize data to have one data point per sample (the mean) for wobbling cells
wob.Length.mean <- data.frame(Mean=colMeans(wob.Length, na.rm = TRUE)) # Creates a data frame with the mean value for each samples
wob.Length.mean

wob.Speed.Mean.mean <- data.frame(Mean=colMeans(wob.Speed.Mean, na.rm = TRUE)) # Creates a data frame with the mean value for each samples
wob.Speed.Mean.mean

wob.Straightness.mean <- data.frame(Mean=colMeans(wob.Straightness, na.rm = TRUE)) # Creates a data frame with the mean value for each samples
wob.Straightness.mean

# Organize data to have one data point per sample (the mean) for migrating cells
mig.Length.mean <- data.frame(Mean= colMeans(mig.Length, na.rm = TRUE)) # Creates a data frame with the mean value for each samples
mig.Length.mean

mig.Speed.Mean.mean <- data.frame(Mean=colMeans(mig.Speed.Mean, na.rm = TRUE)) # Creates a data frame with the mean value for each samples
mig.Speed.Mean.mean

mig.Straightness.mean <- data.frame(Mean=colMeans(mig.Straightness, na.rm = TRUE)) # Creates a data frame with the mean value for each samples
mig.Straightness.mean

# Organize data to have one data point per sample (the mean) for migrating + Long migrating cells
mig2.Length.mean <- data.frame(Mean=colMeans(mig2.Length, na.rm = TRUE)) # Creates a data frame with the mean value for each samples
mig2.Length.mean

mig2.Speed.Mean.mean <- data.frame(Mean=colMeans(mig2.Speed.Mean, na.rm = TRUE)) # Creates a data frame with the mean value for each samples
mig2.Speed.Mean.mean

mig2.Straightness.mean <- data.frame(Mean=colMeans(mig2.Straightness, na.rm = TRUE)) # Creates a data frame with the mean value for each samples
mig2.Straightness.mean

# Organize data to have one data point per sample (the mean) for Long migrating cells
Lmig.Length.mean <- data.frame(Mean=colMeans(Lmig.Length, na.rm = TRUE)) # Creates a data frame with the mean value for each samples
Lmig.Length.mean

Lmig.Speed.Mean.mean <- data.frame(Mean=colMeans(Lmig.Speed.Mean, na.rm = TRUE)) # Creates a data frame with the mean value for each samples
Lmig.Speed.Mean.mean

Lmig.Straightness.mean <- data.frame(Mean=colMeans(Lmig.Straightness, na.rm = TRUE)) # Creates a data frame with the mean value for each samples
Lmig.Straightness.mean

# Write outputs as .csv files
write.csv(all.Length.mean, "all.Length.mean.csv", row.names = FALSE)
write.csv(all.Speed.Mean.mean, "all.Speed.Mean.mean.csv", row.names = FALSE)
write.csv(all.Straightness.mean, "all.Straightness.mean.csv", row.names = FALSE)
write.csv(st.Length.mean, "st.Length.mean.csv", row.names = FALSE)
write.csv(st.Speed.Mean.mean, "st.Speed.Mean.mean.csv", row.names = FALSE)
write.csv(st.Straightness.mean, "st.Straightness.mean.csv", row.names = FALSE)
write.csv(wob.Length.mean, "wob.Length.mean.csv", row.names = FALSE)
write.csv(wob.Speed.Mean.mean, "wob.Speed.Mean.mean.csv", row.names = FALSE)
write.csv(wob.Straightness.mean, "wob.Straightness.mean.csv", row.names = FALSE)
write.csv(mig.Length.mean, "mig.Length.mean.csv", row.names = FALSE)
write.csv(mig.Speed.Mean.mean, "mig.Speed.Mean.mean.csv", row.names = FALSE)
write.csv(mig.Straightness.mean, "mig.Straightness.mean.csv", row.names = FALSE)
write.csv(mig2.Length.mean, "mig2.Length.mean.csv", row.names = FALSE)
write.csv(mig2.Speed.Mean.mean, "mig2.Speed.Mean.mean.csv", row.names = FALSE)
write.csv(mig2.Straightness.mean, "mig2.Straightness.mean.csv", row.names = FALSE)
write.csv(Lmig.Length.mean, "Lmig.Length.mean.csv", row.names = FALSE)
write.csv(Lmig.Speed.Mean.mean, "Lmig.Speed.Mean.mean.csv", row.names = FALSE)
write.csv(Lmig.Straightness.mean, "Lmig.Straightness.mean.csv", row.names = FALSE)

###### Combine proportions files
files.prop <- list.files(wdir, pattern = "*Proportions.csv", full.names = FALSE) # Get proportion files
files.prop
ldf.prop <- sapply(files.prop, read.csv, simplify = FALSE, USE.NAMES = TRUE) # Read proportion files
ldf.prop
Proportions <- do.call("rbind", ldf.prop) # Merge proportion files together into one data frame
Proportions # Check output

write.csv(Proportions, "Proportions.csv", row.names = FALSE) # Write output as a .csv file