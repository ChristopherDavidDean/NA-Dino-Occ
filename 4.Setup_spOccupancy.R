################################################################################
############ OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS #############
################################################################################

# Christopher D. Dean, Lewis A. Jones, Alfio A. Chiarenza, Sinéad Lyster, 
# Alex Farnsworth, Philip D. Mannion, Richard J. Butler.
# 2023
# Script written by Christopher D. Dean

################################################################################
#                     FILE 4: SETUP DATA FOR SPOCCUPANCY                       #
################################################################################

################################################################################
# 1. INITIAL SETUP
################################################################################

##### Load packages #####
library(spOccupancy)
library(stars)
library(ggplot2)
library(abind)
library(stringr)
library(MCMCvis)
library(purrr)
library(sp)

##### Load in Functions #####
source("0.Functions.R") # Import functions from other R file (must be in same working directory)

##### Set values #####
# Set resolution
res <- 0.1
# Set extent
e <- extent(-155, -72, 22.5, 73)
# Set max limit value
max_val <- 20
max_val_on <- TRUE
bin.type <- "scotese"
target <- "Ceratopsidae"
bins <- read.csv("Data/Occurrences/scotesebins.csv")
bins$bin <- bins$code
bin <- NA # Change this if you want a single season model

# Load occurrence dataset
master.occs.binned.targeted <- read.csv(paste("Prepped_data/Occurrence_Data/", 
                                              bin.type, "/", res, "_", bin.type, 
                                              "_occurrence_dataset.csv", sep = ""))
# Set random number for seed
rand <- round(runif(1, min = 0, max = 999)) 

################################################################################
# 2. PREPARING DATA
################################################################################

#########################
##### INITIAL SETUP #####
#########################

# Get grid cells for occurrences
get_grid(master.occs.binned.targeted, res = res, e = e)

# Prep data for 
prepare_for_spOcc(master.occs.binned.targeted.grid, single = FALSE, bin = bin)

# Get site IDs for sorting later
site_IDs <- as.numeric(rownames(eh_list[[1]][[1]]))

# Palaeo-rotate site IDs for each time bin to enable acquiring relevant covariate data
rotated <- p_rotate(res, e, site_IDs, bins, bin = bin)

# Extract palaeo-covariates
extracted_covs <- extract_p(rotated)

# Format occupancy covariates for spOccupancy
occ_covs <- format_occ_covs(p_cov_list = extracted_covs)

# If generating single season data, reorganise occupancy covariates into suitable format
if(is.na(bin) == F){
  occ_covs <- data.frame(bind_rows(!!! occ_covs))
  colnames(occ_covs) <- sub("$*", "", colnames(occ_covs))
}

# Format detection covariates for spOccupancy
det_covs <- organise_det(siteCoords, extracted_covs, 
                         occ_data = master.occs.binned.targeted.grid, 
                         bin = bin)

# Remove any sites without relevant data to ensure model runs
site_remove(eh_list, occ_covs, det_covs, siteCoords, single = FALSE)

# Transpose collected data into array format for spOccupancy
if(is.na(bin) == T){
  transpose_eh(eh_list, target)
}

#########################
##### SPATIAL SETUP #####
#########################

# Define the proj4 string for NAD83
nad83_proj <- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

# Define my spatial object
sp_points <- SpatialPoints(siteCoords)

# Apply CRS to data
proj4string(sp_points) <- CRS("+proj=longlat +datum=WGS84")

sp_points_proj <- sp::spTransform(sp_points, nad83_proj)

# Organise coordinates for spatial models
coords <- sp_points_proj@coords[,1:2]
names(coords) <- c("X", "Y")

# Add site variable to detection covariates
det_covs$Site <- 1:nrow(coords)

# Add year variable to detection covariates
det_covs$Year <- occ_covs$Year

################################################################################
# 3. COMBINE AND SAVE
################################################################################

# Combined data into correct format depending on single or multi season approach
if(is.na(bin) == T){
  sp.data <- Array_prep(target, sp = TRUE)
}else{
  sp.data <- list(y = eh_list[[1]][[1]], 
                occ.covs = occ_covs, 
                det.covs = det_covs,
                coords = coords)
}

# Save data
if(is.na(bin) == T){
  saveRDS(sp.data, 
          file = paste("Prepped_data/spOccupancy/Multi_season/", res, "/", target, 
                       "_multi_", res, ".rds", sep = ""))
}else{
  saveRDS(sp.data, 
          file = paste("Prepped_data/spOccupancy/Single_season/", res, "/", 
                       bin, "/", target, "_single_", res, ".rds", sep = ""))
}

