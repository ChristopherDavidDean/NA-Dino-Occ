################################################################################
############ OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS #############
################################################################################

# Christopher D. Dean, Alfio Alessandro Chiarenza, Jeffrey W. Doser, Alexander
# Farnsworth, Lewis A. Jones, Sinéad Lyster, Charlotte L. Outhwaite, Paul J. 
# Valdes, Richard J. Butler, Philip D. Mannion.
# 2025
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
library(tictoc)
library(car)

##### Set timer #####
tic("Full code")

##### Load in Functions #####
source("0.Functions.R") # Import functions from other R file (must be in same working directory)

##### Set values #####
# Set resolution
res <-  1
# Set extent
e <- extent(-155, -72, 22.5, 73)
# Set max limit value
max_val <- 10
max_val_on <- TRUE
target <- "Tyrannosauridae"
bins <- read.csv("Data/Occurrences/scotesebins.csv")
bins$bin <- bins$code
nomam <- F
bin <- NA # Change this if you want a single season model
form_cells <- "N" # Change this to "Y" for single season models with cells grouped 
                  # by formations

# Load occurrence dataset
if(nomam == F){
  master.occs.binned.targeted <- read.csv(paste("Prepped_data/Occurrence_Data/scotese/scotese_occurrence_dataset.csv", sep = ""))
}else{
  master.occs.binned.targeted <- read.csv(paste("Prepped_data/Occurrence_Data/scotese/scotese_no_mammal_occurrence_dataset.csv", sep = ""))
}

# Set random number for seed
rand <- round(runif(1, min = 0, max = 999)) 

################################################################################
# 2. PREPARING DATA
################################################################################

#########################
##### INITIAL SETUP #####
#########################

# Get grid cells for occurrences
get_grid(master.occs.binned.targeted, res = res, e = e, formCells = form_cells)

# Prep data 
prepare_for_spOcc(master.occs.binned.targeted.grid, single = FALSE, bin = bin)

# Get site IDs for sorting later
if(form_cells == "Y"){
  site_IDs <- as.numeric(gsub("([0-9]+).*$", "\\1", rownames(eh_list[[1]][[1]])))
}else{
  site_IDs <- as.numeric(rownames(eh_list[[1]][[1]]))
}

# Palaeo-rotate site IDs for each time bin to enable acquiring relevant covariate data
rotated <- p_rotate(res, e, site_IDs, bins, bin = bin)

# Extract palaeo-covariates
tic("Extracting palaeo-covariates")
extracted_covs <- extract_p(rotated)
toc()

# Format occupancy covariates for spOccupancy
occ_covs <- format_occ_covs(p_cov_list = extracted_covs)

# If generating single season data, reorganise occupancy covariates into suitable format
if(is.na(bin) == F){
  occ_covs <- lapply(occ_covs, unlist, use.names = F)
  occ_covs <- as.data.frame(bind_rows(!!! occ_covs))
  colnames(occ_covs) <- sub("$*", "", colnames(occ_covs))
}

# Format detection covariates for spOccupancy
det_covs <- organise_det(siteCoords, extracted_covs, 
                         occ_data = master.occs.binned.targeted.grid, 
                         bin = bin)

# Remove any sites without relevant data to ensure model runs
if(is.na(bin) == T){
  site_remove(eh_list, occ_covs, det_covs, siteCoords, single = FALSE)
}

# Transpose collected data into array format for spOccupancy
if(is.na(bin) == T){
  transpose_eh(eh_list, target)
}

# Add road distance data
Distance <- distance_fun(eh_list)
det_covs$Distance <- Distance

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
rownames(coords) <- 1:nrow(coords)
colnames(coords) <- c("X", "Y")

# Add site variable to detection covariates
det_covs$Site <- 1:nrow(coords)

# Add year variable to detection covariates
det_covs$Year <- occ_covs$Year

###################################
##### MULTICOLLINEARITY CHECK #####
###################################

##### OCCUPANCY #####
# select relevant covariates for testing
if(is.na(bin) == T){
  test <- within(occ_covs, rm(Year, site.effect))
}else{
  test <- within(occ_covs, rm(site.effect))
}

# Format and scale covariates if multi-season
if(is.na(bin) == T){
  occ_test <- bind_cols(lapply(names(test), function(nm) { 
    x <- as.data.frame(stack(test[[nm]])[,1])
    x <- scale(x)
    colnames(x) <- c(nm)
    return(x)
  }))
}else{
  occ_test <- test
}

# Make dummy response variable
occ_test$bogus <- scale(sample(100, size = nrow(occ_test), replace = TRUE))[,1]
# Make function to automatically remove high VIF
vif_fun <- function(occ_test, val){
  occ_test <- janitor::clean_names(occ_test)
  while(TRUE) {
    vifs <- car::vif(lm(bogus ~. , data = occ_test))
    if (max(vifs) < val) {
      break
    }
    highest <- c(names((which(vifs == max(vifs)))))
    occ_test <- occ_test[,-which(names(occ_test) %in% highest)]
  }
  return(occ_test)
}
# Run test
occ_test <- vif_fun(occ_test = occ_test, val = 10)
# Get names of relevant covs
occ_names <- colnames(occ_test)
# Remove dummy response variable
occ_names <- occ_names[! occ_names %in% c("bogus")]
# make relevant covs into appropriate format for running spOccupancy models
occ.form <- as.vector(sapply(occ_names, function(x){
  paste("scale(", x, ")", sep = "")
}))
occ.form <- str_replace(occ.form, "scale(year)", "factor(Year)")

##### DETECTION #####
# select relevant covariates for testing
test <- within(det_covs, rm(Distance, land, Site, Year))

det_test <- bind_cols(lapply(names(test), function(nm) {
  if(is.null(ncol(test[[nm]])) == T){
    x <- as.data.frame(rep(test[[nm]], 4))
  }else{
    x <- as.data.frame(stack(test[[nm]])[,1])
  }
  x <- scale(x)
  colnames(x) <- c(nm)
  return(x)
}))

# Make dummy response variable
det_test$bogus <- scale(sample(100, size = nrow(det_test), replace = TRUE))[,1]
# Run test
det_test <- vif_fun(occ_test = det_test, val = 10)
# Get names of relevant covs
det_names <- colnames(det_test)
# Remove dummy response variable
det_names <- det_names[! det_names %in% c("bogus")]
# make relevant covs into appropriate format for running spOccupancy models
det.form <- as.vector(sapply(det_names, function(x){
  paste("scale(", x, ")", sep = "")
}))
det.form <- str_replace(det.form, "mgvf", "MGVF")
det.form <- c(det.form, "factor(land)", "scale(Distance)", "factor(Year)")

################################################################################
# 3. COMBINE AND SAVE
################################################################################

# Combined data into correct format depending on single or multi-season approach
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
  if(nomam == F){
    saveRDS(sp.data,
            file = paste("Prepped_data/spOccupancy/Updated/Multi_season/", res, "/", target, 
                         "_multi_", res, "_", max_val, ".2.rds", sep = ""))
    saveRDS(occ.form, 
            file = paste("Prepped_data/spOccupancy/Updated/Multi_season/", res, "/", target, 
                         "_multi_", res, "_", max_val, ".2.occ.form.rds", sep = ""))
    saveRDS(det.form, 
            file = paste("Prepped_data/spOccupancy/Updated/Multi_season/", res, "/", target, 
                         "_multi_", res, "_", max_val, ".2.det.form.rds", sep = ""))
  }else{
    saveRDS(sp.data,
            file = paste("Prepped_data/spOccupancy/Updated/Multi_season/", res, "/", target, 
                         "_multi_", res, "_", max_val, "_no_mammal.rds", sep = ""))
    saveRDS(occ.form, 
            file = paste("Prepped_data/spOccupancy/Updated/Multi_season/", res, "/", target, 
                         "_multi_", res, "_", max_val, "_no_mammal.occ.form.rds", sep = ""))
    saveRDS(det.form, 
            file = paste("Prepped_data/spOccupancy/Updated/Multi_season/", res, "/", target, 
                         "_multi_", res, "_", max_val, "_no_mammal.det.form.rds", sep = ""))
  }
}else{
  if(form_cells == "Y"){
    saveRDS(sp.data, 
            file = paste("Prepped_data/spOccupancy/Single_season/", res, "/", 
                         bin, "/", target, "_single_", res, "_", max_val, ".formcells.rds", sep = ""))
  }else{
    saveRDS(sp.data, 
            file = paste("Prepped_data/spOccupancy/Single_season/", res, "/", 
                         bin, "/", target, "_single_", res, "_", max_val, ".rds", sep = ""))
    saveRDS(occ.form, 
            file = paste("Prepped_data/spOccupancy/Single_season/", res, "/", 
                         bin, "/", target, "_single_", res, "_", max_val, ".occ.form.rds", sep = ""))
    saveRDS(det.form, 
            file = paste("Prepped_data/spOccupancy/Single_season/", res, "/", 
                         bin, "/", target, "_single_", res, "_", max_val, ".det.form.rds", sep = ""))
  }
}

##### Close timer loop #####
toc()
