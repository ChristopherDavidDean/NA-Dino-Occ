################################################################################
############ OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS #############
################################################################################

# Christopher D. Dean, Alfio Alessandro Chiarenza, Jeffrey W. Doser, Alexander
# Farnsworth, Lewis A. Jones, Sinéad Lyster, Charlotte L. Outhwaite, Richard J. 
# Butler, Philip D. Mannion.
# 2024
# Script written by Christopher D. Dean and Lewis A. Jones

################################################################################
#         FILE 0: SETTING UP COVARIATE DATA FOR OCCUPANCY MODELLING            #
################################################################################

################################################################################
# 1. INITIAL SETUP
################################################################################

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

# Load necessary packages
library(raster)
library(ncdf4)
library(sf)
library(terra)
library(dplyr)
library(beepr)

# Load in Functions
source("0.Functions.R")

# Set values
res <- 0.5
e <- extent(-155, -72, 22.5, 73) # For 0.5 degree
e <- extent(-155, -72, 23, 73) # For 1 degree

# Setup Folders
dir.create(paste0("Prepped_data", sep =""))
dir.create(paste0("Prepped_data/Covariate_Data", sep =""), showWarnings = FALSE)
dir.create(paste0("Prepped_data/Covariate_Data/PalaeoPrecip", sep =""), showWarnings = FALSE)
dir.create(paste0("Prepped_data/Covariate_Data/PalaeoTemp", sep =""), showWarnings = FALSE)
dir.create(paste0("Prepped_data/Covariate_Data/PalaeoDEM", sep =""), showWarnings = FALSE)
dir.create(paste0("Prepped_data/Covariate_Data/PalaeoSed", sep =""), showWarnings = FALSE)
dir.create(paste0("Prepped_data/Covariate_Data/All_data", sep =""), showWarnings = FALSE)
dir.create(paste0("Prepped_data/Covariate_Data/All_data/0.1deg", sep =""), showWarnings = FALSE)
dir.create(paste0("Prepped_data/Covariate_Data/All_data/0.5deg", sep =""), showWarnings = FALSE)
dir.create(paste0("Prepped_data/Covariate_Data/All_data/1deg", sep =""), showWarnings = FALSE)
dir.create(paste0("Prepped_data/Covariate_Data/All_data/0.1deg/Palaeo", sep =""), showWarnings = FALSE)
dir.create(paste0("Prepped_data/Covariate_Data/All_data/0.5deg/Palaeo", sep =""), showWarnings = FALSE)
dir.create(paste0("Prepped_data/Covariate_Data/All_data/1deg/Palaeo", sep =""), showWarnings = FALSE)
dir.create(paste0("Prepped_data/Covariate_Data/DEM", sep =""), showWarnings = FALSE)
dir.create(paste0("Prepped_data/Covariate_Data/landcvi0201", sep =""), showWarnings = FALSE)
dir.create(paste0("Prepped_data/Covariate_Data/MGVF", sep =""), showWarnings = FALSE)
dir.create(paste0("Prepped_data/Covariate_Data/Outcrop", sep =""), showWarnings = FALSE)
dir.create(paste0("Prepped_data/Covariate_Data/Worldclim", sep =""), showWarnings = FALSE)
dir.create(paste0("Prepped_data/Covariate_Data/PalaeoClimate", sep =""), showWarnings = FALSE)
dir.create(paste0("Prepped_data/Covariate_Data/Relief", sep =""), showWarnings = FALSE)
dir.create(paste0("Prepped_data/Covariate_Data/SedFlux", sep =""), showWarnings = FALSE)

################################################################################
# 2. FORMAT DEM
################################################################################

# Load raster
dem <- raster("Data/Covariate_Data/Elevation_GRID/NA_Elevation.asc") 

# Check projection and update. Follow:
# https://gis.stackexchange.com/questions/291256/reprojecting-raster-between-laea-and-lon-lat-alignment-issues
raster::projection(dem) 
raster::projection(dem) <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs" 
plot(dem)

# Create raster for projecting raster, then project and crop to extent
r1 <- raster(res = res) 
newData <- raster::projectRaster(dem, r1) 
newData <- crop(newData, e) 

# Plot
plot(newData) 

# Save
writeRaster(newData, paste("Prepped_data/Covariate_Data/DEM/DEM_", 
                           res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) 
writeRaster(newData, paste("Prepped_data/Covariate_Data/All_data/", 
                           res, "deg/DEM_", res, ".asc", sep = ""), 
            pattern = "ascii", overwrite = TRUE)

# Add Slope and save
newData <- raster::terrain(newData, opt='slope', unit='degrees', neighbors=8)
plot(newData)
writeRaster(newData, paste("Prepped_data/Covariate_Data/DEM/SLOPE_", 
                           res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) 
writeRaster(newData, paste("Prepped_data/Covariate_Data/All_data/", res, "deg/SLOPE_", 
                           res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) 

# Add Relief
res1 <- res/10
r1 <- raster(res = res1)
newData <- raster::projectRaster(dem, r1) 
test1 <- aggregate(newData, 10, fun = max)
test2 <- aggregate(newData, 10, fun = min)
Relief <- (test1 - test2)
Relief <- crop(Relief, e)
plot(Relief)
writeRaster(Relief, paste("Prepped_data/Covariate_Data/Relief/Relief_", 
                          res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) 
writeRaster(Relief, paste("Prepped_data/Covariate_Data/All_data/", 
                          res, "deg/Relief_", res,".asc", sep = ""), pattern = "ascii", 
            overwrite = TRUE) 

################################################################################
# 3. FORMAT LANDCVI
################################################################################

##########################
##### TWO CATEGORIES #####
##########################

# Load raster
landcvi <- raster::raster("Data/Covariate_Data/landcvi0201/landcvi0201.tif")

# Make vector of unwanted landuses
unwanted <- c(1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 22, 24, 255)
for (i in 1:length(unwanted)){
  landcvi[landcvi == unwanted[i]] <- NA
}
landcvi[landcvi == 7] <- 1
landcvi[landcvi == 8] <- 1
landcvi[landcvi == 9] <- 1
landcvi[landcvi == 19] <- 1
landcvi[landcvi == 10] <- 1
plot(landcvi)
raster::projection(landcvi)
raster::projection(landcvi) <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs" 
r1 <- raster(res = res) 
newData <- raster::projectRaster(landcvi, r1) #project raster
newData <- crop(newData, e) #crop data to extent object
plot(newData) #plot data
writeRaster(newData, paste("Prepped_data/Covariate_Data/landcvi0201/LANDCVI_binary_", 
                           res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
writeRaster(newData, paste("Prepped_data/Covariate_Data/All_data/", res, 
                           "deg/LANDCVI_binary_", res, ".asc", sep = ""), 
            pattern = "ascii", overwrite = TRUE) 

###############################
##### MULTIPLE CATEGORIES #####
###############################

# Load raster
landcvi <- raster::raster("Data/Covariate_Data/landcvi0201/landcvi0201.tif")

# Make vectors of landuses
human_altered <- c(1, 2, 3, 4, 5, 6)
open_terrain <- c(7, 8, 9, 10, 19, 17, 20, 23)
forest <- c(11, 12, 13, 14, 15, 18, 21)
other <- c(16, 24, 255, 22)

for (i in 1:length(human_altered)){
  landcvi[landcvi == human_altered[i]] <- 1
}
for (i in 1:length(open_terrain)){
  landcvi[landcvi == open_terrain[i]] <- 2
}
for (i in 1:length(forest)){
  landcvi[landcvi == forest[i]] <- 3
}
for (i in 1:length(other)){
  landcvi[landcvi == other[i]] <- 4
}

plot(landcvi)
raster::projection(landcvi) <-"+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs" 
newData <- raster::projectRaster(landcvi, r1, method="ngb") 
newData <- crop(newData, e) 
plot(newData)

writeRaster(newData, paste("Prepped_data/Covariate_Data/landcvi0201/LANDCVI_multiple_", 
                           res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
writeRaster(newData, paste("Prepped_data/Covariate_Data/All_data/", res, 
                           "deg/LANDCVI_multiple_", res, ".asc", sep = ""), 
            pattern = "ascii", overwrite = TRUE) 

################################################################################
# 4. FORMAT MGVF
################################################################################

MGVF <- raster("Data/Covariate_Data/MGVF/Average MGVF.tif")
plot(MGVF)
raster::projection(MGVF) #no projection alteration needed
newData <- raster::crop(MGVF, e) #crop data to extent object
r1 <- raster(res = res, ext = e) 
newData <- raster::projectRaster(newData, r1) #project raster
plot(newData) #plot data
writeRaster(newData, paste("Prepped_data/Covariate_Data/MGVF/MGVF_", 
                           res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
writeRaster(newData, paste("Prepped_data/Covariate_Data/All_data/",
                           res, "deg/MGVF_", res, ".asc", sep = ""), pattern = "ascii", 
            overwrite = TRUE) 

################################################################################
# 5. FORMAT WORLDCLIM
################################################################################

##### Precipitation ######
setwd("Data/Covariate_Data/Worldclim/wc2.0_30s_prec/") 
wc <- list.files("./")
wc <- stack(wc)
wc <- crop(wc, e) #crop data to extent object
fact <- res * 120
wc2 <- aggregate(wc, fact = fact, fun = mean)
#wc <- resample(wc, r)
#wc <- crop(wc, e)
newData <- wc2
newData <- mean(newData)
plot(newData)
setwd("/Users/christopherdean/MEGA/Research/PROJECTS/DINO_RANGE/Fresh Setup/Code Data and Github/NA-Dino-Occ/NA-Dino-Occ/")
writeRaster(newData, filename = (paste("Prepped_data/Covariate_Data/Worldclim/WC_Prec_",
                                       res, ".asc", sep = "")), pattern = "ascii",
            overwrite = TRUE)
writeRaster(newData, paste("Prepped_data/Covariate_Data/All_data/", res, 
                           "deg/WC_Prec_", res, ".asc", sep = ""), 
            pattern = "ascii", overwrite = TRUE)

##### Temperature ######
setwd("Data/Covariate_Data/Worldclim/wc2.0_30s_tavg/")
wc <- list.files("./")
wc <- stack(wc)
wc <- crop(wc, e) #crop data to extent object
fact <- res * 120
wc2 <- aggregate(wc, fact = fact, fun = mean)
#wc <- resample(wc, r)
#wc <- crop(wc, e)
newData <- mean(wc2)
plot(newData)
setwd("/Users/christopherdean/MEGA/Research/PROJECTS/DINO_RANGE/Fresh Setup/Code Data and Github/NA-Dino-Occ/NA-Dino-Occ/")
writeRaster(newData, paste("Prepped_data/Covariate_Data/Worldclim/WC_Temp_", 
                           res, ".asc", sep = ""), pattern = "ascii", 
            overwrite = TRUE)
writeRaster(newData, paste("Prepped_data/Covariate_Data/All_data/", res, 
                           "deg/WC_Temp_", res, ".asc", sep = ""), 
            pattern = "ascii", overwrite = TRUE)

################################################################################
# 6. FORMAT SCOTESE PALAEO DATA (PRECIP, TEMP, SEDFLUX)
################################################################################

#################
##### SETUP #####
#################

##### Load packages #####
library(raster)
library(ncdf4)
library(stringr)

##### Setup data #####
# Grab months and convert to lowercase
months <- c(tolower(month.abb))

# Set variables
new_extent <- extent(-180, 180, -90, 90)

# Define raster for resampling
r <- raster(res = res, ext = new_extent)

# Load function
format.p.covs <- function(bin){
  if(bin == "teyen" | bin == "teyeo"){
    # Get file names
    files <- list.files(paste("Data/Covariate_Data/Lewis_Climate/Maastrichtian/",
                              bin, "/", sep = ""), full.names = TRUE)
  }else{
    # Get file names
    files <- list.files(paste("Data/Covariate_Data/Lewis_Climate/Campanian/",
                              bin, "/", sep = ""), full.names = TRUE)
  }
  # Extract only monthly variables
  files <- files[sapply(months, function(x){str_which(files, x)})]
  ##### PRECIP ######
  # Load precipitation data
  precip <- stack(files, varname = "precip_mm_srf")
  # Mask data
  #precip <- mask(x = precip, mask = msk, maskvalue = 0)
  # Assign names
  names(precip) <- months
  # Convert from kg/m2/s to mm/day
  precip <- precip * 86400
  # Driest/wettest quarters:
  dry <- precip[[c(12,1,2)]]
  wet <- precip[[c(6,7,8)]]
  ##### TEMP #####
  # Load temperature data
  temp <- stack(files, varname = "temp_mm_1_5m")
  # Mask data
  #temp <- mask(x = temp, mask = msk, maskvalue = 0)
  # Assign names
  names(temp) <- months
  # Convert kelvin to celsius
  temp <- temp - 273.15
  # Hottest/coolest months
  col <- temp[[c(12,1,2)]]
  hot <- temp[[c(6,7,8)]]
  ##### CALCS #####
  col.mean <- calc(x = col, fun = mean)
  hot.mean <- calc(x = hot, fun = mean)
  wet.mean <- calc(x = wet, fun = mean)
  dry.mean <- calc(x = dry, fun = mean)
  ann.SD <- calc(x = temp, fun = sd)
  # Stack data
  stk <- stack(col.mean, hot.mean, wet.mean, dry.mean, ann.SD)
  # Add names
  names(stk) <- c("col_mean", "hot_mean", "wet_mean", "dry_mean", "ann_sd")
  # Rotate data
  stk <- raster::rotate(stk)
  # Resample data (the extent and resolution must be updated to avoid rgdal issue)
  stk <- resample(x = stk, y = r)
  # Define original GCRS
  crs(stk) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  # Plot data
  plot(stk)
  # Save data
  raster::writeRaster(x = stk, format = "raster", 
                      filename = paste0("Prepped_data/Covariate_Data/All_data/", 
                                        res, "deg/Palaeo/", bin, "_", names(stk), 
                                        "_", res, ".asc", sep = ""), 
                      pattern = "ascii", bylayer = TRUE, overwrite = TRUE)
}

###########################
##### PRECIP AND TEMP #####
###########################

##### Run function #####
format.p.covs("teyen")
format.p.covs("teyeo")
format.p.covs("teyep")
format.p.covs("teyeq")

#########################
##### SEDIMENT FLUX #####
#########################

##### CAMPANIAN #####
# teyep #
teyep_sed <- sf::st_read("Data/Covariate_Data/catchment_data_080223/teyep/",
                       layer = "m18_watersheds")
teyep_sed <- teyep_sed %>%
  select(qs_m3yr) 
sf::st_crs(teyep_sed) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
r1 <- raster(res = res)
teyep_sed <- rasterize(teyep_sed, r1)
teyep_sed <- crop(teyep_sed, new_extent)
writeRaster(teyep_sed, paste("Prepped_data/Covariate_Data/PalaeoSed/teyep_sed_", 
                                res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
writeRaster(teyep_sed, paste("Prepped_data/Covariate_Data/All_data/", res, "deg/Palaeo/teyep_sed_", 
                                res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
# teyeq #
teyeq_sed <- sf::st_read("Data/Covariate_Data/catchment_data_080223/teyeq",
                         layer = "m19_watersheds")
teyeq_sed <- teyeq_sed %>%
  select(qs_m3yr) 
sf::st_crs(teyeq_sed) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
r1 <- raster(res = res)
teyeq_sed <- rasterize(teyeq_sed, r1)
teyeq_sed <- crop(teyep_sed, new_extent)
writeRaster(teyeq_sed, paste("Prepped_data/Covariate_Data/PalaeoSed/teyeq_sed_", 
                             res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
writeRaster(teyeq_sed, paste("Prepped_data/Covariate_Data/All_data/", res, "deg/Palaeo/teyeq_sed_", 
                             res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
##### MAASTRICHTIAN #####
# teyen #
teyen_sed <- sf::st_read("Data/Covariate_Data/catchment_data_080223/teyen",
                         layer = "m16_watersheds")
teyen_sed <- teyen_sed %>%
  select(qs_m3yr) 
sf::st_crs(teyen_sed) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
r1 <- raster(res = res)
teyen_sed <- rasterize(teyen_sed, r1)
teyen_sed <- crop(teyen_sed, new_extent)
writeRaster(teyen_sed, paste("Prepped_data/Covariate_Data/PalaeoSed/teyen_sed_", 
                             res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
writeRaster(teyen_sed, paste("Prepped_data/Covariate_Data/All_data/", res, "deg/Palaeo/teyen_sed_", 
                             res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
# teyeo #
teyeo_sed <- sf::st_read("Data/Covariate_Data/catchment_data_080223/teyeo",
                         layer = "m17_watersheds")
teyeo_sed <- teyeo_sed %>%
  select(qs_m3yr) 
sf::st_crs(teyeo_sed) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
r1 <- raster(res = res)
teyeo_sed <- rasterize(teyeo_sed, r1)
extent(teyeo_sed)
teyeo_sed <- crop(teyeo_sed, new_extent)
writeRaster(teyeo_sed, paste("Prepped_data/Covariate_Data/PalaeoSed/teyeo_sed_", 
                             res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
writeRaster(teyeo_sed, paste("Prepped_data/Covariate_Data/All_data/", res, "deg/Palaeo/teyeo_sed_", 
                             res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
beepr::beep(sound = 3)