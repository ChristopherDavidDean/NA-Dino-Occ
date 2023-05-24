################################################################################
############# OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS ############
################################################################################

# Christopher D. Dean, Lewis A. Jones, Alfio A. Chiarenza, Sinéad Lyster, 
# Alex Farnsworth, Philip D. Mannion, Richard J. Butler.
# 2019
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
source("0.Functions_DR.R")

# Set values
res <- 1
e <- extent(-155, -72, 22.5, 73) # For 0.5 and 0.1 degree
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
# 6. FORMAT OUTCROP AREA
################################################################################

#####################
##### CAMPANIAN #####
#####################

##### TERRESTRIAL #####
CampMix <- sf::st_read("Data/Covariate_Data/Outcrop/Updated_shapefiles",
                    layer = "Campanian_Mixed")
CampTer <- sf::st_read("Data/Covariate_Data/Outcrop/Updated_shapefiles",
                    layer = "Campanian_Terrestrial")
CampOut <- rbind(CampMix, CampTer)
sf::st_crs(CampOut) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
r1 <- raster(ext = e, res = res)
newData <- raster::rasterize(CampOut, r1, getCover = TRUE, progress = "window")
writeRaster(newData, paste("Prepped_data/Covariate_Data/Outcrop/COut_", 
                           res, ".asc", sep = ""), pattern = "ascii", 
            overwrite = TRUE)
writeRaster(newData, paste("Prepped_data/Covariate_Data/All_data/", res, 
                           "deg/COut_", res, ".asc", sep = ""), 
            pattern = "ascii", overwrite = TRUE)

CampAll <- sf::st_read("Data/Covariate_Data/Outcrop/Updated_shapefiles",
                       layer = "Campanian_All")

##### ALL #####
sf::st_crs(CampAll) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
r1 <- raster(ext = e, res = res)
newData <- rasterize(CampAll, r1, getCover = TRUE, progress = "window")
plot(newData)
writeRaster(newData, paste("Prepped_data/Covariate_Data/Outcrop/COut_all_", 
                           res, ".asc", sep = ""), pattern = "ascii", 
            overwrite = TRUE)
writeRaster(newData, paste("Prepped_data/Covariate_Data/All_data/", res, 
                           "deg/COut_all_", res, ".asc", sep = ""), 
            pattern = "ascii", overwrite = TRUE)

#########################
##### MAASTRICHTIAN #####
#########################

##### TERRESTRIAL #####
MaasMix <- sf::st_read("Data/Covariate_Data/Outcrop/Updated_shapefiles",
                       layer = "Maastrichtian_Mixed")
MaasTer <- sf::st_read("Data/Covariate_Data/Outcrop/Updated_shapefiles",
                       layer = "Maastrichtian_Terrestrial")
MaasOut <- rbind(MaasMix, MaasTer)
sf::st_crs(MaasOut) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
r1 <- raster(ext = e, res = res)
newData <- raster::rasterize(MaasOut, r1, getCover = TRUE, progress = "window")
writeRaster(newData, paste("Prepped_data/Covariate_Data/Outcrop/MOut_", 
                           res, ".asc", sep = ""), pattern = "ascii", 
            overwrite = TRUE)
writeRaster(newData, paste("Prepped_data/Covariate_Data/All_data/", res, 
                           "deg/MOut_", res, ".asc", sep = ""), 
            pattern = "ascii", overwrite = TRUE)

##### ALL #####
MaasAll <- sf::st_read("Data/Covariate_Data/Outcrop/Updated_shapefiles",
                       layer = "Maastrichtian_All")
sf::st_crs(MaasAll) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
r1 <- raster(ext = e, res = res)
newData <- rasterize(MaasAll, r1, getCover = TRUE, progress = "window")
writeRaster(newData, paste("Prepped_data/Covariate_Data/Outcrop/MOut_all_", 
                           res, ".asc", sep = ""), pattern = "ascii", 
            overwrite = TRUE)
writeRaster(newData, paste("Prepped_data/Covariate_Data/All_data/", res, 
                           "deg/MOut_all_", res, ".asc", sep = ""), 
            pattern = "ascii", overwrite = TRUE)

#######################
##### ALL OUTCROP #####
#######################

CampAll <- sf::st_read("Data/Covariate_Data/Outcrop/Updated_shapefiles",
                       layer = "Campanian_All")
MaasAll <- sf::st_read("Data/Covariate_Data/Outcrop/Updated_shapefiles",
                       layer = "Maastrichtian_All")
CretAll <- rbind(CampAll, MaasAll)
sf::st_crs(CretAll) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
r1 <- raster(ext = e, res = res)
newData <- rasterize(CretAll, r1, getCover = TRUE, progress = "window")
plot(newData)
writeRaster(newData, paste("Prepped_data/Covariate_Data/Outcrop/AllOut_", 
                           res, ".asc", sep = ""), pattern = "ascii", 
            overwrite = TRUE)
writeRaster(newData, paste("Prepped_data/Covariate_Data/All_data/", res, 
                           "deg/AllOut_", res, ".asc", sep = ""), 
            pattern = "ascii", overwrite = TRUE)

################################################################################
# 6. FORMAT SCOTESE PALAEO DATA (DEM, PRECIP, TEMP, SEDFLUX)
################################################################################

###############
##### DEM #####
###############

##### CAMPANIAN #####
# teyep #
teyep_DEM <- raster("Data/Covariate_Data/Scotese_Wright_2018_Maps_1-88_6minX6min_PaleoDEMS_nc/Map18_PALEOMAP_6min_Late_Cretaceous_75Ma.nc",
            varname = "z")
crs(teyep_DEM) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
r <- raster(res = res) #create raster for projecting raster
teyep_DEM <- raster::projectRaster(teyep_DEM, r) #project raster
plot(teyep_DEM)
writeRaster(teyep_DEM, paste("Prepped_data/Covariate_Data/PalaeoDEM/teyep_DEM_", 
                            res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
writeRaster(teyep_DEM, paste("Prepped_data/Covariate_Data/All_data/", res, "deg/Palaeo/teyep_DEM_", 
                             res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
# teyeq #
teyeq_DEM <- raster("Data/Covariate_Data/Scotese_Wright_2018_Maps_1-88_6minX6min_PaleoDEMS_nc/Map19_PALEOMAP_6min_Late_Cretaceous_80Ma.nc",
                    varname = "z")
crs(teyeq_DEM) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
r <- raster(res = res) #create raster for projecting raster
teyeq_DEM <- raster::projectRaster(teyeq_DEM, r) #project raster
plot(teyeq_DEM)
writeRaster(teyeq_DEM, paste("Prepped_data/Covariate_Data/PalaeoDEM/teyeq_DEM_", 
                             res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
writeRaster(teyeq_DEM, paste("Prepped_data/Covariate_Data/All_data/", res, "deg/Palaeo/teyeq_DEM_", 
                             res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
##### MAASTRICHTIAN #####
# teyen #
teyen_DEM <- raster("Data/Covariate_Data/Scotese_Wright_2018_Maps_1-88_6minX6min_PaleoDEMS_nc/Map16_PALEOMAP_6min_KT_Boundary_65Ma.nc",
                    varname = "z")
crs(teyen_DEM) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
r <- raster(res = res) #create raster for projecting raster
teyen_DEM <- raster::projectRaster(teyen_DEM, r) #project raster
plot(teyen_DEM)
writeRaster(teyen_DEM, paste("Prepped_data/Covariate_Data/PalaeoDEM/teyen_DEM_", 
                             res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
writeRaster(teyen_DEM, paste("Prepped_data/Covariate_Data/All_data/", res, "deg/Palaeo/teyen_DEM_", 
                             res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
# teyeo #
teyeo_DEM <- raster("Data/Covariate_Data/Scotese_Wright_2018_Maps_1-88_6minX6min_PaleoDEMS_nc/Map17_PALEOMAP_6min_Late_Cretaceous_70Ma.nc",
                    varname = "z")
crs(teyeo_DEM) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
r <- raster(res = res) #create raster for projecting raster
teyeo_DEM <- raster::projectRaster(teyeo_DEM, r) #project raster
plot(teyeo_DEM)
writeRaster(teyeo_DEM, paste("Prepped_data/Covariate_Data/PalaeoDEM/teyeo_DEM_", 
                             res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
writeRaster(teyeo_DEM, paste("Prepped_data/Covariate_Data/All_data/", res, "deg/Palaeo/teyeo_DEM_", 
                             res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)

#########################
##### PRECIPITATION #####
#########################

##### CAMPANIAN #####
# teyep #
teyep_precip <- raster("Data/Covariate_Data/Lewis_Climate/Results/Campanian/teyep/max_precip.grd")
crs(teyep_precip) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
r <- raster(res = res) #create raster for projecting raster
teyep_precip <- raster::projectRaster(teyep_precip, r) #project raster
plot(teyep_precip)
writeRaster(teyep_precip, paste("Prepped_data/Covariate_Data/PalaeoPrecip/teyep_precip_", 
                             res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
writeRaster(teyep_precip, paste("Prepped_data/Covariate_Data/All_data/", res, "deg/Palaeo/teyep_precip_", 
                             res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
# teyeq #
teyeq_precip <- raster("Data/Covariate_Data/Lewis_Climate/Results/Campanian/teyeq/max_precip.grd")
crs(teyeq_precip) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
r <- raster(res = res) #create raster for projecting raster
teyeq_precip <- raster::projectRaster(teyeq_precip, r) #project raster
plot(teyeq_precip)
writeRaster(teyeq_precip, paste("Prepped_data/Covariate_Data/PalaeoPrecip/teyeq_precip_", 
                                res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
writeRaster(teyeq_precip, paste("Prepped_data/Covariate_Data/All_data/", res, "deg/Palaeo/teyeq_precip_", 
                                res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)

##### MAASTRICHTIAN #####
# teyen #
teyen_precip <- raster("Data/Covariate_Data/Lewis_Climate/Results/Maastrichtian/teyen/max_precip.grd")
crs(teyep_precip) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
r <- raster(res = res) #create raster for projecting raster
teyen_precip <- raster::projectRaster(teyen_precip, r) #project raster
plot(teyen_precip)
writeRaster(teyen_precip, paste("Prepped_data/Covariate_Data/PalaeoPrecip/teyen_precip_", 
                                res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
writeRaster(teyen_precip, paste("Prepped_data/Covariate_Data/All_data/", res, "deg/Palaeo/teyen_precip_", 
                                res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
# teyeo #
teyeo_precip <- raster("Data/Covariate_Data/Lewis_Climate/Results/Maastrichtian/teyeo/max_precip.grd")
crs(teyeo_precip) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
r <- raster(res = res) #create raster for projecting raster
teyeo_precip <- raster::projectRaster(teyeo_precip, r) #project raster
plot(teyeo_precip)
writeRaster(teyeo_precip, paste("Prepped_data/Covariate_Data/PalaeoPrecip/teyeo_precip_", 
                                res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
writeRaster(teyeo_precip, paste("Prepped_data/Covariate_Data/All_data/", res, "deg/Palaeo/teyeo_precip_", 
                                res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
#######################
##### TEMPERATURE #####
#######################

##### CAMPANIAN #####
# teyep #
teyep_temp <- raster("Data/Covariate_Data/Lewis_Climate/Results/Campanian/teyep/max_temp.grd")
crs(teyep_temp) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
r <- raster(res = res) #create raster for projecting raster
teyep_temp <- raster::projectRaster(teyep_temp, r) #project raster
plot(teyep_temp)
writeRaster(teyep_temp, paste("Prepped_data/Covariate_Data/PalaeoTemp/teyep_temp_", 
                                res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
writeRaster(teyep_temp, paste("Prepped_data/Covariate_Data/All_data/", res, "deg/Palaeo/teyep_temp_", 
                                res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
# teyeq #
teyeq_temp <- raster("Data/Covariate_Data/Lewis_Climate/Results/Campanian/teyeq/max_temp.grd")
crs(teyeq_temp) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
r <- raster(res = res) #create raster for projecting raster
teyeq_temp <- raster::projectRaster(teyeq_temp, r) #project raster
plot(teyeq_temp)
writeRaster(teyeq_temp, paste("Prepped_data/Covariate_Data/PalaeoTemp/teyeq_temp_", 
                                res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
writeRaster(teyeq_temp, paste("Prepped_data/Covariate_Data/All_data/", res, "deg/Palaeo/teyeq_temp_", 
                                res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
##### MAASTRICHTIAN #####
# teyen #
teyen_temp <- raster("Data/Covariate_Data/Lewis_Climate/Results/Maastrichtian/teyen/max_temp.grd")
crs(teyep_temp) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
r <- raster(res = res) #create raster for projecting raster
teyen_temp <- raster::projectRaster(teyen_temp, r) #project raster
plot(teyen_temp)
writeRaster(teyen_temp, paste("Prepped_data/Covariate_Data/PalaeoTemp/teyen_temp_", 
                                res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
writeRaster(teyen_temp, paste("Prepped_data/Covariate_Data/All_data/", res, "deg/Palaeo/teyen_temp_", 
                                res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
# teyeo #
teyeo_temp <- raster("Data/Covariate_Data/Lewis_Climate/Results/Maastrichtian/teyeo/max_temp.grd")
crs(teyeo_temp) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
r <- raster(res = res) #create raster for projecting raster
teyeo_temp <- raster::projectRaster(teyeo_temp, r) #project raster
plot(teyeo_temp)
writeRaster(teyeo_temp, paste("Prepped_data/Covariate_Data/PalaeoTemp/teyeo_temp_", 
                                res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
writeRaster(teyeo_temp, paste("Prepped_data/Covariate_Data/All_data/", res, "deg/Palaeo/teyeo_temp_", 
                                res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)

#########################
##### SEDIMENT FLUX #####
#########################

new_extent <- extent(-180, 180, -90, 90)

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

################################################################################
# A1. FORMAT GETECH PALAEO DATA (PRECIP, TEMP)
################################################################################

#########################
##### PRECIPITATION #####
#########################

##### CAMPANIAN #####
CampPrecip <- raster::raster("Data/Covariate_Data/Climate_Data/Camp_Precip.asc")
crs(CampPrecip) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
r <- raster(res = res) #create raster for projecting raster
CampPrecip <- raster::projectRaster(CampPrecip, r) #project raster
plot(CampPrecip)
writeRaster(CampPrecip, paste("Data/Covariate_Data/Formatted/PalaeoClimate/CampPrecip_", 
                              res, ".asc", sep = ""), pattern = "ascii", 
            overwrite = TRUE)
writeRaster(CampPrecip, paste("Data/Covariate_Data/Formatted/All_data/", res, 
                              "deg/PalaeoClimate/CampPrecip_", res, ".asc", sep = ""), 
            pattern = "ascii", overwrite = TRUE)

##### MAASTRICHTIAN #####
MaasPrecip <- raster("Data/Covariate_Data/Climate_Data/Maas_Precip.asc")
crs(MaasPrecip) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
r <- raster(res = res) #create raster for projecting raster
MaasPrecip <- raster::projectRaster(MaasPrecip, r) #project raster
plot(MaasPrecip)
writeRaster(MaasPrecip, paste("Data/Covariate_Data/Formatted/PalaeoClimate/MaasPrecip_", 
                              res, ".asc", sep = ""), pattern = "ascii", 
            overwrite = TRUE)
writeRaster(MaasPrecip, paste("Data/Covariate_Data/Formatted/All_data/", res,
                              "deg/PalaeoClimate/MaasPrecip_", res, ".asc", sep = ""), 
            pattern = "ascii", overwrite = TRUE)

#######################
##### TEMPERATURE #####
#######################

##### CAMPANIAN #####
CampTemp <- raster::raster("Data/Covariate_Data/Climate_Data/Camp_Temp.asc")
crs(CampTemp) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
r <- raster(res = res) #create raster for projecting raster
CampTemp <- raster::projectRaster(CampTemp, r) #project raster
plot(CampTemp)
writeRaster(CampTemp, paste("Data/Covariate_Data/Formatted/PalaeoClimate/CampTemp_", 
                            res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
writeRaster(CampTemp, paste("Data/Covariate_Data/Formatted/All_data/", 
                            res, "deg/PalaeoClimate/CampTemp_", res, ".asc", sep = ""), 
            pattern = "ascii", overwrite = TRUE)

##### MAASTRICHTIAN #####
MaasTemp <- raster("Data/Covariate_Data/Climate_Data/Maas_Temp.asc")
crs(MaasTemp) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
r <- raster(res = res) #create raster for projecting raster
MaasTemp <- projectRaster(MaasTemp, r) #project raster
plot(MaasTemp)
writeRaster(MaasTemp, paste("Data/Covariate_Data/Formatted/PalaeoClimate/MaasTemp_", 
                            res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) 
writeRaster(MaasTemp, paste("Data/Covariate_Data/Formatted/All_data/", 
                            res, "deg/PalaeoClimate/MaasTemp_", res, ".asc", sep = ""), 
            pattern = "ascii", overwrite = TRUE)