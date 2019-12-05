#========================================= COVARIATES FOR OCCUPANCY MODELLING  ===================================================#
#                                                                                                                                 #
#      PRIMARY AUTHOR: LEWIS A. JONES                                                                                             #
#      CO-AUTHOR: CHRISTOPHER D. DEAN                                                                                             #
#                                                                                                                                 #    
#=================================================================================================================================#

#========================================== REQUIRED PACKAGES & SET-UP =================================================

library(raster)
library(ncdf4)
library(ncdf.tools)
library(rgdal)
library(dplyr)
library(beepr)

# Set working directory
setwd("C:/Users/deancd/Documents/RESEARCH/PROJECTS/DINO_RANGE/NA-Dino-Occ/") # Set your working directory

# Load in Functions
source("0.Functions_DR.R")

#=========================================== SETUP FOLDERS =================================================================

dir.create(paste0("Data/Covariate_Data/Formatted", sep =""), showWarnings = FALSE)
dir.create(paste0("Data/Covariate_Data/Formatted/All_data", sep =""), showWarnings = FALSE)
dir.create(paste0("Data/Covariate_Data/Formatted/All_data/0.1deg", sep =""), showWarnings = FALSE)
dir.create(paste0("Data/Covariate_Data/Formatted/All_data/0.5deg", sep =""), showWarnings = FALSE)
dir.create(paste0("Data/Covariate_Data/Formatted/All_data/1deg", sep =""), showWarnings = FALSE)
dir.create(paste0("Data/Covariate_Data/Formatted/Elevation_GRID", sep =""), showWarnings = FALSE)
dir.create(paste0("Data/Covariate_Data/Formatted/landcvi0201", sep =""), showWarnings = FALSE)
dir.create(paste0("Data/Covariate_Data/Formatted/MGVF", sep =""), showWarnings = FALSE)
dir.create(paste0("Data/Covariate_Data/Formatted/Outcrop", sep =""), showWarnings = FALSE)
dir.create(paste0("Data/Covariate_Data/Formatted/Worldclim", sep =""), showWarnings = FALSE)
dir.create(paste0("Data/Covariate_Data/Formatted/PalaeoClimate", sep =""), showWarnings = FALSE)

#=========================================== FORMAT DEM =================================================================

dem <- raster("Data/Covariate_Data/Elevation_GRID/NA_Elevation.asc") #load raster
projection(dem) #check projection
#original projection is wrong and confused R. Update accordingly. https://gis.stackexchange.com/questions/291256/reprojecting-raster-between-laea-and-lon-lat-alignment-issues
projection(dem) <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs" 
r <- raster(res = res) #create raster for projecting raster
newData <- projectRaster(dem, r) #project raster
newData <- crop(newData, e) #crop data to extent object
plot(newData) #plot data
writeRaster(newData, paste("Data/Covariate_Data/Formatted/Elevation_GRID/DEM_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii
writeRaster(newData, paste("Data/Covariate_Data/Formatted/All_data/", res, "deg/DEM_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii

# Add Slope
newData <- terrain(newData, opt='slope', unit='degrees', neighbors=8)
plot(newData)
writeRaster(newData, paste("Data/Covariate_Data/Formatted/Elevation_GRID/SLOPE_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii
writeRaster(newData, paste("Data/Covariate_Data/Formatted/All_data/", res, "deg/SLOPE_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii

#======================================= FORMAT LANDCVI =================================================================

landcvi <- raster("Data/Covariate_Data/landcvi0201/landcvi0201.tif")
unwanted <- c(1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 22, 24, 255)
for (i in 1:length(unwanted)){
  landcvi[landcvi == unwanted[i]] <- NA
}
landcvi[landcvi == 7] <- 1
landcvi[landcvi == 8] <- 1
landcvi[landcvi == 9] <- 1
landcvi[landcvi == 19] <- 1
landcvi[landcvi == 10] <- 1
writeRaster(landcvi, paste("Data/Covariate_Data/landcvi0201/LANDCVI_test", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii

landcvi <- raster("Data/Covariate_Data/landcvi0201/LANDCVI_test1.asc")
projection(landcvi)
projection(landcvi) <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs" 
r <- raster(res = res)#create raster for projecting raster
test1 <- aggregate(landcvi, 100, fun = 'sum')
newData <- projectRaster(test1, r) #project raster
newData <- crop(newData, e) #crop data to extent object
plot(newData) #plot data
writeRaster(newData, paste("Data/Covariate_Data/Formatted/landcvi0201/LANDCVI_selected_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii
writeRaster(newData, paste("Data/Covariate_Data/Formatted/All_data/", res, "deg/LANDCVI_selected_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii

test2 <- raster("Data/Covariate_Data/Formatted/All_data/0.5deg/DEM_0.5.asc")


#================================== FORMAT AVERAGE MGVF =================================================================

MGVF <- raster("Data/Covariate_Data/MGVF/Average MGVF.tif")
plot(MGVF)
projection(MGVF) #no projection alteration needed
r <- raster(res = res) #create raster for projecting raster
newData <- projectRaster(MGVF, r) #project raster
newData <- crop(newData, e) #crop data to extent object
plot(newData) #plot data
writeRaster(newData, paste("Data/Covariate_Data/Formatted/MGVF/MGVF_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii
writeRaster(newData, paste("Data/Covariate_Data/Formatted/All_data/", res, "deg/MGVF_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii

#===================================== FORMAT WORLDCLIM =================================================================

#wind
setwd("Data/Covariate_Data/Worldclim/wc2.0_30s_wind/")
wc <- list.files("./")
wc <- stack(wc)
wc <- crop(wc, e) #crop data to extent object
wc <- resample(wc, r)
wc <- crop(wc, e)
newData <- mean(wc)
plot(newData)
setwd("C:/Users/deancd/Documents/RESEARCH/PROJECTS/DINO_RANGE/NA-Dino-Occ/")
writeRaster(newData, paste("Data/Covariate_Data/Formatted/Worldclim/", "WC_Wind_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii
writeRaster(newData, paste("Data/Covariate_Data/Formatted/All_data/", res, "deg/WC_Wind_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii

#precipitation
setwd("Data/Covariate_Data/Worldclim/wc2.0_30s_prec/")
wc <- list.files("./")
wc <- stack(wc)
wc <- crop(wc, e) #crop data to extent object
wc <- resample(wc, r)
wc <- crop(wc, e)
newData <- wc
newData <- mean(newData)
plot(newData)
setwd("C:/Users/deancd/Documents/RESEARCH/PROJECTS/DINO_RANGE/NA-Dino-Occ/")
writeRaster(newData, filename = (paste("Data/Covariate_Data/Formatted/Worldclim/WC_Prec_", res, ".asc", sep = "")), pattern = "ascii", overwrite = TRUE) #write ascii
writeRaster(newData, paste("Data/Covariate_Data/Formatted/All_data/", res, "deg/WC_Prec_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii

#temperature
setwd("Data/Covariate_Data/Worldclim/wc2.0_30s_tavg/")
wc <- list.files("./")
wc <- stack(wc)
wc <- crop(wc, e) #crop data to extent object
wc <- resample(wc, r)
wc <- crop(wc, e)
newData <- mean(wc)
plot(newData)
setwd("C:/Users/deancd/Documents/RESEARCH/PROJECTS/DINO_RANGE/NA-Dino-Occ/")
writeRaster(newData, paste("Data/Covariate_Data/Formatted/Worldclim/WC_Temp_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii
writeRaster(newData, paste("Data/Covariate_Data/Formatted/All_data/", res, "deg/WC_Temp_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii

#========================================= OUTCROP AREA ===================================================================

# Campanian #
CampOut <- readOGR(dsn = "Data/Covariate_Data/Outcrop/Campanian/Campanian.shp")
CampOut <- spTransform(CampOut,  CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
CampOut <- CampOut[grepl(c("Sedimentary"), CampOut$ROCKTYPE), ]
#plot(CampOut)
e <- floor(e)
r <- raster(ext = e, res = res)
CampOut <- rasterize(CampOut, r, getCover = TRUE)
plot(CampOut)
writeRaster(CampOut, paste("Data/Covariate_Data/Formatted/Outcrop/Camp_out_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii
writeRaster(CampOut, paste("Data/Covariate_Data/Formatted/All_data/", res, "deg/Camp_out_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)

# Maastichtian #
MaasOut <- readOGR(dsn = "Data/Covariate_Data/Outcrop/Maastrichtian/Maastrichtian.shp")
MaasOut <- spTransform(MaasOut,  CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
MaasOut <- MaasOut[grepl(c("Sedimentary"), MaasOut$ROCKTYPE), ]
#plot(MaasOut)
e <- floor(e)
r <- raster(ext = e, res = res)
MaasOut <- rasterize(MaasOut, r, getCover = TRUE)
plot(MaasOut)
writeRaster(MaasOut, paste("Data/Covariate_Data/Formatted/Outcrop/Maas_out_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii
writeRaster(MaasOut, paste("Data/Covariate_Data/Formatted/All_data/", res, "deg/Maas_out_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)

#========================================= PALAEO DATA ================================================================

# Precip
CampPrecip <- raster("Data/Covariate_Data/Climate_Data/Camp_Precip.asc")
projection(CampPrecip) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
r <- raster(res = res) #create raster for projecting raster
CampPrecip <- projectRaster(CampPrecip, r) #project raster
plot(CampPrecip)
writeRaster(CampPrecip, paste("Data/Covariate_Data/Formatted/PalaeoClimate/CampPrecip_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii
writeRaster(CampPrecip, paste("Data/Covariate_Data/Formatted/All_data/", res, "deg/CampPrecip_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)

MaasPrecip <- raster("Data/Covariate_Data/Climate_Data/Maas_Precip.asc")
projection(MaasPrecip) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
r <- raster(res = res) #create raster for projecting raster
MaasPrecip <- projectRaster(MaasPrecip, r) #project raster
plot(MaasPrecip)
writeRaster(MaasPrecip, paste("Data/Covariate_Data/Formatted/PalaeoClimate/MaasPrecip_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii
writeRaster(MaasPrecip, paste("Data/Covariate_Data/Formatted/All_data/", res, "deg/MaasPrecip_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)

# Temp
CampTemp <- raster("Data/Covariate_Data/Climate_Data/Camp_Temp.asc")
projection(CampTemp) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
r <- raster(res = res) #create raster for projecting raster
CampTemp <- projectRaster(CampTemp, r) #project raster
plot(CampTemp)
writeRaster(CampTemp, paste("Data/Covariate_Data/Formatted/PalaeoClimate/CampTemp_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii
writeRaster(CampTemp, paste("Data/Covariate_Data/Formatted/All_data/", res, "deg/CampTemp_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)

MaasTemp <- raster("Data/Covariate_Data/Climate_Data/Maas_Temp.asc")
projection(MaasTemp) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
r <- raster(res = res) #create raster for projecting raster
MaasTemp <- projectRaster(MaasTemp, r) #project raster
plot(MaasTemp)
writeRaster(MaasTemp, paste("Data/Covariate_Data/Formatted/PalaeoClimate/MaasTemp_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii
writeRaster(MaasTemp, paste("Data/Covariate_Data/Formatted/All_data/", res, "deg/MaasTemp_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)
beepr::beep(sound = 3)