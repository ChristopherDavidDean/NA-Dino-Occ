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

# Set working directory
setwd("C:/Users/deancd/Documents/RESEARCH/PROJECTS/DINO_RANGE/NA-Dino-Occ/") # Set your working directory

# Load in Functions
source("0.Functions_DR.R")

# Set extent and res - blocked out for running file direct from 1.Setup_occpancy
#get_extent(camp.occs)
#res <- 0.5

#=========================================== FORMAT DEM =================================================================

dem <- raster("Lewis_Occupancy_data/Data/Elevation_GRID/NA_Elevation.asc") #load raster
projection(dem) #check projection
#original projection is wrong and confused R. Update accordingly. https://gis.stackexchange.com/questions/291256/reprojecting-raster-between-laea-and-lon-lat-alignment-issues
projection(dem) <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs" 
r <- raster(res = res) #create raster for projecting raster
newData <- projectRaster(dem, r) #project raster
newData <- crop(newData, e) #crop data to extent object
plot(newData) #plot data
writeRaster(newData, paste("Lewis_Occupancy_data/Data/Formatted/Elevation_GRID/DEM_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii
writeRaster(newData, paste("Lewis_Occupancy_data/Data/Formatted/All_data/", res, "deg/DEM_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii

# Add Slope
newData <- terrain(newData, opt='slope', unit='degrees', neighbors=8)
plot(newData)
writeRaster(newData, paste("Lewis_Occupancy_data/Data/Formatted/Elevation_GRID/SLOPE_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii
writeRaster(newData, paste("Lewis_Occupancy_data/Data/Formatted/All_data/", res, "deg/SLOPE_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii

#======================================= FORMAT LANDCVI =================================================================

landcvi <- raster("Lewis_Occupancy_data/Data/landcvi0201/landcvi0201.tif")
projection(landcvi)
projection(landcvi) <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs" 
r <- raster(res = res) #create raster for projecting raster
newData <- projectRaster(landcvi, r, method = "ngb") #project raster
newData <- crop(newData, e) #crop data to extent object
plot(newData) #plot data
writeRaster(newData, paste("Lewis_Occupancy_data/Data/Formatted/landcvi0201/LANDCVI_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii
writeRaster(newData, paste("Lewis_Occupancy_data/Data/Formatted/All_data/", res, "deg/LANDCVI_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii

#================================== FORMAT AVERAGE MGVF =================================================================

MGVF <- raster("Lewis_Occupancy_data/Data/MGVF/Average MGVF.tif")
MGVF <- raster("Lewis_Occupancy_data/Data/MGVF/Average MGVF.tif")
plot(MGVF)
projection(MGVF) #no projection alteration needed
r <- raster(res = res) #create raster for projecting raster
newData <- projectRaster(MGVF, r) #project raster
newData <- crop(newData, e) #crop data to extent object
plot(newData) #plot data
writeRaster(newData, paste("Lewis_Occupancy_data/Data/Formatted/MGVF/MGVF_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii
writeRaster(newData, paste("Lewis_Occupancy_data/Data/Formatted/All_data/", res, "deg/MGVF_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii

#================================== FORMAT MODIS GLC =================================================================

GLC <- raster("Lewis_Occupancy_data/Data/MODIS_GLC/GlobalLandCover_tif/LCType.tif")
plot(GLC)
projection(GLC) #no projection alteration needed
r <- raster(res = res) #create raster for projecting raster
newData <- projectRaster(GLC, r, mathod = "ngb") #project raster
newData <- crop(newData, e) #crop data to extent object
plot(newData) #plot data
writeRaster(newData, paste("Lewis_Occupancy_data/Data/Formatted/MODIS_GLC/GLC_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii
writeRaster(newData, paste("Lewis_Occupancy_data/Data/Formatted/All_data/", res, "deg/GLC_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii

#===================================== FORMAT GHCN2  =================================================================

GHCN2 <- list.files("Lewis_Occupancy_data/Data/UoD_GHCN2", pattern = "nc")

for(i in 1:length(GHCN2)){
  name <- tools::file_path_sans_ext(GHCN2[i])
  GHC <- raster(paste("Lewis_Occupancy_data/Data/UoD_GHCN2/", GHCN2[i], sep = ""))
  GHC <- rotate(GHC)
  r <- raster(res = res) #create raster for projecting raster
  newData <- projectRaster(GHC, r) #project raster
  newData <- crop(newData, e) #crop data to extent object
  plot(newData) #plot data
  writeRaster(newData, paste("Lewis_Occupancy_data/Data/Formatted/UoD_GHCN2/", name, "_", res, ".asc", sep = ""), overwrite = TRUE)
  writeRaster(newData, paste("Lewis_Occupancy_data/Data/Formatted/All_data/", res, "deg/", name, "_", res, ".asc", sep = ""), overwrite = TRUE)
}

#===================================== FORMAT WORLDCLIM =================================================================

#wind
setwd("Lewis_Occupancy_data/Data/Worldclim/wc2.0_30s_wind/")
wc <- list.files("./")
wc <- stack(wc)
wc <- crop(wc, e) #crop data to extent object
wc <- resample(wc, r)
wc <- crop(wc, e)
newData <- mean(wc)
plot(newData)
writeRaster(newData, paste("Lewis_Occupancy_data/Data/Formatted/Worldclim/", "WC_Wind_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii
writeRaster(newData, paste("Lewis_Occupancy_data/Data/Formatted/All_data/", res, "deg/WC_Wind_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii

#precipitation
setwd("Lewis_Occupancy_data/Data/Worldclim/wc2.0_30s_prec/")
wc <- list.files("./")
wc <- stack(wc)
wc <- crop(wc, e) #crop data to extent object
wc <- resample(wc, r)
wc <- crop(wc, e)
newData <- wc
newData <- mean(newData)
plot(newData)
writeRaster(newData, paste("Lewis_Occupancy_data/Data/Formatted/Worldclim/WC_Prec_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii
writeRaster(newData, paste("Lewis_Occupancy_data/Data/Formatted/All_data/", res, "deg/WC_Prec_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii

#temperature
setwd("Lewis_Occupancy_data/Data/Worldclim/wc2.0_30s_tavg/")
wc <- list.files("./")
wc <- stack(wc)
wc <- crop(wc, e) #crop data to extent object
wc <- resample(wc, r)
wc <- crop(wc, e)
newData <- mean(wc)
plot(newData)
writeRaster(newData, paste("Lewis_Occupancy_data/Data/Formatted/Worldclim/WC_Temp_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii
writeRaster(newData, paste("Lewis_Occupancy_data/Data/Formatted/All_data/", res, "deg/WC_Temp_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii

#========================================= OUTCROP AREA ===================================================================

# Campanian #
CampOut <- readOGR(dsn = "Lewis_Occupancy_data/Data/Outcrop/Campanian/Campanian.shp")
CampOut <- spTransform(CampOut,  CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
CampOut <- CampOut[grepl(c("Sedimentary"), CampOut$ROCKTYPE), ]
#plot(CampOut)
r <- raster(ext = e, res = res)
CampOut <- rasterize(CampOut, r, getCover = TRUE)
extent(CampOut)
plot(CampOut)
writeRaster(CampOut, paste("Lewis_Occupancy_data/Data/Formatted/Outcrop/Camp_out_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii
writeRaster(CampOut, paste("Lewis_Occupancy_data/Data/Formatted/All_data/", res, "deg/Camp_out_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)

# Maastichtian #
MaasOut <- readOGR(dsn = "Lewis_Occupancy_data/Data/Outcrop/Maastrichtian/Maastrichtian.shp")
MaasOut <- spTransform(MaasOut,  CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
MaasOut <- MaasOut[grepl(c("Sedimentary"), MaasOut$ROCKTYPE), ]
#plot(MaasOut)
r <- raster(ext = e, res = res)
MaasOut <- rasterize(MaasOut, r, getCover = TRUE)
plot(MaasOut)
writeRaster(MaasOut, paste("Lewis_Occupancy_data/Data/Formatted/Outcrop/Maas_out_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii
writeRaster(MaasOut, paste("Lewis_Occupancy_data/Data/Formatted/All_data/", res, "deg/Maas_out_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)

#========================================= PALAEO DATA ================================================================

# Precip
CampPrecip <- raster("Lewis_Occupancy_data/Data/Climate_Data/Camp_Precip.asc")
projection(CampPrecip) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
r <- raster(res = res) #create raster for projecting raster
CampPrecip <- projectRaster(CampPrecip, r) #project raster
plot(CampPrecip)
writeRaster(CampPrecip, paste("Lewis_Occupancy_data/Data/Formatted/PalaeoClimate/CampPrecip_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii
writeRaster(CampPrecip, paste("Lewis_Occupancy_data/Data/Formatted/All_data/", res, "deg/CampPrecip_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)

MaasPrecip <- raster("Lewis_Occupancy_data/Data/Climate_Data/Maas_Precip.asc")
projection(MaasPrecip) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
r <- raster(res = res) #create raster for projecting raster
MaasPrecip <- projectRaster(MaasPrecip, r) #project raster
plot(MaasPrecip)
writeRaster(MaasPrecip, paste("Lewis_Occupancy_data/Data/Formatted/PalaeoClimate/MaasPrecip_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii
writeRaster(MaasPrecip, paste("Lewis_Occupancy_data/Data/Formatted/All_data/", res, "deg/MaasPrecip_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)

# Temp
CampTemp <- raster("Lewis_Occupancy_data/Data/Climate_Data/Camp_Temp.asc")
projection(CampTemp) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
r <- raster(res = res) #create raster for projecting raster
CampTemp <- projectRaster(CampTemp, r) #project raster
plot(CampTemp)
writeRaster(CampTemp, paste("Lewis_Occupancy_data/Data/Formatted/PalaeoClimate/CampTemp_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii
writeRaster(CampTemp, paste("Lewis_Occupancy_data/Data/Formatted/All_data/", res, "deg/CampTemp_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)

MaasTemp <- raster("Lewis_Occupancy_data/Data/Climate_Data/Maas_Temp.asc")
projection(MaasTemp) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
r <- raster(res = res) #create raster for projecting raster
MaasTemp <- projectRaster(MaasTemp, r) #project raster
plot(MaasTemp)
writeRaster(MaasTemp, paste("Lewis_Occupancy_data/Data/Formatted/PalaeoClimate/MaasTemp_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE) #write ascii
writeRaster(MaasTemp, paste("Lewis_Occupancy_data/Data/Formatted/All_data/", res, "deg/MaasTemp_", res, ".asc", sep = ""), pattern = "ascii", overwrite = TRUE)