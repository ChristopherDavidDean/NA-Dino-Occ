## -----------------------------------------------------------------------------
##
## Script name: prepare-climate-data.R
##
## Purpose of script: Prepare climate data for analyses
##
## Author: Dr Lewis Jones, edited by Christopher D. Dean
##
## Last update: 2022-09-16
##

##### Load packages #####
library(raster)
library(ncdf4)
library(stringr)

##### Setup data #####
# Grab months and convert to lowercase
months <- c(tolower(month.abb))

# Set variables
e <- extent(-180, 180, -90, 90)
res <- 0.5
bin <- "teyeq"

# Define raster for resampling
r <- raster(res = res, ext = e)

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
  # Make folders
  dir.create(paste0("Data/Covariate_Data/Lewis_Climate/Results/", sep =""), 
             showWarnings = FALSE)
  if(bin == "teyen" | bin == "teyeo"){
    dir.create(paste0("Data/Covariate_Data/Lewis_Climate/Results/Maastrichtian/", 
                      sep =""), showWarnings = FALSE)
    dir.create(paste0("Data/Covariate_Data/Lewis_Climate/Results/Maastrichtian/", 
                      bin, "/", sep =""), showWarnings = FALSE)
    # Save data
    raster::writeRaster(x = stk, format = "raster",
                        filename = paste0("Data/Covariate_Data/Lewis_Climate/Results/Maastrichtian/",
                                          bin, "/", names(stk), ".grd", sep = ""), 
                        bylayer = TRUE, overwrite = TRUE)
  }else{
    dir.create(paste0("Data/Covariate_Data/Lewis_Climate/Results/Campanian/", 
                      sep =""), showWarnings = FALSE)
    dir.create(paste0("Data/Covariate_Data/Lewis_Climate/Results/Campanian/", 
                      bin, "/", sep =""), showWarnings = FALSE)
    # Save data
    raster::writeRaster(x = stk, format = "raster",
                        filename = paste0("Data/Covariate_Data/Lewis_Climate/Results/Campanian/",
                                          bin, "/", names(stk), ".grd", sep = ""), 
                        bylayer = TRUE, overwrite = TRUE)
  }
}

##### Run function #####
format.p.covs(bin)
