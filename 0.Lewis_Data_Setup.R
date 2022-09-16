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
# Load packages ----------------------------------------------------------------
library(raster)
library(ncdf4)
library(stringr)

# Data preparation -------------------------------------------------------------
# Grab months and convert to lowercase
months <- c(tolower(month.abb))

e <- extent(-180, 180, 0, 90)
# Define raster for resampling
r <- raster(res = res, ext = e)

#================================= Maastrichtian ===============================

#===== TEYEN ====

# Load mask 
msk <- raster(
  "Data/Covariate_Data/Lewis_Climate/Maastrichtian/teyen/teyen.qrparm.mask.nc", 
  varname = "lsm")
# Get file names
files <- list.files("Data/Covariate_Data/Lewis_Climate/Maastrichtian/teyen/", 
                    full.names = TRUE)
# Extract only monthly variables
files <- files[sapply(months, function(x){str_which(files, x)})]
# Load precipitation data
precip <- stack(files, varname = "precip_mm_srf")
# Mask data
precip <- mask(x = precip, mask = msk, maskvalue = 0)
# Assign names
names(precip) <- months
# Convert from kg/m2/s to mm/day
precip <- precip * 86400
# Load temperature data
temp <- stack(files, varname = "temp_mm_srf")
# Mask data
temp <- mask(x = temp, mask = msk, maskvalue = 0)
# Assign names
names(temp) <- months
# Convert kelvin to celsius
temp <- temp - 273.15
# Calculate max and min
max_precip <- calc(x = precip, fun = max)
min_precip <- calc(x = precip, fun = min)
max_temp <- calc(x = temp, fun = max)
min_temp <- calc(x = temp, fun = min)
# Stack data
stk <- stack(max_precip, min_precip, max_temp, min_temp)
# Add names
names(stk) <- c("max_precip", "min_precip", "max_temp", "min_temp")
# Rotate data
stk <- raster::rotate(stk)
# Resample data (the extent and resolution must be updated to avoid rgdal issue)
stk <- resample(x = stk, y = r)
# Define original GCRS
crs(stk) <- gcrs
# Plot data
plot(stk)
# Make folder
dir.create(paste0("Data/Covariate_Data/Lewis_Climate/Results/", sep =""), showWarnings = FALSE)
dir.create(paste0("Data/Covariate_Data/Lewis_Climate/Results/Maastrichtian/", sep =""), showWarnings = FALSE)
dir.create(paste0("Data/Covariate_Data/Lewis_Climate/Results/Maastrichtian/teyen/", sep =""), showWarnings = FALSE)
# Save data
writeRaster(
  x = stk,
  filename = paste0("Data/Covariate_Data/Lewis_Climate/Results/Maastrichtian/teyen/", 
                    names(stk), ".grd"), bylayer = TRUE, overwrite = TRUE)

#===== TEYEO ====

# Load mask 
msk <- raster(
  "Data/Covariate_Data/Lewis_Climate/Maastrichtian/teyeo/teyeo.qrparm.mask.nc", 
  varname = "lsm")
# Get file names
files <- list.files("Data/Covariate_Data/Lewis_Climate/Maastrichtian/teyeo/", 
                    full.names = TRUE)
# Extract only monthly variables
files <- files[sapply(months, function(x){str_which(files, x)})]
# Load precipitation data
precip <- stack(files, varname = "precip_mm_srf")
# Mask data
precip <- mask(x = precip, mask = msk, maskvalue = 0)
# Assign names
names(precip) <- months
# Convert from kg/m2/s to mm/day
precip <- precip * 86400
# Load temperature data
temp <- stack(files, varname = "temp_mm_srf")
# Mask data
temp <- mask(x = temp, mask = msk, maskvalue = 0)
# Assign names
names(temp) <- months
# Convert kelvin to celsius
temp <- temp - 273.15
# Calculate max and min
max_precip <- calc(x = precip, fun = max)
min_precip <- calc(x = precip, fun = min)
max_temp <- calc(x = temp, fun = max)
min_temp <- calc(x = temp, fun = min)
# Stack data
stk <- stack(max_precip, min_precip, max_temp, min_temp)
# Add names
names(stk) <- c("max_precip", "min_precip", "max_temp", "min_temp")
# Rotate data
stk <- raster::rotate(stk)
# Resample data (the extent and resolution must be updated to avoid rgdal issue)
stk <- resample(x = stk, y = r)
# Define original GCRS
crs(stk) <- gcrs
# Plot data
plot(stk)
# Make folder
dir.create(paste0("Data/Covariate_Data/Lewis_Climate/Results/Maastrichtian/teyeo/", sep =""), showWarnings = FALSE)
# Save data
writeRaster(
  x = stk,
  filename = paste0("Data/Covariate_Data/Lewis_Climate/Results/Maastrichtian/teyeo/", 
                    names(stk), ".grd"), bylayer = TRUE, overwrite = TRUE)

#================================= Campanian ===================================

#===== TEYEP====
# Load mask 
msk <- raster("Data/Covariate_Data/Lewis_Climate/Campanian/teyep/teyep.qrparm.mask.nc", 
              varname = "lsm")
# Get file names
files <- list.files("./data/raw-data/climate/Campanian/teyep/",
                    full.names = TRUE)
# Extract only monthly variables
files <- files[sapply(months, function(x){str_which(files, x)})]
# Load precipitation data
precip <- stack(files, varname = "precip_mm_srf")
# Mask data
precip <- mask(x = precip, mask = msk, maskvalue = 0)
# Assign names
names(precip) <- months
# Convert from kg/m2/s to mm/day
precip <- precip * 86400
# Load temperature data
temp <- stack(files, varname = "temp_mm_srf")
# Mask data
temp <- mask(x = temp, mask = msk, maskvalue = 0)
# Assign names
names(temp) <- months
# Convert kelvin to celsius
temp <- temp - 273.15
# Calculate max and min
max_precip <- calc(x = precip, fun = max)
min_precip <- calc(x = precip, fun = min)
max_temp <- calc(x = temp, fun = max)
min_temp <- calc(x = temp, fun = min)
# Stack data
stk <- stack(max_precip, min_precip, max_temp, min_temp)
# Add names
names(stk) <- c("max_precip", "min_precip", "max_temp", "min_temp")
# Rotate data
stk <- raster::rotate(stk)
# Resample data (the extent and resolution must be updated to avoid rgdal issue)
stk <- resample(x = stk, y = r)
# Define original GCRS
crs(stk) <- gcrs
# Plot data
plot(stk)
# Make folders
dir.create(paste0("Data/Covariate_Data/Lewis_Climate/Results/", sep =""), showWarnings = FALSE)
dir.create(paste0("Data/Covariate_Data/Lewis_Climate/Results/Campanian/", sep =""), showWarnings = FALSE)
dir.create(paste0("Data/Covariate_Data/Lewis_Climate/Results/Campanian/teyep/", sep =""), showWarnings = FALSE)
# Save data
writeRaster(
  x = stk,
  filename = paste0("Data/Covariate_Data/Lewis_Climate/Results/Campanian/teyep/", 
                    names(stk), ".grd"), bylayer = TRUE, overwrite = TRUE)
#===== TEYEQ====
# Load mask 
msk <- raster("Data/Covariate_Data/Lewis_Climate/Campanian/teyeq/teyeq.qrparm.mask.nc", 
              varname = "lsm")
# Get file names
files <- list.files("./data/raw-data/climate/Campanian/teyeq/",
                    full.names = TRUE)
# Extract only monthly variables
files <- files[sapply(months, function(x){str_which(files, x)})]
# Load precipitation data
precip <- stack(files, varname = "precip_mm_srf")
# Mask data
precip <- mask(x = precip, mask = msk, maskvalue = 0)
# Assign names
names(precip) <- months
# Convert from kg/m2/s to mm/day
precip <- precip * 86400
# Load temperature data
temp <- stack(files, varname = "temp_mm_srf")
# Mask data
temp <- mask(x = temp, mask = msk, maskvalue = 0)
# Assign names
names(temp) <- months
# Convert kelvin to celsius
temp <- temp - 273.15
# Calculate max and min
max_precip <- calc(x = precip, fun = max)
min_precip <- calc(x = precip, fun = min)
max_temp <- calc(x = temp, fun = max)
min_temp <- calc(x = temp, fun = min)
# Stack data
stk <- stack(max_precip, min_precip, max_temp, min_temp)
# Add names
names(stk) <- c("max_precip", "min_precip", "max_temp", "min_temp")
# Rotate data
stk <- raster::rotate(stk)
# Resample data (the extent and resolution must be updated to avoid rgdal issue)
stk <- resample(x = stk, y = r)
# Define original GCRS
crs(stk) <- gcrs
# Plot data
plot(stk)
# Make folders
dir.create(paste0("Data/Covariate_Data/Lewis_Climate/Results/Campanian/teyeq/", sep =""), showWarnings = FALSE)
# Save data
writeRaster(
  x = stk,
  filename = paste0("Data/Covariate_Data/Lewis_Climate/Results/Campanian/teyeq/", 
                    names(stk), ".grd"), bylayer = TRUE, overwrite = TRUE)
#-------------------------------------------------------------------------------