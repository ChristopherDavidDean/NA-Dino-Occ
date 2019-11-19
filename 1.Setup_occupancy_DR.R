#===============================================================================================================================================
#============================================== OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS ==========================================
#===============================================================================================================================================

# Christopher D. Dean, Lewis A. Jones, Alfio A. Chiarenza, Alex Farnsworth, María V. Jiménez‐Franco, Richard J. Butler.
# 2019
# Script written by Christopher D. Dean and Lewis A. Jones

#============================================== FILE 1: SETTING UP FILES FOR OCCUPANCY MODELLING ===============================================

#============================================== INITIAL SETUP ===============================================

# Broad PBDB download: https://paleobiodb.org/data1.2/occs/list.csv?base_name=tetrapoda&taxon_reso=genus&pres=regular&interval=Campanian,Maastrichtian&cc=NOA&envtype=terr&show=full,strat,lith,env

# Set working directory
setwd("C:/Users/deancd/Documents/RESEARCH/PROJECTS/DINO_RANGE/NA-Dino-Occ/") # Set your working directory

#==== If first time using: ====
# Make vector of package names
packages <- c("beepr", "raster", "dplyr", "lattice", "rasterVis", "sp", "maps", "maptools", "parallel") #list your packages here

# Install packages
ipak(packages)

# Load Required packages
library(beepr)
library(raster)
library(dplyr)
library(lattice)
library(rasterVis)
library(sp)
library(maps)
library(maptools)
library(parallel)

#==== if not: ====
# Load in Functions
source("0.Functions_DR.R") # Import functions from other R file (must be in same working directory)

#=============================================== DATA SETUP ===============================================

# Read in files
camp.occs <- read.csv("Data/Occurrences/Broad_Colls_Data_Campanian.csv", stringsAsFactors = FALSE) # Load in occurrences

# Set working resolution and extent. Note: these values should be the same ones used for the file 1.Setup_occupancy_DR.
res <- 0.5
get_extent(camp.occs)

#==== Visualise grid cells for collections and occurrences ====

camp.colls <- camp.occs %>% # Make unique collections for visualisation
  dplyr::select(collection_no, lat, lng) %>%
  dplyr::distinct()
maas.colls <- maas.occs %>% # Make unique collections for visualisation
  dplyr::select(collection_no, lat, lng) %>%
  dplyr::distinct()

get_grid_im(camp.occs, res, "Collections")
get_grid_im(maas.occs, res, "Collections")
get_grid_im(camp.colls, res, "Collections")
get_grid_im(maas.colls, res, "Collections")

#==== Testing Targets ====
# target_maker(camp.occs, "family", "Ceratopsidae")
# Ceratops <- camp.occs.targeted %>%
#   filter(Target == "Ceratopsidae")
# get_grid_im(Ceratops, 1, "Ceratopsidae Occurrences")


#=============================================== OCCUPANCY SETUP ===============================================

# Set Target taxa
target = c("Ceratopsidae", "Tyrannosauridae", "Hadrosauridae") # set Targets
target_maker(camp.occs, "family", target) # run target_maker

# Check resolution data
res_data(camp.occs.targeted, target, 0.1, 1, 0.1) # Makes list of dataframes, each containing information about various cell resolutions. Can also see Naive occupancy estimates.
Res_results_list$Ceratopsidae
Res_results_list$Tyrannosauridae
Res_results_list$Hadrosauridae

# Prepare data for unmarked - Campanian
all_results_for_unmarked(data = camp.occs.targeted, res = res, ext = e, target = target, subsamp = TRUE)

#=============================================== COVARIATE SETUP ===============================================

# Convert rasters to desired resolution
source("0.Format_data_DR.R") # Import extent and run cleaning/import of raster datasets. Running this will take a while, 
# but will automatically update rasters so that they are of the desired resolution and place them in the appropriate folder
# for running the next steps. 

# Add covariates to Occurrence spreadsheet
get_cov_from_stack(Final, res = res)

# Clean/split data to just relevant covariates
Camp_Covs <- cov_dat[, -grep("Maas_out_", colnames(cov_dat))]
#Maas_Covs <- cov_dat[, -grep("Camp_out_", colnames(cov_dat))]

#===== PALAEO-DATA =====
# Data
CampPrecip <- raster("Lewis_Occupancy_data/Data/Formatted/PalaeoClimate/CampPrecip.asc")
CampTemp <- raster("Lewis_Occupancy_data/Data/Formatted/PalaeoClimate/CampTemp.asc")
#MaasPrecip <- raster("Lewis_Occupancy_data/Data/Formatted/PalaeoClimate/MaasPrecip.asc")
#MaasTemp <- raster("Lewis_Occupancy_data/Data/Formatted/PalaeoClimate/MaasTemp.asc")

# Camp Stack Extract
CampStack <- stack(CampPrecip, CampTemp)
NAcells <- Camp_Covs %>% # Clean for NA values in paleolat/lng - Record cells with NAs
  filter(is.na(paleolat)) %>%
  distinct(cells.1) # Remove cells with NAs
Camp_Covs <- Camp_Covs %>%
  filter(!is.na(paleolat))
xy <- SpatialPointsDataFrame(cbind.data.frame(Camp_Covs$paleolng, Camp_Covs$paleolat), Camp_Covs) # Use palaeo lat/long to extract from palaeo GCMs
Camp_Covs <- extract(CampStack, xy, sp = TRUE, cellnumbers = TRUE) # Extract Palaeo Temp and Rainfall
Camp_Covs <- as.data.frame(Camp_Covs) # Conver to dataframe

# Make Master Camp Spreadsheet
Camp_Cov_Master <- Camp_Covs
# Clean data to just gridcells and associated covariates
Camp_Covs <- Camp_Covs[,c(129:147, 150:152)]
Camp_Covs <- Camp_Covs %>%
  distinct()
Camp_Covs <- aggregate(x = Camp_Covs, by = list(Camp_Covs$cells.1), FUN = "mean")
Camp_Covs <- Camp_Covs[, 2:23]
Camp_Covs[,22] <- Camp_Covs[,22] - 273.15 # Convert from Kelvin to degrees Celsius

# Maas Stack Extract
#MaasStack <- stack(MaasPrecip, MaasTemp)
#NAcells <- cov_dat %>% # Clean for NA values in paleolat/lng - Record cells with NAs
#  filter(is.na(paleolat)) %>%
#  distinct(cells.1)
#cov_dat_complete<- cov_dat %>%
#  filter(!is.na(paleolat))
#xy <- SpatialPointsDataFrame(cbind.data.frame(cov_dat_complete$paleolng, cov_dat_complete$paleolat), cov_dat_complete)
#Maas_Dat <- extract(MaasStack, xy, sp = TRUE, cellnumbers = TRUE)
#Maas_Dat <- as.data.frame(Maas_Dat)