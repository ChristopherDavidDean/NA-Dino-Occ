#===============================================================================================================================================
#============================================== OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS ==========================================
#===============================================================================================================================================

# Christopher D. Dean, Lewis A. Jones, Alfio A. Chiarenza, Alex Farnsworth, María V. Jiménez‐Franco, Richard J. Butler.
# 2019
# Script written by Christopher D. Dean and Lewis A. Jones

#============================================== FILE 1: SETTING UP FILES FOR OCCUPANCY MODELLING ===============================================

#============================================== INITIAL SETUP ===============================================

# Broad PBDB download: http://paleobiodb.org/data1.2/occs/list.csv?datainfo&rowcount&base_name=Tetrapoda&interval=Campanian,Maastrichtian&cc=NOA,^MX&envtype=terr&show=full,strat,lith,env,timebins,timecompare,ref

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Set your working directory

#==== If first time using: ====
# Make vector of package names
packages <- c("beepr", "raster", "dplyr", "lattice", "rasterVis", "sp", "maps", "maptools", "parallel") #list your packages here

# Install packages
ipak(packages)

# Load Required packages
library(beepr)
library(raster)
library(tidyr)
library(plyr)
library(dplyr)
library(lattice)
library(rasterVis)
library(sp)
library(maps)
library(maptools)
library(parallel)
library(tibble)
library(divDyn)

#==== if not: ====
# Load in Functions
source("0.Functions_DR.R") # Import functions from other R file (must be in same working directory)

#==== Load master spreadsheets and organise ====
master.occs <- read.csv("Data/Occurrences/occs_all_stages_tetrapoda_28_06_2022.csv", skip = 20)
master.occs <- master.occs %>%
  filter(max_ma < 85) # Remove older occurrences than we are looking at

formations <- read.csv("Data/Occurrences/Formations.csv") # Load formations
formations <- formations[,1:11] # Remove additional columns
formations <- formations[order(formations$Formation),] # Reorganise formations
formations$forbinning <- 1:nrow(formations) # Provide if for formations
formations$Range <- formations$max_age - formations$min_age # Calculate formation range
formations$Diversity <- 0 # Add in dummy variables to ensure code works (sorry!)
formations$Occurrences <- 0 # Add in dummy variables to ensure code works (sorry!)
colnames(formations)[1] <- "formation" # Change to allow for further analysis

#==== Set working resolution and extent ====
res <- 0.1
#get_extent(master.occs)
e <- extent(-155, -72, 22.5, 73)

#==== Manually setup extent for figures ====
#maxLat <- round_any((max(occurrence_dataset$lat) + 7), 0.5)  #get value for defining extent, and increase by x for visualisation purposes
#minLat <- round_any((min(occurrence_dataset$lat) - 7), 0.5)  #get value for defining extent, and increase by x for visualisation purposes
#maxLng <- round_any((max(occurrence_dataset$lng) + 10), 0.5)  #get value for defining extent, and increase by x for visualisation purposes
#minLng <- round_any((min(occurrence_dataset$lng) - 10), 0.5) #get value for defining extent, and increase by x for visualisation purposes
#e <<- extent(minLng, maxLng, minLat, maxLat)

#==== Convert rasters to desired resolution ====
# source("0.Format_data_DR.R") # Import extent and run cleaning/import of raster datasets. Running this will take a while, 
# but will automatically update rasters so that they are of the desired resolution and place them in the appropriate folder
# for running the next steps. 

#=============================================== TIME BINNING =======================================================

# Standard Bin setup - trim to fit relevant time frame. 
data(stages)
stages <- stages[80:81,] # Set stages to range from Campanian to Maastrichtian

# Attach new age constraints from formations data
master.occs$old_max_ma <- master.occs$max_ma # Move old values to legacy
master.occs$old_min_ma <- master.occs$min_ma
test <- merge(master.occs, formations, by = "formation", all = T) # Merge occurrences with formation ages
test <- test %>%
  filter(is.na(occurrence_no) == F) # Remove non-occurrences
test$max_ma <- test$max_age # Add new ages to ma_ma and min_ma
test$min_ma <- test$min_age
master.occs <- test[,1:140] # Remove extra columns from final dataset.
master.occs <- master.occs %>% 
  filter(is.na(max_ma) == F) %>% # Remove occurrences without a max age (occurrences which can't be matched to a formation)
  filter(max_ma <84.6) # Remove occurrences older than necessary.

#==== STAGE ====
bins <- stages
colnames(bins) <- c("sys","system", "series", "stage", "short", "max_ma", "mid_ma", "min_ma", "dur", "stg", "systemCol", "seriesCol", "col")
bins$bin <- c(1,2)
master.occs$max_ma[master.occs$max_ma > 83.59] <- 83.59
master.occs$min_ma[master.occs$min_m < 66] <- 66

master.occs.stage <- bin_time(master.occs, bins, method = "majority")
         
#==== SUB STAGE ====
bins <- read.csv("Data/Occurrences/substages.csv")
master.occs.stage <- bin_time(master.occs, bins, method = "majority")

#==== FORMATION BINS ====
formations <- formations %>%
  dplyr::filter(min_age < 84) # Limit what formations are included
bin_limits <- c(2.5, max(formations$max_age), 66) # Set user defined bin size - change the first number to vary resolution in graphs.
Scoring_Grid_1(formations = formations, res = 0.01)
#Scoring_Grid_2(formations = formations, res = 0.01)
newBins(score_grid = score_grid, formations = formations, bin_limits = bin_limits, allbins = allbins, stages = stages, smallamalg = TRUE)

bins <- binlist %>% # Remove non-useful bins outside of range
  filter(bottom < 84)
colnames(bins) <- c("bin", "min_ma", "max_ma", "mid_ma") # Rename bin names to match across
bins[1,2] <- 66 # Cap youngest bin at 66 Ma. 

# If using majority method, all occurrences with minimum age of <66 Ma must be capped to 66 Ma, otherwise bin_time() breaks.
master.occs$min_ma[master.occs$min_ma<66] <- 66

# Bin occurrences
master.occs.form <- bin_time(master.occs, bins, method = 'majority') 

 #=============================================== DATA SETUP ===============================================

# Select relevant occurrences for bin. In this instance, Campanian.
bin.name <- "Campanian"
#bin.name <- "Maastrichtian"
bin.occs <- master.occs.stage %>% 
  filter(bin_assignment == 1)

#==== Visualise grid cells for collections and occurrences ====
bin.colls <- bin.occs %>% # Make unique collections for visualisation
  dplyr::select(collection_no, lat, lng) %>%
  dplyr::distinct()

bin.dinos <- bin.occs %>%
  filter(class == "Ornithischia" | class == "Saurischia")

get_grid_im(bin.occs, res, "Occurrences", ext = e)
get_grid_im(bin.dinos, res, "Dinosaur Occurrences", ext = e)
get_grid_im(bin.colls, res, "Collections", ext = e)

#==== Testing Targets ====
# target_maker(bin.occs, "family", "Ceratopsidae")
# Ceratops <- bin.occs.targeted %>%
#   filter(Target == "Ceratopsidae")
# get_grid_im(Ceratops, 1, "Ceratopsidae Occurrences")

#=============================================== OCCUPANCY SETUP ===============================================

# Set Target taxa
target = c("Ceratopsidae", "Tyrannosauridae", "Hadrosauridae") # set Targets
target_maker(bin.occs, "family", target) # run target_maker

# Check resolution data
vect <- c(0.1, 0.5, 1)
res_data(bin.occs.targeted, target, vect = vect, single = TRUE) # Makes list of dataframes, each containing information about various cell resolutions. Can also see Naive occupancy estimates.
Res_results_list$Ceratopsidae
Res_results_list$Tyrannosauridae
Res_results_list$Hadrosauridae

# Prepare data for unmarked
all_results_for_unmarked(data = bin.occs.targeted, name = bin.name, res = res, ext = e, target = target, subsamp = FALSE, single = TRUE)
all_results_for_unmarked(data = bin.occs.targeted, name = bin.name, res = res, ext = e, target = target, subsamp = TRUE, sampval = 10, single = FALSE)

# Testing reduction of naive occupancy with subsampling
stri <- c()
for (t in 1:50){
  all_results_for_unmarked(data = bin.occs.targeted, name = bin.name, res = res, ext = e, target = target, subsamp = TRUE, sampval = 10, single = FALSE)
  stri <- c(stri, comp_num)
} # Mean score of 2.68 occupied sites dropped.

#===== Prepare data for multispecies =====
prepare_for_multispecies(bin.occs.targeted, res, e, level = "species", target)
prepare_for_multispecies(bin.occs.targeted, res, e, level = "genus", target)

#=============================================== VISUAL CHECKS ===============================================

# Independently find all Ceratopsian occurrences and visualise.
Cera <- Final %>%
  filter(family == "Ceratopsidae")# %>%
#  filter(Coll_count == "Non-singleton")
get_grid_im(Cera, res, "Ceratopsian Occurrences", ext = e)

# Use unmarked dataset to find Ceratopsian occupied cells.
Cera_occupied <- SS_unmarked_Ceratopsidae %>% 
  filter_all(any_vars(. %in% c(1)))
cells <- as.numeric(rownames(Cera_occupied))
values <- base::rep(1,length(cells))

# Generate raster to compare against original dataset.
gen_raster(cells, values, res, ext = e)

#=============================================== COVARIATE SETUP ===============================================

# Grab hi resolution covariates (taken directly from collection co-ordinates)
alt_cov_grab(Final)

# Test using gen_raster
gen_raster(hires_cov_dat$cells, hires_cov_dat$mean_prec, res, ext = e)

# Add covariates to Occurrence spreadsheet
get_cov_from_stack(Final, res = res)

# Clean/split data to just relevant covariates
Camp_Covs <- cov_dat[, -grep("Maas_out_", colnames(cov_dat))]

#===== PALAEO-DATA =====
# Data
CampPrecip <- raster("Data/Covariate_Data/Climate_Data/Camp_Precip.asc")
projection(CampPrecip) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
CampTemp <- raster("Data/Covariate_Data/Climate_Data/Camp_Temp.asc")
projection(CampTemp) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Camp Stack Extract
CampStack <- stack(CampPrecip, CampTemp)

NAcells <- Final %>% # Clean for NA values in paleolat/lng - Record cells with NAs
  filter(is.na(paleolat)) %>%
  distinct(cells.1) # Remove cells with NAs
Final2 <- Final %>%
  drop_na("paleolat", "paleolng")

xy <- SpatialPointsDataFrame(cbind.data.frame(Final2$paleolng, Final2$paleolat), Final2) # Use palaeo lat/long to extract from palaeo GCMs
Camp_rain <- raster::extract(CampPrecip, xy, sp = TRUE, cellnumbers = TRUE)
Camp_rain <- as.data.frame(Camp_rain)
Camp_rain <- Camp_rain %>%
  dplyr::group_by(cells) %>%
  dplyr::summarize(mean_Cpre = mean(Camp_Precip, na.rm = TRUE))

Camp_temp <- raster::extract(CampTemp, xy, sp = TRUE, cellnumbers = TRUE) # Extract Palaeo Temp and Rainfall
Camp_temp <- as.data.frame(Camp_temp) # Conver to dataframe
Camp_temp <- Camp_temp %>%
  dplyr::group_by(cells) %>%
  dplyr::summarize(mean_Ctemp = mean(Camp_Temp, na.rm = TRUE))
Camp_temp$mean_Ctemp <- Camp_temp$mean_Ctemp - 273.15 # Convert from Kelvin to degrees Celsius
Camp_Occ_Covs <- left_join(Camp_rain, Camp_temp, by = "cells")

# Visualise values on modern lat/long
gen_raster(Camp_Occ_Covs$cells, Camp_Occ_Covs$mean_Ctemp, res, e)
gen_raster(Camp_Occ_Covs$cells, Camp_Occ_Covs$mean_Cpre, res, e)

# Combine datasets
full_camp_covs <- full_join(hires_cov_dat, Camp_Occ_Covs, by = "cells")

# Save covariates
write.csv(full_camp_covs, file.path(paste("Results/", res, "/CampanianHiResCovariates.csv", sep="")))

# Make Master Camp Spreadsheet
Camp_Cov_Master <- Camp_Covs
Camp_Cov_Master <- Camp_Cov_Master %>%
  select(collection_no, collection_name, accepted_name, identified_rank, family, genus, formation, member, 
         cells.1, Camp_out_1, DEM_1, LANDCVI_selected_1, MGVF_1, 
         SLOPE_1, WC_Prec_1, WC_Temp_1, cells.2, CampPrecip_1, CampTemp_1, colls_per_cell)

write.csv(Camp_Cov_Master, file.path(paste("Results/", res, "/CampanianMasterSheet.csv", sep="")))

# Clean data to just gridcells and associated covariates
Camp_Covs <- Camp_Covs %>%
  select(cells.1, Camp_out_1, DEM_1, LANDCVI_selected_1, MGVF_1, 
  SLOPE_1, WC_Prec_1, WC_Temp_1, colls_per_cell) 
Camp_Covs <- Camp_Covs %>%
  distinct()
Camp_Covs <- aggregate(x = Camp_Covs, by = list(Camp_Covs$cells.1), FUN = "mean", na.rm = TRUE)
Camp_Covs$CampTemp <- Camp_Covs$CampTemp - 273.15 # Convert from Kelvin to degrees Celsius
write.csv(Camp_Covs, file.path(paste("Results/", res, "/CampanianCovariates.csv", sep="")))
