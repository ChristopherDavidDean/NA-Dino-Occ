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
library(rgeos)
library(RColorBrewer)
library(chronosphere)

#==== if not: ====
# Load in Functions
source("0.Functions_DR.R") # Import functions from other R file (must be in same working directory)

#==== Load master spreadsheets and organise ====
master.occs <- read.csv("Data/Occurrences/occs_all_stages_tetrapoda_28_06_2022.csv", skip = 20)

formations <- read.csv("Data/Occurrences/Formations.csv") # Load formations
formations <- formations[,1:11] # Remove additional columns
formations <- formations[order(formations$Formation),] # Reorganise formations
formations$forbinning <- 1:nrow(formations) # Provide ID for formations
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

#=============================================== TIME BINNING =======================================================

# Standard Bin setup - trim to fit relevant time frame. 
data(stages)
stages <- stages[80:81,] # Set stages to range from Campanian to Maastrichtian

# Attach new age constraints from formations data and setup for binning
master.occs$old_max_ma <- master.occs$max_ma # Move old values to legacy
master.occs$old_min_ma <- master.occs$min_ma
master.occs <- merge(master.occs, formations, by = "formation", all = T) # Merge occurrences with formation ages
master.occs <- master.occs %>%
  filter(is.na(occurrence_no) == F) # Remove non-occurrences
master.occs$max_ma <- master.occs$max_age # Add new ages to ma_ma and min_ma
master.occs$min_ma <- master.occs$min_age
master.occs <- master.occs[,1:140] # Remove extra columns from final dataset.
master.occs <- master.occs %>% 
  filter(is.na(max_ma) == F)# %>% # Remove occurrences without a max age (occurrences which can't be matched to a formation)
master.occs$max_ma[master.occs$max_ma > 83.59] <- 83.59
master.occs$min_ma[master.occs$min_m < 66] <- 66 # If using majority method, all occurrences with minimum age of <66 Ma must be capped to 66 Ma, otherwise bin_time() breaks.

#==== STAGE ====
bins <- stages
colnames(bins) <- c("sys","system", "series", "stage", "short", "max_ma", "mid_ma", "min_ma", "dur", "stg", "systemCol", "seriesCol", "col")
bins <- arrange(bins, (max_ma))
bins$bin <- c("S.1","S.2")
master.occs.binned <- bin_time(master.occs, bins, method = "majority")
bins$code <- c("teyeo", "teyeq")
bin.type <- "stage"
         
#==== SUB STAGE ====
bins <- read.csv("Data/Occurrences/substages.csv")
bins$bin <- c("SB.1","SB.2","SB.3","SB.4","SB.5")
master.occs.binned <- bin_time(master.occs, bins, method = "majority")
bin.type <- "substage"

#==== SCOTESE BINS ====
bins <- read.csv("Data/Occurrences/scotesebins.csv")
bins$bin <- c("SC.1","SC.2","SC.3","SC.4")
master.occs.binned <- bin_time(master.occs, bins, method = "majority")
bin.type <- "scotese"

#==== FORMATION BINS ====
bin_limits <- c(2.5, max(formations$max_age), 66) # Set user defined bin size - change the first number to vary resolution in graphs.
Scoring_Grid_1(formations = formations, res = 0.01) # run score grid to work out appropriate bin points
#Scoring_Grid_2(formations = formations, res = 0.01)
newBins(score_grid = score_grid, formations = formations, bin_limits = bin_limits, allbins = allbins, stages = stages, smallamalg = TRUE) # Run new bins to generate bins
bins <- binlist %>% # Remove non-useful bins outside of range, add range for bins
  filter(bottom < 84) %>%
  mutate(range = top-bottom)
colnames(bins) <- c("bin", "min_ma", "max_ma", "mid_ma", "range") # Rename bin names to match across
bins[1,2] <- 66 # Cap youngest bin at 66 Ma. 
master.occs.binned <- bin_time(master.occs, bins, method = 'majority') # Bin occurrences
bin.type <- "formation"

 #=============================================== VISUALISATION AND TESTING ===============================================

#==== Visualise grid cells for collections and occurrences ====
bin.colls <- master.occs %>% # Make unique collections for visualisation
  dplyr::select(collection_no, lat, lng) %>%
  dplyr::distinct()

bin.dinos <- master.occs %>%
  filter(class == "Ornithischia" | class == "Saurischia")

get_grid_im(master.occs, res, "Occurrences", ext = e)
get_grid_im(bin.dinos, res, "Dinosaur Occurrences", ext = e)
get_grid_im(bin.colls, res, "Collections", ext = e)

#==== Testing Targets ====
target_maker(bin.occs, "family", "Ceratopsidae")
Ceratops <- bin.occs.targeted %>%
   filter(Target == "Ceratopsidae")
get_grid_im(Ceratops, 1, "Ceratopsidae Occurrences")

#===== Testing reduction of naive occupancy with subsampling =====
stri <- c()
for (t in 1:50){
  all_results_for_unmarked(data = bin.occs.targeted, name = bin.name, res = res, ext = e, target = target, subsamp = TRUE, sampval = 10, single = FALSE)
  stri <- c(stri, comp_num)
} # Mean score of 2.68 occupied sites dropped.

#===== Testing matchup of cells =====
# Independently find all Hadrosaur occurrences and visualise.
Hadro <- Final %>%
  filter(family == "Hadrosauridae")# %>%
#  filter(Coll_count == "Non-singleton")
get_grid_im(Hadro, res, "Hadrosauridae Occurrences", ext = e)

# Use unmarked dataset to find Ceratopsian occupied cells.
Hadro_occupied <- SS_unmarked_Hadrosauridae %>% 
  filter_all(any_vars(. %in% c(1)))
cells <- as.numeric(rownames(Hadro_occupied))
values <- base::rep(1,length(cells))

# Generate raster to compare against original dataset.
gen_raster(cells, values, res, ext = e)

#=============================================== OCCUPANCY SETUP ===============================================

# Set Target taxa
target = c("Ceratopsidae", "Tyrannosauridae", "Hadrosauridae") # set Targets
target_maker(master.occs.binned, "family", target) # run target_maker

# Palaeorotate binned and targetted data using chronosphere - WARNING: THIS STEP TAKES A LONG TIME!!!
dems <- fetch(dat="paleomap", var="dem") # Fetch Scotese DEMs
demord <- matchtime(dems, master.occs.binned$bin_midpoint) # Match up Scotese DEMs with bin midpoints.
master.occs.binned.targeted$mapage <- names(demord) # Provide DEM names to allow for binned rotations.
newCoords <- reconstruct(master.occs.binned.targeted[, c("lng", "lat")], 
                         age=master.occs.binned.targeted[, "mapage"], enumerate=FALSE, 
                         verbose=FALSE) # Palaeorotate occurrences. THIS TAKES AGES.
colnames(newCoords) <- c("plng", "plat") # Rename co-ordinates
master.occs.binned.targeted <- cbind(master.occs.binned.targeted,newCoords) # Attach p-coords back to dataset.

#===== Get all data necessary for running occupancy models ====
# For each time bin - 
for(t in 1:length(bins$bin)){ 
  
  #===== DATA READY FOR UNMARKED =====
  # Select relevant occurrences for bin. In this instance, Campanian.
  bin.name <- bins$bin[t]
  bin.occs <- master.occs.binned.targeted %>% 
    filter(bin_assignment == bin.name)
  
  # Create appropriate folder
  dir.create(paste0("Results/", bin.type, "/", bin.name, "/", sep =""), showWarnings = FALSE)
  dir.create(paste0("Results/", bin.type, "/", bin.name, "/", res, "/", sep =""), showWarnings = FALSE)
  
  # Record resolution target data
  vect <- c(0.1, 0.5, 1)
  res_data(bin.occs, target, vect = vect, single = TRUE) # Makes list of dataframes, each containing information about various cell resolutions. Can also see Naive occupancy estimates.
  Res_results <- do.call(rbind.data.frame, Res_results_list)
  write.csv(Res_results, file.path(paste("Results/", bin.type, "/", bin.name, "/Targeted_res_stats.csv", sep="")))
  
  # Prepare data for unmarked
  all_results_for_unmarked(data = bin.occs, name = bin.name, res = res, ext = e, target = target, subsamp = TRUE, sampval = 5, single = FALSE)
  write.csv(Final, file.path(paste("Results/", bin.type, "/", bin.name, "/", res, "/", bin.name, "_dataset.csv", sep = "")))
  
  #===== SITE DETECTION COVARIATE DATA =====
  # Grab hi resolution covariates (taken directly from collection co-ordinates)
  alt_cov_grab(Final, res = res, out = T)
  
  #===== SITE OCCUPANCY COVARIATE DATA =====
  if(bin.type == "stage" | bin.type == "scotese"){
    # Load rasters
    wc <- list.files(paste("Data/Covariate_Data/Lewis_Climate/Results/", bins$stage[t], "/", bins$code[t], "/", sep = ""), pattern="grd")
    stacked <- raster::stack(paste("Data/Covariate_Data/Lewis_Climate/Results/", bins$stage[t], "/", bins$code[t], "/", wc, sep =""))
    # Extract info
    demord[t]
    coll42$elev <- extract(dem42, coll42[, c("plng", "plat")])
    
    # NOTE - ADD PALAEOLAT COVARIATE BACK IN FROM DATA HERE
    # Save info
  }
}

#=============================================== OLD PALAEO-DATA ====================================================

#===== GETECH DATA =====
# Data
CampPrecip <- raster("Data/Covariate_Data/Climate_Data/Camp_Precip.asc")
projection(CampPrecip) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
CampTemp <- raster("Data/Covariate_Data/Climate_Data/Camp_Temp.asc")
projection(CampTemp) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Camp Stack Extract
CampStack <- stack(CampPrecip, CampTemp)

NAcells <- Final %>% # Clean for NA values in paleolat/lng - Record cells with NAs
  filter(is.na(paleolat)) %>%
  distinct(cells) # Remove cells with NAs
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
#Camp_Cov_Master <- Camp_Covs
#Camp_Cov_Master <- Camp_Cov_Master %>%
#  select(collection_no, collection_name, accepted_name, identified_rank, family, genus, formation, member, 
#         cells.1, Camp_out_1, DEM_1, LANDCVI_selected_1, MGVF_1, 
#         SLOPE_1, WC_Prec_1, WC_Temp_1, cells.2, CampPrecip_1, CampTemp_1, colls_per_cell)

#write.csv(Camp_Cov_Master, file.path(paste("Results/", res, "/CampanianMasterSheet.csv", sep="")))

# Clean data to just gridcells and associated covariates
#Camp_Covs <- Camp_Covs %>%
#  select(cells.1, Camp_out_1, DEM_1, LANDCVI_selected_1, MGVF_1, 
#  SLOPE_1, WC_Prec_1, WC_Temp_1, colls_per_cell) 
#Camp_Covs <- Camp_Covs %>%
#  distinct()
#Camp_Covs <- aggregate(x = Camp_Covs, by = list(Camp_Covs$cells.1), FUN = "mean", na.rm = TRUE)
#Camp_Covs$CampTemp <- Camp_Covs$CampTemp - 273.15 # Convert from Kelvin to degrees Celsius
#write.csv(Camp_Covs, file.path(paste("Results/", res, "/CampanianCovariates.csv", sep="")))

