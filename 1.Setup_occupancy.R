################################################################################
############ OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS #############
################################################################################

# Christopher D. Dean, Lewis A. Jones, Alfio A. Chiarenza, Sinéad Lyster, 
# Alex Farnsworth, Philip D. Mannion, Richard J. Butler.
# 2023
# Script written by Christopher D. Dean

################################################################################
#               FILE 1: SETTING UP FILES FOR OCCUPANCY MODELLING               #
################################################################################

################################################################################
# 1. INITIAL SETUP
################################################################################

##################################
##### PACKAGES AND FUNCTIONS #####
##################################

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

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
library(palaeoverse)

##### Load in Functions #####
source("0.Functions.R") # Import functions from other R file (must be in same working directory)

###############################
##### DATASETS AND VALUES #####
###############################

# Broad PBDB download
# http://paleobiodb.org/data1.2/occs/list.csv?datainfo&rowcount&base_name=Tetrapoda&interval=Campanian,Maastrichtian&cc=NOA,^MX&envtype=terr&show=full,strat,lith,env,timebins,timecompare,ref

# Load main dataset
master.occs <- read.csv("Data/Occurrences/occs_all_stages_tetrapoda_23_05_2023.csv", skip = 20)
# Load formations
formations <- read.csv("Data/Occurrences/Formations.csv") # Load formations
# Organise formations
formations <- formations[,1:11] # Remove additional columns
formations <- formations[order(formations$Formation),] # Reorganise formations
formations$forbinning <- 1:nrow(formations) # Provide ID for formations
formations$Range <- formations$max_age - formations$min_age # Calculate formation range
formations$Diversity <- 0 # Add in dummy variables to ensure code works (sorry!)
formations$Occurrences <- 0 # Add in dummy variables to ensure code works (sorry!)
colnames(formations)[1] <- "formation" # Change to allow for further analysis

##### Set values #####
# Set resolution
res <- 0.1
# Set extent
e <- extent(-155, -72, 22.5, 73)
# Set max limit value
max_val <- 40
max_val_on <- T

##### Remove occurrences outside bounds #####
master.occs <- master.occs %>%
  filter(lat > 32) %>%
  filter(lat < 60) %>%
  filter(lng > -120) %>%
  filter(lng < -90)

##### Manually setup extent for figures #####
#maxLat <- round_any((max(occurrence_dataset$lat) + 7), 0.5)  #get value for defining extent, and increase by x for visualisation purposes
#minLat <- round_any((min(occurrence_dataset$lat) - 7), 0.5)  #get value for defining extent, and increase by x for visualisation purposes
#maxLng <- round_any((max(occurrence_dataset$lng) + 10), 0.5)  #get value for defining extent, and increase by x for visualisation purposes
#minLng <- round_any((min(occurrence_dataset$lng) - 10), 0.5) #get value for defining extent, and increase by x for visualisation purposes
#e <<- extent(minLng, maxLng, minLat, maxLat)

################################################################################
# 2. TIME BINNING
################################################################################

#################
##### SETUP #####
#################

# Get information on ICS based intervals & trim to fit relevant time frame. 
data(stages)
# Set stages to range from Campanian to Maastrichtian
stages <- stages[80:81,] 
# Set bin to stop messing up later binning
stages[2,8] <- 66 

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
master.occs$min_ma[master.occs$min_m < 66] <- 66.001 # If using majority method, all occurrences with minimum age of <66 Ma must be capped to 66 Ma, otherwise bin_time() breaks.

#################
##### STAGE #####
#################

bins <- stages
colnames(bins) <- c("sys","system", "series", "stage", "short", "max_ma", 
                    "mid_ma", "min_ma", "dur", "stg", "systemCol", "seriesCol", 
                    "col")
bins <- arrange(bins, (max_ma))
bins$bin <- c("teyeo", "teyeq")
master.occs.binned <- bin_time(master.occs, bins, method = "all")
bin.type <- "stage"
         
####################
##### SUBSTAGE #####
####################

bins <- read.csv("Data/Occurrences/substages.csv")
bins$bin <- c("SB.1","SB.2","SB.3","SB.4","SB.5")
master.occs.binned <- bin_time(master.occs, bins, method = "all")
bin.type <- "substage"

###################
##### SCOTESE #####
###################

bins <- read.csv("Data/Occurrences/scotesebins.csv")
bins$bin <- bins$code
master.occs.binned <- bin_time(master.occs, bins, method = "majority")
bin.type <- "scotese"

#####################
##### FORMATION #####
#####################

# Run combined binning function, choosing adjustable window
bin.res <- 2.5
bins <- binning(bin.res, master.occs)
master.occs <- master.occs %>%
  dplyr::filter(max_ma < bins$max_ma[6]) %>%
  dplyr::filter(min_ma > bins$min_ma[1])
master.occs.binned <- bin_time(master.occs, bins, method = "majority")
bin.type <- "formation"

################################################################################
# 3. VISUALISATION AND TESTING
################################################################################

##### Visualise grid cells for collections and occurrences #####
bin.colls <- master.occs %>% # Make unique collections for visualisation
  dplyr::select(collection_no, lat, lng, formation, max_ma, min_ma) %>%
  dplyr::distinct()

bin.dino.coll <- master.occs %>%
  filter(class == "Ornithischia" | class == "Saurischia") %>%
  dplyr::select(collection_no, lat, lng) %>%
  distinct()

get_grid_im(master.occs, res, "Occurrences", ext = e)
get_grid_im(bin.dino.coll, res, "Dinosaur Collections", ext = e)

##### Testing Targets #####
target_maker(bin.occs, "family", "Ceratopsidae")
Ceratops <- bin.occs.targeted %>%
   filter(Target == "Ceratopsidae")
get_grid_im(Ceratops, 1, "Ceratopsidae Occurrences")

##### Testing reduction of naive occupancy with subsampling #####
stri <- c()
for (t in 1:50){
  all_results_for_unmarked(data = bin.occs.targeted, name = bin.name, res = res, 
                           ext = e, target = target, subsamp = TRUE, sampval = 10, 
                           single = FALSE, formCells = "Y")
  stri <- c(stri, comp_num)
} # Mean score of 2.68 occupied sites dropped.

##### Testing matchup of cells #####

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


################################################################################
# 3. OCCUPANCY SETUP AND PREPPED DATA
################################################################################

# Set Target taxa
target = c("Ceratopsidae", "Tyrannosauridae", "Hadrosauridae") # set Targets
target_maker(master.occs.binned, "family", target) # run target_maker
master.occs.binned.targeted <- master.occs.binned.targeted %>%
  dplyr::mutate(Rot_age = case_when(bin_midpoint == 66.75 ~ 65, 
                                    bin_midpoint == 69.75 ~ 70,
                                    bin_midpoint == 74.95 ~ 75,
                                    bin_midpoint == 80.75 ~ 80))

# Palaeo-rotate to get accurate palaeolat/long
master.occs.binned.targeted <- palaeorotate(occdf = master.occs.binned.targeted,
                                            lat = "lat",
                                            lng = "lng",
                                            age = "Rot_age",
                                            method = "point", 
                                            model = c("PALEOMAP")
                                            )

##### Save prepped occurrence data for other approaches #####
write.csv(master.occs.binned.targeted, file.path(paste("Prepped_data/Occurrence_Data/", bin.type, 
                                    "/", res, "_", bin.type, "_occurrence_dataset.csv", sep = "")))

############################################################
##### DATA PREPPED FOR OCCUPANCY IN SPARTA OR UNMARKED #####
############################################################

# Load dataset
master.occs.binned.targeted <- read.csv(paste("Prepped_data/Occurrence_Data/", bin.type, 
                                              "/", res, "_", bin.type, "_occurrence_dataset.csv", sep = ""))

# For each time bin - 
for(t in 1:length(bins$bin)){ 
  
  # Select relevant occurrences for bin. 
  bin.name <- bins$bin[t]
  bin.occs <- master.occs.binned.targeted %>% 
    filter(bin_assignment == bin.name)
  
  # Create appropriate folder
  dir.create(paste0("Prepped_data/Occurrence_Data/", bin.type, "/"), showWarnings = FALSE)
  dir.create(paste0("Prepped_data/Occurrence_Data/", bin.type, "/", bin.name, "/", sep =""), showWarnings = FALSE)
  dir.create(paste0("Prepped_data/Occurrence_Data/", bin.type, "/", bin.name, "/", res, "/", sep =""), 
             showWarnings = FALSE)
  
  # Record resolution target data
  vect <- c(0.1, 0.5, 1)
  res_data(bin.occs, target, vect = vect, single = TRUE) # Makes list of dataframes, each containing information about various cell resolutions. Can also see Naive occupancy estimates.
  Res_results <- do.call(rbind.data.frame, Res_results_list)
  write.csv(Res_results, file.path(paste("Prepped_data/Occurrence_data/", bin.type, "/", bin.name, 
                                         "/Targeted_res_stats.csv", sep="")))
  
  # Prepare data for unmarked
  all_results_for_unmarked(data = bin.occs, name = bin.name, res = res, ext = e, 
                           target = target, single = TRUE, formCells = "N", 
                           max_val_on = max_val_on, max_val = max_val) 
  
  ##### PRECISE AND SURVEY COVARIATE DATA #####
  # Grab hi resolution covariates (taken directly from collection co-ordinates)
  precise_cov(bin.occs, unmarked_Ceratopsidae[[2]], max_val)
  
  ##### SITE COVARIATE DATA: MODERN #####

  wc <- list.files(paste("Prepped_data/Covariate_Data/All_data/", 
                         res, "deg/", sep = ""), 
                   pattern=".asc")
  stacked <- raster::stack(paste("Prepped_data/Covariate_Data/All_data/", 
                                 res, "deg/", wc, 
                                 sep =""))
  bin.occs <- get_cov(bin.occs, stacked)
  
  covs <- unlist(strsplit(wc, ".asc"))
  
  site.covs <- bin.occs %>%
    dplyr::select(siteID, any_of(covs), Coll_count) %>%
    filter(Coll_count == "Non-singleton") %>%
    distinct()

  ##### SITE COVARIATE DATA: PALAEO #####
  if(bin.type == "stage" | bin.type == "scotese"){
    # Load rasters
    wc <- list.files(paste("Prepped_data/Covariate_Data/All_data/", 
                           res, "deg/Palaeo/", sep = ""), 
                     pattern=paste0("^", bins$code[t], ".*", sep = ""))
    stacked <- raster::stack(paste("Prepped_data/Covariate_Data/All_data/", 
                                   res, "deg/Palaeo/", wc, 
                                   sep =""))
    # Set palaeo lat/long
    names(bin.occs)[names(bin.occs) == "p_lng"] <- "plng"
    names(bin.occs)[names(bin.occs) == "p_lat"] <- "plat"

    bin.occs <- get_p_cov(bin.occs, stacked)
  
    pcovs <- c("dry_mean.1", "col_mean.1", "wet_mean.1", "hot_mean.1", "ann_sd.1", paste(bin.name, "_sed_", res, sep = ""))
    
    p.site.covs <- bin.occs %>%
      dplyr::select(siteID, any_of(pcovs), Coll_count) %>%
      dplyr::filter(Coll_count == "Non-singleton") %>%
      dplyr::rename_with(~ sub(paste("^", bin.name, sep = ""), "p", .x), starts_with(bin.name)) %>%
      dplyr::rename_with(~gsub("\\d+", "", .)) %>%
      dplyr:: rename_with(~gsub("\\.", "", .)) %>%
      dplyr::group_by(siteID) %>%
      dplyr:: summarise(dry_mean = mean(dry_mean, na.rm = TRUE), 
                        wet_mean = mean(wet_mean, na.rm = TRUE),
                        col_mean = mean(col_mean, na.rm = TRUE),
                        hot_mean = mean(hot_mean, na.rm = TRUE), 
                        ann_sd = mean(ann_sd, na.rm = TRUE),
                        mean_psed = mean(p_sed_, na.rm = TRUE)
      )
    site.covs <- merge(site.covs, p.site.covs, by = "siteID")
    site.covs <- site.covs %>% select(-Coll_count)
  }
  
  ##### SAVING FILES #####
  write.csv(site.covs, file.path(paste("Prepped_data/Occurrence_Data/", bin.type, 
                                       "/", bin.name, "/", res, "/", 
                                       "site_occupancy_covs_", max_val, ".csv", 
                                       sep = "")))
  write.csv(bin.occs, file.path(paste("Prepped_data/Occurrence_Data/", bin.type, 
                                      "/", bin.name, "/", res,
                                      "/", bin.name, "_dataset.csv", sep = "")))
}
