################################################################################
############ OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS #############
################################################################################

# Christopher D. Dean, Alfio Alessandro Chiarenza, Jeffrey W. Doser, Alexander
# Farnsworth, Lewis A. Jones, Sinéad Lyster, Charlotte L. Outhwaite, Paul J. 
# Valdes, Richard J. Butler, Philip D. Mannion.
# 2025
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
library(parallel)
library(tibble)
library(divDyn)
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
master.occs <- read.csv("Data/Occurrences/occs_tetrapoda_12_06_2024.csv", skip = 22)

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
res <- 0.5
# Set extent
e <- extent(-155, -72, 22.5, 73)
# Set max limit value
max_val <- "None"
max_val_on <- F
nomam <- T

##### Remove occurrences outside bounds #####
master.occs <- master.occs %>%
  filter(lat > 22.5) %>%
  filter(lat < 73) %>%
  filter(lng > -155) %>%
  filter(lng < -72)

###### Remove occurrences outside bounds #####
#master.occs <- master.occs %>%
#  filter(lat > 32) %>%
#  filter(lat < 60) %>%
#  filter(lng > -120) %>%
#  filter(lng < -90)

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
master.occs$min_ma[master.occs$min_ma < 66] <- 66.001 # If using majority method, all occurrences with minimum age of <66 Ma must be capped to 66 Ma, otherwise bin_time() breaks.

##### FILTERING #####
#master.occs <- master.occs %>%
#  dplyr::filter(class == "Ornithischia" | class == "Reptilia")
  
###################
##### SCOTESE #####
###################

bins <- read.csv("Data/Occurrences/scotesebins.csv")
bins$bin <- bins$code
master.occs.binned <- palaeoverse::bin_time(occdf = master.occs, bins = bins, method = "majority")
#master.occs.binned <- bin_time(occdf = master.occs, bins = bins, method = "all")
bin.type <- "scotese"

#####################
##### FORMATION #####
#####################

# Run combined binning function, choosing adjustable window
bin.res <- 2.5
bins <- binning(window = bin.res, occdf = master.occs, formations = formations)
master.occs <- master.occs %>%
  dplyr::filter(max_ma < bins$max_ma[6]) %>%
  dplyr::filter(min_ma > bins$min_ma[1])
bins[6,3] <- bins[7,3]
bins[6,4] <- (bins[6,2] + bins[6,3])/2
bins <- bins[-7, ]

master.occs.binned <- bin_time(occdf = master.occs, bins = bins, method = "majority")
bin.type <- "formation"

##### Checking list length #####
master.occs.binned %>%
  #dplyr::filter(identified_rank == "genus" | identified_rank == "species") %>%
  group_by(collection_no, bin_midpoint) %>%
  summarise(count = n()) %>%
  group_by(bin_midpoint) %>%
  summarise(mean = mean(count))

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
a <- get_grid_im(bin.dino.coll, 1, "Dinosaur Collections", ext = e2)

pdf(paste("Figures/Dino.colls.pdf", sep = ""))
print(a[[1]])
dev.off()

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
                           single = FALSE, formCells = form_cells)
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
# 3. FILTERING TEMPORALLY UNCERTAIN OCCURRENCES
################################################################################

# Setup new formations dataframe
formations2 <- formations
names(formations2)[names(formations2) == "max_age"] <- "max_ma"
names(formations2)[names(formations2) == "min_age"] <- "min_ma"

# Filter to only include formations within temporal bounds
formations2 <- formations2 %>%
  filter(min_ma < 85) %>%
  filter(max_ma > 65)

# Make a grid of all possible time bin/formation combos
if(bin.type == "scotese"){
  comp_df <- expand.grid(A = sort(unique(formations2$formation)), 
                         B = sort(unique(bins$code)))
}
if(bin.type == "formation"){
  comp_df <- expand.grid(A = sort(unique(formations2$formation)), 
                         B = sort(unique(bins$bin)))
}

# Function to calculate percentage overlap between time bins and formations
perc_overlap = function(x.start, x.end, y.start, y.end){
  #  if(x.start == y.start & x.end == y.end){
  #    return(100)
  #  }
  x.len = abs(x.end - x.start)
  # largest start
  max.start = max(c(x.start, y.start))
  min.end = min(c(x.end, y.end))
  overlap = min.end - max.start
  overlap = ifelse(overlap <= 0, 0, overlap)
  perc_overlap = overlap / x.len * 100
  return(perc_overlap)
}

# Function to apply percentage overlap down a list
comp_func <- function(a, b) {
  min_a <- formations2$min_ma[formations2$formation == a]
  max_a <- formations2$max_ma[formations2$formation == a]
  if(bin.type == "scotese"){
    min_b <- bins$min_ma[bins$code == b]
    max_b <- bins$max_ma[bins$code == b]
  }
  if(bin.type == "formation"){
    min_b <- bins$min_ma[bins$bin == b]
    max_b <- bins$max_ma[bins$bin == b]
  }
  return(perc_overlap(min_a, max_a, min_b, max_b))
}

# Apply percentage overlap to dataframe of formations, using time bins
comp_df$overlap <- unlist(Map(comp_func, a = comp_df$A, b = comp_df$B))

# Calculate formations appearance in bins
overlap_in_bins <- comp_df %>%
  arrange(A) %>%
  filter(overlap != 0) %>%
  dplyr::group_by(A) %>%
  dplyr::summarize(sum = sum(overlap), 
                   n_bins = sum(overlap > 0))

# Make a vector of formations appearing in multiple bins
formations_to_remove <- overlap_in_bins %>%
  filter(n_bins > 2)  %>%
  dplyr::pull(A)

# Remove listed formations from master.occs
master.occs.binned <- master.occs.binned %>%
  dplyr::filter(!(formation %in% formations_to_remove))

################################################################################
# 4. TEST: FILTERING MAMMALS AND SQUAMATA
################################################################################

# Optional: remove mammals, squamates and amphibians from analysis to test whether
#           changes in the frequency of small-bodied collections through time is
#           impacting results.
if(nomam == T){
  master.occs.binned <- master.occs.binned %>%
    dplyr::filter(class != "Mammalia") %>%
    dplyr::filter(class != "Amphibia") %>%
    dplyr::filter(order != "Squamata")
}

################################################################################
# 5. OCCUPANCY SETUP AND PREPPED DATA
################################################################################

# Set Target taxa
target = c("Ceratopsidae", "Tyrannosauridae", "Hadrosauridae", "Ankylosauridae") # set Targets
target_maker(data = master.occs.binned, level = "family", target = target) # run target_maker
master.occs.binned.targeted <- master.occs.binned.targeted %>%
  dplyr::mutate(Rot_age = case_when(bin_midpoint == 66.75 ~ 65, 
                                    bin_midpoint == 69.75 ~ 70,
                                    bin_midpoint == 74.95 ~ 75,
                                    bin_midpoint == 80.75 ~ 80))

# Palaeo-rotate to get accurate palaeolat/long
if(bin.type == "scotese"){
  master.occs.binned.targeted <- palaeoverse::palaeorotate(occdf = master.occs.binned.targeted,
                                              lat = "lat",
                                              lng = "lng",
                                              age = "Rot_age",
                                              method = "point", 
                                              model = c("PALEOMAP")
  )
}
if(bin.type == "formation"){
  master.occs.binned.targeted <- palaeorotate(occdf = master.occs.binned.targeted,
                                              lat = "lat",
                                              lng = "lng",
                                              age = "bin_midpoint",
                                              method = "point", 
                                              model = c("PALEOMAP")
  )
}

##### Save prepped occurrence data for other approaches #####
if(nomam == F){
  write.csv(master.occs.binned.targeted, file.path(paste("Prepped_data/Occurrence_Data/", bin.type, 
                                                         "/", bin.type, "_occurrence_dataset.csv", sep = "")))
}else{
  write.csv(master.occs.binned.targeted, file.path(paste("Prepped_data/Occurrence_Data/", bin.type, 
                                                         "/", bin.type, "_no_mammal_occurrence_dataset.csv", sep = "")))
}

# Record resolution target data
vect <- c(0.5, 1)

# Makes list of dataframes, each containing information about various cell resolutions. 
# Can also see Naive occupancy estimates.
for(t in 1:length(bins$bin)){ 
  # Select relevant occurrences for bin. 
  bin.name <- bins$bin[t]
  bin.occs <- master.occs.binned.targeted %>% 
    filter(bin_assignment == bin.name)
  res_data(bin.occs, target, vect = vect, single = TRUE) 
  Res_results <- do.call(rbind.data.frame, Res_results_list)
  if(nomam == F){
    write.csv(Res_results, file.path(paste("Prepped_data/Occurrence_data/", bin.type, "/", 
                                           "/Targeted_res_stats.", bin.name,".csv", sep="")))
  }else{
    write.csv(Res_results, file.path(paste("Prepped_data/Occurrence_data/", bin.type, "/", 
                                           "/Targeted_res_stats.", bin.name,"_no_mammal.csv", sep="")))
  }
}
