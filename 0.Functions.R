################################################################################
############ OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS #############
################################################################################
#
# Christopher D. Dean, Alfio Alessandro Chiarenza, Jeffrey W. Doser, Alexander
# Farnsworth, Lewis A. Jones, Sinéad Lyster, Charlotte L. Outhwaite, Richard J. 
# Butler, Philip D. Mannion.
# 2024
# Script written by Christopher D. Dean and Lewis A. Jones
#
################################################################################
#                 FILE 0: FUNCTIONS FOR OCCUPANCY MODELLING                    #
################################################################################
#                                                                              #
#  PRIMARY AUTHOR: CHRISTOPHER D. DEAN                                         #
#  CO-AUTHOR: LEWIS A. JONES                                                   #
#                                                                              #
#  Selection of functions that work with PBDB data in order to run a variety   #
#  of occupancy models. Functions range from those that visualise occurrences  #
#  in terms of grid cells to those that reorder presence/absence data per      #
#  grid cell so it fits the format of the package unmarked. Information        #
#  regarding each function can be found in the separate sections below.        #
#                                                                              #
################################################################################

################################################################################
# 0. INDEX
################################################################################

# FUNCTION                                                                  LINE
# 1. REQUIRED PACKAGES........................................................75
# 2. GET_EXTENT..............................................................119
# 3. TARGET_MAKER............................................................135
# 4. SITECOORDSFUN...........................................................168
# 5. VIEW CELLS, GEN_RASTER, GET_GRID_IM.....................................185
#--- 5.1. VIEW_CELLS.........................................................191
#--- 5.2. GEN_RASTER.........................................................232
#--- 5.3. GET_GRID_IM........................................................260
# 6. GET_GRID................................................................314
# 7. FORMATION BINNING.......................................................354
#--- 7.1. SCORING_GRID.......................................................364
#--- 7.2. NEWBINS............................................................434
#--- 7.3. BINNING............................................................564
# 8. GET COVARIATE FUNCTIONS.................................................617
#--- 8.1. FIND_COLLECTIONS...................................................623
#--- 8.2. GET_COV............................................................645
#--- 8.3. GET_P_COV..........................................................666
# 9. PREPARE_FOR_RES_DATA & RES_DATA.........................................686
#--- 9.1. PREPARE_FOR_RES_DATA...............................................692
#--- 9.2. RES_DATA...........................................................769
# 10. SPARTA FUNCTIONS.......................................................807
#--- 10.1. GET_OCC...........................................................815
#--- 10.2. GET_DET...........................................................849
#--- 10.3. NAIVE_RES.........................................................896
#--- 10.4. COMB_RES..........................................................944
#--- 10.5. RUN_MODEL........................................................1000
#--- 10.6. PLOT_OCC.........................................................1061
#--- 10.7. PLOT_NAIVE.......................................................1097
# 11. SPOCCUPANCY FUNCTIONS.................................................1134
#--- 11.1. P_ROTATE.........................................................1141
#--- 11.2. FORMAT_OCC_COVS..................................................1177
#--- 11.3. EXTRACT_P........................................................1246
#--- 11.4. SAMPLE_FOR_SPOCC.................................................1283
#--- 11.4. PREPARE_FOR_SPOCC................................................1317
#--- 11.5. ORGANISE_DET.....................................................1412
#--- 11.6. SITE_REMOVE......................................................1523
#--- 11.7. TRANSPOSE_EH.....................................................1570
#--- 11.8. ARRAY_PREP.......................................................1589
#--- 11.9. DISTANCE_FUN.....................................................1617
#--- 11.10. MAKE_TABLE......................................................1642
# 12. SAVE_LATTICE..........................................................1697
# 13. FORESTPLOT2...........................................................1716
# 14. CLEAN_FOR_FIG.........................................................1937

################################################################################
# 1. REQUIRED PACKAGES
################################################################################

library(beepr)
library(rphylopic)
library(tidyr)
library(raster)
library(plyr)
library(dplyr)
library(lattice)
library(rasterVis)
library(sp)
library(maps)
library(maptools)
library(parallel)
library(reshape2)
library(tibble)
library(divDyn)
library(rgeos)
library(RColorBrewer)
library(chronosphere)
library(sparta)
library(gtools)
library(R2jags)
library(palaeoverse)
library(ggplot2)
library(snowfall)
library(ggnewscale)
library(wesanderson)
library(ggthemes)
library(ggpubr)
library(cowplot)
library(magick)
library(rphylopic)
library(RCurl)
library(png)
library(grid)
library(unmarked)
library(MuMIn)
library(dplyr)
library(AICcmodavg)
library(tibble)

################################################################################
# 2. GET_EXTENT
################################################################################

# Setup raster for resolution and extent. Note: these values should be the same 
# ones used for the file 1.Setup_occupancy_DR.

get_extent <- function(data){
  # get values for defining extent, and increase by x for visualisation purposes
  maxLat <- round_any((max(data$lat) + 1), 0.5) 
  minLat <- round_any((min(data$lat) - 1), 0.5)
  maxLng <- round_any((max(data$lng) + 1), 0.5)
  minLng <- round_any((min(data$lng) - 1), 0.5)
  e <<- extent(minLng, maxLng, minLat, maxLat) # build extent object
}

################################################################################
# 3. TARGET_MAKER
################################################################################

# Adds new "Target" column with targeted organisms based on specific requirements
# of the user. Outputs data as named file ending with "targeted".

target_maker <- function (data, level, target){ 
  # Data is entered data. Level is column to search in. target is vector of 
  # chosen organisms.
  
  for (i in 1:length(target)){
    if(i == 1){ # If Target column doesn't exist:
      filtered <- data %>%
        dplyr::mutate(
          Target = dplyr::case_when(
            !!as.name(level) == target[i] ~ target[i])
        )
    } else { # If Target column already exists
      filtered <- filtered %>%
        # stops case_when overwriting existing Target data
        dplyr:: mutate(Target = 
                         dplyr::case_when(!!as.name(level) == target[i] ~ target[i], 
                                          TRUE ~ (as.character(.$Target)) 
          )
        )
    }
  }
  # Name files based on data entered to function
  temp_name <- paste(deparse(substitute(data)),".", "targeted", sep = "") 
  assign(temp_name, filtered, envir = .GlobalEnv)
}

################################################################################
# 4. SITECOORDSFUN
################################################################################

# Uses a raster, extent and list of site IDs to makes a dataframe for latitudinal
# and longitudinal coordinates for each site ID.

siteCoordsFun <- function(res, e, site_IDs){
  #res is the chosen resolution, e is the extent, site_IDs is a vector of grid cell
  # (site) IDs 
  r <- raster(res = res, ext = e)
  siteCoords <- as.data.frame(xyFromCell(r, site_IDs))
  siteCoords$siteID <- site_IDs
  colnames(siteCoords) <- c("lng", "lat", "siteID")
  return(siteCoords)
}

################################################################################
# 5. VIEW CELLS, GEN RASTER & GET_GRID_IM
################################################################################

# Functions to view specific cells of a raster. 

###########################
##### 5.1. VIEW_CELLS #####
###########################

# Function that extracts data from a raster, and makes a new raster with only a 
# chosen vector of cells.

view_cells <- function(chosen_raster, vector_of_cells, res, zero = FALSE){ 
  # Chosen_raster is the raster with cells that you want to view. vector_of_cells 
  # is the vector of cells to view. res is chosen resolution. zero chooses whether 
  # other cells are listed as NAs or 0 values (for levelplot)
  
  # Get relevant info
  total_cells <- ncell(chosen_raster)
  total_layers <- nlayers(chosen_raster)
  raster_extent <- chosen_raster@extent
  
  # Make dataframe of values
  extracted_values <- raster::extract(chosen_raster, vector_of_cells)
  dframe_of_values <- data.frame(vector_of_cells, extracted_values)
  dframe_of_values <- unique(dframe_of_values)
  colnames(dframe_of_values) <- c("Cells", "Vals")
  
  # Make blank dataframe to copy values into
  dframe_of_cells <- data.frame(1:total_cells)
  colnames(dframe_of_cells) <- "Cells"
  
  # Join dataframes
  full_dframe <- left_join(dframe_of_cells, dframe_of_values, by = "Cells")
  
  if (zero == TRUE){
    full_dframe[is.na(full_dframe)] <- 0
  }
  
  # Make and plot raster
  raster_for_values <- raster(res = res, 
                              val = full_dframe$Vals, 
                              ext = raster_extent)
  plot(raster_for_values)
}

###########################
##### 5.2. GEN_RASTER #####
###########################

# Function that makes a raster from scratch, using a vector of cells and associated 
# values. Used for quickly visualising covariate data.

gen_raster <- function(cell_data, value_data, res, ext, zero = FALSE){
  init_raster <- raster(res = res, ext = ext, val = 1)
  total_cells <- ncell(init_raster)
  dframe_of_values <- data.frame(cell_data, value_data) 
  colnames(dframe_of_values) <- c("Cells", "Vals")
 
  # Make blank dataframe to copy values into
  dframe_of_cells <- data.frame(1:total_cells)
  colnames(dframe_of_cells) <- "Cells"
 
  # Join dataframes
  full_dframe <- left_join(dframe_of_cells, dframe_of_values, by = "Cells")
 
  if (zero == TRUE){
     full_dframe[is.na(full_dframe)] <- 0
  }
  raster_for_values <- raster(res = res, val = full_dframe$Vals, ext = ext)
  return(raster_for_values)
}

#
############################
##### 5.3. GET_GRID_IM #####
############################

# Sets raster to dimensions of inputted data ready for visualisation. Is used in 
# vis_grid. 

# Find maps to use as backdrop
countries <- maps::map("world", plot=FALSE, fill = TRUE) 
states <- maps::map("state", plot = FALSE, fill = TRUE)

# Turn maps into spatialpolygons
countries <<- maptools::map2SpatialPolygons(countries, 
                                            IDs = countries$names, 
                                            proj4string = CRS("+proj=longlat")) 
states <<- maptools::map2SpatialPolygons(states, 
                                         IDs = states$names, 
                                         proj4string = CRS("+proj=longlat")) 

get_grid_im <- function(data, res, name, ext){ 
  # Data is first output from combine_data (fossil.colls). Res is chosen 
  # resolution in degrees. name is user inputted string related to data inputted, 
  # for display on graphs. 
  
  xy <- cbind(as.double(data$lng), as.double(data$lat))
  xy <- unique(xy)
  r <- raster::raster(ext = ext, res = res)
  r <- raster::rasterize(xy, r, fun = 'count')
  #r[r > 0] <- 1 # Remove if you want values instead of pure presence/absence.
  
  # find map to use as backdrop
  countries <- maps::map("world", plot=FALSE, fill = TRUE) 
  # Turn map into spatialpolygons
  countries <<- maptools::map2SpatialPolygons(countries, 
                                              IDs = countries$names, 
                                              proj4string = CRS("+proj=longlat")) 
  mapTheme <- rasterVis::rasterTheme(region=brewer.pal(8,"Reds"))
  
  #create levelplot for raster
  (a <- rasterVis::levelplot(r, margin=F, par.settings=mapTheme,  
                             main = paste("Total ", (substitute(name)), 
                                          " per Grid Cell", sep = "")) + 
      # Plots state lines
      latticeExtra::layer(sp.polygons(states, col = "white", fill = NA), under = T)  + 
      # Plots background colour
      latticeExtra::layer(sp.polygons(countries, col = 0, fill = "light grey"), under = T)) 
  (b <- hist(r, breaks = 20,
             main = paste((substitute(name)), " per Grid Cell", sep = ""),
             xlab = "Number of Collections", ylab = "Number of Grid Cells",
             col = "springgreen"))
  r <<- r
  return(list(a, b))
}

################################################################################
# 6. GET_GRID
################################################################################

# Creates a raster of chosen resolution, and attaches associated grid cell IDs 
# to occurrences/collections

get_grid <- function(data, res, e, r = "N", formCells = "N"){ 
  # data is first output from combine_data (fossil.colls). Res is chosen 
  # resolution in degrees
  
  if (class(r) == "character"){
    r <- raster(res = res, val = 1, ext = e) # Value must be added because 
                                             # extract uses values
    r <<- r
  }
  crs <- r@crs
  xy <- SpatialPointsDataFrame(cbind.data.frame(data$lng, data$lat), 
                               data, 
                               proj4string = crs)
  Final <- raster::extract(r, xy, sp = TRUE, cellnumbers = TRUE)
  Final <- as.data.frame(Final)
  singletons <- Final %>%
    dplyr::select(cells, collection_no) %>%
    dplyr::distinct() %>%
    dplyr::group_by(cells) %>%
    dplyr::summarize(Coll_count = n())
  singletons$Coll_count[singletons$Coll_count == 1] <- "Singleton"
  singletons$Coll_count[singletons$Coll_count != "Singleton"] <- "Non-singleton"
  Final <- left_join(Final, singletons, by = "cells")
  if(formCells == "Y"){
    Final$siteID <- paste(Final$cells, Final$formation, sep = "")
  }
  else{
    Final$siteID <- Final$cells
  }
  temp_name <- paste(deparse(substitute(data)), ".grid", sep = "")
  assign(temp_name, Final, envir = .GlobalEnv)
}

################################################################################
# 7. FORMATION BINNING
################################################################################

# Selection of functions for making formation bins. Further details can be found in:
# Dean, C.D., Chiarenza, A.A., Maidment, S.C.R. 2020. Formation binning: a new 
# method for increased temporal resolution in regional studies, applied to the
# Late Cretaceous dinosaur fossil record of North America, Palaeontology, 63(6), 
# 881-901

#############################
##### 7.1. SCORING_GRID #####
#############################

# Takes the dataframe of formations and scores each one at intervals of 0.01 Ma.
# If the formation does not cross the boundary, it gets an automatic maximum 
# binning score of 100. If the formation crosses the boundary, the smallest 
# percentage of the formation that sits either side of that boundary is identified, 
# and the binning score is reduced by that percentage.

Scoring_Grid_1 <- function(formations, res=0.01) { 
  # Requires formation information. Resolution of time lines is set automatically 
  # at 0.01, but can be adjusted.
  
  # Get ages of formations
  max_age <- max(formations$max_age) 
  min_age <- min(formations$min_age) 
  
  # Make 1ma bins in sequence based on max/min ages. 0.0045 added to ensure 
  # formation is never exactly equivalent to a bin.
  allbins <- seq(min_age-1.0045, max_age+1.0045, res) 
  
  # Make a matrix for the scoring 
  score_grid <- matrix(data = NA, nrow = nrow(formations), ncol = length(allbins)) 
  colnames(score_grid) <- allbins 
  rownames(score_grid) <- formations$Formation 
  
  # Set counter and go through each time line
  counter <- 0
  for(i in allbins) { 
    counter <- sum(counter,1) 
    
    # go through each formation 
    for (f in 1:nrow(formations)){ 
      
      # if timeline is between max/min age of a formation (i.e. formation crosses 
      # that line)
      if (i <= formations$max_age[f] && i >= formations$min_age[f]){ 
        # Work out how much of formation is older/younger than timeline
        a <- formations$max_age[f] - i 
        b <- i - formations$min_age[f] 
        
        # Calculate range of formation
        range <- formations$max_age[f] - formations$min_age[f] 
        
        # Work out percentage that sits each side of line, reduce score by that amount.
        if (a > b){
          score_grid[,counter][f] <- (a/range)*100 
        }
        else{ 
          score_grid[,counter][f] <- (b/range)*100 
        }
      }
      
      # Otherwise, just score it 100. 
      else {
        score_grid[,counter][f] = 100 
      }
    }
  }  
  
  # Work out mean score for each time bin and add to grid
  means <- colMeans(score_grid) 
  score_grid <- rbind(score_grid, means) 
  
  # Output score_grid and bins for later use
  score_grid <<- score_grid 
  allbins <<- allbins 
}

########################
##### 7.2. NEWBINS #####
########################

# Looks at the previously generated score_grid and generates appropriate new bins 
# based on those scores. Boundaries are outputted as a list (binlist). If bins 
# are shorter than 0.5 Ma, they are amalgamated into the bins above and below, and 
# a warning is produced.

newBins <- function(score_grid, formations, bin_limits, allbins, stages, 
                    smallamalg = TRUE){ 
  # Takes previously generated score grid, formations, allbins and stages from 
  # DivDyn package. Also require bin_limits, a user made vector of the following: 
  # 1) user chosen time window in which to look to draw bins. Advised to set at 3 Ma. 
  # 2) Hard maximum age of bins
  # 3) Hard minimum age of bins
  
  # Creates broader bins for testing best point to draw a bin
  score_grid<- as.data.frame(score_grid)
  bin_size <- bin_limits[1]
  max_age <- bin_limits[2]
  min_age <- bin_limits[3]
  form_bins <- c(min_age)
  testbin <- seq(min_age, max_age, bin_size) 
  
  # Drawing bins and giving form_bins (vector of bin boundaries)
  for (i in 1:length((testbin)-1)){
    
    # Create a sequence of ages to draw bins within
    seqs <- seq(testbin[i],testbin[i]+bin_size, 1) 
    pasting <- c()
    for (n in 1:bin_size){
      
      # Set up expression to match to score_grid
      pasting <- c(pasting, paste("^",seqs[n],"|", sep = "")) 
    }
    testmatch2 <- paste(pasting, collapse = "")
    testmatch2 <- substr(testmatch2, 1, nchar(testmatch2)-1) 
    
    # Find maximum bin score within this time window, add to vector of bins
    a <- score_grid[grep(testmatch2, names(score_grid))] 
    z <- apply(a,1,which.max)
    form_bins[i+1] <- names(a)[z][nrow(a)]
  }
  form_bins <- as.numeric(unique(form_bins)) 
  
  # If small bin amalgamation is turned on (is on automatically):
  if (smallamalg == TRUE){ 
    
    # Finds all bins which are under 0.5 Ma in length
    range <- (diff(form_bins) < 0.5) 
    range_checker <- c()
    
    # For each bin, if it is under 0.5 Ma in length:
    for (r in 1:length(range)){
      if (range[r] == TRUE){ 
        
        # Find the length of that bin
        difference <- diff(c(form_bins[r], form_bins[r+1])) 
        
        # Throw warning about bin amalgamation
        warning("Original bin ",  r, " removed due to small range: ~", 
                signif(difference, digits = 3), 
                " Ma. The difference in time has been added to the bins above and below.") 
        
        # Add half length of old bin to bin below and above
        form_bins[r] <- form_bins[r]+(difference/2) 
        form_bins[r+1] <- form_bins[r+1]-(difference/2) 
        
        # Record which bin was too small
        range_checker <- c(range_checker, r) 
      }
    }
    
    # If there have been amalgamated bins, remove old amalgamated bins
    if(length(range_checker) > 0){ 
      form_bins <- form_bins[-range_checker] 
    }
  }
  if (smallamalg == FALSE){
    warning("Small bin amalgamation is turned off. Bins may be too short to record occurrences. You are advised to check bins before running further analyses.") 
  }
  form_bins <<- form_bins
  
  # Creating binlist (data.frame of bins and appropriate age info)
  prefix <- "FB."
  suffix <- seq(1:(length(form_bins)-1))
  my_names <- paste(prefix, suffix, sep = "")
  
  # Combine bin data to make dataframe of minimum, maximum and mid point of each new bin
  binlist <- data.frame(bin = my_names, 
                        bottom = as.numeric(form_bins[1:(length(form_bins)-1)]), 
                        top = as.numeric(form_bins[2:(length(form_bins))]))
  binlist$mid <- (binlist$bottom + binlist$top) / 2
  binlist <<- binlist
  
  #Plot new bins using divDyn package
  par(mar = c(4.1, 4.1, 1, 2.1))
  tsplot(stages, boxes=c("short","system"), 
         xlim=c(stages$bottom[1], stages$top[nrow(stages)]),  
         ylim=c(min(colMeans(score_grid), na.rm = TRUE), 100), 
         prop = 0.08, plot.args = list(cex.lab = 2, cex.axis = 2),
         shading=NULL, boxes.col=c("col","systemCol"), labels.args=list(cex=1.8),
         ylab = "Bin Splitting Score") 
  
  # draw new bins as coloured boxes for comparison to traditional bins
  for(n in 1:length(form_bins)){ 
    if(((n %% 2) == 0) == TRUE) next
    else {
      if(n == length(form_bins)){
        if(nrow(binlist) %% 2 == 0){
          next
        }
        else{
          rect(useful_bins[n], 0, useful_bins[n-1], 100, 
               col = rgb(0.89,0.89,0.89,alpha=0.5), border = NA)
        }
      }
      else{
        rect(form_bins[n], 0, form_bins[n+1], 100-0.01, 
             col = rgb(0.89,0.89,0.89,alpha=0.5), border = NA)
      }
    }
  }
  
  # Draw bin splitting score on plot
  lines(allbins, colMeans(score_grid), lwd = 1.75) 
  box(lwd=2)
}

########################
##### 7.3. BINNING #####
########################

# Wrapper for other formation binning functions. Returns dataframe of appropriate bins, 
# and a binned dataset of fossil occurrences assigned to those bins.

binning <- function(window, occdf, formations){
  # Window is the resolution of the formation bins, occdf is the dataframe of 
  # fossil occurrences to be binned.
  
  # Set user defined bin size - change the first number to vary resolution in graphs.
  bin_limits <- c(window, max(formations$max_age), 66) 
  
  # run score grid to work out appropriate bin points
  Scoring_Grid_1(formations = formations, res = 0.01) 
  #Scoring_Grid_2(formations = formations, res = 0.01)
  
  # Run new bins to generate bins
  newBins(score_grid = score_grid, formations = formations, bin_limits = bin_limits, 
          allbins = allbins, stages = stages, smallamalg = TRUE) 
  
  # Remove non-useful bins outside of range, add range for bins
  bins <- binlist %>% 
    filter(bottom < 84) %>%
    mutate(range = top-bottom)
  
  # Rename bin names to match across
  colnames(bins) <- c("bin", "min_ma", "max_ma", "mid_ma", "range") 
  
  # Cap youngest bin at 66 Ma. 
  bins[1,2] <- 66 
  
  # Bin occurrences using palaeoverse
  master.occs.binned <- bin_time(occdf, bins, method = 'majority') 
  bin.type <- "formation"
  
  # Adding midpoint
  master.occs.binned$mid_ma <- (master.occs.binned$max_ma + master.occs.binned$min_ma)/2
  
  # Reorder bins for later
  lookup <- data.frame('currentbins' = sort(unique(master.occs.binned$bin_assignment)), 
                       'newbins' = seq(from = length(unique(master.occs.binned$bin_assignment)), 
                                       to = 1))
  inds <- match(master.occs.binned$bin_assignment, lookup$currentbins)
  master.occs.binned$new_bins[!is.na(inds)] <- lookup$newbins[na.omit(inds)]
  master.occs.binned <<- master.occs.binned
  
  # Select only necessary bins
  bins.present <- sort(unique(master.occs.binned$bin_assignment))
  bins <<- subset(bins, bin %in% bins.present)
}

################################################################################
# 8. GET COVARIATE FUNCTIONS
################################################################################

# Functions to organise covariate data.

#################################
##### 8.1. FIND_COLLECTIONS #####
#################################

# Extracts information about number of collections per cell of chosen data. 

find_collections <- function(data, single = FALSE){ 
  # Data is output from Get_grid. Res is resolution (only necessary for 
  # later functions)
  
  Collections_per_cell <- data %>% # Counting collections per cell
    dplyr::select(collection_no, siteID) %>%
    dplyr::distinct() %>%
    dplyr::group_by(siteID) %>%
    dplyr::summarize(colls_per_cell = n())
  if(single == TRUE){
    # Removing any cells with 1 collection
    Collections_per_cell <- Collections_per_cell[!Collections_per_cell$colls_per_cell == 1,] 
  }
  Collections_per_cell <<- Collections_per_cell
}

########################
##### 8.2. GET_COV #####
########################

# Attaches grid cell IDs from an inputted raster to occurrences/collections.

get_cov <- function(data, raster, colls = TRUE){
  # data is first output from get_grid. Raster is a chosen raster file, which 
  # can be a raster stack. 
  
  xy <- SpatialPointsDataFrame(cbind.data.frame(data$lng, data$lat), data)
  cov_dat <- raster::extract(raster, xy, sp = TRUE, cellnumbers = FALSE)
  cov_dat <- as.data.frame(cov_dat)
  if(colls == TRUE){
    colls <- find_collections(data)
    cov_dat <- merge(cov_dat, colls, by = "siteID")
  }else{
    cov_dat <- cov_dat
  }
}

##########################
##### 8.3. GET_P_COV #####
##########################

# Attaches grid cell IDs from an inputted raster to occurrences/collections.

get_p_cov <- function(data, raster, colls = TRUE){
  # data is first output from get_grid. Raster is a chosen raster file, which 
  # can be a raster stack. 
  xy <- SpatialPointsDataFrame(cbind.data.frame(data$plng, data$plat), data)
  cov_dat <- raster::extract(raster, xy, sp = TRUE, cellnumbers = FALSE)
  cov_dat <- as.data.frame(cov_dat)
  if(colls == TRUE){
    colls <- find_collections(data)
    p_cov_dat <- merge(cov_dat, colls, by = "siteID")
  }else{
    p_cov_dat <- cov_dat
  }
}

################################################################################
# 9. PREPARE_FOR_RES_DATA AND RES DATA
################################################################################

# Functions to test quality of data at different resolutions of grid cells.

#####################################
##### 9.1. PREPARE_FOR_RES_DATA #####
#####################################

# Produces summary of key stats for data at a specified resolution of grid cell. 
# Used in res_data.

prepare_for_res_data <- function(data, target, single = TRUE){ 
  # Data is first output from combine_data (fossil.colls). Target is chosen 
  # taxon group of interest.
  
  # Select appropriate cells
  rel_data <- data %>% 
    dplyr::select(collection_no, siteID, Target) 
  
  # Make anything that's not the target group a 0
  rel_data$Target[rel_data$Target != target] <- 0 
  
  # Make target's a 1
  rel_data$Target[rel_data$Target == target] <- 1 
  
  # Make all NA's a 0
  rel_data$Target[is.na(rel_data$Target)] <- 0 
  rel_data$Target <- as.numeric(rel_data$Target)
  
  # For all collections, give a mean score of presences and absences
  coll_data <- rel_data %>% 
    dplyr::group_by(collection_no) %>%
    dplyr::summarize(mean(Target)) 
  
  # Anything above a 0 has presences, therefore can be counted as 1
  coll_data$`mean(Target)` <- ceiling(coll_data$`mean(Target)`)
  
  #Join with cell IDs, remove old target, clean column name
  joined_data <- dplyr::left_join(coll_data, rel_data) %>% 
    dplyr::select(-Target, Pres.Abs = `mean(Target)`, collection_no) 
  
  # Remove duplicates of remaining collections
  joined_data <- joined_data %>% dplyr::distinct() 
  cells.removed <- NA
  
  if(single == TRUE){
    
    # Create table for removing singleton cells (cells with only one collection)
    id.table <- table(joined_data$siteID) 
    
    # Record how many cells removed
    cells.removed <- sum(id.table == 1) 
    removed <- subset(joined_data, siteID %in% names(id.table[id.table == 1]))
    
    # Remove cells with only one collection
    joined_data <- subset(joined_data, siteID %in% names(id.table[id.table > 1])) 
  }
  
  # Get data for calculating naive occupancy
  prestest<- joined_data %>% 
    dplyr::group_by(siteID) %>%
    dplyr::summarize(ceiling(mean(Pres.Abs)))
  rem.prestest <- removed %>%
    dplyr::group_by(siteID) %>%
    dplyr::summarize(ceiling(mean(Pres.Abs)))
  
  # Generate results
  results <- c(nrow(prestest), # number of cells
               sum(prestest$`ceiling(mean(Pres.Abs))`), # Number of occupied cells
               sum(prestest$`ceiling(mean(Pres.Abs))`)/nrow(prestest)*100, #naive occupancy
               nrow(joined_data), # total number of collections
               mean(table(joined_data$siteID)), # mean number of collections in each cell
               min(table(joined_data$siteID)), # min number of collections in each cell
               max(table(joined_data$siteID)), # max number of collections in each cell
               median(table(joined_data$siteID)), # median number of collections in each cell
               cells.removed, # Number of cells with one collection removed
               sum(rem.prestest$`ceiling(mean(Pres.Abs))`)
               ) 
  results <<- results
}

#########################
##### 9.2. RES_DATA #####
#########################

# Carries out prepare_for_res_data over a sequence of resolutions, and outputs 
# as a data.frame.

res_data <- function(data, target, single = TRUE, vect, formCells = "N"){ 
  # Data is first output from combine_data (fossil.colls). Target is chosen 
  # group to test. Vect is a vector of resolutions to find data for.
  
  s1 <- vect
  Res_results_list <- list()
   for(t in 1:length(target)){
    Res_results <- data.frame(matrix(ncol = 10, nrow = length(s1)))
    colnames(Res_results) <- c("No.Cells", "Occupied.cells", "Naive.occ", 
                               "Total.Colls", "Mean.Colls", "Min.Colls", 
                               "Max.Colls", "Median.Colls", 
                               "No.Singleton.Cells.Removed.",
                               "No.Singleton.Targeted.Cells.Removed")
    row.names(Res_results) <- s1
    for (i in (1:length(s1))){
      test <- get_grid(data, s1[i], formCells = formCells)
      if (single == FALSE){
        prepare_for_res_data(test, target[t], single = FALSE)
      }
      else {
        prepare_for_res_data(test, target[t])
      }
      Res_results[i,]<- results
    }
  Res_results_list[[t]] <- Res_results
   }
names(Res_results_list) <- target
Res_results_list <<- Res_results_list
}


################################################################################
# 10. SPARTA FUNCTIONS
################################################################################

# Selection of functions to generate occupancy probability estimates using the 
# 'sparta' package. Some of these functions have been adapted from functions
# available in the 'sparta' package.

#########################
##### 10.1. GET_OCC #####
#########################

# Modified version of the function plot.occDet() from the package 'sparta',
# available on github. Takes occupancy results from occDetFunc() and provides them
# in a dataframe for plotting at a later date.

get_occ <- function(x, y = NULL, main = x$SPP_NAME, reg_agg = '', ...){
  
  # gets summary output from the BUGS files 
  spp_data <- as.data.frame(x$BUGSoutput$summary)
  
  if(reg_agg != '') reg_agg <- paste0('.r_', reg_agg)
  
  # get rows we are interested in
  ### take psi.fs rows - these are the yearly proportion of occupied cells ###
  spp_data$X <- row.names(spp_data)
  new_data <- spp_data[grepl(paste0("^psi.fs", reg_agg, "\\["),spp_data$X),]
  new_data$year <- (Year = (x$min_year - 1) + 
                      as.numeric(gsub(paste0("psi.fs", reg_agg), "", 
                                      gsub("\\[|\\]","", row.names(new_data)))))
  
  # rename columns, otherwise ggplot doesn't work properly    
  names(new_data) <- gsub("2.5%","quant_025", names(new_data))
  names(new_data) <- gsub("97.5%","quant_975", names(new_data))
  
  # Add rhat T/F column
  new_data$rhat_threshold[new_data$Rhat < 1.1] <- 'Good (<1.1)'
  new_data$rhat_threshold[new_data$Rhat > 1.1] <- 'Bad (>1.1)'
  
  occ_summary <<- new_data
}

#########################
##### 10.2. GET_DET #####
#########################

# Modified version of the function plot_DetectionOverTime() from the package
# 'sparta' available on github. Takes detection results from occDetFunc() and 
# provides them in a dataframe for plotting at a later date.

get_det <- function (model, min.yr = NULL, CI = 95){ 
  
  if ((CI > 100) | (CI <= 0)) 
    stop("Credible intervals must be between 0 and 100")
  CI2q <- function(CI) {
    q <- (1 - CI/100)/2
    return(c(q, 1 - q))
  }
  quant <- CI2q(CI)
  sims_list <- model$BUGSoutput$sims.list
  pDet1 <- sims_list$alpha.p
  if ("beta1" %in% names(sims_list)) {
    pDet1 <- apply(pDet1, 2, function(x) x + 180 * sims_list$beta1[, 
                                                                   1] + 180^2 * sims_list$beta2[, 1])
  }
  if ("LL.p" %in% names(sims_list)) {
    pDet2 <- pDet1 + sims_list$LL.p * log(2)
    pDet4 <- pDet1 + sims_list$LL.p * log(4)
  }
  else if ("dtype2.p" %in% names(sims_list)) {
    pDet2 <- pDet1 + sims_list$dtype2.p[, 1]
    pDet4 <- pDet1 + sims_list$dtype3.p[, 1]
  }
  else {
    pDet2 <- pDet4 <- NA
  }
  pDet <- melt(list(pDet1, pDet2, pDet4))
  names(pDet) <- c("it", "year", "lgt_pDet", "ListLength")
  pDet$ListLength[pDet$ListLength == 3] <- 4
  pDet$pDet <- inv.logit(pDet$lgt_pDet)
  pDet_summary <- ddply(pDet, .(year, ListLength), summarise, 
                        mean_pDet = mean(pDet, na.rm = T), 
                        lower95CI = quantile(pDet, quant[1], na.rm = T), 
                        upper95CI = quantile(pDet, quant[2], na.rm = T))
  if (!is.null(min.yr)) 
    pDet_summary$year <- pDet_summary$year + min.yr - 1
  pDet_summary <<- pDet_summary
}

###########################
##### 10.3. NAIVE_RES #####
###########################

# Function that provides naive occupancy results based on a list of target taxa
# and a dataframe of occurrences.

naive_res <- function(target, data){
  # Target is a vector of named taxonomic targets; data is a dataframe of occurrences.
  # It must have bin_midpoint present to be able to appropriately bin taxa.
  
  # Make empty results table
  results <- data.frame(matrix(ncol = 6, nrow = 0))
  
  # For each target taxon, calculate raw occupancy through time
  for(i in target){
    temp <- data
    
    # Set everything that's not the target to 0
    temp$Target[which(temp$Target != i | is.na(temp$Target))] <- 0
    
    # Set the target to 1 and make column numeric
    temp$Target[which(temp$Target == i)] <- 1
    temp$Target <- as.numeric(temp$Target)
    
    # Summarise data to establish occupied vs non-occupied grid cells per bin
    temp <- temp %>%
      dplyr::select(bin_midpoint, siteID, Target) %>%
      dplyr::distinct() %>%
      dplyr::group_by(bin_midpoint, siteID) %>%
      dplyr::summarise(sites = ceiling(mean(Target))) 
    
    # Make into dataframe and rearrange data
    temp <- as.data.frame(table(temp$bin_midpoint, temp$sites))
    temp <- dcast(temp, Var1 ~ Var2)
    temp$Total <- temp$`0`+temp$`1`
    temp$perc <- temp$`1`/temp$Total
    temp$Target <- i
    
    # Save to results
    results <- rbind(results, temp)
  }
  
  # Make results numeric for plotting
  results$Var1 <- as.numeric(as.character(results$Var1))
  results <<- results
}

##########################
##### 10.4. COMB_RES #####
##########################

# Combines results produced from sparta models. Is run from within 'run_model'.

comb_res <- function(occ, det, naive, target){
  # combining results
  occ.res <- occ %>%
    dplyr::select(mean, quant_025, quant_975, year, rhat_threshold) %>%
    gather(Data, value, c(mean))
  colnames(occ.res)[1] <- "lower95CI"
  colnames(occ.res)[2] <- "upper95CI"
  occ.res$group <- "occ.det"
  occ.res$Data[which(occ.res$Data == "mean")] <- "Mean occupancy"
  
  det.res <- data.frame(matrix(ncol = 3, nrow = 0))
  for (i in unique(det$ListLength)){
    temp <- det %>%
      filter(ListLength == i) %>%
      gather(Data, value, c(mean_pDet)) %>%
      mutate(Data = paste(Data, ListLength, sep = "_")) %>%
      dplyr::select(-ListLength)
    det.res <- rbind(det.res, temp)
  }
  det.res$group <- "occ.det"
  det.res$Data[which(det.res$Data == "mean_pDet_1")] <- "Mean detection (LL1)"
  det.res$Data[which(det.res$Data == "mean_pDet_2")] <- "Mean detection (LL2)"
  det.res$Data[which(det.res$Data == "mean_pDet_4")] <- "Mean detection (LL4)"
  det.res$rhat_threshold <- NA
  
  naive <- naive %>%
    filter(Target == target) %>%
    mutate(year = rev(seq(1:nrow(.))), Data = 'Naive occupancy') %>%
    dplyr::select(year, perc, Data) %>%
    `colnames<-`(c("year", "value", "Data"))
  naive$group <- "Naive"
  naive$rhat_threshold <- NA
  naive$lower95CI <- NA
  naive$upper95CI <- NA
  
  complete <-rbind(naive, det.res, occ.res)
  
  # Add midpoint
  lookup <<- data.frame('currentbins' = sort(unique(master.occs.grid$bin_assignment)), 
                        'newbins' = 
                          seq(from = length(unique(master.occs.grid$bin_assignment)), 
                              to = 1), 
                        'bin_midpoint' = sort(bins$bin_midpoint))
  
  inds <- match(complete$year, lookup$newbins)
  complete$new_bins <- lookup$bin_midpoint[inds]
  
  assign(paste(target, ".res.comb", sep = ""), complete, envir = .GlobalEnv)
}

###########################
##### 10.5. RUN_MODEL #####
###########################

# Function that runs 'sparta' models and generates a table of results for plotting

run_model <- function(data, target){
  # data is an occurrence dataframe; target is a vector of the targeted taxa used 
  # in analysis.
  
  #==== Occupancy model ====
  # run the model with these data for one species
  if(bin.type == "formation"){
    formattedOccData <- sparta::formatOccData(taxa = data$family,
                                              site = data$siteID,
                                              survey = data$mid_ma,
                                              replicate = data$collection_no,
                                              closure_period = data$new_bins)
  }else{
    formattedOccData <- sparta::formatOccData(taxa = data$family,
                                              site = data$siteID,
                                              survey = data$bin_midpoint,
                                              replicate = data$collection_no,
                                              closure_period = data$new_bins)
  }
  
  
  # Initiate the cluster
  sfInit(parallel = TRUE, cpus = 4)
  
  # Export data to the cluster
  sfExport('formattedOccData')
  sfExport('type')
  sfExport('model.list')
  
  # Run the model in parallel
  system.time({
    para_out <- sfClusterApplyLB(target, occ_mod_function)
  })
  
  # Name each element of this output by the species
  for(i in  1:length(para_out)) names(para_out)[i] <- para_out[[i]]$SPP_NAM
  
  # Stop the cluster
  sfStop()
  
  #==== Obtaining results
  # Get occupancy and detection data
  all.results <- det.res <- data.frame(matrix(ncol = 9, nrow = 0))
  for(i in 1:length(target)){
    temp.occ <- get_occ(para_out[[i]])
    temp.det <- get_det(para_out[[i]])
    temp.comb <- comb_res(temp.occ, temp.det, results, target[i])
    temp.comb$Target <- target[i]
    all.results <- rbind(all.results, temp.comb)
  }
  
  all.results <- list(all.results, para_out)
  return(all.results)
}

##########################
##### 10.6. PLOT_OCC #####
##########################

# Function to plot occupancy and detection probability through time, using results 
# from 'sparta' models. 

plot_occ <- function(res.comb){
  # Res.comb is results table produced by run_model for an individual target taxon.
  
  ggplot(data = subset(res.comb, Data == "Mean occupancy"), aes(x = new_bins, 
                                                                y = value)) +
    geom_blank(aes(color = Data), data = res.comb) +
    geom_ribbon(data = res.comb, aes(x = new_bins, ymin = lower95CI, 
                                     ymax = upper95CI, fill = Data), alpha = 0.2) +
    ylab("Probability") + 
    xlab("Time (Ma)") +
    scale_x_reverse() +
    deeptime::coord_geo(dat = list("stages"), 
                        xlim = c((max(res.comb$new_bins)+1), (min(res.comb$new_bins-1))), 
                        ylim = c(0, 1)) +
    geom_line(data = subset(res.comb, Data != "Naive occupancy"), 
              aes(x = new_bins, y = value, color = Data)) +
    scale_color_manual(values=c("#BDD7E7", "#6BAED6", "#3182BD", 
                                "#DE2D26", "#252424")) +
    scale_fill_manual(values=c("#BDD7E7", "#6BAED6", "#3182BD", 
                               "#DE2D26", "#FFFFFF")) +
    new_scale_color() +
    geom_point(data = subset(res.comb, Data == "Mean occupancy"),
               aes(x = new_bins, y = value, color = rhat_threshold), size = 2) +
    scale_color_manual(name = 'Rhat', values = c('Bad (>1.1)' = 'white',
                                                 'Good (<1.1)' = '#DE2D26')) +
    geom_smooth(method=lm) +
    theme_few()
}

############################
##### 10.7. PLOT_NAIVE #####
############################

# Function to plot naive occupancy estimates for individual taxonomic groups.

plot_naive <- function(res.comb, uuid){
  # res.comb is results table produced by run_model for an individual target taxon,
  # uuid is the rphylopic uuid for the taxonomic group, selected using the get_uuid() 
  # function. 
  silhouette_df <- data.frame(x = c(70), y = c(0.88), 
                              Data = c("Naive occupancy"))
  naive <- res.comb %>%
    filter(Data == "Naive occupancy")
  ggplot(data = naive, aes(x = new_bins, y = value)) +
    ylab("Proportion of total sites") + 
    xlab("Time (Ma)") +
    scale_x_reverse() +
    geom_phylopic(data = silhouette_df, aes(x = x, y = y), 
                  uuid = uuid, 
                  size = 0.14, 
                  alpha = 1, 
                  color = "dark grey") +
    deeptime::coord_geo(dat = list("stages"), 
                        xlim = c((max(res.comb$new_bins)+1), (min(res.comb$new_bins-1))), 
                        ylim = c(0, 1)) +
    geom_line(aes(x = new_bins, y = value, color = Data)) +
    scale_color_manual(breaks = c("Naïve Occupancy", "PAO", "Occupancy Probability", "Detection Probability"),
                       values=c("#252424", "#DE2D26", "#DE2D26", "#3182BD")) +
    scale_fill_manual(breaks = c("Naïve Occupancy", "PAO", "Occupancy Probability", "Detection Probability"), 
                      values=c("#FFFFFF", "white","#DE2D26", "#3182BD")) +
    # scale_color_manual(values=c("#252424")) +
    theme_few() +
    geom_smooth(method=lm) #+
  #theme(legend.position="none")
}

################################################################################
# 11. SPOCCUPANCY FUNCTIONS
################################################################################

# Variety of functions to enable running occupancy models using 'spOccupancy' with
# occurrence data from the PBDB and a variety of palaeo and modern covariates.

##########################
##### 11.1. P_ROTATE #####
##########################

# Function to palaeo-rotate site IDs for each time bin to enable acquiring 
# relevant covariate data

p_rotate <- function(res, e, site_IDs, bins, bin = NA){
  # Requires resolution, extent, a vector of site IDs and the relevant time bins.
  r <- raster(res = res, ext = e)
  siteCoords <- as.data.frame(xyFromCell(r, site_IDs))
  siteCoords$siteID <- site_IDs
  colnames(siteCoords) <- c("lng", "lat", "siteID")
  siteCoords <<- siteCoords
  p_rotate_list <- list()
  if(is.na(bin) == FALSE){
    bins <- bins %>%
      dplyr::filter(code == !!bin)
  }
  # For each time bin, rotate occurrences to correct position.
  for(b in 1:nrow(bins)){
    age <- bins$mid_ma[b]
    temp_data <- cbind(siteCoords, age)
    temp_data <- palaeorotate(occdf = temp_data, 
                              age = "age",
                              method = "point", 
                              model = "PALEOMAP"
    )
    names(temp_data)[names(temp_data) == "p_lat"] <- "plat"
    names(temp_data)[names(temp_data) == "p_lng"] <- "plng"
    p_rotate_list[[b]] <- temp_data
    names(p_rotate_list)[[b]] <- bins$code[b]
  }
  return(p_rotate_list)
}

#################################
##### 11.2. FORMAT_OCC_COVS #####
#################################

# Function to format the covariates potentially influencing the occupancy of 
# taxa. These covariates are produced from palaeo-climatic models. 

format_occ_covs <- function(p_cov_list){
  # Requires output from function 'extract_p()' (11.3, below).
  wet <- data.frame(siteNo = rep(1:nrow(p_cov_list[[1]])))
  dry <- data.frame(siteNo = rep(1:nrow(p_cov_list[[1]])))
  col <- data.frame(siteNo = rep(1:nrow(p_cov_list[[1]])))
  hot <- data.frame(siteNo = rep(1:nrow(p_cov_list[[1]])))
  ann <- data.frame(siteNo = rep(1:nrow(p_cov_list[[1]])))
  
  for(t in rev(names(p_cov_list))){
    temp_ann_sd<- p_cov_list[[t]] %>%
      dplyr::select(ann_sd.1)
    colnames(temp_ann_sd) <- t
    ann <- cbind(ann, temp_ann_sd)
    
    temp_wet <- p_cov_list[[t]] %>%
      dplyr::select(wet_mean.1)
    colnames(temp_wet) <- t
    wet <- cbind(wet, temp_wet)
    
    temp_dry <- p_cov_list[[t]] %>%
      dplyr::select(dry_mean.1)
    colnames(temp_dry) <- t
    dry <- cbind(dry, temp_dry)
    
    temp_col <- p_cov_list[[t]] %>%
      dplyr::select(col_mean.1)
    colnames(temp_col) <- t
    col <- cbind(col, temp_col)
    
    temp_hot <- p_cov_list[[t]] %>%
      dplyr::select(hot_mean.1)
    colnames(temp_hot) <- t
    hot <- cbind(hot, temp_hot)
  }
  hot <- hot %>%
    dplyr::select(-siteNo)
  col <- col %>%
    dplyr::select(-siteNo)
  wet <- wet %>%
    dplyr::select(-siteNo)
  dry <- dry %>%
    dplyr::select(-siteNo)
  ann <- ann %>%
    dplyr::select(-siteNo)
  
  if(length(names(p_cov_list)) > 1){
    Year <- data.frame(teyeq = rep(1, length(site_IDs)), 
                       teyep = rep(2, length(site_IDs)), 
                       teyeo = rep(3, length(site_IDs)),
                       teyen = rep(4, length(site_IDs)))
  }
  site.effect <- c(1:length(site_IDs))
  if(length(names(p_cov_list)) > 1){
    occ_covs <- list(wet, dry, hot, col, ann, Year, site.effect)
    names(occ_covs) <- c("wet", "dry", "hot", "col", "ann", "Year", "site.effect")
  }else{
    occ_covs <- list(wet, dry, hot, col, ann, site.effect)
    names(occ_covs) <- c("wet", "dry", "hot", "col", "ann", "site.effect")
  }
  return(occ_covs)
}

###########################
##### 11.3. EXTRACT_P #####
###########################

# Function for extracting the palaeo co-ordinates for sites necessary to calculate 
# relevant covariates. Additionally finds the nearest relevant sed. flux value for
# occurrences plotted within the Western Interior Seaway (which has no associated
# sed. flux values).

extract_p <- function(p_rotate_list){
  # Takes output from 'p_rotate()' function.
  full_p_covs <- list()
  for(t in names(p_rotate_list)){
    wc <- list.files(paste("Prepped_data/Covariate_Data/All_data/", 
                           res, "deg/Palaeo/", sep = ""), 
                     pattern=paste0("^", t, ".*", sep = ""))
    stacked <- raster::stack(paste("Prepped_data/Covariate_Data/All_data/", 
                                   res, "deg/Palaeo/", wc, 
                                   sep =""))
    full_p_covs[[which(!is.na(str_locate(bins$bin, t)), arr.ind=TRUE)[1,1]]] <- get_p_cov(p_rotate_list[[t]], stacked, colls = FALSE)
    colnames(full_p_covs[[which(!is.na(str_locate(bins$bin, t)), arr.ind=TRUE)[1,1]]]) <- gsub(paste(t, "_", sep = ""), "", colnames(full_p_covs[[which(!is.na(str_locate(bins$bin, t)), arr.ind=TRUE)[1,1]]]))
    names(full_p_covs)[[which(!is.na(str_locate(bins$bin, t)), arr.ind=TRUE)[1,1]]] <- t
    
    # Get nearest value for seds
    xy <- p_rotate_list[[t]][,8:9]
    sed <- raster(paste("Prepped_data/Covariate_Data/All_data/", res, "deg/Palaeo/", 
                        t, "_sed.asc", sep = ""))
    sed <- readAll(sed)
    sampled <- apply(X = xy, MARGIN = 1, FUN = function(xy) sed@data@values[which.min(replace(distanceFromPoints(sed, xy), is.na(sed), NA))])
    # Save only first element (stops issue of equidistant values)
    sampled <- sapply(sampled, "[[", 1)
    full_p_covs[[which(!is.na(str_locate(bins$bin, t)), arr.ind=TRUE)[1,1]]]$sed <- sampled
  }
  full_p_covs <- full_p_covs[!sapply(full_p_covs, is.null)]
  full_p_covs <- full_p_covs
}

##################################
##### 11.4. SAMPLE_FOR_SPOCC #####
##################################

# Function that randomly samples collections (visits) for any sites that contain
# over a specified maximum number of collections. Returns a list of two dataframes;
# one containing a site by collections matrix of detection histories, the other 
# a site by collections matrix of collection numbers.

sample_for_spOcc <- function(for_spOcc, max_val){
  # Requires both the encounter history list produced in 'prepare_for_spOcc()', 
  # and a specified maximum number of collections per site.
  dframe <- for_spOcc[[1]]
  colframe <- for_spOcc[[2]]
  up_colframe <- colframe[,1:max_val]
  up_dframe <- dframe[,1:max_val]
  for(n in 1:nrow(dframe)){ # For each row (site)
    if(any(is.na(colframe[n,])) == FALSE){ #If there are no NAs (max no. of collections)
      samples <- sample(1:ncol(colframe), max_val, replace=FALSE) # Sample 10 positions
      up_colframe[n,] <- colframe[n,samples] # Use those positions to subset columns
      up_dframe[n,] <- dframe[n,samples] # Use those positions to subset columns
    }
    else if(which(is.na(colframe[n,]))[1] > (max_val + 1)){ # If it has more collections than max_val
      samples <- sample(1:(which(is.na(colframe[n,]))[1]-1), max_val, replace=FALSE) # Sample 10 positions
      up_colframe[n,] <- colframe[n,samples] # Use those positions to subset columns
      up_dframe[n,] <- dframe[n,samples] # Use those positions to subset columns
    }else{
      up_colframe[n,] <- colframe[n,1:max_val]
      up_dframe[n,] <- dframe[n,1:max_val]
    }
  }
  up_for_spOcc <- list(up_dframe, up_colframe)
}

###################################
##### 11.5. PREPARE_FOR_SPOCC #####
###################################

# Function that creates a list of site by collections detection histories and 
# collections assessed for detection histories for each time bin, in a format 
# that can be easily adapted for use in functions from the package 'spOccupancy'.

prepare_for_spOcc <- function(data, single = TRUE, bin = NA){ 
  # Requires a time and spatially binned occurrence dataframe with relevant targets.
  rel_data <- data %>% 
    dplyr::select(collection_no, siteID, Target, bin_assignment)
  
  if(is.na(bin) == F){
    rel_data <- rel_data %>%
      dplyr::filter(bin_assignment == bin)
  }
  
  eh_list <- list()
  
  for(b in 1:length(unique(rel_data$bin_assignment))){
    temp_rel_data <- rel_data %>%
      filter(bin_assignment == rev(sort(unique(rel_data$bin_assignment)))[b])
    
    # Make target's a 1
    temp_rel_data$Target[temp_rel_data$Target == target] <- 1 
    temp_rel_data$Target[temp_rel_data$Target != 1] <- NA
    
    # Make anything that's not the target group a 0
    temp_rel_data$Target[is.na(temp_rel_data$Target)] <- 0 
    temp_rel_data$Target <- as.numeric(temp_rel_data$Target)
    
    # For all collections, give a mean score of presences and absences
    temp_coll_data <- temp_rel_data %>% 
      dplyr::group_by(collection_no) %>%
      dplyr::summarize(mean(Target)) 
    
    # Anything above a 0 has presences, therefore can be counted as 1
    temp_coll_data$`mean(Target)` <- ceiling(temp_coll_data$`mean(Target)`) 
    
    # Join with cell IDs, remove old target, clean column name
    temp_joined_data <- dplyr::left_join(temp_coll_data, temp_rel_data) %>% 
      dplyr::select(-Target, Pres.Abs = `mean(Target)`, collection_no) 
    
    # Remove duplicates of remaining collections, then create table for removing 
    # singleton siteID (siteID with only one collection) and remove siteID with 
    # only one collection
    temp_joined_data <- temp_joined_data %>% 
      dplyr::distinct() 
    if(single == TRUE){
      id.table <- table(temp_joined_data$siteID) 
      temp_joined_data <- subset(temp_joined_data, siteID %in% names(id.table[id.table > 1])) 
    }
    
    prestest<- temp_joined_data %>%  
      dplyr::group_by(siteID) %>%
      dplyr::summarize(ceiling(mean(Pres.Abs)))
    
    # Making dataframe for unmarked data
    dframe_for_unmarked <- data.frame(matrix(ncol =  max(table(temp_joined_data$siteID)), 
                                             nrow = length(unique(rel_data$siteID)))) 
    
    # Sort siteID so they match output of covariate data
    temp_joined_data <- temp_joined_data %>% 
      dplyr::arrange(siteID, collection_no) 
    
    colnames(dframe_for_unmarked) <- c(sprintf("y.%d", seq(1,max(table(temp_joined_data$siteID)))))
    row.names(dframe_for_unmarked) <- sort(unique(rel_data$siteID))
    colframe <- dframe_for_unmarked
    for (g in 1:length(unique(rel_data$siteID))){
      counter <- 1
      for (r in 1:nrow(temp_joined_data)){
        if (temp_joined_data[r,3] == sort(unique(rel_data$siteID))[g]){
          dframe_for_unmarked[g, counter] <- as.numeric(temp_joined_data[r, 2])
          colframe[g, counter] <- temp_joined_data[r,1]
          counter <- counter + 1
        }
      }
    }
    combined_list <-list(dframe_for_unmarked, colframe)
    eh_list[[b]] <- combined_list
    names(eh_list)[[b]] <- rev(sort(unique(rel_data$bin_assignment)))[b]
  }
  
  set.seed(rand)
  
  eh_list_samp <- list()
  
  for(s in 1:length(eh_list)){
    eh_list_samp[[s]] <- sample_for_spOcc(eh_list[[s]], max_val = max_val)
    names(eh_list_samp)[[s]] <- rev(sort(unique(rel_data$bin_assignment)))[s]
  }
  assign("eh_list", eh_list_samp, envir = .GlobalEnv)
}  

##############################
##### 11.6. ORGANISE_DET #####
##############################

# Function to format various detection covariates for spOccupancy, including modern 
# geographic influences (e.g. rainfall, MGVF, land use), sampling influences (other 
# occurrences, collections per grid cell) and geological influences (outcrop area,
# sediment flux). Formatted in a list, with appropriate individual formatting for 
# input into spOccupancy functions.

organise_det <- function(siteCoords, extracted_covs, occ_data, bin = NA){
  # Requires the geographic site coordinates (including lat/long and site ID), 
  # a dataframe of extracted covariates, the occurrence dataframe and the bins used.
  
  # Modern det
  wc <- list.files(paste("Prepped_data/Covariate_Data/All_data/", 
                         res, "deg/", sep = ""), 
                   pattern=".asc")
  stacked <- raster::stack(paste("Prepped_data/Covariate_Data/All_data/", 
                                 res, "deg/", wc, 
                                 sep =""))
  covs <- get_cov(siteCoords, stacked, colls = FALSE)
  
  names(covs) <- gsub(".[[:digit:]]", "", names(covs))
  
  # Outcrop
  outcrop <- data.frame(teyeq = covs$Out_teyeq, 
                        teyep = covs$Out_teyep,
                        teyeo = covs$Out_teyeo,
                        teyen = covs$Out_teyen
  )
  if(is.na(bin) == F){
    outcrop <- as.vector(outcrop[bin])
    outcrop <- outcrop[[1]]
  }
  
  # Collections
  coll <- occ_data %>% 
    dplyr::select(collection_no, siteID, bin_assignment) %>%
    dplyr::distinct() %>%
    dplyr::group_by(siteID, bin_assignment) %>%
    dplyr::summarise(colls = n())
  coll <- pivot_wider(coll, names_from = bin_assignment, values_from = colls)
  column_index <- c("teyeq", "teyep", "teyeo", "teyen")
  coll <- as.data.frame(coll[, column_index])
  if(is.na(bin) == F){
    coll <- as.vector(coll[bin])
    coll <- coll[[1]]
    coll <- coll[!is.na(coll)]
  }
  
  # Occs
  occur_red <- occ_data %>% # Get partial list
    dplyr::filter(is.na(Target) == T) %>% # Remove target organisms to avoid circularity
    dplyr::group_by(siteID, bin_assignment) %>%
    dplyr::summarise(occur = n()) %>%
    mutate(combID = paste(siteID, bin_assignment, sep = "_"))
  occur_ful <- occ_data %>% # Get full list of sites
    dplyr::group_by(siteID, bin_assignment) %>%
    dplyr::summarise(occur = n()) %>%
    dplyr::select(siteID, bin_assignment) %>%
    mutate(combID = paste(siteID, bin_assignment, sep = "_"))
  
  # Find IDs of sites removed through removing target organisms. (e.g. sites with only target organisms)
  IDs <- setdiff(occur_ful$combID, occur_red$combID)
  missing_sites <- subset(occur_ful, combID %in% IDs)
  missing_sites$occur <- 0 # Set occurrence to 0
  # Combine datasets and sort
  test <- rbind(occur_red, missing_sites)
  test <- test %>% dplyr::select(-combID)
  occur <- test 
  occur <- pivot_wider(occur, names_from = bin_assignment, values_from = occur)
  occur <- occur[order(occur$siteID),]
  
  column_index <- c("teyeq", "teyep", "teyeo", "teyen")
  occur <- as.data.frame(occur[, column_index])
  if(is.na(bin) == F){
    occur <- as.vector(occur[bin])
    occur <- occur[[1]]
    occur <- occur[!is.na(occur)]
  }
  
  # Palaeo det
  sed <- data.frame(siteNo = rep(1:length(site_IDs)))
  
  for(t in rev(names(extracted_covs))){
    names(extracted_covs[[t]]) <- gsub(".[[:digit:]]", "", names(extracted_covs[[t]]))
    temp_sed <- extracted_covs[[t]] %>%
      dplyr::select(sed)
    colnames(temp_sed) <- t
    sed <- cbind(sed, temp_sed)
  }
  
  sed <- sed %>%
    dplyr::select(-siteNo)
  
  if(is.na(bin) == F){
    sed <- sed[[1]]
  }
  
  # Format for spOccupancy
  det_covs <- list(outcrop = outcrop, 
                   sedflux = sed, 
                   land = as.factor(covs$LANDCVI_multiple), 
                   MGVF = covs$MGVF, 
                   rain = covs$WC_Prec, 
                   temp = covs$WC_Temp, 
                   coll = coll, 
                   occur = occur)
}

#############################
##### 11.7. SITE_REMOVE #####
#############################

# Function to remove any sites that have insufficient data to use within the
# occupancy model (e.g. site which have NA values within the occupancy covariates, 
# or as a result of this have only one total visit across all time periods).

site_remove <- function(eh_list, occ_covs, det_covs, siteCoords, single = TRUE){
  # Find all rows with NAs in occupancy covariates
  removal <- sort(unique(c(which(is.na(occ_covs[[1]]), arr.ind=TRUE)[,1],
                           which(is.na(occ_covs[[2]]), arr.ind=TRUE)[,1], 
                           which(is.na(occ_covs[[3]]), arr.ind=TRUE)[,1], 
                           which(is.na(occ_covs[[4]]), arr.ind=TRUE)[,1])))
  if(length(removal) > 0){
    # Remove those rows from other dataframes, arrange back into lists
    occ_covs <- c(lapply(occ_covs[1:4], function(x) {x <- x[-removal, ]}),
                  lapply(occ_covs[5], function(x) {x <- x[-removal]}))
    eh_list <- lapply(eh_list, lapply, function(x) {x <- x[-removal, ]})
    det_covs <- c(lapply(det_covs[c(1:2, 7)], function(x) {x <- x[-removal, ]}), 
                  lapply(det_covs[3:6], function(x) {x <- x[-removal]}))
    siteCoords <- siteCoords[-removal,]
  }
  if(single == TRUE){
    # Now check for any sites without any history (e.g. due to single collection)
    results <- data.frame(row = 1:nrow(eh_list[[1]][[1]]))
    for(i in 1:length(eh_list)){
      results[,i+1] <- apply(eh_list[[i]][[1]], 1, function(x) all(is.na(x)))
    }
    removal <- results %>%
      filter(V2 == "TRUE" & V3 == "TRUE" & V4 == "TRUE" & V5 == "TRUE") %>%
      .$row
    # Remove those rows from other dataframes, arrange back into lists
    occ_covs <<- c(lapply(occ_covs[1:4], function(x) {x <- x[-removal, ]}), 
                   lapply(occ_covs[5], function(x) {x <- x[-removal]}))
    eh_list <<- lapply(eh_list, lapply, function(x) {x <- x[-removal, ]})
    det_covs <<- c(lapply(det_covs[c(1:2, 7)], function(x) {x <- x[-removal, ]}), 
                   lapply(det_covs[3:6], function(x) {x <- x[-removal]}))
    siteCoords <<- siteCoords[-removal,]
  }else{
    occ_covs <<- occ_covs
    eh_list <<- eh_list
    det_covs <<- det_covs
    siteCoords <<- siteCoords
  }
}

#############################
##### 11.8. TRANSPOSE_EH ####
#############################

# Function to reorganise the detection history lists into an array that is 
# suitable for multi-season spOccupancy. 

transpose_eh <- function(eh_list, target){
  # Requires the detection history list, and a character based name that can be 
  # assigned to the resulting array.
  if(class(eh_list[[1]]) == "list"){
    eh_list <- list(eh_list[[1]][[1]], eh_list[[2]][[1]], 
                    eh_list[[3]][[1]], eh_list[[4]][[1]])
  }
  testArray <- abind(eh_list, along = 3)
  newArray <- aperm(testArray, c(1,3,2))
  assign(paste("EH_array_", target, sep = ""), newArray, envir = .GlobalEnv)
}

###########################
##### 11.9. ARRAY_PREP ####
###########################

# Function to arrange the detection history array and various covariates into the
# correct format for spOccupancy. 

Array_prep <- function(data_suffix, sp = FALSE) {
  # Requires the character based name used in 'transpose_eh()' - 11.8 above.
  # Create the full data frame name
  data_name <- paste0("EH_array_", data_suffix)
  
  # Retrieve the data frame using get()
  data <- get(data_name)
  if(sp == FALSE){
    # Perform operations on the data frame
    revi.data <- list(y = data, 
                      occ.covs = occ_covs, 
                      det.covs = det_covs)
  }else{
    revi.data <- list(y = data, 
                      occ.covs = occ_covs, 
                      det.covs = det_covs,
                      coords = coords)
  }
  revi.data <- revi.data
}

##############################
##### 11.10. DISTANCE_FUN ####
##############################

# Function to take the distances to nearest road for each collection and assign 
# them to the correct places within the detection history site by collection matrix.

distance_fun <- function(eh_list){
  # Requires the detection history list
  road <- read.csv("Data/Covariate_Data/Distance_roads.csv")
  road <- road %>%
    dplyr::distinct()
  final.list <- list()
  for(t in names(eh_list)){
    temp_data <- eh_list[[t]][[2]]
    final.list[[t]] <- temp_data %>% 
      dplyr::mutate(across(everything(), .fns = ~ road$distance[match(., road$collection)]))
  }
  if(is.na(bin) == T){
    Distance <- transpose_eh(final.list, "road_distance")
  }else{
    Distance <- final.list[[1]]
  }
}

############################
##### 11.11. MAKE_TABLE ####
############################

# Function to make a table of results from the output of spOccupancy that can be 
# used for generating appropriate figures and tables. 

make_table <- function(out.sp, target, res, ss = F, bin = NA){
  # Requires best model output from spOccupancy, character based designation of
  # target taxon and the chosen spatial resolution.
  
  # Occupancy
  Mean <- colMeans(as.data.frame(out.sp$beta.samples))
  Quant <- apply(as.data.frame(out.sp$beta.samples), 2, quantile, c(0.025, 0.975))
  Rhat <- out.sp$rhat$beta
  occ <- as.data.frame(t(rbind(Mean, Quant, Rhat)))
  occ$Submodel <- "Occupancy"
  occ$Res <- res
  occ$Group <- target
  occ$Covariate <- rownames(occ)
  if(ss == T){
    occ$Bin <- bin
  }
  # Detection
  Mean <- colMeans(as.data.frame(out.sp$alpha.samples))
  Quant <- apply(as.data.frame(out.sp$alpha.samples), 2, quantile, c(0.025, 0.975))
  Rhat <- out.sp$rhat$alpha
  det <- as.data.frame(t(rbind(Mean, Quant, Rhat)))
  det$Submodel <- "Detection"
  det$Res <- res
  det$Group <- target
  det$Covariate <- rownames(det)
  if(ss == T){
    det$Bin <- bin
  }
  # Random Effects
  if(is.null(out.sp$sigma.sq.p.samples) == F){
    Mean <- mean(out.sp$sigma.sq.p.samples)
    Quant <- apply(as.data.frame(out.sp$sigma.sq.p.samples), 2, quantile, c(0.025, 0.975))
    Rhat <- out.sp$rhat$sigma.sq.p
    var <- as.data.frame(t(rbind(Mean, Quant, Rhat)))
    var$Res <- res
    var$Submodel <- "Detection"
    var$Group <- target
    var$Covariate <- "REV: Site"
    if(ss == T){
      var$Bin <- bin
    }  
    comb <- rbind(occ, det, var)
  }else{
    comb <- rbind(occ, det)
  }
  return(comb)
}

################################################################################
# 12. SAVE_LATTICE
################################################################################

# Function for saving results of unmodelled spatial heterogeneity in detection
# probability

save_lattice <- function(p1, ss = F){
  if(ss == T){
    pdf(paste("Results/spOccupancy/single_season/Figures/Unmodelled.det.", target, ".", res, ".", bin, ".pdf", sep = ""))
    print(p1)
    dev.off()
  }else{
    pdf(paste("Figures/X.unmodelled.det.", target, ".", res, ".pdf", sep = ""))
    print(p1)
    dev.off()
  }
}

################################################################################
# 13. FORESTPLOT2
################################################################################

# Function modified from R package forestplotGG to allow the incorporation of 
# Bayesian Credible Intervals instead of Confidence Intervals. 

forestplot2 <- function(df,
                        name = name,
                        estimate = estimate,
                        CI.max = CI.max,
                        CI.min = CI.min,
                        pvalue = NULL,
                        colour = NULL,
                        shape = NULL,
                        logodds = FALSE,
                        psignif = 0.05,
                        ci = 0.95,
                        ...) {
  
  # Input checks
  stopifnot(is.data.frame(df))
  stopifnot(is.logical(logodds))
  
  # TODO: Add some warnings when name, estimate etc are missing from df and user
  # is not defining the name, estimate etc explicitly.
  
  # Quote input
  name <- enquo(name)
  estimate <- enquo(estimate)
  CI.max <- enquo(CI.max)
  CI.min <- enquo(CI.min)
  pvalue <- enquo(pvalue)
  colour <- enquo(colour)
  shape <- enquo(shape)
  
  args <- list(...)
  
  # TODO: Allow color instead of colour. This will do it, but also breaks other
  # options at the end, so fix those before uncommenting this.
  # args <- enquos(...)
  # if (quo_is_null(colour) && "color" %in% names(args)) {
  #   colour <- args$color
  # }
  
  # Adjust data frame variables
  df <-
    df %>%
    # Convert to a factor to preserve order.
    dplyr::mutate(
      !!name := factor(
        !!name,
        levels = !!name %>% unique() %>% rev(),
        ordered = TRUE
      ),
      # Added here to estimate xbreaks for log odds later
      .xmin = !!CI.min,
      .xmax = !!CI.max,
      # Add a logical variable with the info on whether points will be filled.
      # Defaults to TRUE.
      .filled = TRUE,
      # Add a variable with the estimates to be printed on the right side of y-axis
      .label = sprintf("%.2f", !!estimate)
    )
  
  # Exponentiate the estimates and CIs if logodds
  if (logodds) {
    df <-
      df %>%
      mutate(
        .xmin = exp(.data$.xmin),
        .xmax = exp(.data$.xmax),
        !!estimate := exp(!!estimate)
      )
  }
  
  # If pvalue provided, adjust .filled variable
  if (!rlang::quo_is_null(pvalue)) {
    df <-
      df %>%
      dplyr::mutate(.filled = !!pvalue > !!psignif)
  }
  
  # Plot
  g <-
    ggplot2::ggplot(
      df,
      aes(
        x = !!estimate,
        y = !!name
      )
    )
  
  # If logodds, adjust axis scale
  if (logodds) {
    if ("xtickbreaks" %in% names(args)) {
      g <-
        g +
        scale_x_continuous(
          trans = "log10",
          breaks = args$xtickbreaks
        )
    } else {
      g <-
        g +
        scale_x_continuous(
          trans = "log10",
          breaks = scales::log_breaks(n = 7)
        )
    }
  }
  
  g <-
    g +
    # Add custom theme
    theme_forest() +
    # Add Nightingale colour palette
    scale_colour_ng_d() +
    scale_fill_ng_d() +
    # Add striped background
    geom_stripes() +
    # Add vertical line at null point
    geom_vline(
      xintercept = ifelse(test = logodds, yes = 1, no = 0),
      linetype = "solid",
      size = 0.4,
      colour = "black"
    )
  
  g <-
    g +
    # And point+errorbars
    geom_effect(
      ggplot2::aes(
        xmin = .data$.xmin,
        xmax = .data$.xmax,
        colour = !!colour,
        shape = !!shape,
        filled = .data$.filled
      ),
      position = ggstance::position_dodgev(height = 0.5)
    ) +
    # Define the shapes to be used manually
    ggplot2::scale_shape_manual(values = c(21L, 22L, 23L, 24L, 25L)) +
    guides(
      colour = guide_legend(reverse = TRUE),
      shape = guide_legend(reverse = TRUE)
    )
  
  # Limits adjustment
  #
  # # Extend the shorter x-axis side to mirror the longer one
  # xext <-
  #   c(
  #     df[[quo_name(xmin)]],
  #     df[[quo_name(xmax)]]
  #   ) %>%
  #   abs() %>%
  #   max()
  # g <-
  #   g +
  #   ggplot2::expand_limits(x = c(-xext, xext))
  
  # If no groups specified (through either colour or shape), show estimate values.
  # ### Note: I had to switch back to row number as the y aesthetic in order to use
  # ### continuous scale and to be able to add a secondary axis for labels on the
  # ### right. Any other solution was too time consuming.
  # ### I also had to simplify the fi statement below cause when colour and/or
  # ### shape is specified a legend is added even if they are of length 1L and
  # ### messes up right side visuals.
  # if (
  #   (quo_is_null(colour) || length(unique(df[[quo_name(colour)]])) == 1L) &&
  #   (quo_is_null(shape) || length(unique(df[[quo_name(shape)]])) == 1L)
  # ) {
  # if ( quo_is_null(colour) && quo_is_null(shape)){
  #   g <-
  #     g +
  #     geom_text(
  #       aes(label = .label),
  #       x = 1.1 * xext,
  #       hjust = 1
  #     ) +
  #     expand_limits(x = 1.1 * xext)
  # } else {
  #   g <-
  #     g +
  #     ggplot2::scale_y_continuous(
  #       breaks = df %>% pull(.name_order),
  #       labels = df %>% pull(!!name)
  #     )
  # }
  
  # Pass through graphical parameters and define defaults values for some.
  if ("title" %in% names(args)) {
    g <- g + labs(title = args$title)
  }
  if ("subtitle" %in% names(args)) {
    g <- g + labs(subtitle = args$subtitle)
  }
  if ("caption" %in% names(args)) {
    g <- g + labs(caption = args$caption)
  }
  if ("xlab" %in% names(args)) {
    g <- g + labs(x = args$xlab)
  }
  if (!"ylab" %in% names(args)) {
    args$ylab <- ""
  }
  g <- g + labs(y = args$ylab)
  if ("xlim" %in% names(args)) {
    g <- g + coord_cartesian(xlim = args$xlim)
  }
  if ("ylim" %in% names(args)) {
    g <- g + ylim(args$ylim)
  }
  if ("xtickbreaks" %in% names(args) & !logodds) {
    g <- g + scale_x_continuous(breaks = args$xtickbreaks)
  }
  g
}

################################################################################
# 14. CLEAN_FOR_FIG
################################################################################

# Function for cleaning results tables prior to use in figures. 

clean_for_fig <- function(cov.table){
  cov.table[cov.table == "scale(ann)"] <- "Ann"
  cov.table[cov.table == "scale(dry)"] <- "Dry"
  cov.table[cov.table == "scale(col)"] <- "Cold"
  cov.table[cov.table == "scale(wet)"] <- "Wet"
  cov.table[cov.table == "scale(hot)"] <- "Hot"
  cov.table[cov.table == "scale(Distance)"] <- "Distance"
  cov.table[cov.table == "scale(sedflux)"] <- "Sediment Flux" 
  cov.table[cov.table == "scale(rain)"] <- "Rainfall" 
  cov.table[cov.table == "scale(MGVF)"] <- "MGVF"
  cov.table[cov.table == "scale(occur)"] <- "Occurrences"
  cov.table[cov.table == "scale(outcrop)"] <- "Outcrop area"
  cov.table[cov.table == "factor(Year)2"] <- "Bin 2 (75 Ma)"
  cov.table[cov.table == "factor(Year)3"] <- "Bin 3 (69 Ma)"
  cov.table[cov.table == "factor(Year)4"] <- "Bin 4 (66.7 Ma)"
  cov.table[cov.table == "factor(land)2"] <- "Land cover (2)"
  cov.table[cov.table == "factor(land)3"] <- "Land cover (3)"
  cov.table[cov.table == "factor(land)4"] <- "Land cover (4)"
  cov.table[cov.table == "(Intercept)"] <- "Intercept"
  cov.table[cov.table == "scale(coll)"] <- "Collections"
  cov.table[cov.table == "Random Effect Variance; Site"] <- "REV: Site"
  cov.table[cov.table == "teyen"] <- "Bin 4 (66.7 Ma)"
  cov.table[cov.table == "teyeo"] <- "Bin 3 (69 Ma)"
  cov.table[cov.table == "teyep"] <- "Bin 2 (75 Ma)"
  cov.table[cov.table == "teyeq"] <- "Bin 1 (80.8 Ma)"
  return(cov.table)
}
