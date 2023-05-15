################################################################################
#                       FUNCTIONS FOR OCCUPANCY MODELLING                      #
################################################################################
#                                                                              #
#  PRIMARY AUTHOR: CHRISTOPHER D. DEAN                                         #
#  CO-AUTHOR: LEWIS A. JONES                                                   #
#                                                                              #
#  Selection of functions that work with PBDB data in order to run a variety   #
#  of occupancy models. Functions range from those that visualise occurrences  #
#  in terms of grid cells to those that reorder presence/absence data per      #
#  grid cell so it fits the format of the package unmarked. Information        #
#  regarding each function can be found in the seperate sections below.        #
#                                                                              #
################################################################################

#======================= iPAK AND REQUIRED PACKAGES ============================

# function that automatically installs necessary packages that the user is lacking.

#===== iPAK =====
ipak <- function(pkg){ # Function to install packages. Read in character vector 
                       # of any packages required. 
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
} 

#===== REQUIRED PACKAGES =====

library(beepr)
library(rphylopic)
library(tidyr)
library(raster)
library(plyr)
library(dplyr)
library(lattice)
library(rasterVis)
library(sp)
library(rgdal)
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

#============================== GET_EXTENT =====================================

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

#============================== TARGET_MAKER ===================================

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

#========================== VIEW_CELLS & GEN_RASTER ============================

# Functions to view specific cells of a raster. 

#==== view_cells ====

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

#==== gen_raster ====

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
  raster_for_values <<- raster(res = res, val = full_dframe$Vals, ext = ext)
  plot(raster_for_values)
}

#============================== GET_GRID =======================================

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

#=========================== GET_COV FUNCTIONS =================================

# Functions to organise covariate data.

#===== GET_COLLECTIONS =====

# Extracts information about number of collections per cell of chosen data. 

find_collections <- function(data, single = FALSE){ 
  # Data is output from Get_grid. Res is resolution (only neccessary for 
  # later functions)
  
  Collections_per_cell <- data %>% # Counting collections per cell
    dplyr::select(collection_no, cells) %>%
    dplyr::distinct() %>%
    dplyr::group_by(cells) %>%
    dplyr::summarize(colls_per_cell = n())
  if(single == TRUE){
    # Removing any cells with < 1 collection
    Collections_per_cell <- Collections_per_cell[!Collections_per_cell$colls_per_cell == 1,] 
  }
  Collections_per_cell <<- Collections_per_cell
}

#===== GET_COV =====

# Attaches grid cell IDs from an inputted raster to occurrences/collections.

get_cov <- function(data, raster){
  # data is first output from get_grid. Raster is a chosen raster file, which 
  # can be a raster stack. 
  
  xy <- SpatialPointsDataFrame(cbind.data.frame(data$lng, data$lat), data)
  cov_dat <- raster::extract(raster, xy, sp = TRUE, cellnumbers = TRUE)
  cov_dat <- as.data.frame(cov_dat)
  colls <- find_collections(data)
  colnames(colls)[1] <- "cells.1"
  cov_dat <<- merge(cov_dat, colls, by = "cells.1")
}

#===== GET_COV_FROM_STACK =====

# Creates a raster stack of chosen resolution from previously stored raster files, 
# and attaches associated grid cell IDs to occurrences/collections. 

get_cov_from_stack <- function(data, res){ 
  # data is first output from combine_data (fossil.colls) and chosen resolution. 
  
  # Find rasters of appropriate resolution
  grids <- list.files(paste("Data/Covariate_Data/Formatted/All_data/", 
                            res, "deg/", sep = ""), pattern = "asc") 
  # stack those rasters
  raster <- raster::stack(paste0("Data/Covariate_Data/Formatted/All_data/", 
                                 res, "deg/", grids)) 
  get_cov(data, raster)
  CovStack <<- raster
}

#===== HIRES_COV_DAT =====

# Creates dataframe of covariate data associated with relevant grid cells, taken
# from original hi-resolution rasters. Covariate values are created from the mean 
# value of collections within larger grid cells of chosen resolution. 

alt_cov_grab <- function(data, res, out = T){
  # ADD INFO HERE
  
  wc <- list.files("Data/Covariate_Data/Formatted_For_Precise/")
  wc <- wc[!grepl('0.*', wc)] # Remove folders based on resolution
  wc <- stack(paste0("Data/Covariate_Data/Formatted_For_Precise/", wc, sep = ""))
  projection(wc) <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs" 
  xy <- SpatialPointsDataFrame(cbind.data.frame(data$lng, data$lat), data, 
                               proj4string = CRS("+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs"))
  hires_cov_dat <- raster::extract(wc, xy, sp = TRUE, cellnumbers = TRUE)
  hires_cov_dat <- as.data.frame(hires_cov_dat)
  counting_colls <- hires_cov_dat %>%
    dplyr::select(cells, collection_no) %>%
    dplyr::distinct() %>%
    dplyr::group_by(cells) %>%
    dplyr::summarize(Coll_count = n())
  hires_cov_dat <- hires_cov_dat %>%
    dplyr::group_by(cells) %>%
    dplyr::summarize(mean_DEM = mean(DEM, na.rm = TRUE), 
              mean_MGVF = mean(MGVF, na.rm = TRUE),
              mean_prec = mean(prec, na.rm = TRUE),
              mean_temp = mean(temp, na.rm = TRUE),
              max_DEM = max(DEM, na.rm = TRUE),
              min_DEM = min(DEM, na.rm = TRUE),
              relief = (max(DEM, na.rm = TRUE) - min(DEM, na.rm = TRUE)))
  hires_cov_dat <- cbind(hires_cov_dat, counting_colls$Coll_count)

  if (out == T){
    if(bin.type == "formation"){
      stop("Formation bins unsuitable with outcrop covariates. Please set 'out' argument of alt_cov_grab() to F.")
    }
    wc <- list.files(paste0("Data/Covariate_Data/Formatted_For_Precise/", res, sep = ""))
    wc <- stack(paste0("Data/Covariate_Data/Formatted_For_Precise/", res, "/", wc, sep = ""))
    projection(wc) <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs" 
    xy <- SpatialPointsDataFrame(cbind.data.frame(data$lng, data$lat), data, 
                                 proj4string = CRS("+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs"))
    out_dat <- raster::extract(wc, xy, sp = TRUE, cellnumbers = TRUE)
    out_dat <- as.data.frame(out_dat)
    colnames(out_dat) <- sub("_0.*", "", colnames(out_dat))
    if(bin.type == "stage" | bin.type == "subtage" | bin.type == "scotese"){
      if(bins$stage[t] == "Maastrichtian"){
        out_dat <- out_dat %>%
          dplyr::group_by(cells) %>%
          dplyr::summarize(Maas_terr_outcrop = mean(MOut, na.rm = TRUE), 
                           Maas_all_outcrop = mean(MOut_all, na.rm = TRUE), 
                           Cret_outcrop = mean(Out_all, na.rm = TRUE))
        hires_cov_dat <- cbind(hires_cov_dat, out_dat[,2:4])
      }
      if(bins$stage[t] == "Campanian"){
        out_dat <- out_dat %>%
          dplyr::group_by(cells) %>%
          dplyr::summarize(Camp_terr_outcrop = mean(COut, na.rm = TRUE), 
                           Camp_all_outcrop = mean(COut_all, na.rm = TRUE), 
                           Cret_outcrop = mean(Out_all, na.rm = TRUE))
        hires_cov_dat <- cbind(hires_cov_dat, out_dat[,2:4])
      }
      if(bin.type == "subtage" | bin.type == "scotese"){
        warning(paste("Outcrop covariate produced at stage level. Please bear in mind for further analysis."))
      }
    }
  }
  hires_cov_dat <<- hires_cov_dat
  write.csv(hires_cov_dat, 
            file.path(paste("Results/", bin.type, "/", bin.name, "/", res, 
                            "/site_detection_covs.csv", sep="")))
}

#============================== GET_GRID_IM ====================================

#==== Set up background info ====

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

#==== Get_grid_im ====

# Sets raster to dimensions of inputted data ready for visualisation. Is used in 
# vis_grid. 

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
  mapTheme <- rasterVis::rasterTheme(region=brewer.pal(8,"Greens"))
  
  #create levelplot for raster
  print(rasterVis::levelplot(r, margin=F, par.settings=mapTheme,  
                             main = paste("Total ", (substitute(name)), 
                                          " per Grid Cell", sep = "")) + 
    # Plots state lines
    latticeExtra::layer(sp.polygons(states, col = "white", fill = NA), under = T)  + 
    # Plots background colour
    latticeExtra::layer(sp.polygons(countries, col = 0, fill = "light grey"), under = T)) 
  hist(r, breaks = 20,
       main = paste((substitute(name)), " per Grid Cell", sep = ""),
       xlab = "Number of Collections", ylab = "Number of Grid Cells",
       col = "springgreen")
  r <<- r
}

#====================== PREPARE_FOR_RES_DATA AND RES DATA ======================

# Functions to test quality of data at different resolutions of grid cells.

#===== PREPARE_FOR_RES_DATA =====

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

#===== RES_DATA =====

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

#=============== PREPARE_FOR_UNMARKED AND ALL_RESULTS_FOR_UNMARKED =============

# Functions for converting data into the correct format for unmarked (occupancy 
# modelling package). Can be run individually or for multiple Targets and Resolutions. 


#===== PREPARE_FOR_UNMARKED =====

# WHAT DOES THIS FUNCTION DO?

prepare_for_unmarked <- function(data, target, single = TRUE){ 
  # data is output from Get_Grid. target is specified group to examine. 
  
  # Select relevant info
  rel_data <- data %>% 
    dplyr::select(collection_no, siteID, Target)
  
  # Make target's a 1
  rel_data$Target[rel_data$Target == target] <- 1 
  rel_data$Target[rel_data$Target != 1] <- NA
  
  # Make anything that's not the target group a 0
  rel_data$Target[is.na(rel_data$Target)] <- 0 
  rel_data$Target <- as.numeric(rel_data$Target)
  
  # For all collections, give a mean score of presences and absences
  coll_data <- rel_data %>% 
    dplyr::group_by(collection_no) %>%
    dplyr::summarize(mean(Target)) 
  
  # Anything above a 0 has presences, therefore can be counted as 1
  coll_data$`mean(Target)` <- ceiling(coll_data$`mean(Target)`) 
  
  # Join with cell IDs, remove old target, clean column name
  joined_data <- dplyr::left_join(coll_data, rel_data) %>% 
    dplyr::select(-Target, Pres.Abs = `mean(Target)`, collection_no) 
  
  # Remove duplicates of remaining collections, then create table for removing 
  # singleton siteID (siteID with only one collection) and remove siteID with 
  # only one collection
  joined_data <- joined_data %>% dplyr::distinct() 
  if (single == TRUE){
    id.table <- table(joined_data$siteID) 
    joined_data <- subset(joined_data, siteID %in% names(id.table[id.table > 1])) 
  }
  
  # Get data for calculating naive occupancy
  prestest<- joined_data %>%  
    dplyr::group_by(siteID) %>%
    dplyr::summarize(ceiling(mean(Pres.Abs)))
  
  # Generate results
  results <- c(nrow(prestest), # number of siteID
               sum(prestest$`ceiling(mean(Pres.Abs))`)/nrow(prestest)*100, #naive occupancy
               nrow(joined_data), # total number of collections
               mean(table(joined_data$siteID)), # mean number of collections in each cell
               min(table(joined_data$siteID)), # min number of collections in each cell
               max(table(joined_data$siteID)), # max number of collections in each cell
               median(table(joined_data$siteID))) # median number of collections in each cell
  
  # Making dataframe for unmarked data
  dframe_for_unmarked <- data.frame(matrix(ncol = results[6], nrow = results[1])) 

  # Sort siteID so they match output of covariate data
  joined_data <- joined_data %>% 
    dplyr::arrange(siteID, collection_no) 
  colnames(dframe_for_unmarked) <- c(sprintf("y.%d", seq(1,results[6])))
  row.names(dframe_for_unmarked) <- unique(joined_data$siteID)
  test <- unique(joined_data$siteID)
  colframe <- dframe_for_unmarked
  for (g in 1:results[1]){
    counter <- 1
    for (r in 1:nrow(joined_data)){
      if (joined_data[r,3] == test[g]){
        dframe_for_unmarked[g, counter] <- as.numeric(joined_data[r, 2])
        colframe[g, counter] <- joined_data[r,1]
        counter <- counter + 1
      }
    }
  }

  for_unmarked <- list(dframe_for_unmarked, colframe)
  temp_name <- paste("unmarked_", target, sep = "")
  assign(temp_name, for_unmarked, envir = .GlobalEnv)
}

#===== SAMPLE_FOR_UNMARKED =====

sample_for_unmarked <- function(for_unmarked, max_val){
  dframe <- for_unmarked[[1]]
  colframe <- for_unmarked[[2]]
  up_colframe <- colframe[,1:max_val]
  up_dframe <- dframe[,1:max_val]
  for(n in 1:nrow(dframe)){ # For each row (site)
    if(any(is.na(colframe[n,])) == FALSE){
      samples <- sample(1:ncol(colframe), 10, replace=FALSE) # Sample 10 positions
      up_colframe[n,] <- colframe[n,samples] # Use those positions to subset columns
      up_dframe[n,] <- dframe[n,samples] # Use those positions to subset columns
    }
    else if(which(is.na(colframe[n,]))[1] > (max_val + 1)){ # If it has more collections than max_val
      samples <- sample(1:(which(is.na(colframe[n,]))[1]-1), 10, replace=FALSE) # Sample 10 positions
      up_colframe[n,] <- colframe[n,samples] # Use those positions to subset columns
      up_dframe[n,] <- dframe[n,samples] # Use those positions to subset columns
    }else{
      up_colframe[n,] <- colframe[n,1:max_val]
      up_dframe[n,] <- dframe[n,1:max_val]
    }
  }
  up_for_unmarked <- list(up_dframe, up_colframe)
  temp_name <- paste("up_unmarked_", target, sep = "")
  assign(temp_name, up_for_unmarked, envir = .GlobalEnv)
}

#===== ALL_RESULTS_FOR_UNMARKED =====

# loop that takes basic combined data and writes multiple .csv files into Results 
# folder in current directory for chosen grid cells resolutions and targets in 
# correct format for unmarked. Sound rings when function has finished running.

all_results_for_unmarked <- function(data, res, target, ext, name, single = TRUE, 
                                     formCells = "N", max_val_on = TRUE, max_val = 10){ 
  # data is first output from combined_data (fossil.colls). res is vectors of 
  # chosen resolutions. target is vector of chosen targets. 
  
  for (r in 1:length(res)){
    ptm <- proc.time()
    test1 <- get_grid(data, res[r], ext, formCells = formCells)
    for (q in 1:length(target)){
      if(single == FALSE){
        test2 <- prepare_for_unmarked(test1, target[q], single = FALSE)
        if(max_val_on == TRUE){
          test2 <- sample_for_unmarked(test2, max_val)
        }
      }
      else {
        test2 <- prepare_for_unmarked(test1, target[q])
        if(max_val_on == TRUE){
          test2 <- sample_for_unmarked(test2, max_val)
        }
      }
      # Create folders, remove warning if they already exist.
      dir.create(paste0("Results/", bin.type, "/", sep = ""), showWarnings = FALSE) 
      dir.create(paste0("Results/", bin.type, "/", bin.name, "/", sep =""), 
                 showWarnings = FALSE) 
      dir.create(paste0("Results/", bin.type, "/", bin.name, "/", res, "/", 
                        sep = ""), showWarnings = FALSE) 

      if(max_val_on == TRUE){
        temp_name_1 <- paste(name, ".", res[r], ".", target[q], ".dframe.",  max_val, sep = "")
        temp_name_2 <- paste(name, ".", res[r], ".", target[q], ".colframe.", max_val, sep = "")
      }else{
        temp_name_1 <- paste(name, ".", res[r], ".", target[q], ".dframe",  sep = "")
        temp_name_2 <- paste(name, ".", res[r], ".", target[q], ".colframe", sep = "")
      }
      
      write.csv(test2[[1]], file.path(paste("Results/", bin.type, "/", bin.name, "/", 
                                        res, "/", temp_name_1, ".csv", sep="")))
      write.csv(test2[[2]], file.path(paste("Results/", bin.type, "/", bin.name, "/", 
                                       res, "/", temp_name_2, ".csv", sep="")))
    }
    proc.time() - ptm
  }
  beepr::beep(sound = 3)
}

#=========================PREPARE_FOR_MULTISPECIES =============================

# Function to prepare PBDB data for use in multispecies occupancy modelling. 
# Generates data at a chosen taxonomic level.

prepare_for_multispecies <- function(data, res, ext, level = "genus", target, 
                                     formCells = "N"){
  
  # WHAT IS NEEDED HERE??
  
  # Set up potential inputs
  TYPE <- c("species", "genus") 
  
  # if entered level doesn't match potential inputs, throw error.
  if (is.na(pmatch(level, TYPE))){ 
    stop("Invalid level. Choose either species or genus") 
  }
  
  # Run get_grid so that occurrences have cell reference IDs
  gridded <- get_grid(data, res, ext, formCells = formCells) 
  
  # Make dataframe lookup table out of chosen targets. Assigned numbers for targets 
  # in order of entered targets. 
  targets <- data.frame(target, 1:length(target)) 
  
  # rename columns
  colnames(targets) <- c("Target", "Code") 
  
  # If input equals "species" level
  if (level == "species"){ 
    
    # take output from Get_Grid, filter to only include taxa at species rank, 
    # select only name and targeted info and keep distinct rows
    targeted <- gridded %>% 
      dplyr::filter(accepted_rank == "species") %>% 
      dplyr::select(accepted_name, Target) %>% 
      distinct() 
    
    # Change column name
    colnames(targeted)[1] <- "Name" 
  }
  else{ # If input equals "genus" level
    
    # take output from Get_Grid, filter to only include taxa at species rank, 
    # select only name and targeted info and keep distinct rows
    targeted <- gridded %>% 
      dplyr::select(genus, Target) %>% 
      distinct() 
    
    # Change column name
    colnames(targeted)[1] <- "Name" 
  }
  
  # Take targeted info and turn any blank entries into NAs, then remove any NAs 
  # in "Name" column. Should now be left with just useful names.
  targeted <- targeted %>% dplyr::na_if("") %>% 
    drop_na(Name) 
  
  # Attached Target IDs to genera names.
  targeted <- merge(targeted, targets, by = "Target", all.x = TRUE)[,2:3]  
  
  # Sort into alphabetical order
  targeted <- targeted[order(targeted$Name),] 
  
  # Remove any NAs from targeted column (i.e. any organisms/records we don't 
  # want - all targeted taxa should have a number associated with them)
  target_keep <- na.omit(targeted) 
  
  # If input equals "species" level, use reshape to arrange into taxa x site 
  # matrix, then keep only records that were Targeted (i.e. taxa of interest)
  if(level == "species"){ 
    species.site <- reshape2::dcast(gridded, siteID ~ accepted_name, length) 
    species.site.final <- species.site[, target_keep$Name] 
  }
  # If input equals "genera" level, use reshape to arrange into taxa x site 
  # matrix, then keep only records that were Targeted (i.e. taxa of interest)
  else{ 
    species.site <- dcast(gridded, siteID ~ genus, length) 
    species.site.final <- species.site[, target_keep$Name] 
  }
  # Bind cells back to taxa x site matrix
  species.site.final <- cbind(species.site$siteID, species.site.final) 
  
  # Rename just added column 
  colnames(species.site.final)[1] <- "Cell" 
  
  # Save information to folders and global environment
  temp_name <- paste(deparse(substitute(data)), ".", res, ".", level, 
                     ".multispecies", sep = "")
  dir.create(paste0("Results/", res, sep =""), showWarnings = FALSE)
  write.csv(species.site.final, file.path(paste("Results/", res, "/", temp_name, 
                                                ".csv", sep="")))
  species.site.final <<- species.site.final
  target.cov <<- target_keep$Code
}

#=========================== FORMATION BINNING =================================

#==== Scoring_Grid 1 ====

# WHAT DOES THIS DO...

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

#==== Scoring_Grid 2 ====

# Create a scoring grid ignoring formations with length longer than the 3rd quantile. 
# In his way, long ranging formations don't bias the creation of bins, especially 
# when they appear during the same time interval.

Scoring_Grid_2 <- function(formations, res=0.01) { 
  # Requires formation information. Resolution of time lines is set automatically 
  # at 0.01, but can be adjusted.
  
  # Find max/mins ages of formations
  max_age <- max(formations$max_age) 
  min_age <- min(formations$min_age) 
  
  # Make 1ma bins in sequence based on max/min ages. 0.0045 added to ensure 
  # formation is never exactly equivalent to a bin.
  allbins <- seq(min_age-1.0045, max_age+1.0045, res) 
  
  # makes a matrix for the scoring of each time line in terms of how good it is 
  # to be a bin boundary
  score_grid <- matrix(data = NA, nrow = nrow(formations), ncol = length(allbins)) 
  colnames(score_grid) <- allbins 
  rownames(score_grid) <- formations$Formation 
  
  # Set counter and go through each time line
  counter <- 0
  for(i in allbins) {
    counter <- sum(counter,1)
    
    # Go through each formation 
    for (f in 1:nrow(formations)){ 
      
      # If formation range is less than the 3rd Quantile
      if (formations$Range[f] < quantile(formations$Range, 0.75, type = 7)) { 
        
        # if timeline is between max/min age of a formation (i.e. formation 
        # crosses that line)
        if (i <= formations$max_age[f] && i >= formations$min_age[f]){ 
          
          # Work out how much of formation is younger/older than timeline
          a <- formations$max_age[f] - i 
          b <- i - formations$min_age[f] 
          
          # Calculate range of formation
          range <- formations$max_age[f] - formations$min_age[f] 
          
          # Work out percentage that sits each side of line, reduce score by 
          # that amount
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
      
      # If formation range is longer than mean formation range, skip (bin drawing 
      # isn't affected)
      else { 
        next
      }
    }
  }  
  
  # Remove effect of formations longer than mean formation range
  score_grid <- na.omit(score_grid) 
  
  # Work out mean score for each time bin and add to grid
  means <- colMeans(score_grid) 
  score_grid <- rbind(score_grid, means)
  
  # Output results
  score_grid <<- score_grid
  allbins <<- allbins
}
#================================ NEWBINS ======================================

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

#================================= GET_OCC =====================================

# Modified version of the function plot.occDet() from the package 'sparta',
# available on github. Takes occupancy results from occDetFunc() and provides them
# in a dataframe for plotting at a later date.

get_occ <- function(x, y = NULL, main = x$SPP_NAME, reg_agg = '', ...){
  
  # WHAT DOES THIS FUNCTION NEED????
  
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

#================================= GET_DET =====================================

# Modified version of the function plot_DetectionOverTime() from the package
# 'sparta' available on github. Takes detection results from occDetFunc() and 
# provides them in a dataframe for plotting at a later date.

get_det <- function (model, min.yr = NULL, CI = 95){ 
  
  # WHAT DOES THIS FUNCTION NEED???
  
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

#================================= BINNING =====================================

# WHAT DOES THIS FUNCTION DO???

binning <- function(window){

  # WHAT DOES THIS FUNCTION NEED TO RUN?
  
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
  master.occs.binned <- bin_time(master.occs, bins, method = 'majority') 
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

#=============================== NAIVE.RES =====================================

# WHAT DOES THIS FUNCTION DO?

naive.res <- function(target, data){
  
  # WHAT DOES THIS FUNCTION NEED TO RUN?
  
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
      distinct() %>%
      group_by(bin_midpoint, siteID) %>%
      summarise(sites = ceiling(mean(Target))) 
    
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


#================================= COMB.RES ====================================

# WHAT DOES THIS FUNCTION DO?

comb.res <- function(occ, det, naive, target){
  
  # WHAT DOES THIS FUNCTION NEED TO RUN?
  
  #===== combining results =====
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
  lookup <<- data.frame('currentbins' = sort(unique(master.occs.binned$bin_assignment)), 
                       'newbins' = 
                         seq(from = length(unique(master.occs.binned$bin_assignment)), 
                             to = 1), 
                       'bin_midpoint' = bins$mid_ma)
  
  inds <- match(complete$year, lookup$newbins)
  complete$new_bins <- lookup$bin_midpoint[inds]
  
  assign(paste(target, ".res.comb", sep = ""), complete, envir = .GlobalEnv)
}

#================================= RUN.MODEL ===================================

# WHAT DOES THIS FUNCTION DO?

run.model <- function(data, target){
  
  # WHAT DOES THIS FUNCTION NEED TO RUN?
  
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
    temp.comb <- comb.res(temp.occ, temp.det, results, target[i])
    temp.comb$Target <- target[i]
    all.results <- rbind(all.results, temp.comb)
  }
  all.results <<- all.results
}

#================================= PLOT.OCC ====================================

# WHAT DOES THIS FUNCTION DO?

plot.occ <- function(res.comb){
  
  # WHAT DOES THIS FUNCTION NEED TO RUN?
  
  ggplot(data = subset(res.comb, Data == "Mean occupancy"), aes(x = new_bins, 
                                                                y = value)) +
    geom_blank(aes(color = Data), data = res.comb) +
    geom_ribbon(data = res.comb, aes(x = new_bins, ymin = lower95CI, 
                                     ymax = upper95CI, fill = Data), alpha = 0.2) +
    ylab("Probability") + 
    xlab("Time (Ma)") +
    scale_x_reverse() +
    coord_geo(dat = list("stages"), 
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

#=============================== PLOT.NAIVE ====================================

# WHAT DOES THIS FUNCTION DO?

plot.naive <- function(res.comb){
  
  # WHAT DOES THIS FUNCTION NEED TO RUN?
  
  naive <- res.comb %>%
    filter(Data == "Naive occupancy")
  ggplot(data = naive, aes(x = new_bins, y = value)) +
    ylab("Proportion of total sites") + 
    xlab("Time (Ma)") +
    scale_x_reverse() +
    coord_geo(dat = list("stages"), 
              xlim = c((max(res.comb$new_bins)+1), (min(res.comb$new_bins-1))), 
              ylim = c(0, 1)) +
    geom_line(aes(x = new_bins, y = value, color = Data)) +
    scale_color_manual(values=c("#252424")) +
    theme_few() +
    geom_smooth(method=lm) +
    theme(legend.position="none")
}

#============================== OCCURRENCE.PLOT ================================

# WHAT DOES THIS FUNCTION DO?

occurrence.plot <- function(data, target){
  
  # WHAT DOES THIS FUNCTION NEED TO RUN?
  
  occ.comb <- data.frame(matrix(ncol = 0, nrow = 0))
  for(i in 1:length(target)){
    temp <- master.occs.binned %>%
      dplyr::filter(Target == target[i])
    temp <- as.data.frame(table(temp$bin_midpoint))
    temp$Var1 <- as.numeric(as.character(temp$Var1))
    temp$Family <- target[i]
    occ.comb <- bind_rows(occ.comb, temp)
  }
  ggplot(data = occ.comb, aes(x = Var1, y = Freq, color = Family)) +
    ylab("No. of occurrences") + 
    xlab("Time (Ma)") +
    scale_x_reverse() +
    coord_geo(dat = list("stages"), xlim = c((max(occ.comb$Var1)+1), 
                                             (min(occ.comb$Var1-1))), 
              ylim = c(0, (max(occ.comb$Freq)+10))) +
    geom_line(aes(x = Var1, y = Freq)) +
    theme_few() +
    scale_color_viridis(discrete=TRUE) 
}

if(max_val_on == TRUE){
  test2 <- sample_for_unmarked(test2, max_val)
}

################################################################################
################################################################################
################################################################################

# OLD FUNCTIONS


#===== SUBSAMP_FOR_UNMARKED =====

# Takes data prepared for unmarked, and standardizes it down to a total of X site
# visits (collections) for each gridsquare through subsampling.  

SubSamp_for_unmarked <- function(data, target, sampval = 10, trials = 100){ 
  # Data is output from prepare_for_unmarked. Outputs same data, but subsampled 
  # to sampval site visits. sampval by default set to 10.
  
  new_dframe_for_unmarked <- data.frame()
  if(sampval > ncol(data)){
    stop(paste("Sampling value (", sampval, 
               ") is larger than maximum number of collections per grid cell (", 
               ncol(data),"). Please choose a lower sampling value.", sep = ""))
  }
  
  # For each row in unmarked ready data:
  for (n in 1:nrow(data)){ 
    
    # If number of collections in a site is less than set subsample value
    if(rowSums(is.na(data[n,])) > (NCOL(data)-(sampval+1))){ 
      namecols <- colnames(new_dframe_for_unmarked)
      tempdat <- data[n,1:sampval]
      colnames(tempdat) <- namecols
      
      # Ignore, and add all obs to new dataframe
      new_dframe_for_unmarked <- rbind(new_dframe_for_unmarked, tempdat) 
    }
    else{ # otherwise (total collections (observations) is greater than chosen subsampled value):
      
      # Make temporary dataframe
      temp_dframe <- data.frame(1:sampval) 
      
      # For X trials
      for (t in 1:trials){ 
        
        # Take current row of data
        temp_vec <- data[n,]
        
        # remove NA values
        temp_vec <- temp_vec[!is.na(temp_vec)] 
        
        # Sample a chosen sampling value of the data, without replacement
        temp_dframe <- cbind(temp_dframe, sample(temp_vec, sampval, replace = FALSE)) 
      }
      
      # Remove first column (cleaning)
      temp_dframe <- temp_dframe[2:101] 
      
      # Find mean of columns (i.e. mean number of subsampled target occurrences 
      # within X trials)
      num <- round(mean(colSums(temp_dframe))) 
      
      # # Find mean of columns (i.e. mean number of subsampled target occurrences 
      # within X trials)
      samp <- sample(c(rep(1, num), rep(0, (sampval-num)))) 
      
      # Add subsampled results to dataframe
      new_dframe_for_unmarked <- rbind(new_dframe_for_unmarked, samp) 
      
      # Rename row name
      rownames(new_dframe_for_unmarked)[n] <- rownames(data)[n] 
    }
  }
  new_dframe_for_unmarked <- new_dframe_for_unmarked[,1:sampval]
  colnames(new_dframe_for_unmarked) <-  c(sprintf("y.%d", seq(1,sampval)))
  
  # Comparison between subsampled and original datasets
  ori_data <- rowSums(data, na.rm = T)
  new_data <- rowSums(new_dframe_for_unmarked, na.rm = T)
  comp <- as.data.frame(cbind(ori_data, new_data))
  comp[comp > 1] <- 1
  comp$comp <- ifelse(comp$ori_data > comp$new_data, 1, 0)
  warning(paste("Decrease in naive occupancy of ", sum(comp$comp), 
                " sites, equivalent to ", 
                round((sum(comp$comp)/nrow(comp))*100, digits = 2), 
                "%. Information about lost sites can be found in 'comp'.", sep = ""))
  comp_num <<- sum(comp$comp)
  comp <<- comp
  
  # Assign name
  temp_name <- paste("SS_unmarked_", target, sep = "")
  assign(temp_name, new_dframe_for_unmarked, envir = .GlobalEnv)
}
