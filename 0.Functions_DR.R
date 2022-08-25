#=============================================== FUNCTIONS FOR OCCUPANCY MODELLING ===============================================#
#                                                                                                                                 #
#      PRIMARY AUTHOR: CHRISTOPHER D. DEAN                                                                                        #
#      CO-AUTHOR: LEWIS A. JONES                                                                                                  #
#                                                                                                                                 #
#      Selection of functions that work with PBDB data in order to set up occupancy models. Functions range from those that       #
#      visualise occurrences in terms of grid cells to those that reorder presence/absence data per grid cell so it fits          #
#      the format of the package unmarked. Information regarding each function can be found in the seperate sections below.       #
#                                                                                                                                 #   
#=================================================================================================================================#

#=========================================== iPAK AND REQUIRED PACKAGES ================================================

# function that automatically installs necessary packages that the user is lacking.

#===== iPAK =====
ipak <- function(pkg){ # Function to install packages. Read in character vector of any packages required. 
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
} 

#===== REQUIRED PACKAGES =====

# must load in this order before reading data in, otherwise select within dplyr is overwritten.
library(beepr)
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

#=============================================== GET_EXTENT =============================================================

# Setup raster for resolution and extent. Note: these values should be the same ones used for the file 1.Setup_occupancy_DR.

get_extent <- function(data){
  maxLat <- round_any((max(data$lat) + 1), 0.5)  #get value for defining extent, and increase by x for visualisation purposes
  minLat <- round_any((min(data$lat) - 1), 0.5)  #get value for defining extent, and increase by x for visualisation purposes
  maxLng <- round_any((max(data$lng) + 1), 0.5)  #get value for defining extent, and increase by x for visualisation purposes
  minLng <- round_any((min(data$lng) - 1), 0.5) #get value for defining extent, and increase by x for visualisation purposes
  e <<- extent(minLng, maxLng, minLat, maxLat) # build extent object
}

#=============================================== TARGET_MAKER ===========================================================

# Adds new "Target" column with targeted organisms based on specific requirements of the user. Outputs data as named file ending with "targeted".

target_maker <- function (data, level, target){ # Data is entered data. Level is column to search in. target is vector of chosen organisms.
  for (i in 1:length(target)){
    if(i == 1){ # If Target column doesn't exist:
      filtered <- data %>%
        dplyr::mutate(
          Target = dplyr::case_when(
            !!as.name(level) == target[i] ~ target[i])
        )
    } else { # If Target column already exists
      filtered <- filtered %>%
        dplyr:: mutate(
          Target = dplyr::case_when(
            !!as.name(level) == target[i] ~ target[i], 
            TRUE ~ (as.character(.$Target)) # stops case_when overwriting existing Target data
          )
        )
    }
  }
  temp_name <- paste(deparse(substitute(data)),".", "targeted", sep = "") #Name files based on data entered to function
  assign(temp_name, filtered, envir = .GlobalEnv)
}

#=============================================== VIEW_CELLS & GEN_RASTER ================================================

# Functions to view specific cells of a raster. 

#==== view_cells ====

# Function that extracts data from a raster, and makes a new raster with only a chosen vector of cells.

view_cells <- function(chosen_raster, vector_of_cells, res, zero = FALSE){ # Chosen_raster is the raster with cells that you want to view. vector_of_cells is the vector of cells to view. res is chosen resolution. zero chooses whether other cells are listed as NAs or 0 values (for levelplot)
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
  raster_for_values <- raster(res = res, val = full_dframe$Vals, ext = raster_extent)
  plot(raster_for_values)
}

#==== gen_raster ====

# Function that extracts data from a raster from scratch, using a vector of cells and associated values. Used for quickly visualising covariate data.

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
  plot(raster_for_values)
}

#=============================================== BIN_SPLITTER =======================================================

# Takes combined data and splits it into user defined bins based off a vector. 

bin_splitter <- function(bins, fossils){ # takes first output from combined_data and a vector of time bins. 
  for (s in 1:(length(bins)-1)){
    temp_data <- fossils %>%
      dplyr::filter(mid_ma < bins[s] & mid_ma > bins[s + 1])
    temp_name <- paste(deparse(substitute(fossils)), ".bin", s, sep = "") #Name files based on data entered to function
    assign(temp_name, temp_data, envir = .GlobalEnv)
  }
}

#=============================================== GET_GRID ===========================================================

# Creates a raster of chosen resolution, and attaches associated grid cell IDs to occurrences/collections

get_grid <- function(data, res, e, r = "N"){ # data is first output from combine_data (fossil.colls). Res is chosen resolution in degrees
  if (class(r) == "character"){
    r <- raster(res = res, val = 1, ext = e) # Value must be added because extract uses values
    r <<- r
  }
  crs <- r@crs
  xy <- SpatialPointsDataFrame(cbind.data.frame(data$lng, data$lat), data, proj4string = crs)
  Final <- raster::extract(r, xy, sp = TRUE, cellnumbers = TRUE)
  Final <<- as.data.frame(Final)
}

#=============================================== GET_COV FUNCTIONS ===========================================================

# Functions to organise covariate data

#===== GET_COLLECTIONS =====

# Extracts information about number of collections per cell of chosen data. 

find_collections <- function(data, single = FALSE){ # Data is output from Get_grid. Res is resolution (only neccessary for later functions)
  Collections_per_cell <- data %>% # Counting collections per cell
    dplyr::select(collection_name, cells) %>%
    dplyr::distinct() %>%
    dplyr::group_by(cells) %>%
    dplyr::summarize(colls_per_cell = n())
  if(single == TRUE){
    Collections_per_cell <- Collections_per_cell[!Collections_per_cell$colls_per_cell == 1,] # Removing any cells with < 1 collection
  }
  Collections_per_cell <<- Collections_per_cell
}

#===== GET_COV =====

# Attaches grid cell IDs from an inputted raster to occurrences/collections.

get_cov <- function(data, raster){ # data is first output from get_grid. Raster is a chosen raster file, which can be a raster stack. 
  xy <- SpatialPointsDataFrame(cbind.data.frame(data$lng, data$lat), data)
  cov_dat <- raster::extract(raster, xy, sp = TRUE, cellnumbers = TRUE)
  cov_dat <- as.data.frame(cov_dat)
  colls <- find_collections(data)
  colnames(colls)[1] <- "cells.1"
  cov_dat <<- merge(cov_dat, colls, by = "cells.1")
}

#===== GET_COV_FROM_STACK =====

# Creates a raster stack of chosen resolution from previously stored raster files, and attaches associated grid cell IDs to occurrences/collections. 

get_cov_from_stack <- function(data, res){ # data is first output from combine_data (fossil.colls) and chosen resolution. 
  grids <- list.files(paste("Data/Covariate_Data/Formatted/All_data/", res, "deg/", sep = ""), pattern = "asc") # Find rasters of appropriate resolution
  raster <- raster::stack(paste0("Data/Covariate_Data/Formatted/All_data/", res, "deg/", grids)) # stack those rasters
  get_cov(data, raster)
  CovStack <<- raster
}

#===== HIRES_COV_DAT =====

# Creates dataframe of covariate data associated with relevant grid cells, taken from original hi-resolution rasters. Covariate values are created from the mean value of collections within larger grid cells of chosen resolution. 

alt_cov_grab <- function(data){
  wc <- list.files("Data/Covariate_Data/Formatted_For_Precise/")
  wc <- wc[-1]
  wc <- stack(paste0("Data/Covariate_Data/Formatted_For_Precise/", wc, sep = ""))
  projection(wc) <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs" 
  xy <- SpatialPointsDataFrame(cbind.data.frame(data$lng, data$lat), data, proj4string = CRS("+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs"))
  hires_cov_dat <- raster::extract(wc, xy, sp = TRUE, cellnumbers = TRUE)
  hires_cov_dat <- as.data.frame(hires_cov_dat)
  hires_cov_dat <- hires_cov_dat %>%
    dplyr::group_by(cells) %>%
    dplyr::summarize(mean_DEM = mean(DEM, na.rm = TRUE), 
              mean_MGVF = mean(MGVF, na.rm = TRUE),
              mean_prec = mean(prec, na.rm = TRUE),
              mean_temp = mean(temp, na.rm = TRUE))
  hires_cov_dat <<- hires_cov_dat
}

#=============================================== GET_GRID_IM ===========================================================

# Set up background info
countries <- maps::map("world", plot=FALSE, fill = TRUE) # find map to use as backdrop
states <- maps::map("state", plot = FALSE, fill = TRUE)
countries <<- maptools::map2SpatialPolygons(countries, IDs = countries$names, proj4string = CRS("+proj=longlat")) # Turn map into spatialpolygons
states <<- maptools::map2SpatialPolygons(states, IDs = states$names, proj4string = CRS("+proj=longlat")) # Turn map into spatialpolygons

# Sets raster to dimensions of inputted data ready for visualisation. Is used in vis_grid. 

get_grid_im <- function(data, res, name, ext){ # Data is first output from combine_data (fossil.colls). Res is chosen resolution in degrees. name is user inputted string related to data inputted, for display on graphs. 
  xy <- cbind(as.double(data$lng), as.double(data$lat))
  r <- raster::raster(ext = ext, res = res)
  r <- raster::rasterize(xy, r, fun = 'count')
  #r[r > 0] <- 1 # Remove if you want values instead of pure presence/absence.
  countries <- maps::map("world", plot=FALSE, fill = TRUE) # find map to use as backdrop
  countries <<- maptools::map2SpatialPolygons(countries, IDs = countries$names, proj4string = CRS("+proj=longlat")) # Turn map into spatialpolygons
  mapTheme <- rasterVis::rasterTheme(region=brewer.pal(8,"Greens"))
  print(rasterVis::levelplot(r, margin=F, par.settings=mapTheme,  main = paste("Total ", (substitute(name)), " per Grid Cell", sep = "")) + #create levelplot for raster
    latticeExtra::layer(sp.polygons(states, col = "white", fill = NA), under = T)  + # Plots state lines
    latticeExtra::layer(sp.polygons(countries, col = 0, fill = "light grey"), under = T)) # Plots background colour
  hist(r, breaks = 20,
       main = paste((substitute(name)), " per Grid Cell", sep = ""),
       xlab = "Number of Collections", ylab = "Number of Grid Cells",
       col = "springgreen")
  r <<- r
}

#=============================================== PREPARE_FOR_RES_DATA AND RES DATA ===========================================================

# Functions to test quality of data at different resolutions of grid cells.

#===== PREPARE_FOR_RES_DATA =====

# Produces summary of key stats for data at a specified resolution of grid cell. Used in res_data.

prepare_for_res_data <- function(data, target, single = TRUE){ # Data is first output from combine_data (fossil.colls). Target is chosen taxon group of interest
  rel_data <- data %>% 
    dplyr::select(collection_name, cells, Target) # Select appropriate cells
  rel_data$Target[rel_data$Target != target] <- 0 # Make anything that's not the target group a 0
  rel_data$Target[rel_data$Target == target] <- 1 # Make target's a 1
  rel_data$Target[is.na(rel_data$Target)] <- 0 # Make all NA's a 0
  rel_data$Target <- as.numeric(rel_data$Target)
  coll_data <- rel_data %>% # For all collections, give a mean score of presences and absences
    dplyr::group_by(collection_name) %>%
    dplyr::summarize(mean(Target)) 
  coll_data$`mean(Target)` <- ceiling(coll_data$`mean(Target)`) #anything above a 0 has presences, therefore can be counted as 1
  joined_data <- dplyr::left_join(coll_data, rel_data) %>% #Join with cell IDs
    dplyr::select(-Target, Pres.Abs = `mean(Target)`, collection_name) #Remove old target, clean column name
  joined_data <- joined_data %>% dplyr::distinct() # Remove duplicates of remaining collections
  cells.removed <- NA
  if(single == TRUE){
    id.table <- table(joined_data$cells) # Create table for removing singleton cells (cells with only one collection)
    cells.removed <- sum(id.table == 1) # Record how many cells removed
    joined_data <- subset(joined_data, cells %in% names(id.table[id.table > 1])) # Remove cells with only one collection
  }
  prestest<- joined_data %>%  # Get data for calculating naive occupancy
    dplyr::group_by(cells) %>%
    dplyr::summarize(ceiling(mean(Pres.Abs)))
  results <- c(nrow(prestest), # number of cells
               sum(prestest$`ceiling(mean(Pres.Abs))`)/nrow(prestest)*100, #naive occupancy
               nrow(joined_data), # total number of collections
               mean(table(joined_data$cells)), # mean number of collections in each cell
               min(table(joined_data$cells)), # min number of collections in each cell
               max(table(joined_data$cells)), # max number of collections in each cell
               median(table(joined_data$cells)), # median number of collections in each cell
               cells.removed) # Number of cells with one collection removed
  results <<- results
}

#===== RES_DATA =====

# Carries out prepare_for_res_data over a sequence of resolutions, and outputs as a data.frame.

res_data <- function(data, target, single = TRUE, start, fin, int){ # Data is first output from combine_data (fossil.colls). Target is chosen group to test. start is first resolution, fin is last resolution,int is the interval to count between start and fin.
  s1 <- seq(start, fin, int)
  Res_results_list <- list()
   for(t in 1:length(target)){
    Res_results <- data.frame(matrix(ncol = 8, nrow = length(s1)))
    colnames(Res_results) <- c("No.Cells", "Naive.occ", "Total.Colls", 
                               "Mean.Colls", "Min.Colls", 
                               "Max.Colls", "Median.Colls", "No.Singleton.Cells.Removed.")
    row.names(Res_results) <- s1
    for (i in (1:length(s1))){
      test <- get_grid(data, s1[i])
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

#=============================================== PREPARE_FOR_UNMARKED AND ALL_RESULTS_FOR_UNMARKED ===========================================================

# Functions for converting data into the correct format for unmarked (occupancy modelling package). Can be run individually or for multiple Targets and Resolutions. 

#===== PREPARE_FOR_UNMARKED =====

# Converts data generated by Get_grid into the correct format for unmarked. 

prepare_for_unmarked <- function(data, target, single = TRUE){ # data is output from Get_Grid. target is specified group to examine. 
  rel_data <- data %>% 
    dplyr::select(collection_name, cells, Target)
  rel_data$Target[rel_data$Target == target] <- 1 # Make target's a 1
  rel_data$Target[rel_data$Target != 1] <- NA
  rel_data$Target[is.na(rel_data$Target)] <- 0 # Make anything that's not the target group a 0
  rel_data$Target <- as.numeric(rel_data$Target)
  coll_data <- rel_data %>% # For all collections, give a mean score of presences and absences
    dplyr::group_by(collection_name) %>%
    dplyr::summarize(mean(Target)) 
  coll_data$`mean(Target)` <- ceiling(coll_data$`mean(Target)`) #anything above a 0 has presences, therefore can be counted as 1
  joined_data <- dplyr::left_join(coll_data, rel_data) %>% #Join with cell IDs
    dplyr::select(-Target, Pres.Abs = `mean(Target)`, collection_name) #Remove old target, clean column name
  joined_data <- joined_data %>% dplyr::distinct() # Remove duplicates of remaining collections
  if (single == TRUE){
    id.table <- table(joined_data$cells) # Create table for removing singleton cells (cells with only one collection)
    joined_data <- subset(joined_data, cells %in% names(id.table[id.table > 1])) # Remove cells with only one collection
  }
  prestest<- joined_data %>%  # Get data for calculating naive occupancy
    dplyr::group_by(cells) %>%
    dplyr::summarize(ceiling(mean(Pres.Abs)))
  results <- c(nrow(prestest), # number of cells
               sum(prestest$`ceiling(mean(Pres.Abs))`)/nrow(prestest)*100, #naive occupancy
               nrow(joined_data), # total number of collections
               mean(table(joined_data$cells)), # mean number of collections in each cell
               min(table(joined_data$cells)), # min number of collections in each cell
               max(table(joined_data$cells)), # max number of collections in each cell
               median(table(joined_data$cells))) # median number of collections in each cell
  dframe_for_unmarked <- data.frame(matrix(ncol = results[6], nrow = results[1])) # Making dataframe for unmarked data
  joined_data <- joined_data %>% 
    dplyr::arrange(cells, collection_name) # Sort cells so they match output of covariate data
  colnames(dframe_for_unmarked) <- c(sprintf("y.%d", seq(1,results[6])))
  row.names(dframe_for_unmarked) <- unique(joined_data$cells)
  test <- unique(joined_data$cells)
  for (g in 1:results[1]){
    counter <- 1
    for (r in 1:nrow(joined_data)){
      if (joined_data[r,3] == test[g]){
        dframe_for_unmarked[g, counter] <- as.numeric(joined_data[r, 2])
        counter <- counter + 1
      }
    }
  }
  dframe_for_unmarked <- dframe_for_unmarked
  temp_name <- paste("unmarked_", target, sep = "")
  assign(temp_name, dframe_for_unmarked, envir = .GlobalEnv)
}

#===== SUBSAMP_FOR_UNMARKED =====

# Takes data prepared for unmarked, and standardizes it down to a total of five site visits (collections) for each gridsquare through subsampling.  

SubSamp_for_unmarked <- function(data, target, sampval = 10, trials = 100){ # Data is output from prepare_for_unmarked. Outputs same data, but subsampled to sampval site visits. sampval by default set to 10
  new_dframe_for_unmarked <- data.frame()
  for (n in 1:nrow(data)){ # for each row in unmarked ready data
    if(rowSums(is.na(data[n,])) > (NCOL(data)-(sampval+1))){ # If number of collections in a site is less than set subsample value (sampval)
      new_dframe_for_unmarked <- rbind(new_dframe_for_unmarked, data[n,]) # Ignore, and add all obs to new dataframe
    }
    else{ # otherwise (total collections (observations) is greater than chosen subsampled value (sampval)):
      temp_dframe <- data.frame(1:sampval) # Make temporary dataframe
      for (t in 1:trials){ # For X trials
        temp_vec <- data[n,] # Take current row of data
        temp_vec <- temp_vec[!is.na(temp_vec)] # remove NA values
        temp_dframe <- cbind(temp_dframe, sample(temp_vec, sampval, replace = FALSE)) #sample a chosen sampling value of the data, without replacement
      }
      temp_dframe <- temp_dframe[2:101] # Remove first column (cleaning)
      num <- round(mean(colSums(temp_dframe))) # Find mean of columns (i.e. mean number of subsampled target occurrences within X trials)
      samp <- sample(c(rep(1, num), rep(0, (sampval-num)))) # Generate vector of mean observations, with appropriate observations (1) and lack of detections (0)
      new_dframe_for_unmarked <- rbind(new_dframe_for_unmarked, samp) # Add subsampled results to dataframe
      rownames(new_dframe_for_unmarked)[n] <- rownames(data)[n] # Rename row name
    }
  }
  new_dframe_for_unmarked <- new_dframe_for_unmarked[,1:sampval]
  temp_name <- paste("SS_unmarked_", target, sep = "")
  assign(temp_name, new_dframe_for_unmarked, envir = .GlobalEnv)
}

#===== ALL_RESULTS_FOR_UNMARKED =====

# loop that take basic combined data and writes multiple .csv files into Results folder in current directory for chosen grid cells resolutions and targets in correct format for unmarked. Sound rings when function has finished running.

all_results_for_unmarked <- function(data, res, target, ext, subsamp = TRUE, sampval = 10, single = TRUE){ # data is first output from combined_data (fossil.colls). res is vectors of chosen resolutions. target is vector of chosen targets. 
  for (r in 1:length(res)){
    ptm <- proc.time()
    test1 <- get_grid(data, res[r], ext)
    for (t in 1:length(target)){
      if(single == FALSE){
        test2 <- prepare_for_unmarked(test1, target[t], single = FALSE)
      }
      else {
        test2 <- prepare_for_unmarked(test1, target[t])
      }
      if (subsamp==TRUE){
        test2 <- SubSamp_for_unmarked(test2, target[t], sampval = sampval)
      }
      temp_name <- paste(deparse(substitute(data)), ".", res[r], ".", target[t], sep = "")
      dir.create(paste0("Results/", res, sep =""), showWarnings = FALSE) #stops warnings if folder already exists
      if (subsamp == TRUE){
        write.csv(test2, file.path(paste("Results/", res, "/", temp_name, "SS.csv", sep="")))
      }
      else{
        write.csv(test2, file.path(paste("Results/", res, "/", temp_name, ".csv", sep="")))
      }
    }
    proc.time() - ptm
  }
  beepr::beep(sound = 3)
}

#=============================================== PREPARE_FOR_MULTISPECIES ===========================================================

# Function to prepare PBDB data for use in multispecies occupancy modelling. Generates data at a chosen taxonomic level.

prepare_for_multispecies <- function(data, res, ext, level = "genus", target){
  TYPE <- c("species", "genus") # set up potential inputs
  if (is.na(pmatch(level, TYPE))){ # if entered level doesn't match potential inputs
    stop("Invalid level. Choose either species or genus") # stop script, generate error. 
  }
  gridded <- get_grid(data, res, ext) # Run get_grid so that occurrences have cell reference IDs
  targets <- data.frame(target, 1:length(target)) # Make dataframe lookup table out of chosen targets. Assigned numbers for targets in order of entered targets. 
  colnames(targets) <- c("Target", "Code") # rename columns
  if (level == "species"){ # If input equals "species" level
    targetted <- gridded %>% # take output from Get Grid
      dplyr::filter(accepted_rank == "species") %>% # Filter to only include taxa at species rank
      dplyr::select(accepted_name, Target) %>% # Select only name and Targetted info
      distinct() # Keep distinct rows
    colnames(targetted)[1] <- "Name" # Change column name
  }
  else{ # If input equals "genus" level
    targetted <- gridded %>% # take output from Get Grid
      dplyr::select(genus, Target) %>% # Select only name and Targetted info
      distinct() # Keep distinct rows
    colnames(targetted)[1] <- "Name" # Change column name
  }
  targetted <- targetted %>% dplyr::na_if("") %>% #take targetted info and turn any blank entries into NAs
    drop_na(Name) # Removes any NAs in "Name" column. Should now be left with just useful names
  targetted <- merge(targetted, targets, by = "Target", all.x = TRUE)[,2:3] # attached Target IDs to genera names. 
  targetted <- targetted[order(targetted$Name),] # Sort into alphabetical order
  target_keep <- na.omit(targetted) # Remove any NAs from targetted column (i.e. any organisms/records we don't want - all targetted taxa should have a number associated with them)
  if(level == "species"){ # If input equals "species" level
    species.site <- reshape2::dcast(gridded, cells ~ accepted_name, length) # Use reshape to arrange into taxa x site matrix
    species.site.final <- species.site[, target_keep$Name] # Keep only records that were Targetted (i.e. taxa of interest)
  }
  else{ # If input equals "genera" level
    species.site <- dcast(gridded, cells ~ genus, length) # Use reshape to arrange into taxa x site matrix
    species.site.final <- species.site[, target_keep$Name] # Keep only records that were Targetted (i.e. taxa of interest)
  }
  species.site.final <- cbind(species.site$cells, species.site.final) # Bind cells back to taxa x site matrix
  colnames(species.site.final)[1] <- "Cell" # Rename just added column 
  
  temp_name <- paste(deparse(substitute(data)), ".", res, ".", level, ".multispecies", sep = "")
  dir.create(paste0("Results/", res, sep =""), showWarnings = FALSE) #stops warnings if folder already exists
  write.csv(species.site.final, file.path(paste("Results/", res, "/", temp_name, ".csv", sep="")))
  species.site.final <<- species.site.final
  target.cov <<- target_keep$Code
}


#=============================================== BIN_TIME =============================================================

# Bins occurrences into time bins according to one of 5 methods. Taken from the Palaeoverse package.

bin_time <- function(occdf, bins, method = "mid", reps = 100,
                     scale = "GTS2020", return_error = FALSE) {
  #=== Handling errors ===
  if (is.data.frame(occdf) == FALSE) {
    stop("`occdf` should be a dataframe.")
  }
  if (is.data.frame(bins) == FALSE) {
    stop("`bins` should be a dataframe.")
  }
  
  possible_methods <- c("all", "majority", "random", "point", "mid")
  method_match <- charmatch(method, possible_methods)
  
  if (is.na(method_match) == TRUE) {
    # If the user has entered a non-valid term for the "method" argument,
    # generate an error and warn the user.
    stop("Invalid `method`. Choose either:
  'all', 'majority', 'random', 'point', or 'mid'.")
  } else {
    method <- possible_methods[method_match]
  }
  
  if (scale %in% c("GTS2020", "GTS2012") == FALSE) {
    stop("Invalid `scale`. Choose either 'GTS2020' or 'GTS2012'")
  }
  if (is.numeric(reps) == FALSE) {
    stop("Invalid `reps`. Choose an numeric value.")
  }
  if (is.logical(return_error) == FALSE) {
    stop("Invalid `return_error`.
           Choose a logical value (i.e. TRUE or FALSE).")
  }
  if (class(occdf$max_ma) != class(occdf$min_ma)) {
    stop("Invalid occdf$max_ma or occdf$min_ma.
           Columns should be of the same class.")
  }
  
  if (is.numeric(occdf$max_ma) &&
      max(occdf$max_ma) > max(bins$max_ma)) {
    stop("Maximum age of occurrence data surpasses maximum age of bins")
  }
  
  #=== Sorting non-numeric age designations ===
  if (is.character(occdf$max_ma)) {
    # If entered value for max_ma is character rather than numeric:
    
    # which geological timescale to use?
    if (scale == "GTS2020") {
      df <- palaeoverse::GTS2020
    }
    if (scale == "GTS2012") {
      df <- palaeoverse::GTS2012
    }
    
    # Re-name columns to work with rest of function.
    occdf$tmp_bin <- seq_len(nrow(occdf))
    names(occdf)[names(occdf) == "max_ma"] <- "max_interval"
    names(occdf)[names(occdf) == "min_ma"] <- "min_interval"
    
    # Merge dataframes (max ma)
    occdf <- merge(
      x = occdf,
      y = df[, c("interval_name", "max_ma")],
      by.x = "max_interval",
      by.y = "interval_name",
      all.x = TRUE
    )
    
    # Merge dataframes (min ma)
    occdf <- merge(
      x = occdf,
      y = df[, c("interval_name", "min_ma")],
      by.x = "min_interval",
      by.y = "interval_name",
      all.x = TRUE
    )
    
    # Ensure order of dataframe is maintained after merge
    occdf <- occdf[order(occdf$tmp_bin), ]
    
    occdf <- occdf[, -which(colnames(occdf) == "tmp_bin")]
    
    # If not all intervals can be matched, produce error report
    # and message to fix spellings.
    if (any(is.na(occdf$min_ma)) == TRUE ||
        any(is.na(occdf$max_ma)) == TRUE) {
      # Generate error vector
      error_vec <- which(is.na(occdf$min_ma) | is.na(occdf$max_ma))
      # Should an error vector be returned to the user?
      if (return_error == TRUE) {
        return(error_vec)
      } else {
        # return error message
        stop(paste(c(
          "Unable to match interval to numeric value for all occurrences. Available
  intervals names are accessible via GTS2020 and GTS2012. Please check interval
  spelling for the following rows in `occdf` (note: an error vector can be
  returned with the `return_error` argument):",
          capture.output(print(error_vec))
        ),
        collapse = "\n"
        ))
      }
    }
  }
  
  #=== Reporting Info ===
  
  # Make an empty list that's the length of the occurrence dataframe.
  bin_list <- list()
  bin_list <- sapply(seq_len(nrow(occdf)), function(x) NULL)
  
  # For each occurrence, find all the bins that it is present within, and
  # add as elements to that part of the list.
  for (i in seq_len(nrow(bins))) {
    v <-
      which(occdf$max_ma > bins$min_ma[i] &
              occdf$min_ma < bins$max_ma[i])
    for (j in v) {
      bin_list[[j]] <- append(bin_list[[j]], bins$bin[i])
    }
  }
  
  # Generate id column for data (this is for tracking duplicate rows).
  id <- seq_len(nrow(occdf))
  occdf$id <- id
  
  # Generate empty column for recording the number of bins an occurrence
  # appears in, and empty columns for the new bin allocation and midpoint.
  occdf$n_bins <- NA
  occdf$bin_assignment <- NA
  occdf$bin_midpoint <- NA
  
  # Assign number of bins per occurrence.
  occdf$n_bins <- lengths(bin_list)
  
  # Generate midpoint ages of bins
  bins$mid_ma <- (bins$max_ma + bins$min_ma) / 2
  
  #=== Methods ===
  
  #--- Method 1: Midpoint ---
  if (method == "mid") {
    # If no mid point is present for occurrence age range, add one in a
    # new column.
    rmcol <- FALSE
    if (("mid_ma" %in% colnames(occdf)) == FALSE) {
      occdf$mid_ma <- (occdf$max_ma + occdf$min_ma) / 2
      rmcol <- TRUE
    }
    
    # Assign bin based on midpoint age of the age range
    for (i in seq_len(nrow(bins))) {
      v <-
        which(occdf$mid_ma > bins$min_ma[i] &
                occdf$mid_ma < bins$max_ma[i])
      occdf$bin_assignment[v] <- bins$bin[i]
      occdf$bin_midpoint[v] <- bins$mid_ma[i]
    }
    
    # Remove mid_ma for fossil occurrences (if not already present as input)
    if (rmcol == TRUE) {
      occdf <- occdf[, -which(colnames(occdf) == "mid_ma")]
    }
    
    # Return the dataframe and end the function.
    return(occdf)
  }
  
  #--- Method 2: Point estimates ---
  if (method == "point") {
    # make occurrence list for filling with reps
    occ_list <- list()
    occ_list <- sapply(seq_len(nrow(occdf)), function(x) NULL)
    
    # For each occurrence max/min age, make probability distribution and
    # sample from it. Record that with each occurrence.
    for (i in seq_len(nrow(occdf))) {
      #generate occurrence sequence for sampling
      occ_seq <- seq(from = occdf[i, "min_ma"],
                     to = occdf[i, "max_ma"],
                     by = 0.01)
      #if max/min ages are the same replicate age
      if (length(occ_seq) == 1) {
        occ_list[[i]] <- rep(occ_seq, times = reps)
        next
      }else {
        prob <- dunif(occ_seq,
                      max = max(occ_seq),
                      min = min(occ_seq))
        estimates <-
          sample(
            x = occ_seq,
            size = reps,
            replace = TRUE,
            prob = prob
          )
        occ_list[[i]] <- estimates
      }
    }
    
    occdf$point_estimates <- NA
    #drop cols that are not needed
    occdf <- occdf[, -which(colnames(occdf) == "bin_midpoint")]
    
    occ_df_list <- list()
    occ_df_list <- sapply(1:reps, function(x) NULL)
    
    #add point estimates to each dataframe
    for (i in 1:reps) {
      occdf$point_estimates <- do.call(rbind, occ_list)[, i]
      for (j in seq_len(nrow(bins))){
        vec <- which(occdf$point_estimates <= bins$max_ma[j] &
                       occdf$point_estimates >= bins$min_ma[j])
        occdf$bin_assignment[vec] <- bins$bin[j]
      }
      occ_df_list[[i]] <- occdf
    }
    
    #return list of data
    return(occ_df_list)
  }
  
  
  #--- Method 3: All ---
  if (method == "all") {
    # Duplicate rows by number of bins.
    occdf <- occdf[rep(seq_len(dim(occdf)[1]), occdf$n_bins), ]
    
    # Use id to track unique rows and update bin numbers.
    for (i in id) {
      id_vec <- which(occdf$id == i)
      occdf$bin_assignment[id_vec] <- bin_list[[i]]
    }
    # Add bin midpoints to dataframe
    for (i in seq_len(nrow(occdf))) {
      vec <- which(occdf$bin_assignment == bins$bin[i])
      occdf$bin_midpoint[vec] <- bins$mid_ma[i]
    }
    
    rownames(occdf) <- seq_len(nrow(occdf))
    
    # Return the dataframe and end the function.
    return(occdf)
  }
  
  #--- Method 4: Majority ---
  if (method == "majority") {
    # Setup column for calculating overlap of age range with bin
    occdf$overlap_percentage <- NA
    
    # Run across bin list
    for (i in seq_along(bin_list)) {
      # Dataframe of bins occurrence known to occur in
      tmpbin <- bins[bins$bin %in% bin_list[[i]], ]
      
      # Generate sequence of length 10000 for percentage calculations
      occ_seq <-
        seq(occdf[i, "min_ma"], occdf[i, "max_ma"], length.out = 10000)
      
      # Calculate overlap across known bins
      percentage <- vector()
      for (j in seq_len(nrow(tmpbin))) {
        percentage[j] <-
          (length(
            which(occ_seq >= tmpbin$min_ma[j] &
                    occ_seq <= tmpbin$max_ma[j])
          ) / 10000) * 100
      }

      # Assign bins, bin midpoints and overlap percentage
      occdf[i, "bin_assignment"] <-
        tmpbin$bin[which.max(percentage)]
      occdf[i, "bin_midpoint"] <-
        tmpbin$mid_ma[which.max(percentage)]
      occdf[i, "overlap_percentage"] <-
        percentage[which.max(percentage)]
    }
    return(occdf)
  }
  
  #--- Method 5: Random ---
  if (method == "random") {
    # Generate empty lists for populating
    occ_list <- list()
    occ_list <- sapply(seq_len(nrow(occdf)), function(x) NULL)
    occ_df_list <- list()
    occ_df_list <- sapply(seq_len(reps), function(x) NULL)
    
    # Randomly sample from the list of bins that occurrence appears in, and
    # add to the bin column for the occurrence.
    for (i in seq_along(bin_list)) {
      # Dataframe of bins occurrence known to occur in
      tmpbin <- bins[bins$bin %in% bin_list[[i]], ]
      
      # If occurrence only appears in one bin, assign bin
      if (length(bin_list[[i]]) == 1) {
        occ_list[[i]] <- rep(x = bin_list[[i]], times = reps)
        next
      } else {
        # Randomly sample from possible bins
        occ_list[[i]] <- sample(x = tmpbin$bin,
                                size = reps,
                                replace = TRUE)
      }
    }
    
    #add point estimates to each dataframe
    for (i in 1:reps) {
      occdf$bin_assignment <- do.call(rbind, occ_list)[, i]
      occdf$bin_midpoint <- bins$mid_ma[
        sapply(occdf$bin_assignment, function(x) {
          which(bins$bin == x)}, simplify = TRUE)]
      occ_df_list[[i]] <- occdf
    }
    return(occ_df_list)
  }
}

#=============================================== FORMATION BINNING ===========================================================

#==== Scoring_Grid 1 ====

# Create a scoring grid ignoring formations with length longer than the 3rd quantile. In his way, long ranging formations don't 
# bias the creation of bins, especially when they appear during the same time interval.
Scoring_Grid_1 <- function(formations, res=0.01) { # Requires formation information. Resolution of time lines is set automatically at 0.01, but can be adjusted.
  max_age <- max(formations$max_age) #finds max age of all formations
  min_age <- min(formations$min_age) #finds min age of all formations
  allbins <- seq(min_age-1.0045, max_age+1.0045, res) # Makes 1ma bins in sequence based on max/min ages. 0.0045 added to ensure formation is never exactly equivalent to a bin.
  
  score_grid <- matrix(data = NA, nrow = nrow(formations), ncol = length(allbins)) # makes a matrix for the scoring 
  # of each time line in terms of how good it is to be a bin boundary
  colnames(score_grid) <- allbins # All time lines
  rownames(score_grid) <- formations$Formation # All formations
  
  counter <- 0
  for(i in allbins) { # Go through each time line
    counter <- sum(counter,1) 
    for (f in 1:nrow(formations)){ # go through each formation 
      if (i <= formations$max_age[f] && i >= formations$min_age[f]){ # if timeline is between max/min age of a formation (i.e. formation crosses that line)
        a <- formations$max_age[f] - i # Work out how much of formation is older than timeline
        b <- i - formations$min_age[f] # Work out how much of formation is younger than timeline
        range <- formations$max_age[f] - formations$min_age[f] # Calculate range of formation
        if (a > b){
          score_grid[,counter][f] <- (a/range)*100 # Work out percentage that sits each side of line, reduce score by that amount.
        }
        else{ 
          score_grid[,counter][f] <- (b/range)*100 # Work out percentage that sits each side of line, reduce score by that amount.
        }
      }
      else {
        score_grid[,counter][f] = 100 # Otherwise, just score it 100. 
      }
    }
  }  
  means <- colMeans(score_grid) # Work out mean score for each time bin
  score_grid <- rbind(score_grid, means) # add to grid
  score_grid <<- score_grid # Outputs score_grid for later use
  allbins <<- allbins # Outputs bins for later use
}

#==== Scoring_Grid 2 ====

Scoring_Grid_2 <- function(formations, res=0.01) { # Requires formation information. Resolution of time lines is set automatically at 0.01, but can be adjusted.
  max_age <- max(formations$max_age) #finds max age of all formations
  min_age <- min(formations$min_age) #finds min age of all formations
  allbins <- seq(min_age-1.0045, max_age+1.0045, res) # Makes 1ma bins in sequence based on max/min ages. 0.00045 added to ensure formation is never exactly equivalent to a bin.
  
  score_grid <- matrix(data = NA, nrow = nrow(formations), ncol = length(allbins)) # makes a matrix for the scoring of each time line in terms of how good it is to be a bin boundary
  colnames(score_grid) <- allbins # All time lines
  rownames(score_grid) <- formations$Formation # All formations
  
  counter <- 0
  for(i in allbins) { # Go through each time line
    counter <- sum(counter,1)
    for (f in 1:nrow(formations)){ # go through each formation 
      if (formations$Range[f] < quantile(formations$Range, 0.75, type = 7)) { # If formation range is less than the 3rd Quantile
        if (i <= formations$max_age[f] && i >= formations$min_age[f]){ # if timeline is between max/min age of a formation (i.e. formation crosses that line)
          a <- formations$max_age[f] - i # Work out how much of formation is older than timeline
          b <- i - formations$min_age[f] # Work out how much of formation is younger than timeline
          range <- formations$max_age[f] - formations$min_age[f] # Calculate range of formation
          if (a > b){
            score_grid[,counter][f] <- (a/range)*100
          }
          else{ 
            score_grid[,counter][f] <- (b/range)*100 # Work out percentage that sits each side of line, reduce score by that amount.
          }
        }
        else {
          score_grid[,counter][f] = 100 # Otherwise, just score it 100. 
        }
      }
      else { # If formation range is longer than mean formation range, skip (bin drawing isn't affected)
        next
      }
    }
  }  
  score_grid <- na.omit(score_grid) # Remove effect of formations longer than mean formation range
  means <- colMeans(score_grid) # Work out mean score for each time bin
  score_grid <- rbind(score_grid, means) # add to grid
  score_grid <<- score_grid
  allbins <<- allbins
}
#=============================================== NEWBINS ==============================================================

# Looks at the previously generated score_grid and generates appropriate new bins based on those scores. Boundaries are 
# outputted as a list (binlist). If bins are shorter than 0.5 Ma, they are amalgamated into the bins above and below, and 
# a warning is produced.

newBins <- function(score_grid, formations, bin_limits, allbins, stages, smallamalg = TRUE){ # Takes previously generated score grid, formations, 
  # allbins and stages from DivDyn package. Also require bin_limits, a user made vector of the following: 
  # 1) user chosen time window in which to look to draw bins. Advised to be set at 3 Ma. 
  # 2) Hard maximum age of bins
  # 3) Hard minimum age of bins
  
  score_grid<- as.data.frame(score_grid)
  bin_size <- bin_limits[1]
  max_age <- bin_limits[2]
  min_age <- bin_limits[3]
  form_bins <- c(min_age)
  testbin <- seq(min_age, max_age, bin_size) # Creates broader bins for testing best point to draw a bin
  
  # Drawing bins and giving form_bins (vector of bin boundaries)
  for (i in 1:length((testbin)-1)){
    seqs <- seq(testbin[i],testbin[i]+bin_size, 1) # Creates a sequence of ages to draw bins within
    pasting <- c()
    for (n in 1:bin_size){
      pasting <- c(pasting, paste("^",seqs[n],"|", sep = "")) # Sets up expression to match to score_grid
    }
    testmatch2 <- paste(pasting, collapse = "")
    testmatch2 <- substr(testmatch2, 1, nchar(testmatch2)-1) 
    a <- score_grid[grep(testmatch2, names(score_grid))] # Finds all scores within this time window
    z <- apply(a,1,which.max) # Finds maximum bin score within this time window
    form_bins[i+1] <- names(a)[z][nrow(a)] # Adds maximum bin score to a vector of bins
  }
  form_bins <- as.numeric(unique(form_bins)) # Finds all unique bins
  if (smallamalg == TRUE){ # If small bin amalgamation is turned on (is on automatically):
    range <- (diff(form_bins) < 0.5) # Finds all bins which are under 0.5 Ma in length
    range_checker <- c()
    for (r in 1:length(range)){
      if (range[r] == TRUE){ # If a bin is under 0.5 Ma in length
        difference <- diff(c(form_bins[r], form_bins[r+1])) # Find the length of that bin
        warning("Original bin ",  r, " removed due to small range: ~", signif(difference, digits = 3), " Ma. The difference in time has been added to the bins above and below.") # Generate warning about bin amalgamation
        form_bins[r] <- form_bins[r]+(difference/2) # Adds half length of old bin to bin below
        form_bins[r+1] <- form_bins[r+1]-(difference/2) # Adds half length of old bin to bin above
        range_checker <- c(range_checker, r) # Records which bin was too small
      }
    }
    if(length(range_checker) > 0){ # If there have been amalgamated bins
      form_bins <- form_bins[-range_checker] # Remove old amalgamated bins
    }
  }
  if (smallamalg == FALSE){
    warning("Small bin amalgamation is turned off. Bins may be too short to record occurrences. You are advised to check bins before running further analyses.") # Generate warning if small bin amalgamation is turned off. 
  }
  form_bins <<- form_bins
  
  # Creating binlist (data.frame of bins and appropriate age info)
  prefix <- "FB."
  suffix <- seq(1:(length(form_bins)-1))
  my_names <- paste(prefix, suffix, sep = "")
  binlist <- data.frame(bin = my_names, # Combines bin data to make dataframe of minimum, maximum and mid point of each new bin
                        bottom = as.numeric(form_bins[1:(length(form_bins)-1)]), 
                        top = as.numeric(form_bins[2:(length(form_bins))]))
  binlist$mid <- (binlist$bottom + binlist$top) / 2
  binlist <<- binlist
  
  #Plotting new bins
  par(mar = c(4.1, 4.1, 1, 2.1))
  tsplot(stages, boxes=c("short","system"), # Generates plot using DivDyn package
         xlim=c(stages$bottom[1], stages$top[nrow(stages)]),  ylim=c(min(colMeans(score_grid), na.rm = TRUE), 100), 
         prop = 0.08, plot.args = list(cex.lab = 2, cex.axis = 2),
         shading=NULL, boxes.col=c("col","systemCol"), labels.args=list(cex=1.8),
         ylab = "Bin Splitting Score") 
  for(n in 1:length(form_bins)){ # draws new bins as coloured boxes for comparison to traditional bins
    if(((n %% 2) == 0) == TRUE) next
    else {
      if(n == length(form_bins)){
        if(nrow(binlist) %% 2 == 0){
          next
        }
        else{
          rect(useful_bins[n], 0, useful_bins[n-1], 100, col = rgb(0.89,0.89,0.89,alpha=0.5), border = NA)
        }
      }
      else{
        rect(form_bins[n], 0, form_bins[n+1], 100-0.01, col = rgb(0.89,0.89,0.89,alpha=0.5), border = NA)
      }
    }
  }
  lines(allbins, colMeans(score_grid), lwd = 1.75) # draws bin splitting score on plot
  box(lwd=2)
}

#=============================================== FORMATIONGRAPH =======================================================

# Shows what formations look like through time in comparison to Stages and new Bins.
# Requires formations, form_bins from newBins function and stages from DivDyn package

FormationGraph <- function(formations, form_bins, stages, score_grid_2 = FALSE, draw_by = "Lat", Col = "None", legend = TRUE, STAGE = FALSE){ # Requires formations, form_bins from newBins function and stages from DivDyn package. Has a number of user defined options:
  # score_grid_2: When turned to TRUE, long ranging formations will appear in grey.
  # draw_by: Arranges y axis (and formations) in different ways. Takes the arguments "Lat" (arrange by mean latitude of occurrences for formation), "Max_Age" (arrange by maximum age), and "Number" (arrange by arbitrary number).
  # Col: Colours formations by either "Diversity" (taxic diversity per formation), "Occurrences" (number of occurrences per formation), or "None" (No colour)
  # Legend: When TRUE, draws a legend.
  # STAGE: When TRUE, draws stage bins rather than updated bins.
  TYPE <- c("Lat", "Max_Age", "Number")
  if (is.na(pmatch(draw_by, TYPE))){
    stop("Invalid drawing method. Choose either 'Lat', 'Max_Age', or 'Number'.") # If the user has entered a non-valid term for the "draw_by" argument, generate an error and warn the user. 
  }
  TYPE <- c("None", "Diversity", "Occurrences")
  if (is.na(pmatch(Col, TYPE))){
    stop("Invalid colouring method. Choose either 'None', 'Diversity', or 'Occurrences'.") # If the user has entered a non-valid term for the "Col" argument, generate an error and warn the user.
  }
  fp1 <- data.frame(matrix(ncol = 5, nrow = 0)) # Setup dataframe for recording formations with short durations
  fplong <- data.frame(matrix(ncol = 5, nrow = 0)) # Setup dataframe for recording formations with long durations 
  for (f in 1:nrow(formations)){
    if (formations$Range[f] > quantile(formations$Range, 0.75, type = 7)) {
      if (draw_by == "Lat"){ # If user has set draw_by argument to "Lat"
        fplong <- rbind(fplong, c(formations$min_age[f], formations$Mean_Lat[f], formations$max_age[f], formations$Mean_Lat[f], formations$Diversity[f], formations$Occurrences[f])) # Add latitude and bind other relevant information. 
      }
      if (draw_by == "Max_Age"){
        fplong <- rbind(fplong, c(formations$min_age[f], formations$max_age[f], formations$max_age[f], formations$max_age[f], formations$Diversity[f], formations$Occurrences[f])) # Add Max Age and bind other relevant information. 
      }
      if (draw_by == "Number"){
        fplong <- rbind(fplong, c(formations$min_age[f], formations$forbinning[f], formations$max_age[f], formations$forbinning[f], formations$Diversity[f], formations$Occurrences[f])) # Add Number of formation and bind other relevant information. 
      }
    }
    else{ # Otherwise, do the same things as above but add to dataframe of short formations
      if(draw_by == "Lat"){
        fp1 <- rbind(fp1, c(formations$min_age[f], formations$Mean_Lat[f], formations$max_age[f], formations$Mean_Lat[f], formations$Diversity[f], formations$Occurrences[f]))
      }
      if(draw_by == "Max_Age"){
        fp1 <- rbind(fp1, c(formations$min_age[f], formations$max_age[f], formations$max_age[f], formations$max_age[f], formations$Diversity[f], formations$Occurrences[f]))
      }
      if(draw_by == "Number")
        fp1 <- rbind(fp1, c(formations$min_age[f], formations$forbinning[f], formations$max_age[f], formations$forbinning[f], formations$Diversity[f], formations$Occurrences[f]))
    }
  }
  names <- c("x1", "y1", "x2", "y2", "Diversity", "Occurrences")
  colnames(fp1) <- names
  colnames(fplong) <- names
  fp1$Type <- "Short"
  fplong$Type <- "Long"
  fp <- rbind(fp1, fplong) # Bind dataframes
  layout(matrix(1:1,nrow=1), widths=c(1,1), height = c(1,1)) # Set layout
  if(Col == "Diversity"){ #  If colouring by diversity
    graphcol <- rev(sequential_hcl(10, "BluYl"))
    fp$Col <- graphcol[as.numeric(cut(fp$Diversity,breaks = 10))] # Colour taxic diversity of formations according to those colours
    if (legend == TRUE){
      layout(matrix(1:2,nrow=1), widths=c(0.85,0.15), heights = c(1,1)) #  Draw a legend, if set to TRUE
    }
  }
  if(Col == "Occurrences"){ #  If colouring by diversity
    graphcol <- rev(sequential_hcl(10, "BluYl"))
    fp$Col <- graphcol[as.numeric(cut(fp$Occurrences,breaks = 10))] # Colour taxic diversity of formations according to those colours
    if (legend == TRUE){
      layout(matrix(1:2,nrow=1), widths=c(0.85,0.15), heights = c(1,1)) #  Draw a legend, if set to TRUE
    }
  }
  par(mar = c(4.1, 4.1, 1, 1))
  if (draw_by == "Lat"){ # Sets appropriate y axis for mean latitude of occurrences of each formation
    tsplot(stages, boxes=c("short","system"), # Generates plot using Divdyn package
           xlim=c(stages$bottom[1], stages$top[nrow(stages)]),  ylim=range(fp$y1 - 1, fp$y2 + 1), 
           prop = 0.08, plot.args = list(cex.lab = 2, cex.axis = 2),
           shading=NULL, boxes.col=c("col","systemCol"), labels.args=list(cex=1.8),
           ylab = "Formations by mean occurrence latitude")
  }
  if (draw_by == "Max_Age"){ # Sets appropriate y axis for max age of each formation
    tsplot(stages, boxes=c("short","system"), # Generates plot using Divdyn package
           xlim=c(stages$bottom[1], stages$top[nrow(stages)]),  ylim=range(fp$y1 - 1, fp$y2 + 1), 
           prop = 0.08, plot.args = list(cex.lab = 2, cex.axis = 2),
           shading=NULL, boxes.col=c("col","systemCol"), labels.args=list(cex=1.8),
           ylab = "Formations by Maximum Age")
  }
  if (draw_by == "Number"){ # Sets appropriate y axis for aarbitrary number of each formation
    tsplot(stages, boxes=c("short","system"), # Generates plot using Divdyn package
           xlim=c(stages$bottom[1], stages$top[nrow(stages)]),  ylim=range(fp$y1 - 1, fp$y2 + 1), 
           prop = 0.08, plot.args = list(cex.lab = 2, cex.axis = 2),
           shading=NULL, boxes.col=c("col","systemCol"), labels.args=list(cex=1.8),
           ylab = "Formations (Number order)")
  }
  if(STAGE == TRUE){ # If user has set stage to true
    for(n in 1:nrow(stages)){ # 
      if(((n %% 2) == 0) == TRUE) next
      else {
        if(n == 7){
          if(nrow(binlist) %% 2 == 0){
            next
          }
          else{
            rect(stages$top[max(nrow(stages))], 0, stages$bottom[max(nrow(stages))], # Draw grey boxes on graph as stages
                 max(fp$y2+1, na.rm = TRUE), col = rgb(0.89,0.89,0.89,alpha=0.5), border = NA)
          }
        }
        else{
          rect(stages$top[n], 0, stages$bottom[n], 
               max(fp$y2+1, na.rm = TRUE), col = rgb(0.89,0.89,0.89,alpha=0.5), border = NA)
        }
      }
    }
  }
  else{ # Otherwise, use newly generated bins
    for(n in 1:length(form_bins)){ # draws new bins as coloured boxes for comparison to traditional bins
      if(((n %% 2) == 0) == TRUE) next
      else {
        if(n == length(form_bins)){
          if(nrow(binlist) %% 2 == 0){
            next
          }
          else{
            rect(useful_bins[n], 0, useful_bins[n-1], 
                 max(fp$y2+1, na.rm = TRUE), col = rgb(0.89,0.89,0.89,alpha=0.5), border = NA)
          }
        }
        else{
          rect(form_bins[n], 0, form_bins[n+1], 
               max(fp$y2+1, na.rm = TRUE), col = rgb(0.89,0.89,0.89,alpha=0.5), border = NA)
        }
      }
    }
  }
  if (score_grid_2 == TRUE){
    # segments(fplong$x1, fplong$y1, fplong$x2, fplong$y2, lwd = 2, col = "grey") # Plots formations as lines showing their duration
    if (Col == "Diversity" | Col == "Occurrences"){
      segments(fp1$x1, fp1$y1, fp1$x2, fp1$y2, lwd = 5, col = "grey") 
      segments(fp1$x1, fp1$y1, fp1$x2, fp1$y2, lwd = 3, col = fp$Col[fp$Type == "Short"]) 
      rect(fplong$x1, fplong$y1-0.05, fplong$x2, fplong$y2+0.05, angle = 45, col = "lightgrey", density = 40, border = NA, lwd = 5)
      rect(fplong$x1, fplong$y1-0.05, fplong$x2, fplong$y2+0.05, angle = 45, col = fp$Col[fp$Type == "Long"], density = 40, border = NA, lwd = 3)
    }
    else{
      segments(fp1$x1, fp1$y1, fp1$x2, fp1$y2, lwd = 3) 
      rect(fplong$x1, fplong$y1-0.05, fplong$x2, fplong$y2+0.05, angle = 45, col = "grey", density = 40, border = NA, lwd = 3)
    }
  }
  else{
    if (Col == "Diversity" | Col == "Occurrences"){
      segments(fp$x1, fp$y1, fp$x2, fp$y2, lwd = 5, col = "grey") 
      segments(fp$x1, fp$y1, fp$x2, fp$y2, lwd = 3, col = fp$Col) 
    }
    else{
      segments(fp$x1, fp$y1, fp$x2, fp$y2, lwd = 3) 
    }
  }
  box(lwd=2)
  if (Col == "Diversity" && legend == TRUE){
    par(mar = c(4.1, 0.5, 1, 0))
    legend_image <- as.raster(matrix(rev(graphcol), ncol=1))
    plot(c(0,4),c(0,1),type = 'n', axes = F, xlab = '', ylab = '')
    title("Diversity", cex.main = 1.5, line = -9, adj = 0)
    text(x=1.1, y = seq(0.25,0.75,l=5), cex = 1.25, labels = round(seq(0,max(fp$Diversity),l=5)))
    rasterImage(legend_image, 0, 0.25, 0.5, 0.75)
  }
  if (Col == "Occurrences" && legend == TRUE){
    par(mar = c(4.1, 0.5, 1, 0))
    legend_image <- as.raster(matrix(rev(graphcol), ncol=1))
    plot(c(0,4),c(0,1),type = 'n', axes = F, xlab = '', ylab = '')
    title("Occurrences", cex.main = 1.5, line = -9, adj = 0)
    text(x=1.1, y = seq(0.25,0.75,l=5), cex = 1.25, labels = round(seq(0,max(fp$Occurrences),l=5)))
    rasterImage(legend_image, 0, 0.25, 0.5, 0.75)
  }
  p <<- recordPlot()
}
