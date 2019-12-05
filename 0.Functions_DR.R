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
library(raster)
library(plyr)
library(dplyr)
library(lattice)
library(rasterVis)
library(sp)
library(maps)
library(maptools)
library(parallel)
library(plyr)

#=============================================== GET_EXTENT =============================================================

# Setup raster for resolution and extent. Note: these values should be the same ones used for the file 1.Setup_occupancy_DR.
get_extent <- function(data){
  maxLat <- round_any((max(data$lat) + 3), 0.5)  #get value for defining extent, and increase by x for visualisation purposes
  minLat <- round_any((min(data$lat) - 3), 0.5)  #get value for defining extent, and increase by x for visualisation purposes
  maxLng <- round_any((max(data$lng) + 3), 0.5)  #get value for defining extent, and increase by x for visualisation purposes
  minLng <- round_any((min(data$lng) - 3), 0.5) #get value for defining extent, and increase by x for visualisation purposes
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
#=============================================== COMBINE_DATA ===========================================================

# Adds together target occurrence dataset with a broader collections dataset (allowing maximum possible sampling opportunities) if datasets were not downloaded together
combine_data <- function(fossils, Add_fossils){ # fossils is occurrence dataset, Add_fossils is collections dataset
  d1 <- fossils %>%
    dplyr::select(collection_name, lat, lng, Target, Carb_Clast, Lith) # Target specifies what group you want to look at. Can also add other data in here (lithology etc)
  d2 <- Add_fossils %>%
    dplyr::select(collection_name, lat, lng, Target, Carb_Clast, Lith)
  all_colls <- dplyr::bind_rows(d1, d2) 
  just_colls <- all_colls %>%
    dplyr::select(collection_name, lat, lng) %>%
    dplyr::distinct()
  temp_name <- paste(deparse(substitute(fossils)),".", "comb", sep = "") #Name files based on data entered to function
  assign(temp_name, all_colls, envir = .GlobalEnv)
  temp_name <- paste(deparse(substitute(Add_fossils)),".", "comb", sep = "") #Name files based on data entered to function
  assign(temp_name, just_colls, envir = .GlobalEnv)
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
get_grid <- function(data, res, e){ # data is first output from combine_data (fossil.colls). Res is chosen resolution in degrees
  val <- 64800/res
  r <- raster(res = res, val = val, ext = e)
  xy <- SpatialPointsDataFrame(cbind.data.frame(data$lng, data$lat), data)
  Final <- extract(r, xy, sp = TRUE, cellnumbers = TRUE)
  Final <<- as.data.frame(Final)
}

#=============================================== GET_COV FUNCTIONS ===========================================================

# Functions to organise covariate data

#===== GET_COV =====

# Attaches grid cell IDs from an inputted raster to occurrences/collections.
get_cov <- function(data, raster){ # data is first output from combine_data (fossil.colls). Raster is a chosen raster file, which can be a raster stack. 
  xy <- SpatialPointsDataFrame(cbind.data.frame(data$lng, data$lat), data)
  cov_dat <- extract(raster, xy, sp = TRUE, cellnumbers = TRUE)
  cov_dat <<- as.data.frame(cov_dat)
}

#===== GET_COV_FROM_STACK =====

# Creates a raster stack of chosen resolution from previously stored raster files, and attaches associated grid cell IDs to occurrences/collections. 
get_cov_from_stack <- function(data, res){ # data is first output from combine_data (fossil.colls) and chosen resolution. 
  grids <- list.files(paste("Data/Covariate_Data/Formatted/All_data/", res, "deg/", sep = ""), pattern = "asc") # Find rasters of appropriate resolution
  raster <- raster::stack(paste0("Data/Covariate_Data/Formatted/All_data/", res, "deg/", grids)) # stack those rasters
  get_cov(data, raster)
  CovStack <<- raster
}

#=============================================== GET_GRID_IM ===========================================================

# Sets raster to dimensions of inputted data ready for visualisation. Is used in vis_grid. 
get_grid_im <- function(data, res, name, ext){ # Data is first output from combine_data (fossil.colls). Res is chosen resolution in degrees. name is user inputted string related to data inputted, for display on graphs. 
  xy <- cbind(as.double(data$lng), as.double(data$lat))
  r <- raster::raster(ext = ext, res = res)
  r <- raster::rasterize(xy, r, fun = 'count')
  #r[r > 0] <- 1 # Remove if you want values instead of pure presence/absence.
  countries <- maps::map("world", plot=FALSE, fill = TRUE) # find map to use as backdrop
  countries <<- maptools::map2SpatialPolygons(countries, IDs = countries$names, proj4string = CRS("+proj=longlat")) # Turn map into spatialpolygons
  mapTheme <- rasterVis::rasterTheme(region=brewer.pal(8,"Greens"))
  print(rasterVis::levelplot(r, margin=FALSE, par.settings=mapTheme,  main = paste("Total ", (substitute(name)), " per Grid Cell", sep = "")) + #create levelplot for raster
          latticeExtra::layer(sp.polygons(countries, col = 0, fill = "grey")) + #add countries
          rasterVis::levelplot(r, margin=FALSE, par.settings=mapTheme)) #add levelplot again over the top (messy, but unsure how else to set this up!)
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
all_results_for_unmarked <- function(data, res, target, ext, subsamp = TRUE, single = TRUE){ # data is first output from combined_data (fossil.colls). res is vectors of chosen resolutions. target is vector of chosen targets. 
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
        test2 <- SubSamp_for_unmarked(test2, target[t])
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
