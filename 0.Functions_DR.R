#=============================================== FUNCTIONS FOR OCCUPANCY MODELLING ===========================================================

#      Selection of functions that work with PBDB data in order to set up occupancy models. Functions range from those that
#      visualise occurrences in terms of grid cells to those that reorder presence/absence data per grid cell so it fits 
#      the format of the package unmarked. Information regarding each function can be found in the seperate sections below. 

#=============================================== iPAK AND REQUIRED PACKAGES =====================================================

# fucntion that automatically installs necessary packages that the user is lacking.

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
library(dplyr)
library(lattice)
library(rasterVis)
library(sp)
library(maps)
library(maptools)
library(parallel)

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

# Adds together target occurrence dataset with a broader collections dataset (allowing maximum possible sampling opportunities)
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
get_grid <- function(data, res){ # data is first output from combine_data (fossil.colls). Res is chosen resolution in degrees
  val <- 64800/res
  r <- raster(res = res, val = val)
  xy <- SpatialPointsDataFrame(cbind.data.frame(data$lng, data$lat), data)
  Final <- extract(r, xy, sp = TRUE, cellnumbers = TRUE)
  Final <<- as.data.frame(Final)
}

#=============================================== GET_GRID_IM ===========================================================

# Sets raster to dimensions of inputted data ready for visualisation. Is used in vis_grid. 
get_grid_im <- function(data, res, name){ # Data is first output from combine_data (fossil.colls). Res is chosen resolution in degrees. name is user inputted string related to data inputted, for display on graphs. 
  min_lat <- floor(min(data$lat, na.rm = TRUE)/10)*10 
  max_lat <- ceiling(max(data$lat, na.rm = TRUE)/10)*10
  min_lon <- floor(min(data$lng, na.rm = TRUE)/10)*10
  max_lon <- ceiling(max(data$lng, na.rm = TRUE)/10)*10
  Lat_dim <- (max_lat - min_lat)/res
  Lon_dim <- (max_lon - min_lon)/res
  xy <- cbind(as.double(data$lng), as.double(data$lat))
  r <- raster::raster(ncol=Lon_dim, nrow=Lat_dim, xmn=min_lon, xmx=max_lon, ymn=min_lat, ymx=max_lat)
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
       xlab = "Number of Occurrences", ylab = "Number of Grid Cells",
       col = "springgreen")
  r <<- r
}

#=============================================== PREPARE_FOR_RES_DATA AND RES DATA ===========================================================

# Functions to test quality of data at different resolutions of grid cells.

#===== PREPARE_FOR_RES_DATA =====
# Produces summary of key stats for data at a specified resolution of grid cell. Used in res_data.
prepare_for_res_data <- function(data, target){ # Data is first output from combine_data (fossil.colls). Res is chosen resolution in degrees
  rel_data <- data %>% 
    dplyr::select(collection_name, cells, Target) # Select appropriate cells
  rel_data$Target[rel_data$Target != target] <- 0 # Make anything that's not the target group a 0
  rel_data$Target[rel_data$Target == target] <- 1 # Make target's a 1
  coll_data <- rel_data %>% # For all collections, give a mean score of presences and absences
    dplyr::group_by(collection_name) %>%
    dplyr::summarize(mean(Target)) 
  coll_data$`mean(Target)` <- ceiling(coll_data$`mean(Target)`) #anything above a 0 has presences, therefore can be counted as 1
  joined_data <- dplyr::left_join(coll_data, rel_data) %>% #Join with cell IDs
    dplyr::select(-Target, Pres.Abs = `mean(Target)`, collection_name) #Remove old target, clean column name
  joined_data <- joined_data %>% dplyr::distinct() # Remove duplicates of remaining collections
  id.table <- table(joined_data$cells) # Create table for removing singleton cells (cells with only one collection)
  joined_data <- subset(joined_data, cells %in% names(id.table[id.table > 1])) # Remove cells with only one collection
  temp_name <- paste("pres.abs.data.",deparse(substitute(target)), sep = "") #Name files based on data entered to function
  assign(temp_name, joined_data, envir = .GlobalEnv)
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
  results <<- results
}

#===== RES_DATA =====
# Carries out prepare_for_res_data over a sequence of resolutions, and outputs as a data.frame.
res_data <- function(data, target, start, fin, int){ # Data is first output from combine_data (fossil.colls). Target is chosen group to test. start is first resolution, fin is last resolution,int is the interval to count between start and fin.
  s1 <- seq(start, fin, int)
  Res_results <- data.frame(matrix(ncol = 7, nrow = length(s1)))
  colnames(Res_results) <- c("No.Cells", "Naive.occ", "Total.Colls", 
                             "Mean.Colls", "Min.Colls", 
                             "Max.Colls", "Median.Colls")
  row.names(Res_results) <- s1
  for (i in (1:length(s1))){
    test <- get_grid(data, s1[i])
    prepare_for_res_data(test, target)
    Res_results[i,]<- results
  }
  temp_name <- paste("Res_results.",deparse(substitute(target)), sep = "") #Name files based on data entered to function
  assign(temp_name, Res_results, envir = .GlobalEnv)
}

#=============================================== SETUP_CELL_INFO ===========================================================

# Creates information about each grid cell for use as covariate in occupancy modelling. Only works if you are including lith info, and will need to be adjusted for broader use. 
setup_cell_info <- function(data, res){ # Data is output from Get_grid. Res is resolution (only neccessary for later functions)
  Collections_per_cell <- data %>% # Counting collections per cell
    dplyr::select(collection_name, cells) %>%
    dplyr::distinct() %>%
    dplyr::group_by(cells) %>%
    dplyr::summarize(colls_per_cell = n())
  Collections_per_cell <- Collections_per_cell[!Collections_per_cell$colls_per_cell == 1,] # Removing any cells with < 1 collection
  Occurrences_per_cell <- data %>% # Counting occurrences per cell
    dplyr::select(collection_name, cells) %>%
    dplyr::group_by(cells) %>%
    dplyr::summarize(occs_per_cell = n())
  cell_data <- Collections_per_cell %>% left_join(Occurrences_per_cell) # Joining collections and occurrences for cells
  Carb_Perc <- data %>% #calculating percentage of carbonate collections per cell
    dplyr::select(collection_name, Carb_Clast, cells) %>%
    dplyr::group_by(cells) %>%
    dplyr::summarize(perc_carb = mean(Carb_Clast == "Carbonate"))
  cell_data <- dplyr::left_join(cell_data, Carb_Perc) # Joining to cell data
  cell_data$perc_carb <- (round(cell_data$perc_carb, 2)) * 100 # rounding to 2 dp.
  cell_data <<- cell_data
}

#=============================================== PREPARE_FOR_UNMARKED AND ALL_RESULTS_FOR_UNMARKED ===========================================================

# Functions for converting data into the correct format for unmarked (occupancy modelling package). Can be run individually or for multiple Targets and Resolutions. 

#===== PREPARE_FOR_UNMARKED =====
# Converts data generated by Get_grid into the correct format for unmarked. 

prepare_for_unmarked <- function(data, target){ # data is output from Get_Grid. target is specified group to examine. 
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
  id.table <- table(joined_data$cells) # Create table for removing singleton cells (cells with only one collection)
  joined_data <- subset(joined_data, cells %in% names(id.table[id.table > 1])) # Remove cells with only one collection
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
}

#===== SUBSAMP_FOR_UNMARKED =====
# Takes data prepared for unmarked, and standardizes it down to a total of five site visits (collections) for each gridsquare through subsampling.  
SubSamp_for_unmarked <- function(data){ # Data is output from prepare_for_unmarked. Outputs same data, but subsampled to five site visits. 
  new_dframe_for_unmarked <- data.frame()
  for (n in 1:nrow(data)){
    if(rowSums(is.na(data[n,])) > (NCOL(data)-6)){
      new_dframe_for_unmarked <- rbind(new_dframe_for_unmarked, data[n,])
    }
    else{
      temp_dframe <- data.frame(1:5)
      for (t in 1:100){
        temp_vec <- data[n,]
        temp_vec <- temp_vec[!is.na(temp_vec)]
        temp_dframe <- cbind(temp_dframe, sample(temp_vec, 5, replace = FALSE))
      }
      temp_dframe <- temp_dframe[2:101]
      num <- round(mean(colSums(temp_dframe)))
      samp <- sample(c(rep(1, num), rep(0, (5-num))))
      new_dframe_for_unmarked <- rbind(new_dframe_for_unmarked, samp)
      rownames(new_dframe_for_unmarked)[n] <- rownames(data)[n]
    }
  }
  new_dframe_for_unmarked <- new_dframe_for_unmarked[,1:5]
  SS_dframe_for_unmarked <<- new_dframe_for_unmarked
}

#===== PREPARE_COLL_COVS =====
# Creates covariated data using Get_grid, and transforms into the correct format for unmarked. 
prepare_coll_covs <- function(data){ # data is output from Get_Grid. target is specified group to examine. 
  rel_data <- data %>% 
    dplyr::select(collection_name, cells, Lith, Target) %>%
    dplyr::arrange(cells, collection_name) # Sort cells so they match pres/abs targeted data
  rel_data_dist <- rel_data %>%
    dplyr::select(collection_name, cells, Lith) %>%
    dplyr::distinct() # Remove duplicates of remaining collections
  
  if(length(rel_data_dist$collection_name[duplicated(rel_data_dist$collection_name)]) > 0){ # Produces an error if collections contain multiple lithologies
    duplicate_colls <- length(rel_data_dist$collection_name[duplicated(rel_data_dist$collection_name)])
    print(rel_data_dist$collection_name[duplicated(rel_data_dist$collection_name)])                      
    stop(paste("the above ", duplicate_colls, " listed collections have different lithology records. Revise and try again."))
  }
  
  id.table <- table(rel_data_dist$cells) # Create table for removing singleton cells (cells with only one collection)
  rel_data_dist <- subset(rel_data_dist, cells %in% names(id.table[id.table > 1])) # Remove cells with only one collection
  results <- c(length(unique(rel_data_dist$cells)), # number of cells
               nrow(unique(rel_data_dist)), # total number of collections
               mean(table(rel_data_dist$cells)), # mean number of collections in each cell
               min(table(rel_data_dist$cells)), # min number of collections in each cell
               max(table(rel_data_dist$cells)), # max number of collections in each cell
               median(table(rel_data_dist$cells))) # median number of collections in each cell
  
  dframe_for_unmarked <- data.frame(matrix(ncol = results[5], nrow = results[1])) # Making dataframe for unmarked data
  colnames(dframe_for_unmarked) <- c(sprintf("y.%d", seq(1,results[5])))
  row.names(dframe_for_unmarked) <- unique(rel_data_dist$cells)
  test <- unique(rel_data_dist$cells)
  for (g in 1:results[1]){
    counter <- 1
    for (r in 1:nrow(rel_data_dist)){
      if (rel_data_dist[r,2] == test[g]){
        dframe_for_unmarked[g, counter] <- rel_data_dist[r, 3]
        counter <- counter + 1
      }
    }
  }
  site_covs <- dframe_for_unmarked
}

#===== ALL_RESULTS_FOR_UNMARKED =====
# loop that take basic combined data and writes multiple .csv files into Results folder in current directory for chosen grid cells resolutions and targets in correct format for unmarked. Sound rings when function has finished running.
all_results_for_unmarked <- function(data, res, target, subsamp = FALSE){ # data is first output from combined_data (fossil.colls). res is vectors of chosen resolutions. target is vector of chosen targets. 
  for (r in 1:length(res)){
    ptm <- proc.time()
    test1 <- get_grid(data, res[r])
    for (t in 1:length(target)){
      test2 <- prepare_for_unmarked(test1, target[t])
      if (subsamp==TRUE){
        test2 <- SubSamp_for_unmarked(test2)
      }
      temp_name <- paste(deparse(substitute(data)), ".", res[r], ".", target[t], sep = "")
      write.csv(test2, file.path(paste("Results/", temp_name, ".csv", sep="")))
    }
    proc.time() - ptm
  }
  beepr::beep(sound = 3)
}

#===== ALL_COVS_INFO =====
# Creates data for covariates used in occupancy modelling. Includes both site (cell) and observation (collection) specific covariates
all_covs_info <- function(data, res, outcrop){
  for (r in 1:length(res)){
    test <- get_grid(data, res[r])
    cell_inf <- setup_cell_info(test)
    val <- 64800/res[r]
    ra <- raster(res = res[r], val = val)
    system.time(outcrop_perc <- rasterize(outcrop, ra, getCover = TRUE)) # Parallelize rasterize function
    xy <- SpatialPointsDataFrame(cbind.data.frame(test$lng, test$lat), test)
    rel_outcrop_perc <- extract(outcrop_perc, xy, cellnumbers = TRUE)
    new_cell_inf <- rel_outcrop_perc %>%
      as.data.frame() %>%
      distinct() %>%
      inner_join(cell_inf)
    colnames(new_cell_inf)[2] <- "Perc_outcrop_area"
    site_covs <- prepare_coll_covs(test)
    dir.create(paste0("Results"), showWarnings = FALSE) #stops warnings if folder already exists
    write.csv(new_cell_inf, file.path(paste("Results/cellcovs.", res[r], ".csv", sep="")))
    write.csv(site_covs, file.path(paste("Results/sitecovs.", res[r], ".csv", sep="")))
  }
  beepr::beep(sound = 3)
}

# Below code currently not working - need to understand clustering process and interaction with rasterize function.
#
# all_covs_info_clust <- function(data, res, outcrop){
#   no_cores <- detectCores() - 2 # Set number of cores
#   cl <- makeCluster(no_cores) # Initiate cluster
#   for (r in 1:length(res)){
#     test <- get_grid(data, res[r])
#     cell_inf <- setup_cell_info(test)
#     val <- 64800/res[r]
#     ra <- raster(res = res[r], val = val)
#     
#     # Parallel Processing 
#     features <- 1:nrow(outcrop[,]) # Setup features
#     n <- 50
#     parts <- split(features, cut(features, n)) # Split features for parallel processing
#     z <- clusterEvalQ(cl, library("raster")) #setup necessary packages for cluster
#     clusterExport(cl, c("ra", "outcrop", "parts", "n"), envir=environment()) #setup variables for cluster
#     system.time(rParts <- parLapply(cl = cl, X = 1:n, fun = function(x) rasterize(outcrop[parts[[x]],], ra, getCover = TRUE))) # Parallelize rasterize function
#     
#     outcrop_perc <- do.call(merge, rParts) # Merge all raster parts
#     xy <- SpatialPointsDataFrame(cbind.data.frame(test$lng, test$lat), test)
#     rel_outcrop_perc <- extract(outcrop_perc, xy, cellnumbers = TRUE)
#     new_cell_inf <- rel_outcrop_perc %>%
#       as.data.frame() %>%
#       distinct() %>%
#       inner_join(cell_inf)
#     colnames(new_cell_inf)[2] <- "Perc_outcrop_area"
#     site_covs <- prepare_coll_covs(test)
#     dir.create(paste0("Results"), showWarnings = FALSE) #stops warnings if folder already exists
#     write.csv(new_cell_inf, file.path(paste("Results/cellcovs.", res[r], ".csv", sep="")))
#     write.csv(site_covs, file.path(paste("Results/sitecovs.", res[r], ".csv", sep="")))
#   }
#   beepr::beep(sound = 3)
#   stopCluster(cl) # Finish Parallelizing
# }