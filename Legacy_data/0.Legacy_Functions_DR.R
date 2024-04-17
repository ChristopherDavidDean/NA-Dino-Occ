
#================== iPAK AND REQUIRED PACKAGES =================================

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

#====== GET_EXTENT =============================================================

# Setup raster for resolution and extent. Note: these values should be the same ones used for the file 1.Setup_occupancy_DR.
get_extent <- function(data){
  maxLat <- round_any((max(data$lat) + 3), 0.5)  #get value for defining extent, and increase by x for visualisation purposes
  minLat <- round_any((min(data$lat) - 3), 0.5)  #get value for defining extent, and increase by x for visualisation purposes
  maxLng <- round_any((max(data$lng) + 3), 0.5)  #get value for defining extent, and increase by x for visualisation purposes
  minLng <- round_any((min(data$lng) - 3), 0.5) #get value for defining extent, and increase by x for visualisation purposes
  e <<- extent(minLng, maxLng, minLat, maxLat) # build extent object
}

#====== TARGET_MAKER ===========================================================

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
#====== COMBINE_DATA ===========================================================

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

#========== BIN_SPLITTER =======================================================

# Takes combined data and splits it into user defined bins based off a vector. 
bin_splitter <- function(bins, fossils){ # takes first output from combined_data and a vector of time bins. 
  for (s in 1:(length(bins)-1)){
    temp_data <- fossils %>%
      dplyr::filter(mid_ma < bins[s] & mid_ma > bins[s + 1])
    temp_name <- paste(deparse(substitute(fossils)), ".bin", s, sep = "") #Name files based on data entered to function
    assign(temp_name, temp_data, envir = .GlobalEnv)
  }
}

#===== GET_GRID ================================================================

# Creates a raster of chosen resolution, and attaches associated grid cell IDs to occurrences/collections
get_grid <- function(data, res, e){ # data is first output from combine_data (fossil.colls). Res is chosen resolution in degrees
  val <- 64800/res
  r <- raster(res = res, val = val, ext = e)
  xy <- SpatialPointsDataFrame(cbind.data.frame(data$lng, data$lat), data)
  Final <- extract(r, xy, sp = TRUE, cellnumbers = TRUE)
  Final <<- as.data.frame(Final)
}

#======== GET_COV FUNCTIONS ====================================================

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
  grids <- list.files(paste("Lewis_Occupancy_data/Data/Formatted/All_data/", res, "deg/", sep = ""), pattern = "asc") # Find rasters of appropriate resolution
  raster <- raster::stack(paste0("Lewis_Occupancy_data/Data/Formatted/All_data/", res, "deg/", grids)) # stack those rasters
  get_cov(data, raster)
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
       xlab = "Number of Collections", ylab = "Number of Grid Cells",
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
all_results_for_unmarked <- function(data, res, target, ext, subsamp = TRUE){ # data is first output from combined_data (fossil.colls). res is vectors of chosen resolutions. target is vector of chosen targets. 
  for (r in 1:length(res)){
    ptm <- proc.time()
    test1 <- get_grid(data, res[r], ext)
    for (t in 1:length(target)){
      test2 <- prepare_for_unmarked(test1, target[t])
      if (subsamp==TRUE){
        test2 <- SubSamp_for_unmarked(test2, target[t])
      }
      temp_name <- paste(deparse(substitute(data)), ".", res[r], ".", target[t], sep = "")
      dir.create(paste0("Results"), showWarnings = FALSE) #stops warnings if folder already exists
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


###############################
##### PLOT.NAIVE.UNMARKED #####
###############################

plot.naive.unmarked <- function(res.comb){
  naive <- res.comb %>%
    filter(Data == "Naive Occupancy" | Data ==  "PAO")
  ggplot(data = naive, aes(x = new_bins, y = value, color = Data)) +
    ylab("Proportion of total sites") + 
    xlab("Time (Ma)") +
    scale_x_reverse() +
    deeptime::coord_geo(dat = list("stages"), 
                        xlim = c((max(res.comb$new_bins)+1), (min(res.comb$new_bins-1))), 
                        ylim = c(0, 1)) +
    geom_smooth(method=lm, aes(group = Data), colour = "#3182BD", 
                alpha = 0.2, linewidth = 0.75) +
    geom_line(aes(x = new_bins, y = value, color = Data)) +
    scale_color_manual(values=c("#252424", "#DE2D26")) +
    scale_fill_manual(values=c("#252424", "#DE2D26"))  +
    theme_few() +
    theme(legend.position="none")
}

################################################################################
# 9. PREPARE_FOR_UNMARKED, SAMPLE_FOR_UNMARKED AND ALL_RESULTS_FOR_UNMARKED
################################################################################

# Functions for converting data into the correct format for unmarked (occupancy 
# modelling package). Can be run individually or for multiple Targets and Resolutions. 

################################
##### PREPARE_FOR_UNMARKED #####
################################

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

####################################
##### ALL_RESULTS_FOR_UNMARKED #####
####################################

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
    rand <- round(runif(1, min = 0, max = 999)) # Set random number for seed
    for (q in 1:length(target)){
      if(single == FALSE){
        test2 <- prepare_for_unmarked(test1, target[q], single = FALSE)
        if(max_val_on == TRUE){
          set.seed(rand) # set seed to ensure all targets are have same sampling
          test2 <- sample_for_unmarked(test2, max_val)
          temp_name <- paste("unmarked_", target[q], sep = "")
          assign(temp_name, test2, envir = .GlobalEnv)
        }
      }
      else {
        test2 <- prepare_for_unmarked(test1, target[q])
        if(max_val_on == TRUE){
          set.seed(rand) # set seed to ensure all targets have same sampling
          test2 <- sample_for_unmarked(test2, max_val)
          temp_name <- paste("unmarked_", target[q], sep = "")
          assign(temp_name, test2, envir = .GlobalEnv)
        }
      }
      # Create folders, remove warning if they already exist.
      dir.create(paste0("Prepped_data/Occurrence_Data/", bin.type, "/", sep = ""), showWarnings = FALSE) 
      dir.create(paste0("Prepped_data/Occurrence_Data/", bin.type, "/", bin.name, "/", sep =""), 
                 showWarnings = FALSE) 
      dir.create(paste0("Prepped_data/Occurrence_Data/", bin.type, "/", bin.name, "/", res, "/", 
                        sep = ""), showWarnings = FALSE) 
      
      if(max_val_on == TRUE){
        if(form_cells == "Y"){
          temp_name_1 <- paste(name, ".", res[r], ".", target[q], ".dframe.",  
                               max_val, ".formcells", sep = "")
          temp_name_2 <- paste(name, ".", res[r], ".", target[q], ".colframe.", 
                               max_val, ".formcells", sep = "")
        }else{
          temp_name_1 <- paste(name, ".", res[r], ".", target[q], ".dframe.",  max_val, sep = "")
          temp_name_2 <- paste(name, ".", res[r], ".", target[q], ".colframe.", max_val, sep = "")
        }
      }else{
        if(form_cells == "Y"){
          temp_name_1 <- paste(name, ".", res[r], ".", target[q], ".dframe.formcells",  sep = "")
          temp_name_2 <- paste(name, ".", res[r], ".", target[q], ".colframe.formcells", sep = "")
        } else{
          temp_name_1 <- paste(name, ".", res[r], ".", target[q], ".dframe",  sep = "")
          temp_name_2 <- paste(name, ".", res[r], ".", target[q], ".colframe", sep = "")
        }
      }
      
      write.csv(test2[[1]], file.path(paste("Prepped_data/Occurrence_Data/", bin.type, "/", bin.name, "/", 
                                            res, "/", temp_name_1, ".csv", sep="")))
      write.csv(test2[[2]], file.path(paste("Prepped_data/Occurrence_Data/", bin.type, "/", bin.name, "/", 
                                            res, "/", temp_name_2, ".csv", sep="")))
    }
    proc.time() - ptm
    assign("bin.occs", test1, envir = .GlobalEnv)
  }
  beepr::beep(sound = 3)
}

#######################
##### PRECISE_COV #####
#######################

# Creates dataframe of covariate data associated with relevant grid cells, taken
# from original hi-resolution rasters. Covariate values are created from the mean 
# value of collections within larger grid cells of chosen resolution. 

precise_cov <- function(data, samp.data, max_val){
  # ADD INFO HERE
  
  wc <- list.files("Prepped_data/Covariate_Data/Precise/")
  wc <- wc[!grepl('0.*', wc)] # Remove folders based on resolution
  wc <- stack(paste0("Data/Covariate_Data/Formatted_For_Precise/", wc, sep = ""))
  projection(wc) <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs" 
  xy <- SpatialPointsDataFrame(cbind.data.frame(data$lng, data$lat), data, 
                               proj4string = CRS("+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs"))
  hires_cov_dat <- raster::extract(wc, xy, sp = TRUE, cellnumbers = FALSE)
  data <- as.data.frame(hires_cov_dat)
  
  counting_colls <- data %>%
    dplyr::select(siteID, collection_no) %>%
    dplyr::distinct() %>%
    dplyr::group_by(siteID) %>%
    dplyr::summarize(Coll_count = n())
  hires_cov_dat <- data %>%
    dplyr::group_by(siteID) %>%
    dplyr::summarize(mean_DEM = mean(DEM, na.rm = TRUE), 
                     mean_prec = mean(prec, na.rm = TRUE),
                     mean_temp = mean(temp, na.rm = TRUE))
  hires_cov_dat <- cbind(hires_cov_dat, counting_colls$Coll_count)
  
  if(is.numeric(max_val) == T){
    surv.data <- data.frame(siteID = rep(rownames(samp.data), each = max_val),
                            colls = matrix(t(samp.data), ncol=1, nrow=ncol(samp.data)*nrow(samp.data), byrow=F))
  }else{
    surv.data <- data.frame(siteID = rep(rownames(samp.data), each = ncol(samp.data)),
                            colls = matrix(t(samp.data), ncol=1, nrow=ncol(samp.data)*nrow(samp.data), byrow=F))
  }
  
  for(r in 1:nrow(surv.data)){
    tem <- which(surv.data$colls[r] == data$collection_no)
    if(length(tem) == 0){
      surv.data$temp[r] <- NA
      surv.data$prec[r] <- NA
      surv.data$DEM[r] <- NA
    }else{
      tem
      surv.data$temp[r] <- data$temp[[tem[1]]]
      surv.data$prec[r] <- data$prec[[tem[1]]]
      surv.data$DEM[r] <- data$DEM[[tem[1]]]
    }
  }
  write.csv(hires_cov_dat, file.path(paste("Prepped_data/Occurrence_Data/", bin.type, "/", bin.name, "/", 
                                           res, "/precise_mean_covs_", max_val, ".csv", sep="")))
  write.csv(surv.data, file.path(paste("Prepped_data/Occurrence_Data/", bin.type, "/", bin.name, "/", 
                                       res, "/surv_covs_", max_val, ".csv", sep="")))
}


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
  dir.create(paste0("Prepped_data/", res, sep =""), showWarnings = FALSE)
  write.csv(species.site.final, file.path(paste("Prepped_data/", res, "/", temp_name, 
                                                ".csv", sep="")))
  species.site.final <<- species.site.final
  target.cov <<- target_keep$Code
}

#############################
##### PLOT.OCC.UNMARKED #####
#############################

plot.occ.unmarked <- function(res.comb){
  res.comb <- res.comb %>%
    dplyr::filter(Data != "Null.occ.prob") %>%
    dplyr::filter(Data != "Null.det.prob")
  res.comb[res.comb$Data =="PAO",]["lower95CI"] <- NA
  res.comb[res.comb$Data =="PAO",]["upper95CI"] <- NA
  ggplot(data = subset(res.comb, Data == "Occupancy Probability" | Data == "Detection Probability"), aes(x = new_bins, 
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
    geom_line(data = subset(res.comb, Data == "Occupancy Probability" | Data == "Detection Probability"), 
              aes(x = new_bins, y = value, colour = Data)) +
    scale_color_manual(breaks = c("Naïve Occupancy", "PAO", "Occupancy Probability", "Detection Probability"),
                       values=c("#252424", "#DE2D26", "#DE2D26", "#3182BD")) +
    scale_fill_manual(breaks = c("Naïve Occupancy", "PAO", "Occupancy Probability", "Detection Probability"), 
                      values=c("#FFFFFF", "white","#DE2D26", "#3182BD")) +
    # geom_smooth(method=lm) +
    theme_few()
}


################################################################################
# 19. GET.RESULTS
################################################################################

get.results <- function(target){
  results.list <- c()
  for(t in bins$Bin) {
    if(file.exists(paste("Results/Unmarked/", bin.type, "/", t, "/", res, "/", 
                         target, ".combined.results.", res, ".", t, ".", samp_val,".csv", sep ="")) == T){
      results.list <- c(results.list, paste("Results/Unmarked/", bin.type, "/", t, "/", res, "/", 
                                            target, ".combined.results.", res, ".", t,".", samp_val, ".csv", sep =""))
    }
  }
  temp <- do.call(rbind,lapply(results.list,read.csv))
  s.bins <- bins %>%
    dplyr::select(Bin, mid_ma)
  temp <- merge(temp, s.bins)
  temp <- temp %>%
    dplyr::select(-X)
  
  temp[temp$Parameter =="Occ.prob" | temp$Parameter == "Det.prob",]["X2.5."] <- 
    temp[temp$Parameter =="Occ.prob"| temp$Parameter == "Det.prob",]["Estimate"]-
    temp[temp$Parameter =="Occ.prob"| temp$Parameter == "Det.prob",]["SE"]*1.959964
  temp[temp$Parameter =="Occ.prob" | temp$Parameter == "Det.prob",]["X97.5."] <- 
    temp[temp$Parameter =="Occ.prob"| temp$Parameter == "Det.prob",]["Estimate"]+
    temp[temp$Parameter =="Occ.prob"| temp$Parameter == "Det.prob",]["SE"]*1.959964
  
  temp[temp$Parameter =="Occ.prob",]["Parameter"] <- "Occupancy Probability"
  temp[temp$Parameter =="Det.prob",]["Parameter"] <- "Detection Probability"
  
  names(temp)[names(temp) == "X2.5."] <- "lower95CI"
  names(temp)[names(temp) == "X97.5."] <- "upper95CI"
  names(temp)[names(temp) == "Parameter"] <- "Data"
  names(temp)[names(temp) == "Estimate"] <- "value"
  names(temp)[names(temp) == "mid_ma"] <- "new_bins"
  assign(target, temp, envir = .GlobalEnv)
}