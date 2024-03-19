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
  vect <- c(0.5, 1)
  res_data(bin.occs, target, vect = vect, single = TRUE) # Makes list of dataframes, each containing information about various cell resolutions. Can also see Naive occupancy estimates.
  Res_results <- do.call(rbind.data.frame, Res_results_list)
  write.csv(Res_results, file.path(paste("Prepped_data/Occurrence_data/", bin.type, "/", bin.name, 
                                         "/Targeted_res_stats.csv", sep="")))
  
  # Prepare data for unmarked
  all_results_for_unmarked(data = bin.occs, name = bin.name, res = res, ext = e, 
                           target = target, single = TRUE, formCells = form_cells, 
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
