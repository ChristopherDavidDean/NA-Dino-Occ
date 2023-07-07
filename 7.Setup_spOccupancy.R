################################################################################
############ OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS #############
################################################################################

# Christopher D. Dean, Lewis A. Jones, Alfio A. Chiarenza, Sinéad Lyster, 
# Alex Farnsworth, Philip D. Mannion, Richard J. Butler.
# 2023
# Script written by Christopher D. Dean and Lewis A. Jones

################################################################################
#                FILE XXX: DYNAMIC OCCUPANCY WITH SPOCCUPANCY                  #
################################################################################

################################################################################
# 1. INITIAL SETUP
################################################################################

library(spOccupancy)
library(stars)
library(ggplot2)
library(abind)
library(stringr)

##### Load in Functions #####
source("0.Functions.R") # Import functions from other R file (must be in same working directory)

##### Set values #####
# Set resolution
res <- 1
# Set extent
e <- extent(-155, -72, 22.5, 73)
# Set max limit value
max_val <- 30
max_val_on <- TRUE
bin.type <- "scotese"
target <- "Ceratopsidae"
bins <- read.csv("Data/Occurrences/scotesebins.csv")
bins$bin <- bins$code

master.occs.binned.targeted <- read.csv(paste("Prepped_data/Occurrence_Data/", bin.type, 
                                              "/", res, "_", bin.type, "_occurrence_dataset.csv", sep = ""))

rand <- round(runif(1, min = 0, max = 999)) # Set random number for seed

p_rotate <- function(res, e, site_IDs, bins){
  r <- raster(res = res, ext = e)
  siteCoords <- as.data.frame(xyFromCell(r, site_IDs))
  siteCoords$siteID <- site_IDs
  colnames(siteCoords) <- c("lng", "lat", "siteID")
  siteCoords <<- siteCoords
  p_rotate_list <- list()
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

format_occ_covs <- function(p_cov_list){
  temp <- data.frame(siteNo = rep(1:nrow(p_cov_list[[1]])))
  precip <- data.frame(siteNo = rep(1:nrow(p_cov_list[[1]])))
  DEM <- data.frame(siteNo = rep(1:nrow(p_cov_list[[1]])))
  
  for(t in rev(names(p_cov_list))){
    names(p_cov_list[[t]]) <- gsub(".[[:digit:]]", "", names(p_cov_list[[t]]))
    temp_DEM <- p_cov_list[[t]] %>%
      dplyr::select(DEM)
    colnames(temp_DEM) <- t
    DEM <- cbind(DEM, temp_DEM)
    temp_Precip <- p_cov_list[[t]] %>%
      dplyr::select(precip)
    colnames(temp_Precip) <- t
    precip <- cbind(precip, temp_Precip)
    temp_Temp <- p_cov_list[[t]] %>%
      dplyr::select(temp)
    colnames(temp_Temp) <- t
    temp <- cbind(temp, temp_Temp)
  }
  DEM <- DEM %>%
    dplyr::select(-siteNo)
  temp <- temp %>%
    dplyr::select(-siteNo)
  precip <- precip %>%
    dplyr::select(-siteNo)
  
  Year <- data.frame(teyeq = rep(1, length(site_IDs)), 
                    teyep = rep(2, length(site_IDs)), 
                    teyeo = rep(3, length(site_IDs)),
                    teyen = rep(4, length(site_IDs)))
  site.effect <- c(1:length(site_IDs))
  occ_covs <- list(DEM, temp, precip, Year, site.effect)
  names(occ_covs) <- c("DEM", "temp", "precip", "Year", "site.effect")
  occ_covs <<- occ_covs
}

extract_p <- function(p_rotate_list){
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
  }
  full_p_covs <- full_p_covs
}

prepare_for_spOcc <- function(data){ 
  # Select relevant info
  rel_data <- data %>% 
    dplyr::select(collection_no, siteID, Target, bin_assignment)
  
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
    id.table <- table(temp_joined_data$siteID) 
    temp_joined_data <- subset(temp_joined_data, siteID %in% names(id.table[id.table > 1])) 
    
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
    eh_list_samp[[s]] <- sample_for_unmarked(eh_list[[s]], max_val = max_val)
    names(eh_list_samp)[[s]] <- rev(sort(unique(rel_data$bin_assignment)))[s]
  }
  assign("eh_list", eh_list_samp, envir = .GlobalEnv)
}  

organise_det <- function(siteCoords, extracted_covs){
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
  outcrop <- data.frame(teyeq = covs$COut_all, 
                        teyep = covs$COut_all,
                        teyeo = covs$MOut_all,
                        teyen = covs$MOut_all
  )
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
  det_covs <- list(outcrop = outcrop, sedflux = sed, 
                   land = as.factor(covs$LANDCVI_multiple), 
                   MGVF = covs$MGVF, 
                   rain = covs$WC_Prec, 
                   temp = covs$WC_Temp)
}

site_remove <- function(eh_list, occ_covs, det_covs, siteCoords){
  # Find all rows with NAs in occupancy covariates
  removal <- sort(unique(c(which(is.na(occ_covs[[1]]), arr.ind=TRUE)[,1],
                           which(is.na(occ_covs[[2]]), arr.ind=TRUE)[,1], 
                           which(is.na(occ_covs[[3]]), arr.ind=TRUE)[,1], 
                           which(is.na(occ_covs[[4]]), arr.ind=TRUE)[,1])))
  # Remove those rows from other dataframes, arrange back into lists
  occ_covs <- c(lapply(occ_covs[1:4], function(x) {x <- x[-removal, ]}),
                lapply(occ_covs[5], function(x) {x <- x[-removal]}))
  eh_list <- lapply(eh_list, lapply, function(x) {x <- x[-removal, ]})
  det_covs <- c(lapply(det_covs[1:2], function(x) {x <- x[-removal, ]}), 
                 lapply(det_covs[3:6], function(x) {x <- x[-removal]}))
  siteCoords <- siteCoords[-removal,]
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
  det_covs <<- c(lapply(det_covs[1:2], function(x) {x <- x[-removal, ]}), 
                lapply(det_covs[3:6], function(x) {x <- x[-removal]}))
  siteCoords <<- siteCoords[-removal,]
}

transpose_eh <- function(eh_list, target){
  eh_list_adjusted <- list(eh_list[[1]][[1]], eh_list[[2]][[1]], 
                           eh_list[[3]][[1]], eh_list[[4]][[1]])
  testArray <- abind(eh_list_adjusted, along = 3)
  newArray <- aperm(testArray, c(1,3,2))
  assign(paste("EH_array_", target, sep = ""), newArray, envir = .GlobalEnv)
}


##### RUNNING CODE #####

get_grid(master.occs.binned.targeted, res = res, e = e)

prepare_for_spOcc(master.occs.binned.targeted.grid)

site_IDs <- sort(unique(master.occs.binned.targeted.grid$siteID))

rotated <- p_rotate(res, e, site_IDs, bins)

extracted_covs <- extract_p(rotated)

occ_covs <- format_occ_covs(extracted_covs)

det_covs <- organise_det(siteCoords, extracted_covs)

site_remove(eh_list, occ_covs, det_covs, siteCoords)

transpose_eh(eh_list, target)

revi.data <- list(y = EH_array_Ceratopsidae, 
                  occ.covs = occ_covs, 
                  det.covs = det_covs)

# Quick plot
raw.occ.prob <- apply(revi.data$y, 2, mean, na.rm = TRUE)
plot(1:4, raw.occ.prob, pch = 16, 
     xlab = 'Year', ylab = 'Raw Occurrence Proportion', 
     cex = 1.5, frame = FALSE, ylim = c(0, 1))

revi.occ.formula <- ~ scale(temp) + scale(precip) + scale(DEM) + (1 | site.effect)
revi.det.formula <- ~ scale(outcrop) + scale(MGVF) + scale(rain) + scale(temp)

z.inits <- apply(revi.data$y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
revi.inits <- list(beta = 0,         # occurrence coefficients
                   alpha = 0,        # detection coefficients
                   sigma.sq.psi = 1, # occurrence random effect variances
                   z = z.inits)      # latent occurrence values

revi.priors <- list(beta.normal = list(mean = 0, var = 2.72), 
                      alpha.normal = list(mean = 0, var = 2.72), 
                      sigma.sq.psi.ig = list(a = 0.1, b = 0.1))
ar1 <- TRUE
n.chains <- 5
n.batch <- 200
batch.length <- 25
n.samples <- n.batch * batch.length 
n.burn <- 2000
n.thin <- 12

out <- tPGOcc(occ.formula = revi.occ.formula, 
              det.formula = revi.det.formula, 
              data = revi.data, 
              n.batch = n.batch, 
              batch.length = batch.length,
              inits = revi.inits,
              priors = revi.priors,
              ar1 = ar1,
              n.burn = n.burn, 
              n.thin = n.thin, 
              n.chains = n.chains, 
              n.report = 50)
summary(out)
waicOcc(out)

ppcOut <- ppcOcc(out, fit.stat = 'freeman-tukey', group = 1)
summary(ppcOut)






##### SPATIAL #####

# Load the sp package
library(sp)

# Define the proj4 string for NAD83
nad83_proj <- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

# Define my spatial object
sp_points <- SpatialPoints(siteCoords)

proj4string(sp_points) <- CRS("+proj=longlat +datum=WGS84")

sp_points_proj <- sp::spTransform(sp_points, nad83_proj)

coords <- sp_points_proj@coords[,1:2]
names(coords) <- c("X", "Y")


revi.data <- list(y = EH_array_Ceratopsidae, 
                  occ.covs = occ_covs, 
                  det.covs = det_covs, 
                  coords = coords
)

# Quick plot
raw.occ.prob <- apply(revi.data$y, 2, mean, na.rm = TRUE)
plot(1:4, raw.occ.prob, pch = 16, 
     xlab = 'Year', ylab = 'Raw Occurrence Proportion', 
     cex = 1.5, frame = FALSE, ylim = c(0, 1))

revi.sp.occ.formula <- ~ scale(temp) + scale(precip) + scale(DEM)
revi.sp.det.formula <- ~ scale(outcrop) + scale(MGVF) + scale(rain) + scale(temp)

z.inits <- apply(revi.data$y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
# Pair-wise distance between all sites
dist.hbef <- dist(revi.data$coords)
revi.sp.inits <- list(beta = 0, alpha = 0, z = z.inits,
                      sigma.sq = 1, phi = 3 / mean(dist.hbef), 
                      sigma.sq.t = 1.5, rho = 0.2)

revi.sp.priors <- list(beta.normal = list(mean = 0, var = 2.72), 
                       alpha.normal = list(mean = 0, var = 2.72), 
                       sigma.sq.t.ig = c(2, 0.5), 
                       rho.unif = c(-1, 1),
                       sigma.sq.ig = c(2, 1), 
                       phi.unif = c(3 / max(dist.hbef), 3 / min(dist.hbef)))
cov.model <- 'exponential'
n.neighbors <- 5
ar1 <- TRUE
n.batch <- 600
batch.length <- 25
# Total number of samples
n.batch * batch.length

n.burn <- 10000
n.thin <- 20 

# Approx. run time: ~ 2.5 min
out.sp <- stPGOcc(occ.formula = revi.sp.occ.formula, 
                  det.formula = revi.sp.det.formula, 
                  data = revi.data, 
                  inits = revi.sp.inits, 
                  priors = revi.sp.priors, 
                  cov.model = cov.model, 
                  n.neighbors = n.neighbors,
                  n.batch = n.batch, 
                  batch.length = batch.length, 
                  verbose = TRUE, 
                  ar1 = ar1,
                  n.report = 200,
                  n.burn = n.burn, 
                  n.thin = n.thin, 
                  n.chains = 3) 
summary(out.sp)
waicOcc(out.sp)

ppcOut <- ppcOcc(out.sp, fit.stat = 'freeman-tukey', group = 1)
summary(ppcOut)


