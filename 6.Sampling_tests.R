################################################################################
############ OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS #############
################################################################################

# Christopher D. Dean, Lewis A. Jones, Alfio A. Chiarenza, Sinéad Lyster, 
# Alex Farnsworth, Philip D. Mannion, Richard J. Butler.
# 2023
# Script written by Christopher D. Dean and Lewis A. Jones

################################################################################
#             FILE 6: TESTING IMPACT OF RANDOM SAMPLING ON RESULTS             #
################################################################################

################################################################################
# 1. INITIAL SETUP
################################################################################

library(ggplot2)
library(stringr)
library(tidyverse)

##### Load in Functions #####
source("0.Functions.R") # Import functions from other R file (must be in same working directory)

# Set resolution
res <- 0.5

# Set extent
e <- extent(-155, -72, 22.5, 73)

max_val <- 30

bin.name <- "teyen"

# Set bin type
bin.type <- "scotese"

# Load main dataset
master.occs <- read.csv(paste("Prepped_data/Occurrence_Data/", bin.type, "/", 
                              res, "_", bin.type, "_occurrence_dataset.csv", sep = ""))
# Load formations
formations <- read.csv("Data/Occurrences/Formations.csv") # Load formations

# Set bins 
bins <- master.occs %>%
  dplyr::select(bin_assignment, bin_midpoint) %>%
  distinct()

# Create folders for storing test results
dir.create(paste0("Tests"), showWarnings = FALSE)
dir.create(paste0("Tests/", res, "/", sep =""), showWarnings = FALSE)

# Set bin to analyse
bin.occs <- master.occs %>% 
  filter(bin_assignment == bin.name)

# Set targets
target <- "Ceratopsidae" 

################################################################################
# 2. RUNNING TEST
################################################################################
reps <- 100

#set target
best.model.all <- c()

combined.res.all <- data.frame(Parameter = c(), 
                               Estimate = c(), 
                               '2.5%' = c(), 
                               '97.5%' = c(), 
                               SE = c(), 
                               Target = c(), 
                               Rep = c())
temp.model.res.all <- data.frame(name = c(), 
                             value = c(), 
                             '2.5 %' = c(), 
                             '97.5 %' = c(),
                             Rep = c())
occ_gof_all <- c()

# Standard error function
standard_error <- function(x) sd(x) / sqrt(length(x))

for(f in 1:reps){
  tryCatch({
    # Prepare data for unmarked
    print(f)
    all_results_for_unmarked(data = bin.occs, name = bin.name, res = res, ext = e, 
                             target = target, single = TRUE, formCells = "N", 
                             max_val_on = TRUE, max_val = max_val)
    
    data <- unmarked_Ceratopsidae[[1]]
    # Get covariates
    get.all.covs(data.grid)
    
    ##### NOW LOOP THROUGH TARGETS #####
  
    # Turn into matrix
    y <- as.matrix(data)
    
    ####################################
    ##### SETUP DATA FOR MODELLING #####
    ####################################
    
    ##### Unstandardised values of covariates ####
    # Detection
    DEMs.orig <- sitecov[,"mean_DEM"]      
    RAIN.orig <- sitecov[,"mean_prec"]
    TEMP.orig <- sitecov[,"mean_temp"]
    MGVF.orig <- sitecov[,"MGVF"]
    OUTs.orig <- sitecov[,"MOut"]
    OUTa.orig <- sitecov[,"AllOut"]
    LAND.orig <- sitecov[,"LANDCVI_multiple"]
    SEDf.orig <- sitecov[,"mean_psed"]
    
    # Occupancy
    Ptem.orig <- sitecov[,"mean_ptemp"]
    Ppre.orig <- sitecov[,"mean_pprec"]
    PDEM.orig <- sitecov[,"mean_pDEM"]
    
    ##### Scaling and Mode Setup #####
    siteCovs <- data.frame(DEMs = DEMs.orig, OUTs = OUTs.orig,
                           RAIN = RAIN.orig, MGVF = MGVF.orig, 
                           TEMP = TEMP.orig, LAND = LAND.orig, 
                           OUTa = OUTa.orig, SEDf = SEDf.orig,
                           Ppre = Ppre.orig, Ptem = Ptem.orig, 
                           PDEM = PDEM.orig)
    
    # Make Model
    umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs)
    
    # Scale covariates within model
    umf@siteCovs$DEMs <- scale(umf@siteCovs$DEMs)
    umf@siteCovs$OUTs <- scale(umf@siteCovs$OUTs)
    umf@siteCovs$OUTa <- scale(umf@siteCovs$OUTa)
    umf@siteCovs$RAIN <- scale(umf@siteCovs$RAIN)
    umf@siteCovs$MGVF <- scale(umf@siteCovs$MGVF)
    umf@siteCovs$SEDf <- scale(umf@siteCovs$SEDf)
    umf@siteCovs$Ppre <- scale(umf@siteCovs$Ppre)
    umf@siteCovs$Ptem <- scale(umf@siteCovs$Ptem)
    umf@siteCovs$PDEM <- scale(umf@siteCovs$PDEM)
    
    ################################################################################
    # 2. RUNNING MODELS
    ################################################################################
    
    ################################
    ##### NULL MODEL (NO COVS) #####
    ################################
    
    # Get summary of overall occupancy model
    summary(umf) 
    
    # Get naive occupancy
    naive.occ <- y %>% 
      as.data.frame() %>%
      filter(if_any(where(is.numeric), ~ .x > 0)) %>%
      nrow()
    naive.occ <- data.frame(Parameter = "Naive Occupancy", 
                            Estimate = naive.occ/nrow(y))
    
    # Make model without covariates
    summary(fm1 <- occu(~1 ~1, data=umf))
    
    # Get estimates for occupancy and detection from basic model
    print(occ.null <- backTransform(fm1, "state")) 
    print(det.null <- backTransform(fm1, "det")) 
    
    # Compile results
    null.res <- rbind(c(Estimate = mean(occ.null@estimate), 
                        confint(occ.null, level = 0.95)),
                      c(Estimate = mean(det.null@estimate), 
                        confint(det.null, level = 0.95)))
    colnames(null.res) <- c("Estimate", "2.5%", "97.5%")
    null.res <- as.data.frame(null.res)
    null.res$Parameter <- c("Null.occ.prob", "Null.det.prob")
    
    ###############################################
    ##### COVARIATE MODEL AND MODEL SELECTION #####
    ###############################################
    
    # Fit full model
    full <- occu(formula =  ~ OUTa + MGVF + TEMP + LAND + RAIN + SEDf # det
                 ~ Ppre + Ptem + PDEM, # occ
                 data = umf)
    
    # Use dredge to automatically carry out model selection
    (modelList <- dredge(full, rank = "AIC")) 
    
    # If no clear best fit model, combine models for averaged model
    occu_dredge_95 <- get.models(modelList, subset = cumsum(weight) <= 0.95)
    mod.av <- model.avg(occu_dredge_95, fit = TRUE, rank = "AICc")
    
    # Report best model and top 5 models
    best.model <- occu_dredge_95[[1]]
    best.model
    
    # Run MacKenzie and Bailey Goodness-of-fit test (WARNING: MIGHT TAKE A WHILE)
   # system.time(occ_gof <- mb.gof.test(full, nsim = 100, plot.hist = T))
   # occ_gof$p.value
    
    # Examine the effect of covariates from averaged model
    (temp.model.res <- coef(mod.av) %>% 
        enframe())
    av.CI <- as.data.frame(confint(mod.av, type='det', method = 'normal'))
    av.CI$name <- rownames(av.CI)
    temp.model.res <- merge(temp.model.res, av.CI, by = "name")
    temp.model.res$Rep <- f
    
    
    ##### MODEL STATS #####
    # Proportion of area occupied
    re <- unmarked::ranef(best.model)
    EBUP <- bup(re, stat="mean")
    CI <- confint(re, level=0.95)
    SE <- standard_error(re@post)
    (PAO <- rbind(PAO = c(Estimate = sum(EBUP), colSums(CI)) / nrow(y)))
    PAO <- as.data.frame(PAO)
    PAO$Parameter <- "PAO"
    
    # Estimate of Detection Prob. per site (model averaged)
    det.prob <- unmarked::predict(mod.av, type="det") # Predict detection for sites/vists
    (det.prob <- rbind(Det.prob = c(Estimate = mean(det.prob$fit), 
                                    SE = mean(det.prob$se.fit))))
    det.prob <- as.data.frame(det.prob)
    det.prob$Parameter <- "Det.prob"  
    
    # Estimate of Occupancy Prob. per site (model averaged)
    occ.prob <- unmarked::predict(mod.av, type="state") # Predict occupancy for sites/vists
    (occ.prob <- rbind(occ.prob = c(Estimate = mean(occ.prob$fit), 
                                    SE = mean(occ.prob$se.fit))))
    occ.prob <- as.data.frame(occ.prob)
    occ.prob$Parameter <- "Occ.prob"
    occ.det.prob <- merge(occ.prob, det.prob, all = T, sort = T)
    
    # Combine results
    comb <- merge(occ.det.prob, PAO, all = T, sort = F)
    comb <- comb[order(comb$'2.5%'),]
    # Null results and model stats
    null.res <- merge(naive.occ, null.res, all = T, sort = F)
    combined.res <- merge(null.res, comb, all = T, sort = F)
    combined.res$Rep <- f
    
    
    # Saving results
  # best.model.all <- c(best.model.all, best.model)
    combined.res.all <- rbind(combined.res.all, combined.res)
    temp.model.res.all <- rbind(temp.model.res.all, temp.model.res)
  # occ_gof_all <- c(occ_gof_all, occ_gof$p.value)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

write.csv(temp.model.res.all, paste("Tests/", res, "/", bin.name, ".", target, 
                         ".", max_val, ".param.csv", sep = ""))
write.csv(combined.res.all, paste("Tests/", res, "/", bin.name, ".", target, 
                                    ".", max_val, ".ests.csv", sep = ""))

################################################################################
# 3. INDIVIDUAL ANALYSIS
################################################################################

#################
##### SETUP #####
#################

# Select resolution and bin.name
res <- 0.5
bin.name <- "teyen"

# Load comparative data and rename
temp <- list.files(path = paste("Tests/", res, sep = ""), pattern=paste("^", bin.name, sep = ""), full.names =  TRUE)
Tests <- lapply(temp, read.csv)
names(Tests) <- gsub(paste("Tests/", res, "/", bin.name, ".", sep = ""), "", temp)
names(Tests) <- gsub("\\.csv$", "", names(Tests))

# Subset types of test and add sample value to new column
Ests <- Tests[grep('ests', names(Tests))]
Param <- Tests[grep('param', names(Tests))]
Ests <- map2(Ests, names(Ests), ~ mutate(.x, Samp = gsub(".*?([0-9]+).*", "\\1", .y)))
Param <- map2(Param, names(Param), ~ mutate(.x, Samp = gsub(".*?([0-9]+).*", "\\1", .y)))

# Combine into dataframe for analysis
Ests <- do.call("rbind", Ests)
Param <- do.call("rbind", Param)

#####################
##### ESTIMATES #####
#####################

# Reorder factors
Ests$Samp <- factor(Ests$Samp, levels = c("5", "10", "20", "30"))

# Plot results
ggplot(Ests, aes(x=Parameter, y=Estimate, fill = Samp)) + 
  geom_boxplot()

######################
##### PARAMETERS #####
######################

# Reorder factors
Param$Samp <- factor(Param$Samp, levels = c("5", "10", "20", "30"))

# Visual assessment
Param %>%
  ggplot(aes(x=name, y=value, fill = Samp)) + 
  geom_boxplot()

# Identify parameters which were statistically significant
Param <- Param %>% 
  rowwise() %>% 
  mutate(Signif = as.numeric(!between(0,X2.5..,X97.5..))) 

# Generate table of answers
table(Param$Signif, Param$Samp)

################################################################################
# 3. COMBINED ANALYSIS
################################################################################

#################
##### SETUP #####
#################

# Select resolution and bin.name
bin.name <- "teyen"

# Load comparative data and rename
temp <- list.files(path = "Tests/0.5", pattern=paste("^", bin.name, sep = ""), full.names =  TRUE)
Tests <- lapply(temp, read.csv)
names(Tests) <- gsub(paste("Tests/0.5/", bin.name, ".", sep = ""), "", temp)
names(Tests) <- gsub("\\.csv$", "", names(Tests))

# Subset types of test and add sample value to new column
Ests <- Tests[grep('ests', names(Tests))]
Param <- Tests[grep('param', names(Tests))]
Ests <- map2(Ests, names(Ests), ~ mutate(.x, Samp = gsub(".*?([0-9]+).*", "\\1", .y)))
Param <- map2(Param, names(Param), ~ mutate(.x, Samp = gsub(".*?([0-9]+).*", "\\1", .y)))

# Combine into dataframe for analysis
Ests <- do.call("rbind", Ests)
Param <- do.call("rbind", Param)

Ests$Res <- 0.5
Param$Res <- 0.5

# Load comparative data and rename
temp <- list.files(path = "Tests/1", pattern=paste("^", bin.name, sep = ""), full.names =  TRUE)
Tests <- lapply(temp, read.csv)
names(Tests) <- gsub(paste("Tests/1/", bin.name, ".", sep = ""), "", temp)
names(Tests) <- gsub("\\.csv$", "", names(Tests))

# Subset types of test and add sample value to new column
Ests2 <- Tests[grep('ests', names(Tests))]
Param2 <- Tests[grep('param', names(Tests))]
Ests2 <- map2(Ests2, names(Ests2), ~ mutate(.x, Samp = gsub(".*?([0-9]+).*", "\\1", .y)))
Param2 <- map2(Param2, names(Param2), ~ mutate(.x, Samp = gsub(".*?([0-9]+).*", "\\1", .y)))

# Combine into dataframe for analysis
Ests2 <- do.call("rbind", Ests2)
Param2 <- do.call("rbind", Param2)

Ests2$Res <- 1
Param2$Res <- 1

# Combine all
Ests <- rbind(Ests, Ests2)
Param <- rbind(Param, Param2)

#####################
##### ESTIMATES #####
#####################

# Reorder factors
Ests$Res <- factor(Ests$Res, levels = c("0.5", "1"))

# Plot results
ggplot(Ests, aes(x=Parameter, y=Estimate, fill = Res)) + 
  geom_boxplot()

Ests %>%
  ggplot(aes(x=Parameter, y=Estimate, fill = Res)) +
  geom_boxplot() +
  facet_wrap(~Samp, scales = 'free')

######################
##### PARAMETERS #####
######################

# Reorder factors
Param$Res <- factor(Param$Res, levels = c("0.5", "1"))

# Visual assessment
Param %>%
  ggplot(aes(x=name, y=value, fill = Res)) + 
  geom_boxplot()

# Identify parameters which were statistically significant
Param <- Param %>% 
  rowwise() %>% 
  mutate(Signif = as.numeric(!between(0,X2.5..,X97.5..))) 

# Generate table of answers
table(Param$Signif, Param$Res, Param$Samp)

################################################################################
# A1. EDITS TO FUNCTIONS
################################################################################
precise_cov <- function(data, samp.data){
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
  
  surv.data <- data.frame(siteID = rep(rownames(samp.data), each = 10),
                          colls = matrix(t(samp.data), ncol=1, nrow=ncol(samp.data)*nrow(samp.data), byrow=F))
  
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
  surv.data <<- surv.data
  hires_cov_dat <<- hires_cov_dat
}


get.all.covs <- function(bin.occs){
  
  precise_cov(bin.occs, unmarked_Ceratopsidae[[2]])
  
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
  # Load rasters
  wc <- list.files(paste("Prepped_data/Covariate_Data/All_data/", 
                         res, "deg/Palaeo/", sep = ""), 
                   pattern=paste0("^", bins$code[t], ".*", sep = ""))
  stacked <- raster::stack(paste("Prepped_data/Covariate_Data/All_data/", 
                                 res, "deg/Palaeo/", wc, 
                                 sep =""))
  # Set palaeo lat/long
  names(bin.occs)[names(bin.occs) == "p_lng_PALEOMAP"] <- "plng"
  names(bin.occs)[names(bin.occs) == "p_lat_PALEOMAP"] <- "plat"
  
  bin.occs <- get_p_cov(bin.occs, stacked)
  
  pcovs <- unlist(strsplit(wc, ".asc"))
  
  p.site.covs <- bin.occs %>%
    dplyr::select(siteID, any_of(pcovs), Coll_count) %>%
    dplyr::filter(Coll_count == "Non-singleton") %>%
    dplyr::rename_with(~ sub(paste("^", bin.name, sep = ""), "p", .x), starts_with(bin.name)) %>%
    dplyr::rename_with(~gsub("\\d+", "", .)) %>%
    dplyr:: rename_with(~gsub("\\.", "", .)) %>%
    dplyr::group_by(siteID) %>%
    dplyr:: summarise(mean_ptemp = mean(p_temp_, na.rm = TRUE), 
                      mean_pprec = mean(p_precip_, na.rm = TRUE),
                      mean_pDEM = mean(p_DEM_, na.rm = TRUE),
                      mean_psed = mean(p_sed_, na.rm = TRUE)
    )
  site.covs <- merge(site.covs, p.site.covs, by = "siteID")
  site.covs <- site.covs %>% select(-Coll_count)
  
  # Setup survey covariates
  survcov <- surv.data %>%
    select(colls, temp, prec, DEM)
  names(survcov)[names(survcov) == "DEM"] <- "DEM_surv"
  
  # Setup site covariates
  precise.sitecov <- hires_cov_dat %>%
    select(siteID, mean_DEM, mean_prec, mean_temp)
  sitecov <- merge(precise.sitecov, site.covs, by = "siteID")
  sitecov <- select(sitecov, -c(siteID))
  
  # Remove all digits at end of the column names
  colnames(sitecov) <- sub("_\\d.*", "", colnames(sitecov))
  
  # Set categorical variables
  sitecov$LANDCVI_multiple <- as.character(sitecov$LANDCVI_multiple)
  sitecov$LANDCVI_binary <- as.character(sitecov$LANDCVI_binary)
  
  sitecov <<- sitecov
  survcov <<- survcov
}

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
          set.seed(rand) # set seed to ensure all targets are have same sampling
          test2 <- sample_for_unmarked(test2, max_val)
          temp_name <- paste("unmarked_", target[q], sep = "")
          assign(temp_name, test2, envir = .GlobalEnv)
        }
      }
      # Create folders, remove warning if they already exist.
    }
    proc.time() - ptm
  }
}

