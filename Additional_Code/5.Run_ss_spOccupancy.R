################################################################################
############ OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS #############
################################################################################

# Christopher D. Dean, Alfio Alessandro Chiarenza, Jeffrey W. Doser, Alexander
# Farnsworth, Lewis A. Jones, Sinéad Lyster, Charlotte L. Outhwaite, Paul J. 
# Valdes, Richard J. Butler, Philip D. Mannion.
# 2024
# Script written by Christopher D. Dean

################################################################################
#                    FILE 5: RUNNING SINGLE SEASON OCCUPANCY                   #
################################################################################

################################################################################
# 1. INITIAL SETUP
################################################################################

library(spOccupancy)
library(stars)
library(ggplot2)
library(abind)
library(stringr)
library(MCMCvis)
library(purrr)
library(sp)

##### Load in Functions #####
source("0.Functions.R") # Import functions from other R file (must be in same working directory)

##### Set values #####
# Set resolution
res1 <- c(0.5, 1)
# Set extent
e <- extent(-155, -72, 22.5, 73)
bin.type <- "scotese"
target1 <- c("Hadrosauridae", "Ceratopsidae", "Tyrannosauridae")
bins <- read.csv("Data/Occurrences/scotesebins.csv")
bins$bin <- bins$code
form_cells <- "N"
bin1 <- c("teyen", "teyeo", "teyep", "teyeq")

################################################################################
# 2. START LOOP TO GENERATE RESULTS
################################################################################

for(b in bin1){
  bin <- b
  for(r in res1){
    res <- r
    for(t in target1){
      target <- t
      
      if(form_cells == "Y"){
        # Load occurrence dataset
        sp.data <- readRDS(file = paste("Prepped_data/spOccupancy/Single_season/", res, "/", 
                                        bin, "/", target, "_single_", res, ".formcells.rds", sep = ""))
      }else{
        # Load occurrence dataset
        sp.data <- readRDS(file = paste("Prepped_data/spOccupancy/Single_season/", res, "/", 
                                        bin, "/", target, "_single_", res, ".rds", sep = ""))
        # Load occ. covariates adjusted for multicollinearity
        occ.form.1 <- readRDS(file = paste("Prepped_data/spOccupancy/Single_season/", res, "/", 
                                        bin, "/", target, "_single_", res, ".occ.form.rds", sep = ""))
        # Load det. covariates adjusted for multicollinearity
        det.form.1 <- readRDS(file = paste("Prepped_data/spOccupancy/Single_season/", res, "/", 
                                        bin, "/", target, "_single_", res, ".det.form.rds", sep = ""))
        det.form.1 <- det.form.1[!(det.form.1 %in% c("factor(Year)", "scale(temp)"))]
      }
      
      # Adjustments to saved data (whoops...)
      sp.data[[1]] <- as.matrix(sp.data[[1]])
      
      
      # Pair-wise distances between all sites
      dist.hbef <- dist(sp.data$coords)
      # Exponential covariance model
      cov.model <- "exponential"
      # Specify list of inits
      oven.inits <- list(alpha = 0, 
                         beta = 0, 
                         z = apply(sp.data$y, 1, max, na.rm = TRUE), 
                         sigma.sq = 2, 
                         phi = 3 / mean(dist.hbef), 
                         w = rep(0, nrow(sp.data$y)))
      batch.length <- 25
      n.batch <- 800
      n.burn <- 2000
      n.thin <- 10
      n.chains <- 3
      oven.tuning <- list(phi = 1)
      # accept.rate = 0.43 by default, so we do not specify it.
      min.dist <- min(dist.hbef)
      max.dist <- max(dist.hbef)
      oven.priors <- list(beta.normal = list(mean = 0, var = 2.72), 
                          alpha.normal = list(mean = 0, var = 2.72), 
                          sigma.sq.ig = c(2, 1), 
                          phi.unif = c(3/max.dist, 3/min.dist))
      n.omp.threads <- 1
      verbose <- TRUE
      n.report <- 100 # Report progress at every 100th batch.
      
      ################################################################################
      # 3. MODEL SELECTION
      ################################################################################
      
      ################################
      ##### DETECTION COVARIATES #####
      ################################
      
      # Set total potential covariates
      occ.form <- unlist(lapply(1:length(occ.form.1), 
                                function(x) combn(occ.form.1, x, simplify = FALSE)), 
                         recursive = FALSE)
      det.form <- unlist(lapply(1:length(det.form.1), 
                                function(x) combn(det.form.1, x, simplify = FALSE)), 
                         recursive = FALSE)
      occ.form <- append(occ.form, "1")
      det.form <- append(det.form, "1")
      
      det.wrapper <- function(det.form){
        revi.sp.det.formula <- formula(paste("~", paste(det.form, collapse = " + ")))
        revi.sp.occ.formula <- formula(paste("~", paste(occ.form.1, collapse = " + ")))
        print(paste("Running model ", det.form, " out of ", length(det.form), sep = ""))
        out.sp <- spOccupancy::spPGOcc(occ.formula = revi.sp.occ.formula, 
                                       det.formula = revi.sp.det.formula, 
                                       data = sp.data, 
                                       inits = oven.inits, 
                                       n.batch = n.batch, 
                                       batch.length = batch.length, 
                                       priors = oven.priors, 
                                       cov.model = cov.model, 
                                       NNGP = TRUE, 
                                       n.neighbors = 5,
                                       tuning = oven.tuning, 
                                       n.report = n.report, 
                                       n.burn = n.burn, 
                                       n.thin = n.thin, 
                                       n.chains = n.chains)
        det_waic <- t(as.data.frame(spOccupancy::waicOcc(out.sp)))
        rownames(det_waic) <- c(paste("det ~", revi.sp.det.formula)[2])
        return(det_waic)
      }
      
      # Initiate the cluster
      sfInit(parallel = TRUE, cpus = 4)
      
      # Export data to the cluster
      sfExport('sp.data')
      sfExport('occ.form.1')
      sfExport('oven.inits')
      sfExport('n.batch')
      sfExport('batch.length')
      sfExport('oven.priors')
      sfExport('cov.model')
      sfExport('oven.tuning')
      sfExport('n.report')
      sfExport('n.burn')
      sfExport('n.thin')
      sfExport('n.chains')
      sfExport('dist.hbef')
      
      # Run the model in parallel
      system.time({
        det.waic <- sfClusterApplyLB(det.form, det.wrapper)
      })
      
      # Stop the cluster
      sfStop()
      
      # Unlist
      det.waic <- as.data.frame(do.call(rbind, det.waic))
      det.waic$model <- rownames(det.waic)
      
      # Set det.
      best.det <- as.formula(paste(gsub("det", "", det.waic[which.min(det.waic$WAIC), 4]), "+ (1 | Site)", sep = " "))
      
      ################################
      ##### OCCUPANCY COVARIATES #####
      ################################
      
      # Make function for snowfall
      occ.wrapper <- function(occ.form){
        revi.sp.occ.formula <- formula(paste("~", paste(occ.form, collapse = " + ")))
        revi.sp.det.formula <- best.det
        out.sp <- spOccupancy::spPGOcc(occ.formula = revi.sp.occ.formula, 
                                       det.formula = revi.sp.det.formula, 
                                       data = sp.data, 
                                       inits = oven.inits, 
                                       n.batch = n.batch, 
                                       batch.length = batch.length, 
                                       priors = oven.priors, 
                                       cov.model = cov.model, 
                                       NNGP = TRUE, 
                                       n.neighbors = 5,
                                       tuning = oven.tuning, 
                                       n.report = n.report, 
                                       n.burn = n.burn, 
                                       n.thin = n.thin, 
                                       n.chains = n.chains)
        occ_waic <- t(as.data.frame(spOccupancy::waicOcc(out.sp)))
        rownames(occ_waic) <- c(paste("occ ~", revi.sp.occ.formula)[2])
        return(occ_waic)
      }
      
      # Initiate the cluster
      sfInit(parallel = TRUE, cpus = 4)
      
      # Export data to the cluster
      sfExport('sp.data')
      sfExport('best.det')
      sfExport('oven.inits')
      sfExport('n.batch')
      sfExport('batch.length')
      sfExport('oven.priors')
      sfExport('cov.model')
      sfExport('oven.tuning')
      sfExport('n.report')
      sfExport('n.burn')
      sfExport('n.thin')
      sfExport('n.chains')
      sfExport('dist.hbef')
      
      # Run the model in parallel
      system.time({
        occ.waic <- sfClusterApplyLB(occ.form, occ.wrapper)
      })
      
      # Stop the cluster
      sfStop()
      
      # Unlist
      occ.waic <- as.data.frame(do.call(rbind, occ.waic))
      occ.waic$model <- rownames(occ.waic)
      
      best.occ <- as.formula(paste(gsub("occ", "", occ.waic[which.min(occ.waic$WAIC), 4]), sep = " "))
      
      ################################################################################
      # 3. BEST MODEL AND OTHER INFO
      ################################################################################
      
      # Covariates
      revi.sp.det.formula <- best.det
      revi.sp.occ.formula <- best.occ
      
      # Best Model
      out.sp <- spOccupancy::spPGOcc(occ.formula = revi.sp.occ.formula, 
                                     det.formula = revi.sp.det.formula, 
                                     data = sp.data, 
                                     inits = oven.inits, 
                                     n.batch = n.batch, 
                                     batch.length = batch.length, 
                                     priors = oven.priors, 
                                     cov.model = cov.model, 
                                     NNGP = TRUE, 
                                     n.neighbors = 5,
                                     tuning = oven.tuning, 
                                     n.report = n.report, 
                                     n.burn = n.burn, 
                                     n.thin = n.thin, 
                                     n.chains = n.chains)
      
      # Summary
      summary(out.sp)
      
      # Goodness-of-fit across space
      GOFspace <- spOccupancy::ppcOcc(out.sp, fit.stat = 'freeman-tukey', group = 1)
      summary(GOFspace)
      # Goodness-of-fit across replicates
      GOFrep <- spOccupancy::ppcOcc(out.sp, fit.stat = 'freeman-tukey', group = 2)
      summary(GOFrep)
      
      # Occupancy covariates
      MCMCvis::MCMCplot(out.sp$beta.samples, ref_ovl = TRUE, ci = c(50, 95))
      # Detection covariates
      MCMCvis::MCMCplot(out.sp$alpha.samples, ref_ovl = TRUE, ci = c(50, 95))
      
      
      saveRDS(out.sp, file = paste("Results/spOccupancy/single_season/",
                                   res, "/", bin, ".", target, ".best.model.rds", 
                                   sep = ""))
      saveRDS(GOFspace, file = paste("Results/spOccupancy/single_season/",
                                     res, "/", bin, ".", target, ".gof.space.rds", 
                                     sep = ""))
      saveRDS(GOFrep, file = paste("Results/spOccupancy/single_season/",
                                   res, "/", bin, ".", target, ".gof.rep.rds", 
                                   sep = ""))
      
      ################################################################################
      # 5. PLOTTING SITE-LEVEL DETECTION PROBABILITY RANDOM EFFECTS
      ################################################################################
      #if(is.null(out.sp$alpha.star.samples) == T){
      #  next
      #}
      ## Find random effects sizes
      #alpha.star.means <- apply(out.sp$alpha.star.samples, 2, mean)
      #
      ## Create siteIDs (individual cell numbers used)
      #siteIDs <- as.numeric(rownames(sp.data$y))
      #
      ## Find coordinates of each site
      #siteCoords <- siteCoordsFun(res = res, e = e, siteIDs)
      #
      ## Generate a raster using random effects
      #raster_for_values <- gen_raster(siteCoords$siteID, alpha.star.means, res = res, ext = e)
      #
      ## Find map to use as backdrop
      #countries <- maps::map("world", plot=FALSE, fill = TRUE) 
      ## Turn map into spatialpolygons
      #countries <<- maptools::map2SpatialPolygons(countries, 
      #                                            IDs = countries$names, 
      #                                            proj4string = CRS("+proj=longlat")) 
      #mapTheme <- rasterVis::rasterTheme(region=brewer.pal(8,"Reds"))
      #
      #(p1 <- rasterVis::levelplot(raster_for_values, margin=T, par.settings=mapTheme) + 
      #    # Plots state lines
      #    latticeExtra::layer(sp.polygons(states, col = "white", fill = NA), under = T)  + 
      #    # Plots background colour
      #    latticeExtra::layer(sp.polygons(countries, col = 0, fill = "light grey"), under = T))
      ##save_lattice(p1, ss = T)
    }
  }
}
