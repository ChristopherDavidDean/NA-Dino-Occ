################################################################################
############ OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS #############
################################################################################

# Christopher D. Dean, Alfio Alessandro Chiarenza, Jeffrey W. Doser, Alexander
# Farnsworth, Lewis A. Jones, Sinéad Lyster, Charlotte L. Outhwaite, Paul J. 
# Valdes, Richard J. Butler, Philip D. Mannion.
# 2025
# Script written by Christopher D. Dean

################################################################################
#                 FILE 6: RUNNING MULTI-SEASON OCCUPANCY LOOP                 #
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
library(tictoc)

##### Load in Functions #####
source("0.Functions.R") # Import functions from other R file (must be in same working directory)

##### Set values #####
# Set resolution
res1 <- c(0.5, 1)
# Set extent
e <- extent(-155, -72, 22.5, 73)
# Set max limit value
max_val_on <- TRUE
max_val1 <- c(10, 40)
# Set bins
bin.type <- "scotese"
bins <- read.csv("Data/Occurrences/scotesebins.csv")
bins$bin <- bins$code
# Set target taxa
target1 <- c("Ankylosauridae", "Ceratopsidae", "Hadrosauridae", "Tyrannosauridae")
# Set small taxa removed dataset
nomam <- F

################################################################################
# 2. START LOOP
################################################################################
for(m in max_val1){
  max_val <- m
  for(r in res1){
    res <- r
    for(t in target1){
      target <- t
      # Load occurrence dataset
      if(nomam == F){
        sp.data <- readRDS(file = paste("Prepped_data/spOccupancy/Updated/Multi_season/", res, "/", 
                                        target, "_multi_", res, "_", max_val, ".rds", sep = ""))
        occ.form.1 <- readRDS(file = paste("Prepped_data/spOccupancy/Updated/Multi_season/", res, "/", 
                                           target, "_multi_", res, "_", max_val, ".occ.form.rds", sep = ""))
        
        det.form.1 <- readRDS(file = paste("Prepped_data/spOccupancy/Updated/Multi_season/", res, "/", 
                                           target, "_multi_", res, "_", max_val, ".det.form.rds", sep = ""))
      }else{
        sp.data <- readRDS(file = paste("Prepped_data/spOccupancy/Updated/Multi_season/", res, "/", 
                                        target, "_multi_", res, "_", max_val, "_no_mammal.rds", sep = ""))
        occ.form.1 <- readRDS(file = paste("Prepped_data/spOccupancy/Updated/Multi_season/", res, "/", 
                                           target, "_multi_", res, "_", max_val, "_no_mammal.occ.form.rds", sep = ""))
        
        det.form.1 <- readRDS(file = paste("Prepped_data/spOccupancy/Updated/Multi_season/", res, "/", 
                                           target, "_multi_", res, "_", max_val, "_no_mammal.det.form.rds", sep = ""))
      }

      det.form.1 <- det.form.1[!(det.form.1 %in% c("scale(temp)"))]
      
      z.inits <- apply(sp.data$y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
      
      # Pair-wise distance between all sites
      dist.hbef <- dist(sp.data$coords)
      
      revi.sp.inits <- list(beta = 0, alpha = 0, z = z.inits,
                            sigma.sq = 1, phi = 3 / mean(dist.hbef), 
                            sigma.sq.t = 1.5, rho = 0.2)
      
      revi.sp.priors <- list(beta.normal = list(mean = 0, var = 2.72), 
                             alpha.normal = list(mean = 0, var = 2.72), 
                             sigma.sq.t.ig = c(2, 0.5), 
                             rho.unif = c(-1, 1),
                             sigma.sq.ig = c(2, 1), 
                             phi.unif = c(3 / max(dist.hbef), 3 / min(dist.hbef))) # can be fine tuned
      cov.model <- 'exponential'
      n.neighbors <- 15
      ar1 <- TRUE
      n.batch <- 800
      batch.length <- 25
      n.burn <- 10000
      n.thin <- 10 
      
      # Total number of samples
      n.batch * batch.length
      
      # Set total potential covariates
      occ.form <- unlist(lapply(1:length(occ.form.1), 
                                function(x) combn(occ.form.1, x, simplify = FALSE)), 
                         recursive = FALSE)
      det.form <- unlist(lapply(1:length(det.form.1), 
                                function(x) combn(det.form.1, x, simplify = FALSE)), 
                         recursive = FALSE)
      occ.form <- append(occ.form, "1")
      det.form <- append(det.form, "1")
      
      ################################################################################
      # 3. MODEL SELECTION
      ################################################################################
      
      ################################
      ##### DETECTION COVARIATES #####
      ################################
      
      # Make function for snowfall
      det.wrapper <- function(det.form){
        tictoc::tic("Detection covariate loop")
        revi.sp.det.formula <- formula(paste("~", paste(det.form, collapse = " + ")))
        revi.sp.occ.formula <- formula(paste("~", paste(occ.form.1, collapse = " + ")))
        print(paste("Running model ", det.form, " out of ", length(det.form), sep = ""))
        out.sp <- spOccupancy::stPGOcc(occ.formula = revi.sp.occ.formula, 
                                       det.formula = revi.sp.det.formula, 
                                       data = sp.data, 
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
        det_waic <- t(as.data.frame(spOccupancy::waicOcc(out.sp)))
        rownames(det_waic) <- c(paste("det ~", revi.sp.det.formula)[2])
        tictoc::toc()
        return(det_waic)
      }
      
      # Initiate the cluster
      sfInit(parallel = TRUE, cpus = 4)
      
      # Export data to the cluster
      sfExport('sp.data')
      sfExport('z.inits')
      sfExport('dist.hbef')
      sfExport('revi.sp.inits')
      sfExport('revi.sp.priors')
      sfExport('cov.model')
      sfExport('n.neighbors')
      sfExport('ar1')
      sfExport('n.batch')
      sfExport('batch.length')
      sfExport('n.burn')
      sfExport('n.thin')
      sfExport('occ.form.1')
      
      # Run the model in parallel
      system.time({
        det.waic <- sfClusterApplyLB(det.form, det.wrapper)
      })
      
      # Stop the cluster
      sfStop()
      
      det.waic2 <- det.waic
      
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
        out.sp <- spOccupancy::stPGOcc(occ.formula = revi.sp.occ.formula, 
                                       det.formula = revi.sp.det.formula, 
                                       data = sp.data, 
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
        occ_waic <- t(as.data.frame(spOccupancy::waicOcc(out.sp)))
        rownames(occ_waic) <- c(paste("occ ~", revi.sp.occ.formula)[2])
        return(occ_waic)
      }
      
      # Initiate the cluster
      sfInit(parallel = TRUE, cpus = 4)
      
      # Export data to the cluster
      sfExport('sp.data')
      sfExport('z.inits')
      sfExport('dist.hbef')
      sfExport('revi.sp.inits')
      sfExport('revi.sp.priors')
      sfExport('cov.model')
      sfExport('n.neighbors')
      sfExport('ar1')
      sfExport('n.batch')
      sfExport('batch.length')
      sfExport('n.burn')
      sfExport('n.thin')
      sfExport('best.det')
      
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
      # 4. RESULTS FOR BEST MODEL
      ################################################################################
      
      # Running best model
      out.sp <- stPGOcc(occ.formula = best.occ, 
                        det.formula = best.det, 
                        data = sp.data, 
                        inits = revi.sp.inits, 
                        priors = revi.sp.priors, 
                        cov.model = cov.model, 
                        n.neighbors = 15,
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
      
      # Goodness-of-fit across space
      GOFspace <- ppcOcc(out.sp, fit.stat = 'freeman-tukey', group = 1)
      summary(GOFspace)
      # Goodness-of-fit across replicates
      GOFrep <- ppcOcc(out.sp, fit.stat = 'freeman-tukey', group = 2)
      summary(GOFrep)
      
      # Occupancy
      MCMCplot(out.sp$beta.samples, ref_ovl = TRUE, ci = c(50, 95))
      # Detection
      MCMCplot(out.sp$alpha.samples, ref_ovl = TRUE, ci = c(50, 95))
      
      occ.waic$submodel <- "occupancy"
      det.waic$submodel <- "detection"
      
      all.waic <- rbind(occ.waic, det.waic)
      
      if(nomam == F){
        saveRDS(make_table(out.sp, target, res = res), file = paste("Results/spOccupancy/NEW/",
                                                                    res, "/", target, "_", max_val,".2.cov.table.rds", 
                                                                    sep = ""))
        saveRDS(out.sp, file = paste("Results/spOccupancy/NEW/",
                                     res, "/", target,"_", max_val, ".best.model.rds", 
                                     sep = ""))
        saveRDS(GOFspace, file = paste("Results/spOccupancy/NEW/",
                                       res, "/", target, "_", max_val,".gof.space.rds", 
                                       sep = ""))
        saveRDS(GOFrep, file = paste("Results/spOccupancy/NEW/",
                                     res, "/", target, "_", max_val,".gof.rep.rds", 
                                     sep = ""))
        saveRDS(all.waic, file = paste("Results/spOccupancy/NEW/",
                                       res, "/", target, "_", max_val,".all.waic.rds", 
                                       sep = ""))
      }else{
        saveRDS(make_table(out.sp, target, res = res), file = paste("Results/spOccupancy/NEW/",
                                                                    res, "/", target, "_", max_val,"_no_mammal.cov.table.rds", 
                                                                    sep = ""))
        saveRDS(out.sp, file = paste("Results/spOccupancy/NEW/",
                                     res, "/", target,"_", max_val, "_no_mammal.best.model.rds", 
                                     sep = ""))
        saveRDS(GOFspace, file = paste("Results/spOccupancy/NEW/",
                                       res, "/", target, "_", max_val,"_no_mammal.gof.space.rds", 
                                       sep = ""))
        saveRDS(GOFrep, file = paste("Results/spOccupancy/NEW/",
                                     res, "/", target, "_", max_val,"_no_mammal.gof.rep.rds", 
                                     sep = ""))
        saveRDS(all.waic, file = paste("Results/spOccupancy/NEW/",
                                       res, "/", target, "_", max_val,"_no_mammal.all.waic.rds", 
                                       sep = ""))
      }

      ################################################################################
      # 5. PLOTTING SITE-LEVEL DETECTION PROBABILITY RANDOM EFFECTS
      ################################################################################
      #
      ## Find random effects sizes
      #alpha.star.means <- apply(out.sp$alpha.star.samples, 2, mean)
      #
      ## Create siteIDs (individual cell numbers used)
      #siteIDs <- as.numeric(rownames(sp.data$y[,1,]))
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
      #(p2 <- rasterVis::levelplot(raster_for_values, margin=T, par.settings=mapTheme) + 
      #    # Plots state lines
      #    latticeExtra::layer(sp.polygons(states, col = "white", fill = NA), under = T)  + 
      #    # Plots background colour
      #    latticeExtra::layer(sp.polygons(countries, col = 0, fill = "light grey"), under = T))
      #
      #save_lattice(p2)
    }
  }
}
  