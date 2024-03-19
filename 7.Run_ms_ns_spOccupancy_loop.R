################################################################################
############ OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS #############
################################################################################

# Christopher D. Dean, Lewis A. Jones, Alfio A. Chiarenza, Sinéad Lyster, 
# Alex Farnsworth, Philip D. Mannion, Richard J. Butler.
# 2023
# Script written by Christopher D. Dean

################################################################################
#          FILE 7: RUNNING NON-SPATIAL MULTI-SEASON OCCUPANCY LOOP             #
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
res1 <- c(0.5,1)
# Set extent
e <- extent(-155, -72, 22.5, 73)
# Set max limit value
max_val_on <- TRUE
bin.type <- "scotese"
target1 <- c("Tyrannosauridae", "Hadrosauridae", "Ceratopsidae")
bins <- read.csv("Data/Occurrences/scotesebins.csv")
bins$bin <- bins$code

################################################################################
# 2. START LOOP
################################################################################
for(r in res1){
  res <- r
  max_val <- 40
  for(t in target1){
    target <- t
    # Load occurrence dataset
    sp.data <- readRDS(file = paste("Prepped_data/spOccupancy/Multi_season/", res, "/", 
                                    target, "_multi_", res, ".rds", sep = ""))
    
    # Specify model parameters
    z.inits <- apply(sp.data$y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
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
    
    # Set total potential covariates
    occ.form <- c("scale(ann)", "scale(hot)", "scale(col)", "scale(wet)", "scale(dry)")
    det.form <- c("scale(outcrop)", "scale(MGVF)", "scale(rain)", "factor(Year)", 
                  "factor(land)", "scale(occur)", "Scale(coll)", "scale(sedflux)", "scale(Distance)")
    occ.form <- unlist(lapply(1:length(occ.form), 
                              function(x) combn(occ.form, x, simplify = FALSE)), 
                       recursive = FALSE)
    det.form <- unlist(lapply(1:length(det.form), 
                              function(x) combn(det.form, x, simplify = FALSE)), 
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
      tictoc::tic("Detection Covariate")
      revi.det.formula <- formula(paste("~", paste(det.form, collapse = " + ")))
      revi.occ.formula <- ~ scale(col) + scale(ann) + scale(wet) + scale(dry) + scale(hot)
      print(paste("Running model ", det.form, " out of ", length(det.form), sep = ""))
      out <- spOccupancy::tPGOcc(occ.formula = revi.occ.formula, 
                    det.formula = revi.det.formula, 
                    data = sp.data, 
                    n.batch = n.batch, 
                    batch.length = batch.length,
                    inits = revi.inits,
                    priors = revi.priors,
                    ar1 = ar1,
                    n.burn = n.burn, 
                    n.thin = n.thin, 
                    n.chains = n.chains, 
                    n.report = 50)
      det_waic <- t(as.data.frame(spOccupancy::waicOcc(out)))
      rownames(det_waic) <- c(paste("det ~", revi.det.formula)[2])
      tictoc::toc()
      return(det_waic)
    }
    
    # Initiate the cluster
    sfInit(parallel = TRUE, cpus = 4)
    
    # Export data to the cluster
    sfExport('sp.data')
    sfExport('z.inits')
    sfExport('revi.inits')
    sfExport('revi.priors')
    sfExport('n.burn')
    sfExport('n.thin')
    sfExport('batch.length')
    sfExport('n.chains') 
    sfExport('ar1')
    sfExport('n.batch')
    sfExport('n.samples')

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
      revi.occ.formula <- formula(paste("~", paste(occ.form, collapse = " + ")))
      revi.det.formula <- best.det
      out <- spOccupancy::tPGOcc(occ.formula = revi.occ.formula, 
                                 det.formula = revi.det.formula, 
                                 data = sp.data, 
                                 n.batch = n.batch, 
                                 batch.length = batch.length,
                                 inits = revi.inits,
                                 priors = revi.priors,
                                 ar1 = ar1,
                                 n.burn = n.burn, 
                                 n.thin = n.thin, 
                                 n.chains = n.chains, 
                                 n.report = 50)
      occ_waic <- t(as.data.frame(spOccupancy::waicOcc(out)))
      rownames(occ_waic) <- c(paste("occ ~", revi.occ.formula)[2])
      return(occ_waic)
    }
    
    # Initiate the cluster
    sfInit(parallel = TRUE, cpus = 4)
    
    # Export data to the cluster
    sfExport('sp.data')
    sfExport('z.inits')
    sfExport('revi.inits')
    sfExport('revi.priors')
    sfExport('n.burn')
    sfExport('n.thin')
    sfExport('batch.length')
    sfExport('n.chains') 
    sfExport('ar1')
    sfExport('n.batch')
    sfExport('n.samples')
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
    out <- spOccupancy::tPGOcc(occ.formula = best.occ, 
                               det.formula = best.det, 
                               data = sp.data, 
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
    
    # Goodness-of-fit across space
    GOFspace <- ppcOcc(out, fit.stat = 'freeman-tukey', group = 1)
    summary(GOFspace)
    # Goodness-of-fit across replicates
    GOFrep <- ppcOcc(out, fit.stat = 'freeman-tukey', group = 2)
    summary(GOFrep)
    
    # Occupancy
    #MCMCplot(out$beta.samples, ref_ovl = TRUE, ci = c(50, 95))
    # Detection
    #MCMCplot(out$alpha.samples, ref_ovl = TRUE, ci = c(50, 95))
    
    occ.waic$submodel <- "occupancy"
    det.waic$submodel <- "detection"
    
    all.waic <- rbind(occ.waic, det.waic)
    
    saveRDS(make.table(out, target, res), file = paste("Results/spOccupancy/Non.spatial/",
                                                     res, "/", target, ".cov.table.rds", 
                                                     sep = ""))
    saveRDS(out, file = paste("Results/spOccupancy/Non.spatial/",
                                 res, "/", target, ".best.model.rds", 
                                 sep = ""))
    saveRDS(GOFspace, file = paste("Results/spOccupancy/Non.spatial/",
                                   res, "/", target, ".gof.space.rds", 
                                   sep = ""))
    saveRDS(GOFrep, file = paste("Results/spOccupancy/Non.spatial/",
                                 res, "/", target, ".gof.rep.rds", 
                                 sep = ""))
    saveRDS(all.waic, file = paste("Results/spOccupancy/Non.spatial/",
                                   res, "/", target, ".all.waic.rds", 
                                   sep = ""))
    
    ################################################################################
    # 5. PLOTTING SITE-LEVEL DETECTION PROBABILITY RANDOM EFFECTS
    ################################################################################
    
    ## Find random effects sizes
    #alpha.star.means <- apply(out$alpha.star.samples, 2, mean)
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
    #save.lattice(p2)
  }
}
