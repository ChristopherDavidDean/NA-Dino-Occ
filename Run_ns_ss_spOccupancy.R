################################################################################
############ OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS #############
################################################################################

# Christopher D. Dean, Alfio Alessandro Chiarenza, Jeffrey W. Doser, Alexander
# Farnsworth, Lewis A. Jones, Sinéad Lyster, Charlotte L. Outhwaite, Richard J. 
# Butler, Philip D. Mannion.
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
form_cells <- "Y"
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
      sp.data <- readRDS(file = paste("Prepped_data/spOccupancy/Single_season/", res, "/", 
                                        bin, "/", target, "_single_", res, ".formcells.rds", sep = ""))
      # Adjustments to saved data (whoops...)
      sp.data[[1]] <- as.matrix(sp.data[[1]])
      
      # Set parameters
      oven.inits <- list(alpha = c(0, 0, 0, 0), 
                         beta = c(0, 0, 0), 
                         z = apply(sp.data$y, 1, max, na.rm = TRUE))

      oven.inits <- list(alpha = 0, 
                         beta = 0, 
                         z = apply(sp.data$y, 1, max, na.rm = TRUE))
      
      oven.priors <- list(alpha.normal = list(mean = 0, var = 2.72), 
                          beta.normal = list(mean = 0, var = 2.72))
      
      n.samples <- 5000
      n.burn <- 3000
      n.thin <- 2
      n.chains <- 3
      
      ################################################################################
      # 3. MODEL SELECTION
      ################################################################################
      
      ################################
      ##### DETECTION COVARIATES #####
      ################################
      
      # Set total potential covariates
      occ.form <- c("scale(ann)", "scale(hot)", "scale(col)", "scale(wet)", "scale(dry)")
      det.form <- c("scale(outcrop)", "scale(MGVF)", "scale(rain)", "factor(land)", 
                    "scale(occur)", "scale(sedflux)", "scale(Distance)", "scale(coll)")
      occ.form <- unlist(lapply(1:length(occ.form), 
                                function(x) combn(occ.form, x, simplify = FALSE)), 
                         recursive = FALSE)
      det.form <- unlist(lapply(1:length(det.form), 
                                function(x) combn(det.form, x, simplify = FALSE)), 
                         recursive = FALSE)
      occ.form <- append(occ.form, "1")
      det.form <- append(det.form, "1")
      
      det.wrapper <- function(det.form){
        revi.sp.det.formula <- formula(paste("~", paste(det.form, collapse = " + ")))
        revi.sp.occ.formula <- ~ scale(col) + scale(ann) + scale(wet) + scale(dry) + scale(hot)
        print(paste("Running model ", det.form, " out of ", length(det.form), sep = ""))
        out <- spOccupancy::PGOcc(occ.formula = revi.sp.occ.formula, 
                                       det.formula = revi.sp.det.formula, 
                                       data = sp.data, 
                                       inits = oven.inits, 
                                       n.samples = n.samples, 
                                       priors = oven.priors, 
                                       n.omp.threads = 1, 
                                       verbose = TRUE, 
                                       n.report = 1000, 
                                       n.burn = n.burn, 
                                       n.thin = n.thin, 
                                       n.chains = n.chains)

        det_waic <- t(as.data.frame(spOccupancy::waicOcc(out)))
        rownames(det_waic) <- c(paste("det ~", revi.sp.det.formula)[2])
        return(det_waic)
      }
      
      # Initiate the cluster
      sfInit(parallel = TRUE, cpus = 4)
      
      # Export data to the cluster
      sfExport('sp.data')
      sfExport('oven.inits')
      sfExport('oven.priors')
      sfExport('n.samples')
      sfExport('n.burn')
      sfExport('n.thin')
      sfExport('n.chains')
      
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
        out <- spOccupancy::PGOcc(occ.formula = revi.sp.occ.formula, 
                                       det.formula = revi.sp.det.formula, 
                                       data = sp.data, 
                                     inits = oven.inits, 
                                     n.samples = n.samples, 
                                     priors = oven.priors, 
                                     n.omp.threads = 1, 
                                     verbose = TRUE, 
                                     n.report = 1000, 
                                     n.burn = n.burn, 
                                     n.thin = n.thin, 
                                     n.chains = n.chains)
        
        occ_waic <- t(as.data.frame(spOccupancy::waicOcc(out)))
        rownames(occ_waic) <- c(paste("occ ~", revi.sp.occ.formula)[2])
        return(occ_waic)
      }
      
      # Initiate the cluster
      sfInit(parallel = TRUE, cpus = 4)
      
      # Export data to the cluster
      sfExport('sp.data')
      sfExport('best.det')
      sfExport('oven.inits')
      sfExport('oven.priors')
      sfExport('n.samples')
      sfExport('n.burn')
      sfExport('n.thin')
      sfExport('n.chains')
      
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
      
      revi.sp.det.formula <- ~ scale(Distance) + scale(outcrop) + (1 | Site)
      revi.sp.occ.formula <- ~ scale(col) + scale(ann) + scale(wet) + scale(dry) + scale(hot)
      
      # Covariates
      revi.sp.det.formula <- best.det
      revi.sp.occ.formula <- best.occ
      
      # Best Model
      out <- spOccupancy::PGOcc(occ.formula = revi.sp.occ.formula, 
                                   det.formula = revi.sp.det.formula, 
                                   data = sp.data, 
                                   inits = oven.inits, 
                                   n.samples = n.samples, 
                                   priors = oven.priors, 
                                   n.omp.threads = 1, 
                                   verbose = TRUE, 
                                   n.report = 1000, 
                                   n.burn = n.burn, 
                                   n.thin = n.thin, 
                                   n.chains = n.chains)
      
      # Summary
      summary(out)
      
      # Goodness-of-fit across space
      GOFspace <- spOccupancy::ppcOcc(out, fit.stat = 'freeman-tukey', group = 1)
      summary(GOFspace)
      # Goodness-of-fit across replicates
      GOFrep <- spOccupancy::ppcOcc(out, fit.stat = 'freeman-tukey', group = 2)
      summary(GOFrep)
      
      # Occupancy covariates
      MCMCvis::MCMCplot(out$beta.samples, ref_ovl = TRUE, ci = c(50, 95))
      # Detection covariates
      MCMCvis::MCMCplot(out$alpha.samples, ref_ovl = TRUE, ci = c(50, 95))
      
      
      saveRDS(out, file = paste("Results/spOccupancy/single_season/",
                                   res, "/", bin, ".", target, ".formcells.best.model.rds", 
                                   sep = ""))
      saveRDS(GOFspace, file = paste("Results/spOccupancy/single_season/",
                                     res, "/", bin, ".", target, ".formcells.gof.space.rds", 
                                     sep = ""))
      saveRDS(GOFrep, file = paste("Results/spOccupancy/single_season/",
                                   res, "/", bin, ".", target, ".formcells.gof.rep.rds", 
                                   sep = ""))
      
    }
  }
}
