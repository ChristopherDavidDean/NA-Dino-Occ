################################################################################
############ OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS #############
################################################################################

# Christopher D. Dean, Lewis A. Jones, Alfio A. Chiarenza, Sinéad Lyster, 
# Alex Farnsworth, Philip D. Mannion, Richard J. Butler.
# 2023
# Script written by Christopher D. Dean

################################################################################
#                    FILE 6: RUNNING MULTI-SEASON OCCUPANCY                    #
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
res <- 0.5
# Set extent
e <- extent(-155, -72, 22.5, 73)
# Set max limit value
max_val <- 40
max_val_on <- TRUE
bin.type <- "scotese"
target <- "Tyrannosauridae"
bins <- read.csv("Data/Occurrences/scotesebins.csv")
bins$bin <- bins$code

# Load occurrence dataset
sp.data <- readRDS(file = paste("Prepped_data/spOccupancy/Multi_season/", res, "/", 
                    target, "_multi_", res, ".rds", sep = ""))

################################################################################
# 2. SETUP MODELLING
################################################################################

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
occ.form <- c("scale(ann)", "scale(hot)", "scale(col)", "scale(wet)", "scale(dry)")
det.form <- c("scale(outcrop)", "scale(MGVF)", "scale(rain)", "factor(Year)", "factor(land)")
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
  revi.sp.det.formula <- formula(paste("~", paste(det.form, collapse = " + ")))
  revi.sp.occ.formula <- ~ scale(col) + scale(ann) + scale(wet) + scale(dry) + scale(hot)
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
  rownames(occ_waic) <- c(paste("occ ~", revi.sp.det.formula)[2])
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

best.occ <- as.formula(paste(gsub("occ", "", occ_waic[which.min(occ_waic$WAIC), 4]), sep = " "))

#######################
##### BEST MODELS #####
#######################

# 0.1
# Hadrosaur Best Model

# Ceratopsidae Best model

# Tyrannosauridae Best model

  
# 0.5
# Hadrosaur Best Model
#revi.sp.occ.formula <- ~ scale(col) + scale(wet) + scale(dry)
#revi.sp.det.formula <- ~ scale(rain) + factor(Year) + (1 | Site)

# Ceratopsidae Best model
#revi.sp.occ.formula <- ~ scale(ann)
#revi.sp.det.formula <- ~ scale(outcrop) + factor(land) + factor(Year) + (1 | Site)

# Tyrannosauridae Best model


# 1 
# Hadrosaur Best Model

# Ceratopsidae Best model
 
# Tyrannosauridae Best model
 
  
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

saveRDS(make.table(out.sp, target), file = paste("Results/spOccupancy/",
                                       res, "/", target, ".cov.table.rds", 
                                       sep = ""))
saveRDS(out.sp, file = paste("Results/spOccupancy/",
                                                 res, "/", target, ".best.model.rds", 
                                                 sep = ""))
saveRDS(GOFspace, file = paste("Results/spOccupancy/",
                                                 res, "/", target, ".gof.space.rds", 
                                                 sep = ""))
saveRDS(GOFrep, file = paste("Results/spOccupancy/",
                                                 res, "/", target, ".gof.rep.rds", 
                                                 sep = ""))

################################################################################
# 5. PLOTTING SITE-LEVEL DETECTION PROBABILITY RANDOM EFFECTS
################################################################################

# Find random effects sizes
alpha.star.means <- apply(out.sp$alpha.star.samples, 2, mean)

# Create siteIDs (individual cell numbers used)
siteIDs <- as.numeric(rownames(sp.data$y[,1,]))

# Find coordinates of each site
siteCoords <- siteCoordsFun(res = res, e = e, siteIDs)

# Generate a raster using random effects
gen_raster(siteCoords$siteID, alpha.star.means, res = res, ext = e)

# Find map to use as backdrop
countries <- maps::map("world", plot=FALSE, fill = TRUE) 
# Turn map into spatialpolygons
countries <<- maptools::map2SpatialPolygons(countries, 
                                            IDs = countries$names, 
                                            proj4string = CRS("+proj=longlat")) 
mapTheme <- rasterVis::rasterTheme(region=brewer.pal(8,"Reds"))

(p2 <- rasterVis::levelplot(raster_for_values, margin=T, par.settings=mapTheme) + 
        # Plots state lines
        latticeExtra::layer(sp.polygons(states, col = "white", fill = NA), under = T)  + 
        # Plots background colour
        latticeExtra::layer(sp.polygons(countries, col = 0, fill = "light grey"), under = T))

save.lattice(p2)

################################################################################
# A1. ADDITIONAL - NON-SPATIAL MODEL
################################################################################

# Prep the array for the individual target chosen
sp.data <- Array_prep(target, sp = FALSE)

# Quick plot to check occupancy
raw.occ.prob <- apply(sp.data$y, 2, mean, na.rm = TRUE)
plot(1:4, raw.occ.prob, pch = 16, 
     xlab = 'Year', ylab = 'Raw Occurrence Proportion', 
     cex = 1.5, frame = FALSE, ylim = c(0, 1))

# Specify formula
revi.occ.formula <- ~ scale(ann) + scale(dry)
revi.det.formula <- ~ scale(outcrop) + factor(land) 

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

# Run model
out <- tPGOcc(occ.formula = revi.occ.formula, 
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
# Model report/waic
summary(out)
waicOcc(out)

# Goodness-of-fit test
ppcOut <- ppcOcc(out, fit.stat = 'freeman-tukey', group = 1)
summary(ppcOut)
ppcOut <- ppcOcc(out, fit.stat = 'freeman-tukey', group = 2)
summary(ppcOut)

# Occupancy covariates
MCMCplot(out$beta.samples, ref_ovl = TRUE, ci = c(50, 95))
# Detection covariates
MCMCplot(out$alpha.samples, ref_ovl = TRUE, ci = c(50, 95))