################################################################################
############ OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS #############
################################################################################

# Christopher D. Dean, Lewis A. Jones, Alfio A. Chiarenza, Sinéad Lyster, 
# Alex Farnsworth, Philip D. Mannion, Richard J. Butler.
# 2023
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
res <- 0.5
# Set extent
e <- extent(-155, -72, 22.5, 73)
bin.type <- "scotese"
target <- "Hadrosauridae"
bins <- read.csv("Data/Occurrences/scotesebins.csv")
bins$bin <- bins$code
bin <- "teyen" # Change this if you want a single season model

# Load occurrence dataset
sp.data <- readRDS(file = paste("Prepped_data/spOccupancy/Single_season/", res, "/", 
                                bin, "/", target, "_single_", res, ".rds", sep = ""))

################################################################################
# 2. SETUP MODELLING
################################################################################

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

# Adjustments to saved data (whoops...)
sp.data[[1]] <- as.matrix(sp.data[[1]])

###########################
##### MODEL SELECTION #####
###########################

# Set total potential covariates
occ.form <- c("scale(ann)", "scale(hot)", "scale(col)", "scale(wet)", "scale(dry)")
det.form <- c("scale(outcrop)", "scale(MGVF)", "scale(rain)", "factor(land)")
occ.form <- unlist(lapply(1:length(occ.form), 
                          function(x) combn(occ.form, x, simplify = FALSE)), 
                   recursive = FALSE)
det.form <- unlist(lapply(1:length(det.form), 
                          function(x) combn(det.form, x, simplify = FALSE)), 
                   recursive = FALSE)
occ.form <- append(occ.form, "1")
det.form <- append(det.form, "1")

# Make det WAIC dataframe
det_waic <- data.frame(model = NA, elpd = NA, dP = NA, WAIC = NA)

# Hold Occ, find det
for(d in 1:length(det.form)){
  revi.sp.det.formula <- formula(paste("~", paste(det.form[[d]], collapse = " + ")))
  revi.sp.occ.formula <- ~ scale(col) + scale(ann) + scale(wet) + scale(dry) + scale(hot)
  # Approx. run time: ~ 2.5 min
  print(paste("Running model ", d, " out of ", length(det.form), sep = ""))
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
  
  det_waic <- rbind(det_waic, c(paste("det ~", revi.sp.det.formula)[2], spOccupancy::waicOcc(out.sp)))
}

# Set det.
best.det <- ~ factor(outcrop) + (1 | Site)

# Hold det, find occ
occ_waic <- data.frame(model = NA, elpd = NA, dP = NA, WAIC = NA)

# Hold det, find occ
for(d in 1:length(occ.form)){
  revi.sp.occ.formula <- formula(paste("~", paste(occ.form[[d]], collapse = " + ")))
  revi.sp.det.formula <- best.det
  # Approx. run time: ~ 2.5 min
  print(paste("Running model ", d, " out of ", length(occ.form), sep = ""))
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
  occ_waic <- rbind(occ_waic, c(paste(" ~", revi.sp.occ.formula)[2], spOccupancy::waicOcc(out.sp)))
}

################################################################################
# 3. BEST MODEL AND OTHER INFO
################################################################################

# Covariates
revi.sp.det.formula <- ~ scale(outcrop) + (1 | Site)
revi.sp.occ.formula <- ~ scale(hot) + scale(ann) 

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

# Goodness-of-fit tests
ppc.sp.out <- spOccupancy::ppcOcc(out.sp, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.sp.out)
ppc.sp.out <- spOccupancy::ppcOcc(out.sp, fit.stat = 'freeman-tukey', group = 2)
summary(ppc.sp.out)

# Occupancy covariates
MCMCvis::MCMCplot(out.sp$beta.samples, ref_ovl = TRUE, ci = c(50, 95))
# Detection covariates
MCMCvis::MCMCplot(out.sp$alpha.samples, ref_ovl = TRUE, ci = c(50, 95))

