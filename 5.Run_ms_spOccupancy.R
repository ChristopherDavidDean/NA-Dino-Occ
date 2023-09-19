################################################################################
############ OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS #############
################################################################################

# Christopher D. Dean, Lewis A. Jones, Alfio A. Chiarenza, Sinéad Lyster, 
# Alex Farnsworth, Philip D. Mannion, Richard J. Butler.
# 2023
# Script written by Christopher D. Dean and Lewis A. Jones

################################################################################
#                    FILE 5: RUNNING MULTI-SEASON OCCUPANCY                    #
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
max_val <- 30
max_val_on <- TRUE
bin.type <- "scotese"
target <- "Hadrosauridae"
bins <- read.csv("Data/Occurrences/scotesebins.csv")
bins$bin <- bins$code

# Load occurrence dataset
sp.data <- readRDS(file = paste("Prepped_data/spOccupancy/Multi_season/", res, "/", 
                    target, "_multi_", res, ".rds", sep = ""))

################################################################################
# 2. SETUP MODELLING
################################################################################

sp.data$det.covs$Year <- sp.data$occ.covs$Year

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
# Total number of samples
n.batch * batch.length
n.burn <- 10000
n.thin <- 10 

###########################
##### MODEL SELECTION #####
###########################

# Set total potential covariates
occ.form <- c("scale(ann)", "scale(hot)", "scale(col)", "scale(wet)", "scale(dry)")
det.form <- c("scale(outcrop)", "scale(MGVF)", "scale(rain)", "factor(Year)", "factor(land)")
occ.form <- unlist(lapply(1:length(occ.form), 
                          function(x) combn(occ.form, x, simplify = FALSE)), 
                   recursive = FALSE)
det.form <- unlist(lapply(1:length(det.form), 
                          function(x) combn(det.form, x, simplify = FALSE)), 
                   recursive = FALSE)

# Make det WAIC dataframe
det_waic <- data.frame(model = NA, elpd = NA, dP = NA, WAIC = NA)

# Hold Occ, find det
for(d in 1:length(det.form)){
  revi.sp.det.formula <- formula(paste("~", paste(det.form[[d]], collapse = " + ")))
  revi.sp.occ.formula <- ~ 1
  # Approx. run time: ~ 2.5 min
  print(paste("Running model ", d, " out of ", length(det.form), sep = ""))
  out.sp <- stPGOcc(occ.formula = revi.sp.occ.formula, 
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
  det_waic <- rbind(det_waic, c(paste("det ~", revi.sp.det.formula)[2], waicOcc(out.sp)))
}

# Set det.
best.det <- ~ scale(MGVF) + factor(Year) - 1 + (1 | Site)

# Hold det, find occ
occ_waic <- data.frame(model = NA, elpd = NA, dP = NA, WAIC = NA)

# Hold det, find occ
for(d in 1:length(occ.form)){
  revi.sp.occ.formula <- formula(paste("~", paste(occ.form[[d]], collapse = " + ")))
  revi.sp.det.formula <- best.det
  # Approx. run time: ~ 2.5 min
  print(paste("Running model ", d, " out of ", length(occ.form), sep = ""))
  out.sp <- stPGOcc(occ.formula = revi.sp.occ.formula, 
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
  occ_waic <- rbind(occ_waic, c(paste(" ~", revi.sp.occ.formula)[2], waicOcc(out.sp)))
}

#######################
##### BEST MODELS #####
#######################

##### 0.1 #####
# Hadrosaur Best Model

# Ceratopsidae Best model

# Tyrannosauridae Best model


##### 0.5 #####
# Hadrosaur Best Model
revi.sp.occ.formula <- ~ scale(dry) 
revi.sp.det.formula <- ~ scale(MGVF) + factor(Year) - 1 + (1 | Site)
# Ceratopsidae Best model
revi.sp.occ.formula <- ~ scale(hot) + scale(wet) + scale(dry) 
revi.sp.det.formula <- ~ scale(outcrop) + factor(land) - 1 + factor(Year) - 1 + (1 | Site)
# Tyrannosauridae Best model


##### 1 #####
# Hadrosaur Best Model

# Ceratopsidae Best model

# Tyrannosauridae Best model

##################################
##### RESULTS FOR BEST MODEL #####
##################################

# Running best model
out.sp <- stPGOcc(occ.formula = revi.sp.occ.formula, 
                  det.formula = revi.sp.det.formula, 
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
ppcOut <- ppcOcc(out.sp, fit.stat = 'freeman-tukey', group = 1)
summary(ppcOut)
# Goodness-of-fit across replicates
ppcOut <- ppcOcc(out.sp, fit.stat = 'freeman-tukey', group = 2)
summary(ppcOut)

# Occupancy
MCMCplot(out.sp$beta.samples, ref_ovl = TRUE, ci = c(50, 95))
# Detection
MCMCplot(out.sp$alpha.samples, ref_ovl = TRUE, ci = c(50, 95))





###### NON SPATIAL ######
# Prep the array for the individual target chosen
sp.data <- Array_prep(target, sp = FALSE)

# Quick plot to check occupancy
raw.occ.prob <- apply(sp.data$y, 2, mean, na.rm = TRUE)
plot(1:4, raw.occ.prob, pch = 16, 
     xlab = 'Year', ylab = 'Raw Occurrence Proportion', 
     cex = 1.5, frame = FALSE, ylim = c(0, 1))

# Specify formula
revi.occ.formula <- ~ scale(col) + scale(hot) + (1 | site.effect)

revi.det.formula <- ~ scale(outcrop)

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

