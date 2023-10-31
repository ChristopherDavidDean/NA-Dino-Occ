################################################################################
############ OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS #############
################################################################################

# Christopher D. Dean, Lewis A. Jones, Alfio A. Chiarenza, Sinéad Lyster, 
# Alex Farnsworth, Philip D. Mannion, Richard J. Butler.
# 2023
# Script written by Christopher D. Dean

################################################################################
#               FILE 3: RUNNING OCCUPANCY MODELS IN RPRESENCE                #
################################################################################

################################################################################
# 1. INITIAL SETUP
################################################################################

#######################################
##### PACKAGES AND VARIABLE SETUP #####
#######################################

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

# Load Packages
library(RPresence)
library(dplyr)
library(AICcmodavg)

# Quick filters for loading data
bin.type <- "scotese"
res <- 1
bin <- "teyen"
target <- "Ceratopsidae"

###################################
##### LOAD AND SORT VARIABLES #####
###################################

##### OCCUPANCY DATA #####
# Load data
data <- read.csv(paste("Prepped_data/Occurrence_Data/", bin.type, "/", bin, "/", 
                       res, "/", bin, ".", res, ".", target, ".dframe.10.csv", 
                       sep = "")) 
# Setup check
data.check <- data %>%
  select(X) %>%
  distinct()

##### SITE COVARIATES #####
# Load site covariates
sitecov <- read.csv(paste("Prepped_data/Occurrence_Data/", bin.type, "/", bin, 
                          "/", res, "/site_occupancy_covs.csv", sep = "")) 
# Load precise site covariates
precise.sitecov <- read.csv(paste("Prepped_data/Occurrence_Data/", bin.type, "/", bin, 
                          "/", res, "/precise_mean_covs.csv", sep = "")) 
precise.sitecov <- precise.sitecov %>% 
  filter(counting_colls.Coll_count > 1)

# Setup checks
site.check <- sitecov %>%
  select(siteID) %>%
  distinct()
precise.check <- precise.sitecov %>%
  select(siteID) %>%
  distinct()

##### SURVEY COVARIATES #####
# Load survey covariates
survcov <- read.csv(paste("Prepped_data/Occurrence_Data/", bin.type, "/", bin, 
                          "/", res, "/surv_covs.csv", sep = "")) 
surv.check <- survcov %>%
  select(siteID) %>%
  distinct()

##### QUICK CHECKS #####
# Do all siteIDs match? Are they all in the same order?
identical(precise.check$siteID, surv.check$siteID)
identical(precise.check$siteID, site.check$siteID)
identical(data.check$X, precise.check$siteID)

##### FINAL SETUP AFTER CHECKS #####
# Setup Encounter histories for RPresence
eh <- data %>%
  `rownames<-`(.[,1]) %>% 
  select(-X)

# Setup survey covariates
survcov <- survcov %>%
  select(colls, temp, prec, DEM)
names(survcov)[names(survcov) == "DEM"] <- "DEM_surv"

# Setup site covariates
precise.sitecov <- precise.sitecov %>%
  select(siteID, mean_DEM, mean_prec, mean_temp)
sitecov <- merge(precise.sitecov, sitecov, by = "siteID")
sitecov <- select(sitecov, -c(siteID, X))

# Remove all digits at end of the column names
colnames(sitecov) <- sub("_\\d.*", "", colnames(sitecov))

# Set categorical variables
sitecov$LANDCVI_multiple <- as.character(sitecov$LANDCVI_multiple)
sitecov$LANDCVI_binary <- as.character(sitecov$LANDCVI_binary)

##### CREATE PAO #####
# Set name
temp_name <- paste("Single Season (",  bin.type, ", ", bin, ") ", 
                   target, ", ", res, " degrees resolution", sep = "")
# Make PAO
pao <- createPao(data = eh, 
                 unitcov = sitecov, 
                 survcov = survcov,
                 unitnames = rownames(eh), 
                 title = temp_name)

################################################################################
# 2. RUNNING MODELS
################################################################################

##### Goodness of Fit Test on full model #####
full.mod <- occMod(
  data = pao, model = list(
    psi ~ mean_ptemp + mean_pprec + mean_pDEM, 
    p ~ mean_psed + LANDCVI_multiple + AllOut + MGVF + WC_Prec + WC_Temp + mean_prec +
      mean_temp + mean_DEM + MOut + temp + prec), 
  modfitboot = 5000,
  type = 'so'
)
full.mod$gof

###########################
##### MODEL SELECTION #####
###########################

##### OCCUPANCY #####
m0 <- occMod(
  data = pao, model = list(
    psi ~ 1, 
    p ~ 1), 
  type = 'so'
)
m1 <- occMod(
  data = pao, model = list(
    psi ~ mean_ptemp, 
    p ~ 1), 
  type = 'so'
)
m2 <- occMod(
  data = pao, model = list(
    psi ~ mean_pprec, 
    p ~ 1), 
  type = 'so'
)
m3 <- occMod(
  data = pao, model = list(
    psi ~ mean_pDEM, 
    p ~ 1), 
  type = 'so'
)
m4 <- occMod(
  data = pao, model = list(
    psi ~ mean_ptemp + mean_pprec, 
    p ~ 1), 
  type = 'so'
)
m5 <- occMod(
  data = pao, model = list(
    psi ~ mean_ptemp + mean_pDEM, 
    p ~ 1), 
  type = 'so'
)
m6 <- occMod(
  data = pao, model = list(
    psi ~ mean_pprec + mean_pDEM, 
    p ~ 1), 
  type = 'so'
)
m7 <- occMod(
  data = pao, model = list(
    psi ~ mean_ptemp + mean_pprec + mean_pDEM, 
    p ~ 1), 
  type = 'so'
)
m8 <- occMod(
  data = pao, model = list(
    psi ~ mean_pprec +  I(mean_pprec^2), 
    p ~ 1), 
  type = 'so'
)
m9 <- occMod(
  data = pao, model = list(
    psi ~ mean_ptemp +  I(mean_ptemp^2), 
    p ~ 1), 
  type = 'so'
)

# Make/show AICc table
aic <- createAicTable(
  list(m0, m1, m2, m3, m4, m5, m6, m7, m8, m9),
  use.aicc = TRUE
)
aic$table

# Best Model
bm.o <- m2

##### DETECTION #####

m0 <- occMod(
  data = pao, 
  model = list(
    psi ~ 1, 
    p ~ 1),
  type = 'so'
)
m1 <- occMod(
  data = pao, 
  model = list(
    psi ~ mean_pprec, 
    p ~ LANDCVI_multiple),
  type = 'so'
)
m2 <- occMod(
  data = pao, 
  model = list(
    psi ~ mean_pprec, 
    p ~ AllOut),
  type = 'so'
)
m3 <- occMod(
  data = pao, 
  model = list(
    psi ~ mean_pprec, 
    p ~ MGVF),
  type = 'so'
)
m4 <- occMod(
  data = pao, 
  model = list(
    psi ~ mean_pprec, 
    p ~ WC_Prec),
  type = 'so'
)
m5 <- occMod(
  data = pao, 
  model = list(
    psi ~ mean_pprec, 
    p ~ WC_Temp),
  type = 'so'
)
m6 <- occMod(
  data = pao, 
  model = list(
    psi ~ mean_pprec, 
    p ~ mean_psed),
  type = 'so'
)
m7 <- occMod(
  data = pao, 
  model = list(
    psi ~ mean_pprec, 
    p ~ LANDCVI_multiple + AllOut),
  type = 'so'
)
m8 <- occMod(
  data = pao, 
  model = list(
    psi ~ mean_pprec, 
    p ~ AllOut + WC_Prec),
  type = 'so'
)
m9 <- occMod(
  data = pao, 
  model = list(
    psi ~ mean_pprec, 
    p ~ mean_psed + AllOut),
  type = 'so'
)
m10 <- occMod(
  data = pao, 
  model = list(
    psi ~ 1, 
    p ~ mean_psed + LANDCVI_multiple + mean_prec + mean_psed),
  type = 'so'
)


m10

# Make/show AICc table
aic <- createAicTable(
  list(m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10),
  use.aicc = TRUE
)
aic$table

# Best Model
bm.od <- m10

##### Check optimisation #####
summary(bm.od)

##### Get beta coefficients #####
coef(object = bm.od, 
     param = 'psi', 
     prob = 0.95)
coef(object = bm.od, 
     param = 'p', 
     prob = 0.95)

##### Transform back to probability #####

