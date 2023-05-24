################################################################################
############ OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS #############
################################################################################

# Christopher D. Dean, Lewis A. Jones, Alfio A. Chiarenza, Sinéad Lyster, 
# Alex Farnsworth, Philip D. Mannion, Richard J. Butler.
# 2019
# Script written by Christopher D. Dean

################################################################################
#               FILE 2.5: RUNNING OCCUPANCY MODELS IN RPRESENCE                #
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

##### Model Setup #####
test <- occMod(
  data = pao, model = list(
    psi ~ mean_ptemp + mean_pprec + mean_pDEM, 
    p ~ mean_psed + LANDCVI_multiple + AllOut + MGVF + WC_Prec), 
  type = 'so'
)

##### Goodness of Fit Test on full model #####

##### Check optimisation #####
summary(test)

##### Get beta coefficients #####
coef(object = test, 
     param = 'psi', 
     prob = 0.95)
coef(object = test, 
     param = 'p', 
     prob = 0.95)

##### Transform back to probability #####

