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

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

##### Load Packages #####
library(RPresence)
library(dplyr)
library(AICcmodavg)

##### Load Data ######

# Quick filters for loading data
bin.type <- "stage"
res <- 1
bin <- "S.1"
stage <- "Maas" # Remember to change this along with bin
target <- "Ceratopsidae"
ss <- 5

# Load occupancy data
data <- read.csv(paste("Results/", bin.type, "/", bin, "/", res, "/", bin, ".", 
                       res, ".", target, "SS.", ss ,".csv", sep = "")) 
# Load site occupancy covariates
site.occ <- read.csv(paste("Results/", bin.type, "/", bin, "/", res, "/", 
                           "site_occupancy_covs.csv", sep = ""))
# Load site detection covariates
site.det <- read.csv(paste("Results/", bin.type, "/", bin, "/", res, "/", 
                           "site_detection_covs.csv", sep = ""))
# Combine
site.covs <- cbind(site.occ[,2:ncol(site.occ)], site.det[,3:ncol(site.det)])
site.covs <- site.covs %>%
  `rownames<-`(.[,1]) %>% 
  select(-cells)

# Setup for RPresence
eh <- data %>%
  `rownames<-`(.[,1]) %>% 
  select(-X)

#===== Create Pao =====
temp_name <- paste("Single Season (",  bin.type, ", ", bin, ") ", 
                   target, ", ", res, " degrees resolution", sep = "")
pao <- createPao(data = eh, 
                        unitcov = site.covs, 
                        unitnames = rownames(eh), 
                        title = temp_name)


################################################################################
# 2. RUNNING MODELS
################################################################################

##### Model Setup #####
test <- occMod(
  data = test_pao, 
  model = list(psi ~ mean_max_ptemp + mean_max_pprec, p ~ mean_prec), 
  type = 'so'
)

##### Goodness of Fit Test on full model #####


##### Check optimisation #####
summary(test)

##### Get beta coefficients #####
coef(object = test, 
     param = 'p', 
     prob = 0.95)
coef(object = test, 
     param = 'psi', 
     prob = 0.95)

##### Transform back to probability #####

