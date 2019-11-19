#===============================================================================================================================================
#============================================== OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS ==========================================
#===============================================================================================================================================

# Christopher D. Dean, Lewis A. Jones, Alfio A. Chiarenza, Alex Farnsworth, María V. Jiménez‐Franco, Richard J. Butler.
# 2019
# Script written by Christopher D. Dean and Lewis A. Jones

# Broad PBDB download: https://paleobiodb.org/data1.2/occs/list.csv?base_name=tetrapoda&taxon_reso=genus&pres=regular&interval=Campanian,Maastrichtian&cc=NOA&envtype=terr&show=full,strat,lith,env

#=============================================== INITIAL SETUP ===============================================

# Set working directory
setwd("C:/Users/deancd/Documents/RESEARCH/PROJECTS/DINO_RANGE/NA-Dino-Occ/") # Set your working directory

#==== If first time using: ====
# Make vector of package names
packages <- c("beepr", "raster", "dplyr", "lattice", "rasterVis", "sp", "maps", "maptools", "parallel") #list your packages here

# Install packages
ipak(packages)

# Load Required packages
library(beepr)
library(raster)
library(dplyr)
library(lattice)
library(rasterVis)
library(sp)
library(maps)
library(maptools)
library(parallel)

#==== if not: ====
# Load in Functions
source("0.Functions_DR.R") # Import functions from other R file (must be in same working directory)

#=============================================== DATA SETUP ===============================================

# Read in files
camp.occs <- read.csv("Data/Occurrences/Broad_Colls_Data_Campanian.csv", stringsAsFactors = FALSE) # Load in occurrences

# Set working resolution and extent. Note: these values should be the same ones used for the file 1.Setup_occupancy_DR.
res <- 0.5
get_extent(camp.occs)

#==== Visualise grid cells for collections and occurrences ====

camp.colls <- camp.occs %>% # Make unique collections for visualisation
  dplyr::select(collection_no, lat, lng) %>%
  dplyr::distinct()
maas.colls <- maas.occs %>% # Make unique collections for visualisation
  dplyr::select(collection_no, lat, lng) %>%
  dplyr::distinct()

get_grid_im(camp.occs, res, "Collections")
get_grid_im(maas.occs, res, "Collections")
get_grid_im(camp.colls, res, "Collections")
get_grid_im(maas.colls, res, "Collections")

#==== Testing Targets ====
# target_maker(camp.occs, "family", "Ceratopsidae")
# Ceratops <- camp.occs.targeted %>%
#   filter(Target == "Ceratopsidae")
# get_grid_im(Ceratops, 1, "Ceratopsidae Occurrences")


#=============================================== OCCUPANCY SETUP ===============================================

# Set Target taxa
target = c("Ceratopsidae", "Tyrannosauridae", "Hadrosauridae") # set Targets
target_maker(camp.occs, "family", target) # run target_maker

# Check resolution data
res_data(camp.occs.targeted, target, 0.1, 1, 0.1) # Makes list of dataframes, each containing information about various cell resolutions. Can also see Naive occupancy estimates.
Res_results_list$Ceratopsidae
Res_results_list$Tyrannosauridae
Res_results_list$Hadrosauridae

# Prepare data for unmarked - Campanian
all_results_for_unmarked(data = camp.occs.targeted, res = res, ext = e, target = target, subsamp = TRUE)

#=============================================== COVARIATE SETUP ===============================================

# Convert rasters to desired resolution
source("0.Format_data_DR.R") # Import extent and run cleaning/import of raster datasets. Running this will take a while, 
# but will automatically update rasters so that they are of the desired resolution and place them in the appropriate folder
# for running the next steps. 

# Add covariates to Occurrence spreadsheet
get_cov_from_stack(Final, res = res)

# Clean/split data to just relevant covariates
Camp_Covs <- cov_dat[, -grep("Maas_out_", colnames(cov_dat))]
#Maas_Covs <- cov_dat[, -grep("Camp_out_", colnames(cov_dat))]

#===== PALAEO-DATA =====
# Data
CampPrecip <- raster("Lewis_Occupancy_data/Data/Formatted/PalaeoClimate/CampPrecip.asc")
CampTemp <- raster("Lewis_Occupancy_data/Data/Formatted/PalaeoClimate/CampTemp.asc")
#MaasPrecip <- raster("Lewis_Occupancy_data/Data/Formatted/PalaeoClimate/MaasPrecip.asc")
#MaasTemp <- raster("Lewis_Occupancy_data/Data/Formatted/PalaeoClimate/MaasTemp.asc")

# Camp Stack Extract
CampStack <- stack(CampPrecip, CampTemp)
NAcells <- Camp_Covs %>% # Clean for NA values in paleolat/lng - Record cells with NAs
  filter(is.na(paleolat)) %>%
  distinct(cells.1) # Remove cells with NAs
Camp_Covs <- Camp_Covs %>%
  filter(!is.na(paleolat))
xy <- SpatialPointsDataFrame(cbind.data.frame(Camp_Covs$paleolng, Camp_Covs$paleolat), Camp_Covs) # Use palaeo lat/long to extract from palaeo GCMs
Camp_Covs <- extract(CampStack, xy, sp = TRUE, cellnumbers = TRUE) # Extract Palaeo Temp and Rainfall
Camp_Covs <- as.data.frame(Camp_Covs) # Conver to dataframe

# Make Master Camp Spreadsheet
Camp_Cov_Master <- Camp_Covs
# Clean data to just gridcells and associated covariates
Camp_Covs <- Camp_Covs[,c(129:147, 150:152)]
Camp_Covs <- Camp_Covs %>%
  distinct()
Camp_Covs <- aggregate(x = Camp_Covs, by = list(Camp_Covs$cells.1), FUN = "mean")
Camp_Covs <- Camp_Covs[, 2:23]
Camp_Covs[,22] <- Camp_Covs[,22] - 273.15 # Convert from Kelvin to degrees Celsius

# Maas Stack Extract
#MaasStack <- stack(MaasPrecip, MaasTemp)
#NAcells <- cov_dat %>% # Clean for NA values in paleolat/lng - Record cells with NAs
#  filter(is.na(paleolat)) %>%
#  distinct(cells.1)
#cov_dat_complete<- cov_dat %>%
#  filter(!is.na(paleolat))
#xy <- SpatialPointsDataFrame(cbind.data.frame(cov_dat_complete$paleolng, cov_dat_complete$paleolat), cov_dat_complete)
#Maas_Dat <- extract(MaasStack, xy, sp = TRUE, cellnumbers = TRUE)
#Maas_Dat <- as.data.frame(Maas_Dat)

#=========================================================== OCCUPANCY =================================================================

library("unmarked")

data <- read.csv("Results/Subsampled/maas.occs.targeted.1.Hadrosauridae.csv") # import data from previous step. Result here is just a test case to show.
# Turn into matrix
y <- as.matrix(data[,2:6])

umf <- unmarkedFrameOccu(y = y)
summary(umf)
summary(fm1 <- occu(~1 ~1, data=umf))

print(occ.null <- backTransform(fm1, "state")) 
print(det.null <- backTransform(fm1, "det")) 


lith.orig <- read.csv("Results/sitecovs.1.csv")
site.covs <- read.csv("Results/cellcovs.1.csv")
site.covs <- site.covs %>%
  dplyr::arrange(cells)
  

OA.orig <- site.covs[,"Perc_outcrop_area"]      # Unstandardised, original values of covariates
Carb.orig <- site.covs[,"perc_carb"]
CPC.orig <- site.covs[,"colls_per_cell"]
lith <- as.matrix(lith.orig)

# Overview of Covariates
covs <- cbind(OA.orig, Carb.orig, CPC.orig)
par(mfrow = c(3,1))
for(i in 1:3){
  hist(covs[,i], breaks = 50, col = "grey", main = colnames(covs)[i])
}
pairs(cbind(OA.orig, Carb.orig, CPC.orig))

# Standardize covariates and mean-impute OA and Carbation
# Compute means and standard deviations
(means <- c(apply(cbind(OA.orig, Carb.orig, CPC.orig), 2, mean)))
(sds <- c(apply(cbind(OA.orig, Carb.orig, CPC.orig), 2, sd)))

# Scale covariates
OA <- (OA.orig - means[1]) / sds[1]
Carb <- (Carb.orig - means[2]) / sds[2]
CPC <- (CPC.orig - means[3]) / sds[3]

# Turn into matrix
y <- as.matrix(data[,2:231]) # set second number as number of sites (i.e. number of columns)

library(unmarked)
umf <- unmarkedFrameOccu(y = y, siteCovs = data.frame(OA = OA, Carb = Carb, CPC = CPC))
summary(umf)

# Fit null model
# Fit a series of models for detection first and do model selection
summary(fm1 <- occu(~1 ~1, data=umf))
summary(fm2 <- occu(~CPC ~1, data=umf))
summary(fm3 <- occu(~CPC+I(CPC^2) ~1, data=umf))
summary(fm4 <- occu(~Carb ~1, data=umf))
summary(fm5 <- occu(~OA+CPC ~1, data=umf))
summary(fm6 <- occu(~OA+I(OA^2)+Carb ~1, data=umf))
summary(fm7 <- occu(~Carb+I(Carb^2) ~1, data=umf))
summary(fm8 <- occu(~OA+Carb+I(Carb^2) ~1, data=umf))
summary(fm9 <- occu(~OA+I(OA^2)+Carb+I(Carb^2) ~1, data=umf))


# Put the fitted models in a "fitList" and rank them by AIC
fms <- fitList("p(.)psi(.)"                     = fm1,
               "p(OA)psi(.)"                      = fm2,
               "p(OA+OA2)psi(.)"                = fm3,
               "p(Carb)psi(.)"                       = fm4,
               "p(OA+Carb)psi(.)"                  = fm5,
               "p(OA+OA2+Carb)psi(.)"            = fm6,
               "p(Carb+Carb2)psi(.)"                  = fm7,
               "p(OA+Carb+Carb2)psi(.)"             = fm8,
               "p(OA+OA2+Carb+carb2)psi(.)"           = fm9)
(ms <- modSel(fms))

library(AICcmodavg)
system.time(gof.boot <- mb.gof.test(fm10, nsim = 1000))
gof.boot


# Create new covariates for prediction ('prediction covs')
pred.OA <- seq(0, 100,,131)    # New covs for prediction
pred.Carb <- seq(0, 100,,131)
p.OA <- (pred.OA - means[1]) / sds[1] # Standardise them like actual covs
p.Carb <- (pred.Carb - means[2]) / sds[2]

# Obtain predictions
newData <- data.frame(Outcrop=p.OA, Carb=0)
pred.occ.OA <- predict(fm9, type="det", newdata=newData, appendData=TRUE)
newData <- data.frame(Outcrop=0, Carb=p.Carb)
pred.occ.Carb <- predict(fm9, type="det", newdata=newData, appendData=TRUE)

# Plot predictions against unstandardized 'prediction covs'
par(mfrow = c(2,1), mar = c(5,5,2,3), cex.lab = 1.2)
plot(pred.occ.OA[[1]] ~ pred.OA, type = "l", lwd = 3, col = "blue", ylim = c(0,1), las = 1, ylab = "Pred. detection prob.", xlab = "Outcrop Area (%)", frame = F)
matlines(pred.OA, pred.occ.OA[,3:4], lty = 1, lwd = 1, col = "grey")
plot(pred.occ.Carb[[1]] ~ pred.Carb, type = "l", lwd = 3, col = "blue", ylim = c(0,1), las = 1, ylab = "Pred. detection prob.", xlab = "Carbonate Collections (%)", frame = F)
matlines(pred.Carb, pred.occ.Carb[,3:4], lty = 1, lwd = 1, col = "grey")


#Look at the occupancy and detection estimates on probability scale
#Note that 'backTransform' will only work without intermediate steps if inquiring about the null model.
print(occ.null <- backTransform(fm1, type="state")) # Occupancy estimate - 0.387
print(det.null <- backTransform(fm1, type="det")) # detection estimate - 0.108
