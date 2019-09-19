#===============================================================================================================================================
#============================================== OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS ==========================================
#===============================================================================================================================================

# Christopher D. Dean 
# 2019

# Broad PBDB download: https://paleobiodb.org/data1.2/occs/list.csv?base_name=tetrapoda&taxon_reso=genus&pres=regular&interval=Campanian,Maastrichtian&cc=NOA&envtype=terr&show=full,strat,lith,env

#=========================================================== DATA SETUP ======================================================================

# Set working directory
setwd("C:/Users/deancd/Documents/RESEARCH/PROJECTS/DINO_RANGE/NA-Dino-Occ/") # Set your working directory

# Load in Functions
source("0.Functions_DR.R") # Import functions from other R file (must be in same working directory)

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

#=============================================== DATA SETUP ===============================================

# Read in files
ct.occs <- read.csv("Data/Occurrences/Non-Rotated/cen_tur.csv", stringsAsFactors = FALSE) # Load in occurrences
ct.colls <- read.csv("Data/Occurrences/Non-Rotated/cen_tur_broad.csv", stringsAsFactors = FALSE) # load in all collections
ct.outcrop <- readShapePoly("Data/Geology/CenTurComplete.shp") #Get outcrop areas as polygons

# Combine datasets
combine_data(ct.occs, ct.colls)

# Visualise grid cells
get_grid_im(ct.occs.comb, 0.5, "Occurrences")
get_grid_im(ct.colls.comb, 0.5, "Collections")

just_arag_bivs <- ct.occs.comb %>%
  filter(Target == "Aragonitic Bivalves")
just_amms <- ct.occs.comb %>%
  filter(Target == "Ammonites")
just_calc_bivs <- ct.occs.comb %>%
  filter(Target == "Calcitic Bivalves")

get_grid_im(just_amms, 0.5, "Ammonite Occurrences")
get_grid_im(just_arag_bivs, 0.5, "Aragonitic Bivalve Occurrences")
get_grid_im(just_calc_bivs, 0.5, "Calcitic Bivalve Occurrences")

# Setup results
target = c("Ammonites", "Aragonitic Bivalves", "Calcitic Bivalves")
res <- c(0.5, 1)

# Run results
all_results_for_unmarked(ct.occs.comb, res, target, subsamp = TRUE)

# Get cell covariates
all_covs_info(ct.occs.comb, res, ct.outcrop)

#=============================================== LOWER CAMPANIAN ===============================================

# Read in files
lc.occs <- read.csv("Data/Occurrences/Non-Rotated/lower_camp_2.csv", stringsAsFactors = FALSE) # Load in occurrences
lc.colls <- read.csv("Data/Occurrences/Non-Rotated/lower_camp_Broad.csv", stringsAsFactors = FALSE) # load in all collections
lc.outcrop <- readShapePoly("Data/Geology/CampComplete.shp") #Get outcrop areas as polygons

# Combine datasets
combine_data(lc.occs, lc.colls)

# Visualise grid cells
get_grid_im(lc.occs.comb, 0.5, "Occurrences")
get_grid_im(lc.colls.comb, 0.5, "Collections")

just_arag_bivs <- lc.occs.comb %>%
  filter(Target == "Aragonitic Bivalves")
just_amms <- lc.occs.comb %>%
  filter(Target == "Ammonites")
just_calc_bivs <- lc.occs.comb %>%
  filter(Target == "Calcitic Bivalves")

get_grid_im(just_amms, 0.5, "Ammonite Occurrences")
get_grid_im(just_arag_bivs, 0.5, "Aragonitic Bivalve Occurrences")
get_grid_im(just_calc_bivs, 0.5, "Calcitic Bivalve Occurrences")

# Setup results
target = c("Ammonites", "Aragonitic Bivalves", "Calcitic Bivalves")
res <- c(0.1, 0.5, 1)

# Run results
all_results_for_unmarked(lc.occs.comb, res, target)

res <- 1

# Get cell covariates
all_covs_info(lc.occs.comb, res, lc.outcrop)
all_covs_info(ct.occs.comb, res, ct.outcrop)

#=========================================================== OCCUPANCY =================================================================

library("unmarked")

data <- read.csv("Results/ct.occs.comb.1.Ammonites.csv") # import data from previous step. Result here is just a test case to show.
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
