#===============================================================================================================================================
#============================================== OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS ==========================================
#===============================================================================================================================================

# Christopher D. Dean, Lewis A. Jones, Alfio A. Chiarenza, Alex Farnsworth, María V. Jiménez‐Franco, Richard J. Butler.
# 2019
# Script written by Christopher D. Dean and Lewis A. Jones

#============================================== FILE 2: RUNNING OCCUPANCY MODELLING WITH UNMARKED ==============================================

#============================================== INITIAL SETUP ===============================================

library(unmarked)
library(MuMIn)
library(dplyr)

# Load data
data <- read.csv("Results/0.5/camp.occs.targeted.0.5.CeratopsidaeSS.csv") # import data from previous step. Result here is just a test case to show.

# Covariates
site.covs <- read.csv("Results/0.5/CampanianCovariates.csv")
site.covs <- site.covs %>%
  dplyr::arrange(cells.1)

# Ensure data is the same between Observations and Covariates
data$X
site.covs$cells.1

data <- data[-121,] # Remove based on comparison between cells in data and sitecovs

# Check covariates for missing data and remove cells from both datasets
site.covs # DEMs - cell 7697
data <- data[-32,]
site.covs <- site.covs[-32,]

# Turn into matrix
y <- as.matrix(data[,2:ncol(data)])

DEMs.orig <- site.covs[,"DEM_0.5"]      # Unstandardised, original values of covariates
#LaCo.orig <- site.covs[,"LANDCVI_selected_1"]
OUTc.orig <- site.covs[,"Camp_out_0.5"]
RAIN.orig <- site.covs[,"WC_Prec_0.5"]
#TEMP.orig <- site.covs[,"WC_Temp_1"]
MGVF.orig <- site.covs[,"MGVF_0.5"]
COLL.orig <- site.covs[,"colls_per_cell"]
CaRa.orig <- site.covs[,"CampPrecip_0.5"]
CaTe.orig <- site.covs[,"CampTemp"]

# Overview of Covariates
covs <- cbind(DEMs.orig, OUTc.orig, RAIN.orig, MGVF.orig, COLL.orig, CaRa.orig, CaTe.orig)
par(mfrow = c(3,3))
for(i in 1:length(covs)){
  hist(covs[,i], breaks = 50, col = "grey", main = colnames(covs)[i])
}
pairs(cbind(DEMs.orig, OUTc.orig, RAIN.orig, MGVF.orig, COLL.orig, CaRa.orig, CaTe.orig))

# Standardize covariates and mean-impute OA and Carbation
# Compute means and standard deviations
(means <- c(apply(cbind(DEMs.orig, OUTc.orig, RAIN.orig, MGVF.orig, COLL.orig, CaRa.orig, CaTe.orig), 2, mean, na.rm = TRUE)))
(sds <- c(apply(cbind(DEMs.orig, OUTc.orig, RAIN.orig, MGVF.orig, COLL.orig, CaRa.orig, CaTe.orig), 2, sd, na.rm = TRUE)))

# Scale covariates
DEMs <- (DEMs.orig - means[1]) / sds[1] #
#LaCo<- (LaCo.orig - means[2]) / sds[2] 
OUTc <- (OUTc.orig - means[2]) / sds[2] # 
RAIN <- (RAIN.orig - means[3]) / sds[3] #
MGVF <- (MGVF.orig - means[4]) / sds[4] #
COLL <- (COLL.orig - means[5]) / sds[5]
CaRa <- (CaRa.orig - means[6]) / sds[6]
CaTe <- (CaTe.orig - means[7]) / sds[7]

#============================================== OCCUPANCY MODELLING =========================================

# Simple quick test
umf <- unmarkedFrameOccu(y = y)
summary(umf)
summary(fm1 <- occu(~1 ~1, data=umf))

data$X
site.covs$cells.1

print(occ.null <- backTransform(fm1, "state")) 
print(det.null <- backTransform(fm1, "det")) 


# Turn into matrix
umf <- unmarkedFrameOccu(y = y, siteCovs = data.frame(DEMs = DEMs, OUTc = OUTc,
                                                      RAIN = RAIN, MGVF = MGVF, 
                                                      COLL = COLL, CaRa = CaRa, 
                                                      CaTe = CaTe))

# Alternate scaling measure
siteCovs <- data.frame(DEMs = DEMs.orig, OUTc = OUTc.orig,
                       RAIN = RAIN.orig, MGVF = MGVF.orig, COLL = COLL.orig, 
                       CaRa = CaRa.orig, CaTe = CaTe.orig)
umf@siteCovs$DEMs <- scale(umf@siteCovs$DEMs)
umf@siteCovs$OUTc <- scale(umf@siteCovs$OUTc)
umf@siteCovs$RAIN <- scale(umf@siteCovs$RAIN)
umf@siteCovs$MGVF <- scale(umf@siteCovs$MGVF)
umf@siteCovs$CaRa <- scale(umf@siteCovs$CaRa)
umf@siteCovs$CaTe <- scale(umf@siteCovs$CaTe)

# umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs)
summary(umf)

# Fit null model
# Fit a series of models for detection first and do model selection

### OCCUPANCY ###
summary(fm1 <- occu(~1 ~1, data=umf))
summary(fm2 <- occu(~1 ~CaRa, data=umf))
summary(fm3 <- occu(~1 ~CaTe, data=umf))
summary(fm4 <- occu(~1 ~CaRa + CaTe, data=umf))
fms <- fitList("p(.)psi(.)"                   = fm1,
               "p(.)psi(CaRa)"                = fm2,
               "p(.)psi(CaTe)"                = fm3,
               "p(.)psi(CaRa + CaTe)"         = fm4
)
(ms <- modSel(fms))

summary(fm1 <- occu(~1 ~CaRa + I(CaRa^2), data=umf))
summary(fm2 <- occu(~1 ~CaRa + I(CaRa^2) + I(CaRa^3), data=umf))
summary(fm3 <- occu(~1 ~CaTe + I(CaTe^2), data=umf))
summary(fm4 <- occu(~1 ~CaTe + I(CaTe^2) + I(CaTe^3), data=umf))
summary(fm5 <- occu(~1 ~CaRa + CaTe + I(CaTe^2), data=umf))
summary(fm6 <- occu(~1 ~CaRa + CaTe + I(CaRa^2), data=umf))
summary(fm7 <- occu(~1 ~CaRa + CaTe + I(CaRa^2) + I(CaTe^2), data=umf))
fms <- fitList("p(.)psi(CaRa+2CaRa)"            = fm1,
               "p(.)psi(CaRa+2CaRa+3CaRa)"      = fm2,
               "p(.)psi(CaTe+2CaTe)"            = fm3, # BEST
               "p(.)psi(CaTe+2CaTe+3CaTe)"      = fm4,
               "p(.)psi(CaRa+CaTe+2CaTw)"       = fm5,
               "p(.)psi(CaRa+CaTe+2CaRa)"       = fm6,
               "p(.)psi(CaRa+CaTe+2CaRa+2CaTe)" = fm7
)
(ms <- modSel(fms))


### DETECTION ###
summary(fm1 <- occu(~1 ~CaTe + I(CaTe^2), data=umf))

summary(fm2 <- occu(~DEMs ~CaTe + I(CaTe^2), data=umf))
summary(fm3 <- occu(~OUTc ~CaTe + I(CaTe^2), data=umf))
summary(fm4 <- occu(~RAIN ~CaTe + I(CaTe^2), data=umf))
summary(fm5 <- occu(~MGVF ~CaTe + I(CaTe^2), data=umf))
summary(fm6 <- occu(~COLL ~CaTe + I(CaTe^2), data=umf))

summary(fm7 <- occu(~DEMs + OUTc ~CaTe + I(CaTe^2), data=umf))
summary(fm8 <- occu(~DEMs + RAIN ~CaTe + I(CaTe^2), data=umf))
summary(fm9 <- occu(~DEMs + MGVF ~CaTe + I(CaTe^2), data=umf))
summary(fm10 <- occu(~DEMs + COLL ~CaTe + I(CaTe^2), data=umf))
summary(fm11 <- occu(~OUTc + RAIN ~CaTe + I(CaTe^2), data=umf))
summary(fm12 <- occu(~OUTc + MGVF ~CaTe + I(CaTe^2), data=umf))
summary(fm13 <- occu(~OUTc + COLL ~CaTe + I(CaTe^2), data=umf))
summary(fm14 <- occu(~RAIN + MGVF ~CaTe + I(CaTe^2), data=umf))
summary(fm15 <- occu(~RAIN + COLL ~CaTe + I(CaTe^2), data=umf))
summary(fm16 <- occu(~MGVF + COLL ~CaTe + I(CaTe^2), data=umf))

summary(fm17 <- occu(~DEMs + OUTc + RAIN ~CaTe + I(CaTe^2), data=umf))
summary(fm18 <- occu(~DEMs + OUTc + MGVF ~CaTe + I(CaTe^2), data=umf))
summary(fm19 <- occu(~DEMs + OUTc + COLL ~CaTe + I(CaTe^2), data=umf))
summary(fm20 <- occu(~DEMs + RAIN + MGVF ~CaTe + I(CaTe^2), data=umf))
summary(fm21 <- occu(~DEMs + RAIN + COLL ~CaTe + I(CaTe^2), data=umf))
summary(fm22 <- occu(~DEMs + MGVF + COLL ~CaTe + I(CaTe^2), data=umf))
summary(fm23 <- occu(~OUTc + RAIN + MGVF ~CaTe + I(CaTe^2), data=umf))
summary(fm24 <- occu(~OUTc + RAIN + COLL ~CaTe + I(CaTe^2), data=umf))
summary(fm25 <- occu(~OUTc + MGVF + COLL ~CaTe + I(CaTe^2), data=umf))
summary(fm26 <- occu(~RAIN + MGVF + COLL ~CaTe + I(CaTe^2), data=umf))

summary(fm27 <- occu(~DEMs + OUTc + RAIN + MGVF ~CaTe + I(CaTe^2), data=umf))
summary(fm28 <- occu(~DEMs + OUTc + RAIN + COLL ~CaTe + I(CaTe^2), data=umf))
summary(fm29 <- occu(~DEMs + OUTc + MGVF + COLL ~CaTe + I(CaTe^2), data=umf))
summary(fm30 <- occu(~DEMs + RAIN + MGVF + COLL ~CaTe + I(CaTe^2), data=umf))
summary(fm31 <- occu(~RAIN + MGVF + COLL + OUTc ~CaTe + I(CaTe^2), data=umf))

summary(fm32 <- occu(~DEMs + OUTc + RAIN + MGVF + COLL ~CaTe + I(CaTe^2), data=umf))

# Put the fitted models in a "fitList" and rank them by AIC
fms <- fitList("p(.)psi(CaTe+CaTe2)"                        = fm1,
               "p(DEMs)psi(CaTe+CaTe2)"                     = fm2,
               "p(OUTc)psi(CaTe+CaTe2)"                     = fm3, # BEST
               "p(RAIN)psi(CaTe+CaTe2)"                     = fm4,
               "p(MGVF)psi(CaTe+CaTe2)"                     = fm5,
               "p(COLL)psi(CaTe+CaTe2)"                     = fm6,
               "p(DEMs+OUTc)psi(CaTe+CaTe2)"                = fm7, # BEST
               "p(DEMS+RAIN)psi(CaTe+CaTe2)"                = fm8,
               "p(DEMS+MGVF)psi(CaTe+CaTe2)"                = fm9, 
               "p(DEMS+COLL)psi(CaTe+CaTe2)"                = fm10,
               "p(OUTc+RAIN)psi(CaTe+CaTe2)"                = fm11,
               "p(OUTc+MGVF)psi(CaTe+CaTe2)"                = fm12, # BEST
               "p(OUTc+COLL)psi(CaTe+CaTe2)"                = fm13,
               "p(RAIN+MGVF)psi(CaTe+CaTe2)"                = fm14,
               "p(RAIN+COLL)psi(CaTe+CaTe2)"                = fm15,
               "p(MGVF+COLL)psi(CaTe+CaTe2)"                = fm16,
               "p(DEMS+OUTc+RAIN)psi(CaTe+CaTe2)"                = fm17,
               "p(DEMS+OUTc+MGVF)psi(CaTe+CaTe2)"                = fm18, # BEST
               "p(DEMS+OUTc+COLL)psi(CaTe+CaTe2)"                = fm19,
               "p(DEMS+RAIN+MGVF)psi(CaTe+CaTe2)"                = fm20,
               "p(DEMS+RAIN+COLL)psi(CaTe+CaTe2)"                = fm21,
               "p(DEMS+MGVF+COLL)psi(CaTe+CaTe2)"                = fm22,
               "p(OUTc+RAIN+MGVF)psi(CaTe+CaTe2)"                = fm23, # BEST
               "p(OUTc+RAIN+COLL)psi(CaTe+CaTe2)"                = fm24,
               "p(OUTc+MGVF+COLL)psi(CaTe+CaTe2)"                = fm25,
               "p(RAIN+MGVF+COLL)psi(CaTe+CaTe2)"                = fm26,
               "p(DEMS+OUTc+RAIN+MGVF)psi(CaTe+CaTe2)"                = fm27,
               "p(DEMS+OUTc+RAIN+COLL)psi(CaTe+CaTe2)"                = fm28,
               "p(DEMS+OUTc+MGVF+COLL)psi(CaTe+CaTe2)"                = fm29,
               "p(DEMS+RAIN+MGVF+COLL)psi(CaTe+CaTe2)"                = fm30,
               "p(RAIN+MGVF+COLL+OUTc)psi(CaTe+CaTe2)"                = fm31,
               "p(DEMS+OUTc+RAIN+MGVF+COLL)psi(CaTe+CaTe2)"                = fm32
               )
(ms <- modSel(fms))



### Proportion of Area Occupied ###
AICbest <- occu(formula = ~ OUTc+MGVF # Best result here
                          ~ CaTe+I(CaTe^2),
                data = umf)
re <- ranef(AICbest)
EBUP <- bup(re, stat="mean")
CI <- confint(re, level=0.9)
rbind(PAO = c(Estimate = sum(EBUP), colSums(CI)) / nrow(y))


### Auto Select - Not using for Now ###
# Using MuMIn
full <- occu(formula =  ~ OUTc +  MGVF + DEMs + COLL
             ~ CaRa + CaTe,
             data = umf)
(modelList <- dredge(full, rank = "AIC"))

occu_dredge_95 <- get.models(modelList, subset = cumsum(weight) <= 0.95)
oc_avg <- model.avg(occu_dredge_95, fit = TRUE)
t(oc_avg$coefficients)


#### PREDICTION ####

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
