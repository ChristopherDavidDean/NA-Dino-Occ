#===============================================================================================================================================
#============================================== OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS ==========================================
#===============================================================================================================================================

# Christopher D. Dean, Lewis A. Jones, Alfio A. Chiarenza, Alex Farnsworth, María V. Jiménez‐Franco, Richard J. Butler.
# 2019
# Script written by Christopher D. Dean and Lewis A. Jones

#============================================== FILE 2: RUNNING OCCUPANCY MODELLING WITH UNMARKED ==============================================

#============================================== INITIAL SETUP ===============================================

library("unmarked")

data <- read.csv("Results/1/camp.occs.targeted.1.CeratopsidaeSS.csv") # import data from previous step. Result here is just a test case to show.
data <- data[-73,]

# Turn into matrix
y <- as.matrix(data[,2:11])

#============================================== OCCUPANCY MODELLING =========================================

# Simple quick test
umf <- unmarkedFrameOccu(y = y)
summary(umf)
summary(fm1 <- occu(~1 ~1, data=umf))

data$X
site.covs$cells.1

print(occ.null <- backTransform(fm1, "state")) 
print(det.null <- backTransform(fm1, "det")) 

# Covariates
site.covs <- read.csv("Results/1/CampanianCovariates.csv")
site.covs <- site.covs %>%
  dplyr::arrange(cells.1)

DEMs.orig <- site.covs[,"DEM_1"]      # Unstandardised, original values of covariates
LaCo.orig <- site.covs[,"LANDCVI_selected_1"]
OUTc.orig <- site.covs[,"Camp_out_1"]
RAIN.orig <- site.covs[,"WC_Prec_1"]
TEMP.orig <- site.covs[,"WC_Temp_1"]
MGVF.orig <- site.covs[,"MGVF_1"]
CaRa.orig <- site.covs[,"CampPrecip_1"]
CaTe.orig <- site.covs[,"CampTemp"]

# Overview of Covariates
covs <- cbind(DEMs.orig, LaCo.orig, OUTc.orig, RAIN.orig, MGVF.orig, CaRa.orig, CaTe.orig)
par(mfrow = c(3,3))
for(i in 1:7){
  hist(covs[,i], breaks = 50, col = "grey", main = colnames(covs)[i])
}
pairs(cbind(DEMs.orig, LaCo.orig, OUTc.orig, RAIN.orig, MGVF.orig, CaRa.orig, CaTe.orig))

# Standardize covariates and mean-impute OA and Carbation
# Compute means and standard deviations
(means <- c(apply(cbind(DEMs.orig, LaCo.orig, OUTc.orig, RAIN.orig, MGVF.orig, CaRa.orig, CaTe.orig), 2, mean, na.rm = TRUE)))
(sds <- c(apply(cbind(DEMs.orig, LaCo.orig, OUTc.orig, RAIN.orig, MGVF.orig, CaRa.orig, CaTe.orig), 2, sd, na.rm = TRUE)))

# Scale covariates
DEMs<- (DEMs.orig - means[1]) / sds[1] #
LaCo<- (LaCo.orig - means[2]) / sds[2] 
OUTc<- (OUTc.orig - means[3]) / sds[3] # 
RAIN<- (RAIN.orig - means[4]) / sds[4] #
MGVF<- (MGVF.orig - means[5]) / sds[5] #
CaRa<- (CaRa.orig - means[6]) / sds[6]
CaTe<- (CaTe.orig - means[7]) / sds[7]


# Turn into matrix
umf <- unmarkedFrameOccu(y = y, siteCovs = data.frame(DEMs = DEMs, LaCo = LaCo, OUTc = OUTc,
                                                      RAIN = RAIN, MGVF = MGVF, CaRa = CaRa, 
                                                      CaTe = CaTe))
summary(umf)

# Fit null model
# Fit a series of models for detection first and do model selection
summary(fm1 <- occu(~1 ~1, data=umf))
summary(fm2 <- occu(~1 ~CaRa, data=umf))
summary(fm3 <- occu(~1 ~CaTe, data=umf))
summary(fm4 <- occu(~1 ~CaRa + CaTe, data=umf))
fms <- fitList("p(.)psi(.)"                   = fm1,
               "p(CaRa)psi(.)"                = fm2,
               "p(CaTe)psi(.)"                = fm3,
               "p(CaRa + CaTe)psi(.)"         = fm4
)
(ms <- modSel(fms))


# WITH DEMS
summary(fm1 <- occu(~1 ~1, data=umf))
summary(fm2 <- occu(~DEMs ~1, data=umf))
summary(fm3 <- occu(~OUTc ~1, data=umf))
summary(fm4 <- occu(~RAIN ~1, data=umf))
summary(fm5 <- occu(~MGVF ~1, data=umf))
summary(fm6 <- occu(~DEMs + OUTc ~1, data=umf))
summary(fm7 <- occu(~DEMs + RAIN ~1, data=umf))
summary(fm8 <- occu(~DEMs + MGVF ~1, data=umf))
summary(fm9 <- occu(~OUTc + RAIN ~1, data=umf))
summary(fm10 <- occu(~OUTc + MGVF ~1, data=umf))
summary(fm11 <- occu(~RAIN + MGVF~1, data=umf))
summary(fm12 <- occu(~DEMs + OUTc + RAIN ~1, data=umf))
summary(fm13 <- occu(~DEMs + OUTc + MGVF ~1, data=umf))
summary(fm14 <- occu(~OUTc + RAIN + MGVF ~1, data=umf))
summary(fm15 <- occu(~DEMs + RAIN + MGVF ~1, data=umf))
summary(fm16 <- occu(~DEMs + + OUTc + RAIN + MGVF ~1, data=umf))

# Put the fitted models in a "fitList" and rank them by AIC
fms <- fitList("p(.)psi(.)"                   = fm1,
               "p(DEMs)psi(.)"                = fm2,
               "p(OUTc)psi(.)"                = fm3,
               "p(RAIN)psi(.)"                = fm4,
               "p(MGVF)psi(.)"                = fm5,
               "p(DEMs+OUTc)psi(.)"           = fm6,
               "p(DEMs+RAIN)psi(.)"           = fm7,
               "p(DEMS+MGVF)psi(.)"           = fm8,
               "p(OUTc+RAIN)psi(.)"           = fm9, 
               "p(OUTc+MGVF)psi(.)"           = fm10,
               "p(RAIN+MGVF)psi(.)"           = fm11,
               "p(DEMs+OUTc+RAIN)psi(.)"      = fm12,
               "p(DEMs+OUTc+MGVF)psi(.)"      = fm13,
               "p(OUTc+RAIN+MGVF)psi(.)"      = fm14,
               "p(DEMs+RAIN+MGVF)psi(.)"      = fm15,
               "p(DEMs+OUTc+RAIN+MGVF)psi(.)" = fm16
               )
(ms <- modSel(fms))

# WITHOUT DEMS
summary(fm1 <- occu(~1 ~1, data=umf))
summary(fm2 <- occu(~OUTc ~1, data=umf))
summary(fm3 <- occu(~RAIN ~1, data=umf))
summary(fm4 <- occu(~MGVF ~1, data=umf))
summary(fm5 <- occu(~OUTc + RAIN ~1, data=umf))
summary(fm6 <- occu(~OUTc + MGVF ~1, data=umf))
summary(fm7 <- occu(~RAIN + MGVF ~1, data=umf))
summary(fm8 <- occu(~OUTc + RAIN + MGVF ~1, data=umf))

# Put the fitted models in a "fitList" and rank them by AIC
fms <- fitList("p(.)psi(.)"                   = fm1,
               "p(OUTc)psi(.)"                = fm2,
               "p(RAIN)psi(.)"                = fm3,
               "p(MGVF)psi(.)"                = fm4,
               "p(OUTc+RAIN)psi(.)"           = fm5, 
               "p(OUTc+MGVF)psi(.)"           = fm6,
               "p(RAIN+MGVF)psi(.)"           = fm7,
               "p(OUTc+RAIN+MGVF)psi(.)"      = fm8
)
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
