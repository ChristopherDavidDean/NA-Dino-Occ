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
data <- read.csv("Results/0.5/Standard_Occ/Subsampled_20/Singleton/camp.occs.targeted.0.5.HadrosauridaeSS.csv") # import data from previous step. Result here is just a test case to show.


# Covariates
site.covs <- read.csv("Results/0.5/CampanianHiResCovariates.csv") # Remember to change resolution!!!
#site.covs <- site.covs %>%
#  dplyr::arrange(cells)

# Ensure data is the same between Observations and Covariates
as.numeric(rownames(data))
site.covs$cells


# Check covariates for missing data and remove cells from both datasets
site.covs <- site.covs[-c(2, 170, 281, 284, 285, 287),] # 0.1
data <- data[-c(2, 170, 281, 284, 285, 287),] # 0.1

site.covs <- site.covs[-c(2, 79, 121),] # 0.5
data <- data[-c(2, 79, 121),] # 0.5

# Turn into matrix
y <- as.matrix(data[,2:ncol(data)])

DEMs.orig <- site.covs[,"mean_DEM"]      # Unstandardised, original values of covariates
RAIN.orig <- site.covs[,"mean_prec"]
TEMP.orig <- site.covs[,"mean_temp"]
MGVF.orig <- site.covs[,"mean_MGVF"]
Cpre.orig <- site.covs[,"mean_Cpre"]
Ctem.orig <- site.covs[,"mean_Ctem"]

# Overview of Covariates
covs <- cbind(DEMs.orig, RAIN.orig, MGVF.orig, TEMP.orig, Cpre.orig, Ctem.orig)
par(mfrow = c(1,1))
for(i in 1:length(covs)){
  hist(covs[,i], breaks = 50, col = "grey", main = colnames(covs)[i])
}
pairs(cbind(DEMs.orig, RAIN.orig, MGVF.orig, TEMP.orig, Cpre.orig, Ctem.orig))


# Compute means and standard deviations
(means <- c(apply(cbind(DEMs.orig, RAIN.orig, MGVF.orig, TEMP.orig, Cpre.orig, Ctem.orig), 2, mean, na.rm = TRUE)))
(sds <- c(apply(cbind(DEMs.orig, RAIN.orig, MGVF.orig, TEMP.orig, Cpre.orig, Ctem.orig), 2, sd, na.rm = TRUE)))

# Scale covariates
DEMs <- (DEMs.orig - means[1]) / sds[1] #
RAIN <- (RAIN.orig - means[2]) / sds[2] #
MGVF <- (MGVF.orig - means[3]) / sds[3] #
TEMP <- (TEMP.orig - means[4]) / sds[4]
Cpre <- (Cpre.orig - means[5]) / sds[5]
Ctem <- (Ctem.orig - means[6]) / sds[6]

#============================================== OCCUPANCY MODELLING =========================================

#===== Simple quick test =====
umf <- unmarkedFrameOccu(y = y)
summary(umf)
summary(fm1 <- occu(~1 ~1, data=umf))

print(occ.null <- backTransform(fm1, "state")) 
print(det.null <- backTransform(fm1, "det")) 

#==================================================================================================

# Turn into matrix
umf <- unmarkedFrameOccu(y = y, siteCovs = data.frame(DEMs = DEMs, 
                                                      RAIN = RAIN, 
                                                      MGVF = MGVF, 
                                                      TEMP = TEMP,
                                                      Cpre = Cpre,
                                                      Ctem = Ctem))

# Alternate scaling measure
#siteCovs <- data.frame(DEMs = DEMs.orig, OUTc = OUTc.orig,
#                      RAIN = RAIN.orig, MGVF = MGVF.orig, TEMP = TEMP.orig, 
#                       Cpre = Cpre.orig, Ctem = Ctem.orig)
#umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs)
#umf@siteCovs$DEMs <- scale(umf@siteCovs$DEMs)
#umf@siteCovs$OUTc <- scale(umf@siteCovs$OUTc)
#umf@siteCovs$RAIN <- scale(umf@siteCovs$RAIN)
#umf@siteCovs$MGVF <- scale(umf@siteCovs$MGVF)
#umf@siteCovs$Cpre <- scale(umf@siteCovs$Cpre)
#umf@siteCovs$Ctem <- scale(umf@siteCovs$Ctem)

summary(umf)

# Fit null model
# Fit a series of models for detection first and do model selection

### DETECTION ###
summary(fm1 <- occu(~1 ~1, data=umf))

summary(fm2 <- occu(~DEMs ~1, data=umf))
summary(fm3 <- occu(~RAIN ~1, data=umf))
summary(fm4 <- occu(~MGVF ~1, data=umf))
summary(fm5 <- occu(~TEMP ~1, data=umf))

summary(fm6 <- occu(~DEMs + RAIN ~1, data=umf))
summary(fm7 <- occu(~DEMs + MGVF ~1, data=umf))
summary(fm8 <- occu(~DEMs + TEMP ~1, data=umf))
summary(fm9 <- occu(~RAIN + MGVF ~1, data=umf))
summary(fm10 <- occu(~RAIN + TEMP ~1, data=umf)) #
summary(fm11 <- occu(~MGVF + TEMP ~1, data=umf))

summary(fm12 <- occu(~DEMs + RAIN + MGVF ~1, data=umf))
summary(fm13 <- occu(~DEMs + RAIN + TEMP ~1, data=umf)) #
summary(fm14 <- occu(~DEMs + MGVF + TEMP ~1, data=umf))
summary(fm15 <- occu(~RAIN + MGVF + TEMP ~1, data=umf)) #

summary(fm16 <- occu(~DEMs + RAIN + MGVF + TEMP ~1, data=umf)) #

# Put the fitted models in a "fitList" and rank them by AIC
fms <- fitList("1"                        = fm1,
               "2"                     = fm2,
               "3"                     = fm3, # BEST
               "4"                     = fm4,
               "5"                     = fm5,
               "6"                     = fm6,
               "7"                = fm7, # BEST
               "8"                = fm8,
               "9"                = fm9, 
               "10"                = fm10,
               "11"                = fm11,
               "12"                = fm12, # BEST
               "13"                = fm13,
               "14"                = fm14,
               "15"                = fm15,
               "16"                = fm16
               )
(ms <- modSel(fms))

# Further tests
summary(fm1 <- occu(~RAIN + TEMP ~1, data=umf))
summary(fm2 <- occu(~RAIN + I(RAIN^2) + TEMP ~1, data=umf))
summary(fm3 <- occu(~RAIN + I(RAIN^2) + I(RAIN^3) + TEMP ~1, data=umf))
summary(fm4 <- occu(~RAIN + TEMP + I(TEMP^2) ~1, data=umf))
summary(fm5 <- occu(~RAIN + TEMP + I(TEMP^2) +I(TEMP^3) ~1, data=umf))
summary(fm6 <- occu(~RAIN + I(RAIN^2) + TEMP + I(TEMP^2) ~1, data=umf))
summary(fm7 <- occu(~RAIN + I(RAIN^2) + I(RAIN^3) + TEMP + I(TEMP^2) ~1, data=umf))
summary(fm8 <- occu(~RAIN + I(RAIN^2) + TEMP + I(TEMP^2) + I(TEMP^3) ~1, data=umf))
summary(fm9 <- occu(~RAIN + I(RAIN^2) + I(RAIN^3) + TEMP + I(TEMP^2) + I(TEMP^3) ~1, data=umf))

fms <- fitList("p(RAIN+TEMP)psi(.)"                                  = fm1, # BEST
               "p(RAIN+RAIN2+TEMP)psi(.)"                            = fm2,
               "p(RAIN+RAIN2+RAIN3+TEMP)psi(.)"                      = fm3, 
               "p(RAIN+TEMP+TEMP2)psi(.)"                            = fm4,
               "p(RAIN+TEMP+TEMP2+TEMP3)psi(.)"                      = fm5,
               "p(RAIN+RAIN2+TEMP+TEMP2)psi(.)"                      = fm6,
               "p(RAIN+RAIN2+RAIN3+TEMP+TEMP2)psi(.)"                = fm7, 
               "p(RAIN+RAIN2+TEMP+TEMP2+TEMP3)psi(.)"                = fm8,
               "p(RAIN+RAIN2+RAIN3+TEMP+TEMP2+TEMP3)psi(.)"          = fm9
)
(ms <- modSel(fms))

best.model <- fm3
confint(best.model, type = "det")


### OCCUPANCY ###
summary(fm1 <- occu(~RAIN + I(RAIN^2) + I(RAIN^3) + TEMP + I(TEMP^2) + I(TEMP^3) ~1, data=umf))
summary(fm2 <- occu(~RAIN + I(RAIN^2) + I(RAIN^3) + TEMP + I(TEMP^2) + I(TEMP^3) ~Ctem, data=umf))
summary(fm3 <- occu(~RAIN + I(RAIN^2) + I(RAIN^3) + TEMP + I(TEMP^2) + I(TEMP^3) ~Cpre, data=umf))
summary(fm4 <- occu(~RAIN + I(RAIN^2) + I(RAIN^3) + TEMP + I(TEMP^2) + I(TEMP^3) ~Cpre + Ctem, data=umf))
fms <- fitList("p(RAIN+RAIN2+RAIN3+TEMP+TEMP2+TEMP3)psi(.)"                   = fm1,
               "p(RAIN+RAIN2+RAIN3+TEMP+TEMP2+TEMP3)psi(Cpre)"                = fm2,
               "p(RAIN+RAIN2+RAIN3+TEMP+TEMP2+TEMP3)psi(Ctem)"                = fm3,
               "p(RAIN+RAIN2+RAIN3+TEMP+TEMP2+TEMP3)psi(Cpre + Ctem)"         = fm4
)

summary(fm1 <- occu(~DEMs + I(DEMs^2) + I(DEMs^3) + TEMP ~1, data=umf))
summary(fm2 <- occu(~DEMs + I(DEMs^2) + I(DEMs^3) + TEMP ~Ctem, data=umf))
summary(fm3 <- occu(~DEMs + I(DEMs^2) + I(DEMs^3) + TEMP ~Cpre, data=umf))
summary(fm4 <- occu(~DEMs + I(DEMs^2) + I(DEMs^3) + TEMP ~Cpre + Ctem, data=umf))
fms <- fitList("p(RAIN+RAIN2+RAIN3+TEMP+TEMP2+TEMP3)psi(.)"                   = fm1,
               "p(RAIN+RAIN2+RAIN3+TEMP+TEMP2+TEMP3)psi(Cpre)"                = fm2,
               "p(RAIN+RAIN2+RAIN3+TEMP+TEMP2+TEMP3)psi(Ctem)"                = fm3,
               "p(RAIN+RAIN2+RAIN3+TEMP+TEMP2+TEMP3)psi(Cpre + Ctem)"         = fm4
)

summary(fm1 <- occu(~RAIN + MGVF + TEMP ~1, data=umf))
summary(fm2 <- occu(~RAIN + MGVF + TEMP ~Ctem, data=umf))
summary(fm3 <- occu(~RAIN + MGVF + TEMP ~Cpre, data=umf))
summary(fm4 <- occu(~RAIN + MGVF + TEMP ~Cpre + Ctem, data=umf))
fms <- fitList("p(RAIN+RAIN2+RAIN3+TEMP+TEMP2+TEMP3)psi(.)"                   = fm1,
               "p(RAIN+RAIN2+RAIN3+TEMP+TEMP2+TEMP3)psi(Cpre)"                = fm2,
               "p(RAIN+RAIN2+RAIN3+TEMP+TEMP2+TEMP3)psi(Ctem)"                = fm3,
               "p(RAIN+RAIN2+RAIN3+TEMP+TEMP2+TEMP3)psi(Cpre + Ctem)"         = fm4
)

(ms <- modSel(fms))

summary(fm1 <- occu(~RAIN + I(RAIN^2) + I(RAIN^3) + TEMP + I(TEMP^2) + I(TEMP^3) ~Cpre + I(Cpre^2), data=umf)) # Null best
summary(fm2 <- occu(~RAIN + I(RAIN^2) + I(RAIN^3) + TEMP + I(TEMP^2) + I(TEMP^3) ~Cpre + I(Cpre^2) + I(Cpre^3), data=umf))
summary(fm3 <- occu(~RAIN + I(RAIN^2) + I(RAIN^3) + TEMP + I(TEMP^2) + I(TEMP^3) ~Ctem + I(Ctem^2), data=umf))
summary(fm4 <- occu(~RAIN + I(RAIN^2) + I(RAIN^3) + TEMP + I(TEMP^2) + I(TEMP^3) ~Ctem + I(Ctem^2) + I(Ctem^3), data=umf))
summary(fm5 <- occu(~RAIN + I(RAIN^2) + I(RAIN^3) + TEMP + I(TEMP^2) + I(TEMP^3) ~Cpre + Ctem + I(Ctem^2), data=umf))
summary(fm6 <- occu(~RAIN + I(RAIN^2) + I(RAIN^3) + TEMP + I(TEMP^2) + I(TEMP^3) ~Cpre + Ctem + I(Cpre^2), data=umf))
summary(fm7 <- occu(~RAIN + I(RAIN^2) + I(RAIN^3) + TEMP + I(TEMP^2) + I(TEMP^3) ~Cpre + Ctem + I(Cpre^2) + I(Ctem^2), data=umf))
summary(fm8 <- occu(~RAIN + I(RAIN^2) + I(RAIN^3) + TEMP + I(TEMP^2) + I(TEMP^3) ~1, data=umf))

fms <- fitList("p(RAIN+RAIN2+RAIN3+TEMP+TEMP2+TEMP3)psi(Cpre+2Cpre)"            = fm1,
               "p(RAIN+RAIN2+RAIN3+TEMP+TEMP2+TEMP3)psi(Cpre+2Cpre+3Cpre)"      = fm2,
               "p(RAIN+RAIN2+RAIN3+TEMP+TEMP2+TEMP3)psi(Ctem+2Ctem)"            = fm3, # BEST
               "p(RAIN+RAIN2+RAIN3+TEMP+TEMP2+TEMP3)psi(Ctem+2Ctem+3Ctem)"      = fm4,
               "p(RAIN+RAIN2+RAIN3+TEMP+TEMP2+TEMP3)psi(Cpre+Ctem+2CaTw)"       = fm5,
               "p(RAIN+RAIN2+RAIN3+TEMP+TEMP2+TEMP3)psi(Cpre+Ctem+2Cpre)"       = fm6,
               "p(RAIN+RAIN2+RAIN3+TEMP+TEMP2+TEMP3)psi(Cpre+Ctem+2Cpre+2Ctem)" = fm7, 
               "Null"                           = fm8
)
(ms <- modSel(fms))


summary(fm1 <- occu(~DEMs + I(DEMs^2) + I(DEMs^3) + TEMP ~Cpre + I(Cpre^2), data=umf)) # Null best
summary(fm2 <- occu(~DEMs + I(DEMs^2) + I(DEMs^3) + TEMP ~Cpre + I(Cpre^2) + I(Cpre^3), data=umf))
summary(fm3 <- occu(~DEMs + I(DEMs^2) + I(DEMs^3) + TEMP ~Ctem + I(Ctem^2), data=umf))
summary(fm4 <- occu(~DEMs + I(DEMs^2) + I(DEMs^3) + TEMP ~Ctem + I(Ctem^2) + I(Ctem^3), data=umf))
summary(fm5 <- occu(~DEMs + I(DEMs^2) + I(DEMs^3) + TEMP ~Cpre + Ctem + I(Ctem^2), data=umf))
summary(fm6 <- occu(~DEMs + I(DEMs^2) + I(DEMs^3) + TEMP ~Cpre + Ctem + I(Cpre^2), data=umf))
summary(fm7 <- occu(~DEMs + I(DEMs^2) + I(DEMs^3) + TEMP ~Cpre + Ctem + I(Cpre^2) + I(Ctem^2), data=umf))
summary(fm8 <- occu(~DEMs + I(DEMs^2) + I(DEMs^3) + TEMP ~1, data=umf))

summary(fm1 <- occu(~RAIN + MGVF + TEMP ~Cpre + I(Cpre^2), data=umf)) # Null best
summary(fm2 <- occu(~RAIN + MGVF + TEMP ~Cpre + I(Cpre^2) + I(Cpre^3), data=umf))
summary(fm3 <- occu(~RAIN + MGVF + TEMP ~Ctem + I(Ctem^2), data=umf))
summary(fm4 <- occu(~RAIN + MGVF + TEMP ~Ctem + I(Ctem^2) + I(Ctem^3), data=umf))
summary(fm5 <- occu(~RAIN + MGVF + TEMP ~Cpre + Ctem + I(Ctem^2), data=umf))
summary(fm6 <- occu(~RAIN + MGVF + TEMP ~Cpre + Ctem + I(Cpre^2), data=umf))
summary(fm7 <- occu(~RAIN + MGVF + TEMP ~Cpre + Ctem + I(Cpre^2) + I(Ctem^2), data=umf))
summary(fm8 <- occu(~RAIN + MGVF + TEMP ~1, data=umf))

fms <- fitList("p(RAIN+RAIN2+RAIN3+TEMP+TEMP2+TEMP3)psi(Cpre+2Cpre)"            = fm1,
               "p(RAIN+RAIN2+RAIN3+TEMP+TEMP2+TEMP3)psi(Cpre+2Cpre+3Cpre)"      = fm2,
               "p(RAIN+RAIN2+RAIN3+TEMP+TEMP2+TEMP3)psi(Ctem+2Ctem)"            = fm3, # BEST
               "p(RAIN+RAIN2+RAIN3+TEMP+TEMP2+TEMP3)psi(Ctem+2Ctem+3Ctem)"      = fm4,
               "p(RAIN+RAIN2+RAIN3+TEMP+TEMP2+TEMP3)psi(Cpre+Ctem+2CaTw)"       = fm5,
               "p(RAIN+RAIN2+RAIN3+TEMP+TEMP2+TEMP3)psi(Cpre+Ctem+2Cpre)"       = fm6,
               "p(RAIN+RAIN2+RAIN3+TEMP+TEMP2+TEMP3)psi(Cpre+Ctem+2Cpre+2Ctem)" = fm7, 
               "Null"                           = fm8
)
(ms <- modSel(fms))



best.model <- fm2

confint(best.model, type = "state")



### Proportion of Area Occupied ###
re <- ranef(best.model)
EBUP <- bup(re, stat="mean")
CI <- confint(re, level=0.9)
rbind(PAO = c(Estimate = sum(EBUP), colSums(CI)) / nrow(y))


### Auto Select - Not using for Now ###
# Using MuMIn
full <- occu(formula =  ~ OUTc +  MGVF + DEMs + TEMP
             ~ Cpre + Ctem,
             data = umf)
(modelList <- dredge(full, rank = "AIC"))

occu_dredge_95 <- get.models(modelList, subset = cumsum(weight) <= 0.95)
oc_avg <- model.avg(occu_dredge_95, fit = TRUE)
t(oc_avg$coefficients)

library(AICcmodavg)
system.time(gof.boot <- mb.gof.test(best.model, nsim = 1000))
gof.boot


#### PREDICTION ####



# Create new covariates for prediction ('prediction covs')
pred.TEMP <- seq(min(TEMP.orig), max(TEMP.orig),,100) # New covs for prediction
pred.RAIN <- seq(min(RAIN.orig), max(RAIN.orig),,100)
pred.DEMs <- seq(min(DEMs.orig), max(DEMs.orig),,100)
pred.MGVF <- seq(min(MGVF.orig), max(MGVF.orig),,100)
p.TEMP <- (pred.TEMP - means[4]) / sds[4] # Standardise them like actual covs
p.RAIN <- (pred.RAIN - means[2]) / sds[2]
p.DEMs <- (pred.DEMs - means[1]) / sds[1]
p.MGVF <- (pred.MGVF - means[3]) / sds[3]

pred.Cpre <- seq(min(Cpre.orig), max(Cpre.orig),,100) # New covs for prediction
p.Cpre <- (pred.Cpre - means[5]) / sds[5]
pred.Ctem <- seq(min(Ctem.orig), max(Ctem.orig),,100) # New covs for prediction
p.Ctem <- (pred.Ctem - means[6]) / sds[6]


# Obtain predictions
newData <- data.frame(TEMP=p.TEMP, RAIN = 0)
pred.det.TEMP <- predict(best.model, type="det", newdata=newData, appendData=TRUE)
newData <- data.frame(TEMP=0, RAIN = p.RAIN)
pred.det.RAIN <- predict(best.model, type="det", newdata=newData, appendData=TRUE)
newData <- data.frame(TEMP=0, RAIN = 0, MGVF = p.MGVF)
pred.det.MGVF <- predict(best.model, type="det", newdata=newData, appendData=TRUE)

newData <- data.frame(Cpre = p.Cpre)
pred.occ.Cpre <- predict(best.model, type="state", newdata=newData, appendData=TRUE)
newData <- data.frame(Ctem = p.Ctem)
pred.occ.Ctem <- predict(best.model, type="state", newdata=newData, appendData=TRUE)

# Plot predictions against unstandardized 'prediction covs'
par(mfrow = c(3,1), mar = c(5,5,2,3), cex.lab = 1.2)
plot(pred.det.TEMP[[1]] ~ pred.TEMP, type = "l", lwd = 3, col = "blue", ylim = c(0,1), las = 1, ylab = "Pred. detection prob.", xlab = "Temperature", frame = F)
matlines(pred.TEMP, pred.det.TEMP[,3:4], lty = 1, lwd = 1, col = "grey")
plot(pred.det.RAIN[[1]] ~ pred.RAIN, type = "l", lwd = 3, col = "blue", ylim = c(0,1), las = 1, ylab = "Pred. detection prob.", xlab = "Precipitation", frame = F)
matlines(pred.RAIN, pred.det.RAIN[,3:4], lty = 1, lwd = 1, col = "grey")
plot(pred.det.MGVF[[1]] ~ pred.MGVF, type = "l", lwd = 3, col = "blue", ylim = c(0,1), las = 1, ylab = "Pred. detection prob.", xlab = "MGVF", frame = F)
matlines(pred.MGVF, pred.det.MGVF[,3:4], lty = 1, lwd = 1, col = "grey")

plot(pred.occ.Cpre[[1]] ~ pred.Cpre, type = "l", lwd = 3, col = "blue", ylim = c(0,1), las = 1, ylab = "Pred. occupancy prob.", xlab = "Palaeo Precip.", frame = F)
matlines(pred.Cpre, pred.occ.Cpre[,3:4], lty = 1, lwd = 1, col = "grey")
plot(pred.occ.Ctem[[1]] ~ pred.Ctem, type = "l", lwd = 3, col = "blue", ylim = c(0,1), las = 1, ylab = "Pred. occupancy prob.", xlab = "Palaeo Temp.", frame = F)
matlines(pred.Ctem, pred.occ.Ctem[,3:4], lty = 1, lwd = 1, col = "grey")

#Look at the occupancy and detection estimates on probability scale
#Note that 'backTransform' will only work without intermediate steps if inquiring about the null model.
print(occ.null <- backTransform(best.model, type="state")) # Occupancy estimate - 0.387
print(det.null <- backTransform(best.model, type="det")) # detection estimate - 0.108




gen_raster <- function(cell_data, value_data, res, ext, zero = FALSE){
  init_raster <- raster(res = res, ext = ext, val = 1)
  total_cells <- ncell(init_raster)
  dframe_of_values <- data.frame(cell_data, value_data) 
  colnames(dframe_of_values) <- c("Cells", "Vals")
  
  # Make blank dataframe to copy values into
  dframe_of_cells <- data.frame(1:total_cells)
  colnames(dframe_of_cells) <- "Cells"
  
  # Join dataframes
  full_dframe <- left_join(dframe_of_cells, dframe_of_values, by = "Cells")
  
  if (zero == TRUE){
    full_dframe[is.na(full_dframe)] <- 0
  }
  mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
  raster_for_values <- raster(res = res, val = full_dframe$Vals, ext = ext)
  plot(raster_for_values, col = mapPalette(100), axes = F, box = F, main = "Detection Probability")
}

library(tidyr)
par(mfrow = c(1,1), mar = c(5,5,2,3), cex.lab = 1.2)
chosen_raster <- rain_ras



selected_covs <- site.covs %>%
  dplyr::select(mean_prec, mean_temp)
colnames(selected_covs) <- c("RAIN", "TEMP")






# Need to take: raster, covariates


# First
rain_ras <- raster("Data/Covariate_Data/Formatted/All_data/0.5deg/WC_Prec_0.5.asc")
temp_ras <- raster("Data/Covariate_Data/Formatted/All_data/0.5deg/WC_Temp_0.5.asc")
chosen_raster <- stack(rain_ras, temp_ras)
names(chosen_raster) <- c("RAIN", "TEMP")



# Scale covariates
means <- colMeans(selected_covs)
sds <- selected_covs %>% summarise_each(~(sd(., na.rm=TRUE)))

means[1]


scaled_covs <- apply(ras, 2, scale)

# Useful Info
total_cells <- ncell(chosen_raster) # Number of cells
raster_extent <- chosen_raster@extent # Extent
ras_res <- res(chosen_raster) # Resolution

# Make Dataframe
ras.d.frame <- as.data.frame(chosen_raster) # Get dataframe from raster
ras.d.frame <- data.frame(1:total_cells, ras.d.frame)
colnames(ras.d.frame)[1] <- c("Cells")

ras.occ.cells <- ras.d.frame %>%
  drop_na() 

test2$

  scale(ras.occ.cells[2])
  
  
test <- (ras.occ.cells$RAIN - means[1])

test / as.numeric(sds[1])


# vector of names 
newData <- data.frame(Cells = ras.occ.cells$Cells, RAIN = (ras.occ.cells$RAIN - means[1]) / as.numeric(sds[1]), TEMP = (ras.occ.cells$TEMP - means[2])/as.numeric(sds[2]))
predCH <- predict(best.model, type="det", newdata=newData)

predCH$Cells <- ras.occ.cells$Cells

dframe_of_cells <- data.frame(1:total_cells)
colnames(dframe_of_cells) <- "Cells"

# Join dataframes
full_dframe <- left_join(dframe_of_cells, predCH, by = "Cells")

raster_for_values <- raster(res = res, val = full_dframe$Predicted, ext = raster_extent)
plot(raster_for_values)





ras_vals <- getValues(chosen_raster)





# Make dataframe of values
extracted_values <- raster::extract(chosen_raster, vector_of_cells)
dframe_of_values <- data.frame(vector_of_cells, extracted_values)
colnames(dframe_of_values) <- c("Cells", "Vals")

# Make blank dataframe to copy values into
dframe_of_cells <- data.frame(1:total_cells)
colnames(dframe_of_cells) <- "Cells"

# Join dataframes
full_dframe <- left_join(dframe_of_cells, dframe_of_values, by = "Cells")

if (zero == TRUE){
  full_dframe[is.na(full_dframe)] <- 0
}

# Make and plot raster
raster_for_values <- raster(res = res, val = full_dframe$Vals, ext = raster_extent)
plot(raster_for_values)























test2 <- as.data.frame(test2, xy = TRUE)

# Load the Swiss landscape data from unmarked
data(Switzerland)             # Load Swiss landscape data in unmarked
CH <- Switzerland

# Get predictions of occupancy prob for each 1km2 quadrat of Switzerland
newData <- data.frame(TEMP = (test2$WC_Temp_0.5 - means[4])/sds[4], RAIN = (test2$WC_Prec_0.5 - means[2])/sds[2])
predCH <- predict(best.model, type="det", newdata=newData)

# Prepare Swiss coordinates and produce map
library(raster)
library(rgdal)

# Define new data frame with coordinates and outcome to be plotted
PARAM <- data.frame(x = test2$x, y = test2$y, z = predCH$Predicted)
r1 <- rasterFromXYZ(PARAM)     # convert into raster object

# Mask quadrats with elevation greater than 2250
elev <- rasterFromXYZ(cbind(CH$x, CH$y, CH$elevation))
elev[elev > 2250] <- NA
r1 <- mask(r1, elev)

# Plot species distribution map (Fig. 10-14 left)
par(mfrow = c(1,2), mar = c(1,2,2,5))
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
plot(r1, col = mapPalette(100), axes = F, box = F, main = "Detection Probability")

mask_oc <- raster("Data/Covariate_Data/Formatted_For_Precise/COut_0.5.asc")
projection(mask_oc) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 " 
mask_oc <- crop(mask_oc, e)
plot(mask_oc)

raster_for_values

mask_oc[(mask_oc[]) == 0] <- NA


r1 <- mask(raster_for_values, mask_oc)
mapTheme <- rasterVis::rasterTheme(region=brewer.pal(9,"YlOrRd"))
print(rasterVis::levelplot(r1, margin=F, par.settings=mapTheme,  main = "Detection Probabiliy of Campanian Hadrosaurs") + #create levelplot for raster
        latticeExtra::layer(sp.polygons(states, col = "white", fill = NA), under = T)  + # Plots state lines
        latticeExtra::layer(sp.polygons(countries, col = 0, fill = "light grey"), under = T)) # Plots background colour



pred.matrix2 <- array(NA, dim = c(100, 100)) # Define arrays
for(i in 1:100){
  for(j in 1:100){
    newData2 <- data.frame(TEMP=p.TEMP[i], RAIN=p.RAIN[j])        # For detection
    pred <- predict(best.model, type="det", newdata=newData2)
    pred.matrix2[i, j] <- pred$Predicted
  }
}

image(x=pred.TEMP, y=p.RAIN, z=pred.matrix2, col = mapPalette(100), axes = FALSE, xlab = "Temperature", ylab = "Rainfall")
contour(x=pred.TEMP, y=p.RAIN, z=pred.matrix2, add = TRUE, lwd = 1.5, col = "blue", labcex = 1.3)
axis(1, at = seq(min(pred.TEMP), max(pred.TEMP), by = 5))
axis(2, at = seq(min(pred.RAIN), max(pred.RAIN), by = 20))
box()
title(main = "Hadrosaur detection prob.", font.main = 1)
matpoints(as.matrix(data[, 2:21]), pch="+", cex=1)


