################################################################################
############ OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS #############
################################################################################

# Christopher D. Dean, Lewis A. Jones, Alfio A. Chiarenza, Sinéad Lyster, 
# Alex Farnsworth, Philip D. Mannion, Richard J. Butler.
# 2023
# Script written by Christopher D. Dean and Lewis A. Jones

################################################################################
#                 FILE 2: RUNNING OCCUPANCY MODELS IN UNMARKED                 #
################################################################################

################################################################################
# 1. INITIAL SETUP
################################################################################

#######################################
##### PACKAGES AND VARIABLE SETUP #####
#######################################

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Set your working directory

# Load Packages 
library(unmarked)
library(MuMIn)
library(dplyr)
library(AICcmodavg)
library(tibble)

# Quick filters for loading data
bin.type <- "scotese"
res <- 0.5
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
y <- data %>%
  `rownames<-`(.[,1]) %>% 
  select(-X)
# Turn into matrix
y <- as.matrix(y)

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

# Standard error function
standard_error <- function(x) sd(x) / sqrt(length(x))

####################################
##### SETUP DATA FOR MODELLING #####
####################################

##### Unstandardised values of covariates ####
# Detection
DEMs.orig <- sitecov[,"mean_DEM"]      
RAIN.orig <- sitecov[,"mean_prec"]
TEMP.orig <- sitecov[,"mean_temp"]
MGVF.orig <- sitecov[,"MGVF"]
OUTs.orig <- sitecov[,"MOut"]
OUTa.orig <- sitecov[,"AllOut"]
LAND.orig <- sitecov[,"LANDCVI_multiple"]
SEDf.orig <- sitecov[,"mean_psed"]

# Occupancy
Ptem.orig <- sitecov[,"mean_ptemp"]
Ppre.orig <- sitecov[,"mean_pprec"]
PDEM.orig <- sitecov[,"mean_pDEM"]

##### Overview of Covariates #####
#covs <- cbind(DEMs.orig, RAIN.orig, TEMP.orig, MGVF.orig, RELI.orig, OUTs.orig, 
#              OUTa.orig, Plat.orig, Ptem.orig, Ppre.orig, PDEM.orig)
#par(mfrow = c(1,1))
#for(i in 1:ncol(covs)){
#  hist(covs[,i], breaks = 50, col = "grey", main = colnames(covs)[i])
#}
#pairs(cbind(DEMs.orig, RAIN.orig, TEMP.orig, MGVF.orig, RELI.orig, OUTs.orig, 
#            OUTa.orig, Plat.orig, Ptem.orig, Ppre.orig, PDEM.orig))

##### Scaling and Mode Setup #####
siteCovs <- data.frame(DEMs = DEMs.orig, OUTs = OUTs.orig,
                       RAIN = RAIN.orig, MGVF = MGVF.orig, 
                       TEMP = TEMP.orig, LAND = LAND.orig, 
                       OUTa = OUTa.orig, SEDf = SEDf.orig,
                       Ppre = Ppre.orig, Ptem = Ptem.orig, 
                       PDEM = PDEM.orig)

# Make Model
umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs)

# Scale covariates within model
umf@siteCovs$DEMs <- scale(umf@siteCovs$DEMs)
umf@siteCovs$OUTs <- scale(umf@siteCovs$OUTs)
umf@siteCovs$OUTa <- scale(umf@siteCovs$OUTa)
umf@siteCovs$RAIN <- scale(umf@siteCovs$RAIN)
umf@siteCovs$MGVF <- scale(umf@siteCovs$MGVF)
umf@siteCovs$SEDf<- scale(umf@siteCovs$SEDf)
umf@siteCovs$Ppre <- scale(umf@siteCovs$Ppre)
umf@siteCovs$Ptem <- scale(umf@siteCovs$Ptem)
umf@siteCovs$PDEM <- scale(umf@siteCovs$PDEM)

##### Old scaling method #####
## Compute means and standard deviations
#(means <- c(apply(cbind(DEMs.orig, RAIN.orig, TEMP.orig, MGVF.orig, RELI.orig, 
# OUTs.orig, OUTa.orig, Plat.orig, Ptem.orig, Ppre.orig, PDEM.orig), 2, mean, na.rm = TRUE)))
#(sds <- c(apply(cbind(DEMs.orig, RAIN.orig, TEMP.orig, MGVF.orig, RELI.orig, 
# OUTs.orig, OUTa.orig, Plat.orig, Ptem.orig, Ppre.orig, PDEM.orig), 2, sd, na.rm = TRUE)))
#
## Scale covariates
#DEMs <- (DEMs.orig - means[1]) / sds[1] #
#RAIN <- (RAIN.orig - means[2]) / sds[2] #
#TEMP <- (TEMP.orig - means[3]) / sds[3]
#MGVF <- (MGVF.orig - means[4]) / sds[4] #
#RELI <- (RELI.orig - means[5]) / sds[5]
#OUTs <- (OUTs.orig - means[6]) / sds[6]
#OUTa <- (OUTa.orig - means[7]) / sds[7]
#Plat <- (Plat.orig - means[8]) / sds[8]
#Ptem <- (Ptem.orig - means[9]) / sds[9]
#Ppre <- (Ppre.orig - means[10]) / sds[10]
#PDEM <- (PDEM.orig - means[11]) / sds[11]


################################################################################
# 2. RUNNING MODELS
################################################################################

################################
##### NULL MODEL (NO COVS) #####
################################

# Get summary of overall occupancy model
summary(umf) 

# Make model without covariates
summary(fm1 <- occu(~1 ~1, data=umf))

# Get estimates for occupancy and detection from basic model
print(occ.null <- backTransform(fm1, "state")) 
print(det.null <- backTransform(fm1, "det")) 

# Compile results
null.res <- rbind(Null.occ.prob = c(Estimate = mean(occ.null@estimate), 
                                    confint(occ.null, level = 0.9)),
                  Null.det.prob = c(Estimate = mean(det.null@estimate), 
                                    confint(det.null, level = 0.9)))
colnames(null.res) <- c("Estimate", "5%", "95%")

###############################################
##### COVARIATE MODEL AND MODEL SELECTION #####
###############################################

# Fit full model
full <- occu(formula =  ~ OUTa + MGVF + TEMP + LAND + RAIN + SEDf # det
                        ~ Ppre + Ptem +  PDEM, # occ
             data = umf)

# Use dredge to automatically carry out model selection
(modelList <- dredge(full, rank = "AIC")) 

# If no clear best fit model, combine models for averaged model
occu_dredge_95 <- get.models(modelList, subset = cumsum(weight) <= 0.95)
oc_avg <- model.avg(occu_dredge_95, fit = TRUE)

# Report best model and top 5 models
best.model <- occu_dredge_95[[1]]
best.model
top.5 <- occu_dredge_95[1:5]
top.5

# Run MacKenzie and Bailey Goodness-of-fit test (WARNING: MIGHT TAKE A WHILE)
system.time(occ_gof <- mb.gof.test(best.model, nsim = 1000, plot.hist = FALSE))
occ_gof$p.value

# Examine the effect of covariates from averaged model
(temp.model.res <- coef(oc_avg) %>% 
  enframe())

# Examine beta co-efficents for best model
confint(best.model, type='det', method = 'normal')
confint(best.model, type='state', method = 'normal')

##### MODEL STATS #####
# Proportion of area occupied
re <- ranef(best.model)
EBUP <- bup(re, stat="mean")
CI <- confint(re, level=0.9)
SE <- standard_error(re@post)
(PAO <- rbind(PAO = c(Estimate = sum(EBUP), colSums(CI)) / nrow(y)))

# Estimate of Detection Prob. per site
det.prob <- predict(best.model, type="det") # Predict detection for sites/vists
det.prob <- det.prob[seq(1, nrow(det.prob), (ncol(data)-1)), ] # Remove duplicate vists to only get prediction per site
(det.prob <- rbind(Det.prob = c(Estimate = mean(det.prob$Predicted), 
                                SE = mean(det.prob$SE))))

# Estimate of Occupancy Prob. per site
occ.prob <- predict(best.model, type="state") # Predict occupancy for sites/vists
occ.prob <- occ.prob[seq(1, nrow(occ.prob), (ncol(data)-1)), ] # Remove duplicate vists to only get prediction per site
(occ.prob <- rbind(occ.prob = c(Estimate = mean(occ.prob$Predicted), 
                                SE = mean(occ.prob$SE))))

# Combine results
comb <- merge(det.prob, PAO, all = T)
comb <- comb[order(comb$'5%'),]
rownames(comb) <- c("PAO", "Det.prob")

################################################################################
# 3. SAVING RESULTS
################################################################################

# Make results directory
dir.create("Results/Unmarked/", showWarnings = FALSE)
dir.create(paste0("Results/Unmarked/", bin.type, "/", sep =""), showWarnings = FALSE)
dir.create(paste0("Results/Unmarked/", bin.type, "/", bin, "/", sep = ""), 
           showWarnings = FALSE)
dir.create(paste0("Results/Unmarked/", bin.type, "/", bin, "/", res, "/",
                  sep =""), showWarnings = FALSE)

# Null results and model stats
combined.res <- merge(null.res, comb, all = T, sort = F)
combined.res$Bin <- bin
rownames(combined.res) <- c("null.occ.prob", "null.det.prob", "PAO", "mean.det.prob")
write.csv(combined.res, paste0("Results/Unmarked/", bin.type, "/", bin, "/", res, "/", 
                               target, ".combined.results.csv", sep =""))
# Average model table
write.csv(temp.model.res, paste0("Results/Unmarked/", bin.type, "/", bin, "/", res, "/", 
                               target, ".model.avergaed.covariates.csv", sep =""))
# Full model list
write.csv(modelList, paste0("Results/Unmarked/", bin.type, "/", bin, "/", res, "/", 
                                 target, ".full.mod.table.csv", sep =""))

################################################################################
# 4. PREDICTION
################################################################################

# Create new covariates for prediction ('prediction covs')
pred.DEMs <- seq(min(DEMs.orig), max(DEMs.orig),,100) # New covs for prediction
pred.RAIN <- seq(min(RAIN.orig), max(RAIN.orig),,100)
pred.TEMP <- seq(min(TEMP.orig), max(TEMP.orig),,100) 
pred.MGVF <- seq(min(MGVF.orig), max(MGVF.orig),,100)
pred.RELI <- seq(min(RELI.orig), max(RELI.orig),,100)
pred.OUTs <- seq(min(OUTs.orig), max(OUTs.orig),,100)
pred.OUTa <- seq(min(OUTa.orig), max(OUTa.orig),,100)

p.DEMs <- (pred.DEMs - means[1]) / sds[1] # Standardise them like actual covs
p.RAIN <- (pred.RAIN - means[2]) / sds[2]
p.TEMP <- (pred.TEMP - means[3]) / sds[3] 
p.MGVF <- (pred.MGVF - means[4]) / sds[4]
p.RELI <- (pred.RELI - means[5]) / sds[5]
p.OUTs <- (pred.OUTs - means[6]) / sds[6]
p.OUTa <- (pred.OUTa - means[7]) / sds[7]

pred.Plat <- seq(min(Plat.orig), max(Plat.orig),,100) # New covs for prediction
pred.Ptem <- seq(min(Ptem.orig), max(Ptem.orig),,100) # New covs for prediction
pred.Ppre <- seq(min(Ppre.orig), max(Ppre.orig),,100) # New covs for prediction
pred.PDEM <- seq(min(PDEM.orig), max(PDEM.orig),,100) # New covs for prediction

p.Plat <- (pred.Plat - means[8]) / sds[8]
p.Ptem <- (pred.Ptem - means[9]) / sds[9]
p.Ppre <- (pred.Ppre - means[10]) / sds[10]
p.PDEM <- (pred.PDEM = means[11]) / sds[11]

# Obtain predictions
newData <- data.frame(RELI=p.RELI, TEMP = 0)
pred.det.RELI <- predict(best.model, type="det", newdata=newData, appendData=TRUE)
newData <- data.frame(RELI=0, TEMP = p.TEMP)
pred.det.TEMP <- predict(best.model, type="det", newdata=newData, appendData=TRUE)

newData <- data.frame(Ppre=p.Ppre)
pred.occ.Ppre <- predict(best.model, type="state", newdata=newData, appendData=TRUE)

newData <- data.frame(Plat=p.Plat, PDEM = 0, Ppre = 0)
pred.occ.Plat <- predict(best.model, type="state", newdata=newData, appendData=TRUE)
newData <- data.frame(Plat = 0, PDEM = p.PDEM, Ppre = 0)
pred.occ.PDEM <- predict(best.model, type="state", newdata=newData, appendData=TRUE)
newData <- data.frame(Plat = 0, PDEM = 0, Ppre = p.Ppre)
pred.occ.Ppre <- predict(best.model, type="state", newdata=newData, appendData=TRUE)

# Plot predictions against unstandardized 'prediction covs'
par(mfrow = c(3,1), mar = c(5,5,2,3), cex.lab = 1.2)
plot(pred.det.RELI[[1]] ~ pred.RELI, type = "l", lwd = 3, col = "blue", ylim = c(0,1), las = 1, ylab = "Pred. detection prob.", xlab = "Relief", frame = F)
matlines(pred.RELI, pred.det.RELI[,3:4], lty = 1, lwd = 1, col = "grey")
plot(pred.det.TEMP[[1]] ~ pred.TEMP, type = "l", lwd = 3, col = "blue", ylim = c(0,1), las = 1, ylab = "Pred. detection prob.", xlab = "Temp", frame = F)
matlines(pred.TEMP, pred.det.TEMP[,3:4], lty = 1, lwd = 1, col = "grey")

plot(pred.occ.Plat[[1]] ~ pred.Plat, type = "l", lwd = 3, col = "blue", ylim = c(0,1), las = 1, ylab = "Pred. occupancy prob.", xlab = "Palaeo Latitude", frame = F)
matlines(pred.Plat, pred.occ.Plat[,3:4], lty = 1, lwd = 1, col = "grey")
plot(pred.occ.PDEM[[1]] ~ pred.PDEM, type = "l", lwd = 3, col = "blue", ylim = c(0,1), las = 1, ylab = "Pred. occupancy prob.", xlab = "Palaeo Temp.", frame = F)
matlines(pred.PDEM, pred.occ.PDEM[,3:4], lty = 1, lwd = 1, col = "grey")
plot(pred.occ.Ppre[[1]] ~ pred.Ppre, type = "l", lwd = 3, col = "blue", ylim = c(0,1), las = 1, ylab = "Pred. occupancy prob.", xlab = "Palaeo Precip.", frame = F)
matlines(pred.Ppre, pred.occ.Ppre[,3:4], lty = 1, lwd = 1, col = "grey")

#Look at the occupancy and detection estimates on probability scale
#Note that 'backTransform' will only work without intermediate steps if inquiring about the null model.

lc <- linearComb(best.model, c(1, 50, 10), type="det") # Intercept=1, bosq=2

backTransform(lc)

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



#===== Old manual model selection =====
#### DETECTION ###
#summary(fm0 <- occu(~1 ~1, data=umf))
#
#summary(fm1 <- occu(~RELI ~1, data=umf))
#summary(fm2 <- occu(~RAIN ~1, data=umf))
#summary(fm3 <- occu(~TEMP ~1, data=umf))
#summary(fm4 <- occu(~MGVF ~1, data=umf))
#summary(fm5 <- occu(~OUTa ~1, data=umf))
#
#summary(fm6 <- occu(~RELI + RAIN ~1, data=umf))
#summary(fm7 <- occu(~RELI + TEMP ~1, data=umf))
#summary(fm8 <- occu(~RELI + MGVF ~1, data=umf))
#summary(fm9 <- occu(~RELI + OUTa ~1, data=umf))
#summary(fm10 <- occu(~RAIN + TEMP ~1, data=umf))
#summary(fm11 <- occu(~RAIN + MGVF ~1, data=umf))
#summary(fm12 <- occu(~RAIN + OUTa ~1, data=umf))
#summary(fm13 <- occu(~TEMP + MGVF ~1, data=umf))
#summary(fm14 <- occu(~TEMP + OUTa ~1, data=umf))
#summary(fm15 <- occu(~MGVF + OUTa ~1, data=umf))
#
#summary(fm16 <- occu(~RELI + RAIN + TEMP ~1, data=umf))
#summary(fm17 <- occu(~RELI + RAIN + MGVF ~1, data=umf))
#summary(fm18 <- occu(~RELI + RAIN + OUTa ~1, data=umf))
#summary(fm19 <- occu(~RELI + TEMP + MGVF ~1, data=umf))
#summary(fm20 <- occu(~RELI + TEMP + OUTa ~1, data=umf))
#summary(fm21 <- occu(~RELI + MGVF + OUTa ~1, data=umf))
#summary(fm22 <- occu(~RAIN + TEMP + MGVF ~1, data=umf))
#summary(fm23 <- occu(~RAIN + TEMP + OUTa ~1, data=umf))
#summary(fm24 <- occu(~RAIN + MGVF + OUTa ~1, data=umf))
#summary(fm25 <- occu(~TEMP + MGVF + OUTa ~1, data=umf))
#
#summary(fm26 <- occu(~RELI + RAIN + TEMP + MGVF ~1, data=umf))
#summary(fm27 <- occu(~RELI + RAIN + TEMP + OUTa ~1, data=umf))
#summary(fm28 <- occu(~RELI + RAIN + MGVF + OUTa ~1, data=umf))
#summary(fm29 <- occu(~RELI + TEMP + MGVF + OUTa ~1, data=umf))
#summary(fm30 <- occu(~RAIN + TEMP + MGVF + OUTa ~1, data=umf))
#
#summary(fm31 <- occu(~RELI + RAIN + TEMP + MGVF + OUTa  ~1, data=umf)) #
#
## Put the fitted models in a "fitList" and rank them by AIC
#fms <- fitList("0" = fm0,
#               "1" = fm1,
#               "2" = fm2, 
#               "3" = fm3,
#               "4" = fm4, 
#               "5" = fm5,
#               "6" = fm6, 
#               "7" = fm7,
#               "8" = fm8,
#               "9" = fm9,
#               "10" = fm10,
#               "11" = fm11, 
#               "12" = fm12,
#               "13" = fm13,
#               "14" = fm14,
#               "15" = fm15,
#               "16" = fm16,
#               "17" = fm17,
#               "18" = fm18,
#               "19" = fm19,
#               "20" = fm20,
#               "21" = fm21,
#               "22" = fm22,
#               "23" = fm23,
#               "24" = fm24,
#               "25" = fm25,
#               "26" = fm26,
#               "27" = fm27,
#               "28" = fm28,
#               "29" = fm29,
#               "30" = fm30, 
#               "31" = fm31
#               )
#(ms <- modSel(fms))
#
## Further tests
#summary(fm1 <- occu(~MGVF + TEMP ~1, data=umf))
#summary(fm2 <- occu(~MGVF + I(MGVF^2) + TEMP ~1, data=umf))
#summary(fm3 <- occu(~MGVF + I(MGVF^2) + I(MGVF^3) + TEMP ~1, data=umf))
#summary(fm4 <- occu(~MGVF + TEMP + I(TEMP^2) ~1, data=umf))
#summary(fm5 <- occu(~MGVF + TEMP + I(TEMP^2) +I(TEMP^3) ~1, data=umf))
#summary(fm6 <- occu(~MGVF + I(MGVF^2) + TEMP + I(TEMP^2) ~1, data=umf))
#summary(fm7 <- occu(~MGVF + I(MGVF^2) + I(MGVF^3) + TEMP + I(TEMP^2) ~1, data=umf))
#summary(fm8 <- occu(~MGVF + I(MGVF^2) + TEMP + I(TEMP^2) + I(TEMP^3) ~1, data=umf))
#summary(fm9 <- occu(~MGVF + I(MGVF^2) + I(MGVF^3) + TEMP + I(TEMP^2) + I(TEMP^3) ~1, data=umf))
#
#fms <- fitList("p(MGVF+TEMP)psi(.)"                                  = fm1, # BEST
#               "p(MGVF+MGVF2+TEMP)psi(.)"                            = fm2,
#               "p(MGVF+MGVF2+MGVF3+TEMP)psi(.)"                      = fm3, 
#               "p(MGVF+TEMP+TEMP2)psi(.)"                            = fm4,
#               "p(MGVF+TEMP+TEMP2+TEMP3)psi(.)"                      = fm5,
#               "p(MGVF+MGVF2+TEMP+TEMP2)psi(.)"                      = fm6,
#               "p(MGVF+MGVF2+MGVF3+TEMP+TEMP2)psi(.)"                = fm7, 
#               "p(MGVF+MGVF2+TEMP+TEMP2+TEMP3)psi(.)"                = fm8,
#               "p(MGVF+MGVF2+MGVF3+TEMP+TEMP2+TEMP3)psi(.)"          = fm9
#)
#(ms <- modSel(fms))
#
## Further tests
#summary(fm1 <- occu(~MGVF + RAIN ~1, data=umf))
#summary(fm2 <- occu(~MGVF + I(MGVF^2) + RAIN ~1, data=umf))
#summary(fm3 <- occu(~MGVF + I(MGVF^2) + I(MGVF^3) + RAIN ~1, data=umf))
#summary(fm4 <- occu(~MGVF + RAIN + I(RAIN^2) ~1, data=umf))
#summary(fm5 <- occu(~MGVF + RAIN + I(RAIN^2) +I(RAIN^3) ~1, data=umf))
#summary(fm6 <- occu(~MGVF + I(MGVF^2) + RAIN + I(RAIN^2) ~1, data=umf))
#summary(fm7 <- occu(~MGVF + I(MGVF^2) + I(MGVF^3) + RAIN + I(RAIN^2) ~1, data=umf))
#summary(fm8 <- occu(~MGVF + I(MGVF^2) + RAIN + I(RAIN^2) + I(RAIN^3) ~1, data=umf))
#summary(fm9 <- occu(~MGVF + I(MGVF^2) + I(MGVF^3) + RAIN + I(RAIN^2) + I(RAIN^3) ~1, data=umf))
#
#fms <- fitList("p(MGVF+RAIN)psi(.)"                                  = fm1, # BEST
#               "p(MGVF+MGVF2+RAIN)psi(.)"                            = fm2,
#               "p(MGVF+MGVF2+MGVF3+RAIN)psi(.)"                      = fm3, 
#               "p(MGVF+RAIN+RAIN2)psi(.)"                            = fm4,
#               "p(MGVF+RAIN+RAIN2+RAIN3)psi(.)"                      = fm5,
#               "p(MGVF+MGVF2+RAIN+RAIN2)psi(.)"                      = fm6,
#               "p(MGVF+MGVF2+MGVF3+RAIN+RAIN2)psi(.)"                = fm7, 
#               "p(MGVF+MGVF2+RAIN+RAIN2+RAIN3)psi(.)"                = fm8,
#               "p(MGVF+MGVF2+MGVF3+RAIN+RAIN2+RAIN3)psi(.)"          = fm9
#)
#(ms <- modSel(fms))
#
#best.model <- fm7
#confint(best.model, type = "det")
#
#### OCCUPANCY ###
#
#summary(fm1 <- occu(~RELI + TEMP ~1, data=umf))
#
#summary(fm2 <- occu(~RELI + TEMP ~Plat, data=umf))
#summary(fm3 <- occu(~RELI + TEMP ~Ptem, data=umf))
#summary(fm4 <- occu(~RELI + TEMP ~Ppre, data=umf))
#summary(fm5 <- occu(~RELI + TEMP ~PDEM, data=umf))
#
#summary(fm5 <- occu(~RELI + TEMP ~Plat + Ptem, data=umf))
#summary(fm6 <- occu(~RELI + TEMP ~Plat + Ppre, data=umf))
#summary(fm7 <- occu(~RELI + TEMP ~Plat + PDEM, data=umf))
#summary(fm8 <- occu(~RELI + TEMP ~Ptem + Ppre, data=umf))
#summary(fm9 <- occu(~RELI + TEMP ~Ptem + PDEM, data=umf))
#summary(fm10 <- occu(~RELI + TEMP ~Ppre + PDEM, data=umf))
#
#summary(fm11 <- occu(~RELI + TEMP ~Plat + Ptem + Ppre, data=umf))
#summary(fm12 <- occu(~RELI + TEMP ~Plat + Ptem + PDEM, data=umf))
#summary(fm13 <- occu(~RELI + TEMP ~Plat + Ppre + PDEM, data=umf))
#summary(fm14 <- occu(~RELI + TEMP ~Ptem + Ppre + PDEM, data=umf))
#
#summary(fm15 <- occu(~RELI + TEMP ~Plat + Ptem + Ppre + PDEM, data=umf))
#
#fms <- fitList("1" = fm1,
#               "2" = fm2,
#               "3" = fm3,
#               "4" = fm4,
#               "5" = fm5,
#               "6" = fm6,
#               "7" = fm7,
#               "8" = fm8, 
#               "9" = fm9, 
#               "10" = fm10,
#               "11" = fm11,
#               "12" = fm12,
#               "13" = fm13,
#               "14" = fm14, 
#               "15" = fm15
#)
#
#(ms <- modSel(fms))
#
#summary(fm1 <- occu(~MGVF + TEMP + RAIN ~plat + I(plat^2), data=umf)) # Null best
#summary(fm2 <- occu(~MGVF + TEMP + RAIN ~plat + I(plat^2) + I(plat^3), data=umf))
#summary(fm3 <- occu(~MGVF + TEMP + RAIN ~Ctem + I(Ctem^2), data=umf))
#summary(fm4 <- occu(~MGVF + TEMP + RAIN ~Ctem + I(Ctem^2) + I(Ctem^3), data=umf))
#summary(fm5 <- occu(~MGVF + TEMP + RAIN ~plat + Ctem + I(Ctem^2), data=umf))
#summary(fm6 <- occu(~MGVF + TEMP + RAIN ~plat + Ctem + I(plat^2), data=umf))
#summary(fm7 <- occu(~MGVF + TEMP + RAIN ~plat + Ctem + I(plat^2) + I(Ctem^2), data=umf))
#summary(fm8 <- occu(~MGVF + TEMP + RAIN ~1, data=umf))
#
#
#summary(fm1 <- occu(~MGVF + I(MGVF^2) + RAIN + I(RAIN^2) + I(RAIN^3) ~Cpre + I(Cpre^2), data=umf)) # Null best
#summary(fm2 <- occu(~MGVF + I(MGVF^2) + RAIN + I(RAIN^2) + I(RAIN^3) ~Cpre + I(Cpre^2) + I(Cpre^3), data=umf))
#summary(fm3 <- occu(~MGVF + I(MGVF^2) + RAIN + I(RAIN^2) + I(RAIN^3) ~plat + I(plat^2), data=umf))
#summary(fm4 <- occu(~MGVF + I(MGVF^2) + RAIN + I(RAIN^2) + I(RAIN^3) ~plat + I(plat^2) + I(plat^3), data=umf))
#summary(fm5 <- occu(~MGVF + I(MGVF^2) + RAIN + I(RAIN^2) + I(RAIN^3) ~Cpre + plat + I(plat^2), data=umf))
#summary(fm6 <- occu(~MGVF + I(MGVF^2) + RAIN + I(RAIN^2) + I(RAIN^3) ~Cpre + plat + I(Cpre^2), data=umf))
#summary(fm7 <- occu(~MGVF + I(MGVF^2) + RAIN + I(RAIN^2) + I(RAIN^3) ~Cpre + plat + I(Cpre^2) + I(plat^2), data=umf))
#summary(fm8 <- occu(~MGVF + I(MGVF^2) + RAIN + I(RAIN^2) + I(RAIN^3) ~1, data=umf))
#
#fms <- fitList("p(MGVF+MGVF2+RAIN+RAIN2+RAIN3)psi(Cpre+2Cpre)"            = fm1,
#               "p(MGVF+MGVF2+RAIN+RAIN2+RAIN3)psi(Cpre+2Cpre+3Cpre)"      = fm2,
#               "p(MGVF+MGVF2+RAIN+RAIN2+RAIN3)psi(Ctem+2Ctem)"            = fm3, # BEST
#               "p(MGVF+MGVF2+RAIN+RAIN2+RAIN3)psi(Ctem+2Ctem+3Ctem)"      = fm4,
#               "p(MGVF+MGVF2+RAIN+RAIN2+RAIN3)psi(Cpre+Ctem+2CaTw)"       = fm5,
#               "p(MGVF+MGVF2+RAIN+RAIN2+RAIN3)psi(Cpre+Ctem+2Cpre)"       = fm6,
#               "p(MGVF+MGVF2+RAIN+RAIN2+RAIN3)psi(Cpre+Ctem+2Cpre+2Ctem)" = fm7, 
#               "Null"                           = fm8
#)
#(ms <- modSel(fms))

#best.model <- fm4

#confint(best.model, type = "state")
