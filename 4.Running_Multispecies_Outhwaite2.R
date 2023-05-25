################################################################################
############ OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS #############
################################################################################

# Christopher D. Dean, Lewis A. Jones, Alfio A. Chiarenza, Sinéad Lyster, 
# Alex Farnsworth, Philip D. Mannion, Richard J. Butler.
# 2023
# Script written by Christopher D. Dean

################################################################################
#               FILE 4: RUNNING OCCUPANCY MODELS IN SPARTA.                    #
################################################################################

################################################################################
# 1. INITIAL SETUP
################################################################################

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Set your working directory

# Load Packages 
#install_github('BiologicalRecordsCentre/sparta')
#install.packages("R2jags")
library(sparta)
library(rphylopic)
library(gtools)
library(R2jags)
library(palaeoverse)
library(ggplot2)
library(snowfall)
library(ggnewscale)
library(divDyn)
data(stages)
stage <- stages
library(deeptime)
library(wesanderson)
library(ggthemes)
library(ggpubr)
library(cowplot)
library(magick)
library(RCurl)
library(png)
library(grid)
library(viridis)

# Get Functions
source('0.Functions_DR.R')

# Load data 
res <- 0.5
master.occs <- read.csv(paste("master.occs.", res, ".csv", sep = ""))

# Load formations
formations <- read.csv("Data/Occurrences/Formations.csv") # Load formations

# Get information on ICS based intervals, then trim to fit relevant time frame. 
stages <- stage[80:81,] # Set stages to range from Campanian to Maastrichtian

################################################################################
# 2. TIME BINNING
################################################################################

#####################
##### FORMATION #####
#####################

# Organise formations
formations <- formations[,1:11] # Remove additional columns
formations <- formations[order(formations$Formation),] # Reorganise formations
formations$forbinning <- 1:nrow(formations) # Provide ID for formations
formations$Range <- formations$max_age - formations$min_age # Calculate formation range
formations$Diversity <- 0 # Add in dummy variables to ensure code works (sorry!)
formations$Occurrences <- 0 # Add in dummy variables to ensure code works (sorry!)
colnames(formations)[1] <- "formation" # Change to allow for further analysis

# Run combined binning function, choosing adjustable window
bin.res <- 2.5
binning(bin.res)
bin.type <- "formation"

####################
##### SUBSTAGE #####
####################

# Load bins, bin accordingly
bins <- read.csv("Data/Occurrences/substages.csv")
bins$bin <- c("SB.1","SB.2","SB.3","SB.4","SB.5")
master.occs.binned <- bin_time(master.occs, bins, method = "all")
bin.type <- "substage"

# Adding midpoint
master.occs.binned$mid_ma <- (master.occs.binned$max_ma + master.occs.binned$min_ma)/2

# Reorder bins for later
lookup <- data.frame('currentbins' = sort(unique(master.occs.binned$bin_assignment)), 
                     'newbins' = seq(from = length(unique(master.occs.binned$bin_assignment)), 
                                     to = 1))
inds <- match(master.occs.binned$bin_assignment, lookup$currentbins)
master.occs.binned$new_bins[!is.na(inds)] <- lookup$newbins[na.omit(inds)]

###################
##### SCOTESE #####
###################

# Load bins, bin accordingly
bins <- read.csv("Data/Occurrences/scotesebins.csv")
bins$bin <- c("SC.1","SC.2","SC.3","SC.4")
master.occs.binned <- bin_time(master.occs, bins, method = "majority")
bin.type <- "scotese"

# Adding midpoint
master.occs.binned$mid_ma <- (master.occs.binned$max_ma + master.occs.binned$min_ma)/2

# Reorder bins for later
lookup <- data.frame('currentbins' = sort(unique(master.occs.binned$bin_assignment)), 
                     'newbins' = seq(from = length(unique(master.occs.binned$bin_assignment)), 
                                     to = 1))
inds <- match(master.occs.binned$bin_assignment, lookup$currentbins)
master.occs.binned$new_bins[!is.na(inds)] <- lookup$newbins[na.omit(inds)]

################################################################################
# 3. NAIVE OCCUPANCY
################################################################################

# Make a list of target taxa
target = c("Hadrosauridae", "Ceratopsidae", "Tyrannosauridae")

# Make results table for naive occupancy
naive.res(target, master.occs.binned)

################################################################################
# 4. RUNNING SPARTA
################################################################################

#==== Run models =====
# Specify parameters within a function that takes a species name and runs the model
occ_mod_function <- function(taxa_name){
  occ_out <- sparta::occDetFunc(taxa_name = taxa_name,
                                n_iterations = 40000,
                                nyr = 2,
                                burnin = 20000, 
                                modeltype = c('sparta', 'catlistlength'),
                                occDetdata = formattedOccData$occDetdata,
                                spp_vis = formattedOccData$spp_vis,
                                write_results = FALSE)  
} 

# Run occupancy model and combine results
run.model(master.occs.binned, target)

################################################################################
# 5. PLOTTING/SAVING RESULTS
################################################################################

#################
##### PLOTS #####
#################

# Plotting number of occurrences
occurrence.plot(master.occs.binned, target)

# Plotting occupancy (naive and modelled)
cera <- all.results %>%
  filter(Target == "Ceratopsidae")
tyran <- all.results %>%
  filter(Target == "Tyrannosauridae")
hadro <- all.results %>%
  filter(Target == "Hadrosauridae")

# Plot modelled results
a <- plot.occ(hadro)
b <- plot.naive(hadro)
c <- plot.occ(tyran)
d <- plot.naive(tyran)
e <- plot.occ(cera)
f <- plot.naive(cera)

# Arrange
p <- ggarrange(a, c, e, b, d, f,
          nrow = 2, ncol = 3,
          align='h', labels=c('A', 'B', 'C',
                              'D', 'E', 'F'),
          legend = "bottom",
          common.legend = T)

if(bin.type =="formation"){
  pdf(paste("Results/Outhwaite/", bin.type, "/", bin.res, 
            "/Plot_", res, ".pdf", sep = ""), width = 11.458, height = 7.292)
} else{
  pdf(paste("Results/Outhwaite/", bin.type, "/Plot_", 
            res, ".pdf", sep = ""), width = 11.458, height = 7.292)
}

# Add phylopics
cowplot::ggdraw() +  
  cowplot::draw_plot(p) +
  cowplot::draw_image("https://images.phylopic.org/images/aeeb30a8-afdc-4e7e-9bcc-574cb290a1f6/raster/1536x575.png?build=140", 
             x = 0.435, y = -0.04, scale = 0.1) +
  cowplot::draw_image("https://images.phylopic.org/images/f3808e65-a95f-4df5-95a0-5f5b46a221f2/raster/1536x505.png?build=140", 
           x = 0.09, y = -0.04, scale = 0.12) +
  cowplot::draw_image("https://images.phylopic.org/images/72be89b9-3f2b-4dc3-b485-e74a5f8b1fbc/raster/1536x512.png?build=140", 
             x = -0.23, y = -0.04, scale = 0.1) 

dev.off()

########################
##### SAVE RESULTS #####
########################

if(bin.type =="formation"){
  write.csv(results, paste("Results/Outhwaite/", bin.type, "/", 
                           bin.res, "/naive.results", res, ".csv", sep = ""))
  write.csv(all.results, paste("Results/Outhwaite/", bin.type,
                               "/", bin.res, "/results", res, ".csv", sep = ""))
}else{
  write.csv(results, paste("Results/Outhwaite/", bin.type,
                           "/naive.results", res, ".csv", sep = ""))
  write.csv(all.results, paste("Results/Outhwaite/", bin.type, 
                               "/results", res, ".csv", sep = ""))
}
