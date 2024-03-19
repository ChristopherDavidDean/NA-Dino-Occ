#==============================================================================#
#============= OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS ==========#
#==============================================================================#

# Christopher D. Dean, Lewis A. Jones, Alfio A. Chiarenza, Sinéad Lyster, 
# Alex Farnsworth, Philip D. Mannion, Richard J. Butler.
# 2019
# Script written by Christopher D. Dean

#==============================================================================#
#=============== FILE 4: RUNNING OCCUPANCY MODELS IN SPARTA ===================#
#==============================================================================#

#============================ INITIAL SETUP ====================================

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Set your working directory

#===== Load Packages =====
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

#==== Get Functions ====
source('0.Functions_DR.R')

#==== Load data =====
master.occs <- read.csv("master.occs.0.5.csv")

# Load formations
formations <- read.csv("Data/Occurrences/Formations.csv") # Load formations
# Organise formations
formations <- formations[,1:11] # Remove additional columns
formations <- formations[order(formations$Formation),] # Reorganise formations
formations$forbinning <- 1:nrow(formations) # Provide ID for formations
formations$Range <- formations$max_age - formations$min_age # Calculate formation range
formations$Diversity <- 0 # Add in dummy variables to ensure code works (sorry!)
formations$Occurrences <- 0 # Add in dummy variables to ensure code works (sorry!)
colnames(formations)[1] <- "formation" # Change to allow for further analysis

# Get information on ICS based intervals, then trim to fit relevant time frame. 
stages <- stage[80:81,] # Set stages to range from Campanian to Maastrichtian

#============================ FORMATION BINS ===================================

# Run combined binning function, choosing adjustable window
binning(2)

#============================ NAIVE OCCUPANCY ==================================

# Make a list of target taxa
target = c("Hadrosauridae", "Ceratopsidae", "Tyrannosauridae")

# Make results table for naive occupancy
naive.res(target, master.occs.binned)

#============================ RUNNING MODEL ====================================

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

#============================ PLOTTING RESULTS =================================

cera <- all.results %>%
  filter(Target == "Ceratopsidae")
tyran <- all.results %>%
  filter(Target == "Tyrannosauridae")
hadro <- all.results %>%
  filter(Target == "Hadrosauridae")


# Get phylopics



library(cowplot)
library(magick)

library(ggplot2)
library(rphylopic)
library(RCurl)
library(png)
library(grid)

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

ggdraw() +  
  draw_plot(p) +
  draw_image("https://images.phylopic.org/images/72be89b9-3f2b-4dc3-b485-e74a5f8b1fbc/raster/1536x512.png?build=140", 
             x = 0.435, y = -0.04, scale = 0.1) +
  draw_image("https://images.phylopic.org/images/f3808e65-a95f-4df5-95a0-5f5b46a221f2/raster/1536x505.png?build=140", 
           x = 0.09, y = -0.04, scale = 0.12) +
  draw_image("https://images.phylopic.org/images/aeeb30a8-afdc-4e7e-9bcc-574cb290a1f6/raster/1536x575.png?build=140", 
             x = -0.23, y = -0.04, scale = 0.1) 

colnames(results)[1] <- "new_bins"
results$group <- "Naive"
test <- merge(results, all.results, all = T)

write.csv()  
