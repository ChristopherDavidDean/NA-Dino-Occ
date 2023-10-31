################################################################################
############ OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS #############
################################################################################

# Christopher D. Dean, Lewis A. Jones, Alfio A. Chiarenza, Sinéad Lyster, 
# Alex Farnsworth, Philip D. Mannion, Richard J. Butler.
# 2023
# Script written by Christopher D. Dean

################################################################################
#               FILE 2: RUNNING OCCUPANCY MODELS IN SPARTA.                    #
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

##### Load in Functions #####
source("0.Functions.R") # Import functions from other R file (must be in same working directory)

# Set resolution
res <- 1

# Set extent
e <- extent(-155, -72, 22.5, 73)

# Set bin type
bin.type <- "formation"

# Load main dataset
master.occs <- read.csv(paste("Prepped_data/Occurrence_Data/", bin.type, "/", 
                              res, "_", bin.type, "_occurrence_dataset.csv", sep = ""))
# Load formations
formations <- read.csv("Data/Occurrences/Formations.csv") # Load formations

# Set bins 
bins <- master.occs %>%
  dplyr::select(bin_assignment, bin_midpoint) %>%
  distinct()

get_grid(master.occs, res = res, e = e)

# Adding midpoint
master.occs.grid$mid_ma <- (master.occs.grid$max_ma + master.occs.grid$min_ma)/2

lookup <- data.frame('currentbins' = sort(unique(master.occs.grid$bin_assignment)), 
                     'newbins' = seq(from = length(unique(master.occs.grid$bin_assignment)), 
                                     to = 1))
inds <- match(master.occs.grid$bin_assignment, lookup$currentbins)
master.occs.grid$new_bins[!is.na(inds)] <- lookup$newbins[na.omit(inds)]

################################################################################
# 3. NAIVE OCCUPANCY
################################################################################

# Make a list of target taxa
target = c("Hadrosauridae", "Ceratopsidae", "Tyrannosauridae")

# Make results table for naive occupancy
naive.res(target, master.occs.grid)

################################################################################
# 4. RUNNING SPARTA
################################################################################

# c('sparta', 'catlistlength')
# c('ranwalk', 'halfcauchy', 'catlistlength')

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
all.results <- run.model(master.occs.grid, target)

################################################################################
# 5. PLOTTING/SAVING RESULTS
################################################################################

#################
##### PLOTS #####
#################

# Plotting number of occurrences
occurrence.plot(master.occs.grid, target)

# Plotting occupancy (naive and modelled)
cera <- all.results[[1]] %>%
  filter(Target == "Ceratopsidae")
tyran <- all.results[[1]] %>%
  filter(Target == "Tyrannosauridae")
hadro <- all.results[[1]] %>%
  filter(Target == "Hadrosauridae")

# Plot modelled results
a <- plot.occ(hadro)
b <- plot.naive(hadro)
c <- plot.occ(tyran)
d <- plot.naive(tyran)
e <- plot.occ(cera)
f <- plot.naive(cera)

# Arrange
p <- ggarrange(b, d, f, a, c, e,
          nrow = 2, ncol = 3,
          align='h', labels=c('A', 'B', 'C',
                              'D', 'E', 'F'),
          legend = "bottom",
          common.legend = T)

#if(bin.type =="formation"){
#  pdf(paste("Results/Outhwaite/", bin.type, "/", bin.res, 
#            "/Plot_", res, ".pdf", sep = ""), width = 11.458, height = 7.292)
#} else{
#  pdf(paste("Results/Outhwaite/", bin.type, "/Plot_", 
#            res, ".pdf", sep = ""), width = 11.458, height = 7.292)
#}

# Add phylopics
(p1 <- cowplot::ggdraw() +  
  cowplot::draw_plot(p) +
  cowplot::draw_image("https://images.phylopic.org/images/aeeb30a8-afdc-4e7e-9bcc-574cb290a1f6/raster/1536x575.png?build=140", 
             x = 0.435, y = 0.44, scale = 0.1) +
  cowplot::draw_image("https://images.phylopic.org/images/f3808e65-a95f-4df5-95a0-5f5b46a221f2/raster/1536x505.png?build=140", 
           x = 0.09, y = 0.44, scale = 0.12) +
  cowplot::draw_image("https://images.phylopic.org/images/72be89b9-3f2b-4dc3-b485-e74a5f8b1fbc/raster/1536x512.png?build=140", 
             x = -0.23, y = 0.44, scale = 0.1))

#dev.off()

# Choose type of model
type <- "s"

# Save figure
ggsave(paste("Figures/1.Sparta.", bin.type, ".", res, ".", type, ".png", sep = ""), plot = p1, 
       device = "png", type = "cairo")

########################
##### SAVE RESULTS #####
########################

if(bin.type =="formation"){
  write.csv(results, paste("Results/Outhwaite/", bin.type, "/", 
                           bin.res, "/naive.results.", res, ".csv", sep = ""))
  write.csv(all.results[[1]], paste("Results/Outhwaite/", bin.type,
                               "/", bin.res, "/results.", res, ".csv", sep = ""))
  saveRDS(all.results[[2]], file = paste("Results/Outhwaite/Posterior_checks/",
                                         bin.type, "/", res, ".models.rds", 
                                         sep = ""))
}else{
  write.csv(results, paste("Results/Outhwaite/", bin.type,
                           "/naive.results", res, ".csv", sep = ""))
  write.csv(all.results[[1]], paste("Results/Outhwaite/", bin.type, 
                               "/results", res, ".csv", sep = ""))
  saveRDS(all.results[[2]], file = paste("Results/Outhwaite/Posterior_checks/",
                                         bin.type, "/", res, ".models.rds", 
                                         sep = ""))
}
