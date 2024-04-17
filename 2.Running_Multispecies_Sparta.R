################################################################################
############ OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS #############
################################################################################

# Christopher D. Dean, Alfio Alessandro Chiarenza, Jeffrey W. Doser, Alexander
# Farnsworth, Lewis A. Jones, Sinéad Lyster, Charlotte L. Outhwaite, Richard J. 
# Butler, Philip D. Mannion.
# 2024
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

# Set formation based cells or not
form_cell <- "N"

# Set model type; 1 for sparta, 2 for half cauchy random walk.
type <- 2

# Load main dataset
master.occs <- read.csv(paste("Prepped_data/Occurrence_Data/", bin.type, "/", 
                              res, "_", bin.type, "_occurrence_dataset.csv", sep = ""))
# Load formations
formations <- read.csv("Data/Occurrences/Formations.csv") # Load formations

# Set bins 
bins <- master.occs %>%
  dplyr::select(bin_assignment, bin_midpoint) %>%
  distinct()

# Make grid
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
naive_res(target, master.occs.grid)

################################################################################
# 4. RUNNING SPARTA
################################################################################

model.list <- list(c('sparta', 'catlistlength'), 
                   c('ranwalk', 'halfcauchy', 'catlistlength'))

#==== Run models =====
# Specify parameters within a function that takes a species name and runs the model
occ_mod_function <- function(taxa_name){
  occ_out <- sparta::occDetFunc(taxa_name = taxa_name,
                                n_iterations = 40000,
                                nyr = 2,
                                burnin = 20000, 
                                modeltype = model.list[[type]],
                                occDetdata = formattedOccData$occDetdata,
                                spp_vis = formattedOccData$spp_vis,
                                write_results = FALSE)  
} 

# Run occupancy model and combine results
all.results <- run_model(master.occs.grid, target)

################################################################################
# 5. PLOTTING/SAVING RESULTS
################################################################################

#################
##### PLOTS #####
#################

# Plotting occupancy (naive and modelled)
cera <- all.results[[1]] %>%
  filter(Target == "Ceratopsidae")
tyran <- all.results[[1]] %>%
  filter(Target == "Tyrannosauridae")
hadro <- all.results[[1]] %>%
  filter(Target == "Hadrosauridae")

c.uuid <- get_uuid(name = "Ceratopsidae", n = 4)[[4]]
t.uuid <- get_uuid(name = "Tyrannosauridae", n = 4)[[4]]
h.uuid <- get_uuid(name = "Edmontosaurus", n = 3)[[3]]

# Plot modelled results
a <- plot_occ(hadro)
b <- plot_naive(hadro, h.uuid)
c <- plot_occ(tyran)
d <- plot_naive(tyran, t.uuid)
e <- plot_occ(cera)
f <- plot_naive(cera, c.uuid)

# Arrange
(p <- ggarrange(f, b, d, e, a, c, 
          nrow = 2, ncol = 3,
          align='h', labels=c('A', 'B', 'C',
                              'D', 'E', 'F'),
          legend = "bottom",
          common.legend = T))

#if(bin.type =="formation"){
#  pdf(paste("Results/Outhwaite/", bin.type, "/", bin.res, 
#            "/Plot_", res, ".pdf", sep = ""), width = 11.458, height = 7.292)
#} else{
#  pdf(paste("Results/Outhwaite/", bin.type, "/Plot_", 
#            res, ".pdf", sep = ""), width = 11.458, height = 7.292)
#}
#dev.off()

# Choose type of model
if(type == 1){
  type <- "sp"
}
if(type == 2){
  type <- "hc.rw"
}



# Save figure
ggsave(paste("Figures/2.Sparta.", bin.type, ".", res, ".", type, ".png", sep = ""), plot = p, 
       device = "png")
ggsave(paste("Figures/2.Sparta.", bin.type, ".", res, ".", type, ".pdf", sep = ""), plot = p, 
       device = "pdf")

########################
##### SAVE RESULTS #####
########################

DIC.res <- data.frame(Taxon = target, 
           DIC = c(all.results[[2]]$Hadrosauridae$BUGSoutput$DIC,
                   all.results[[2]]$Ceratopsidae$BUGSoutput$DIC, 
                   all.results[[2]]$Tyrannosauridae$BUGSoutput$DIC), 
           Res = rep(res, 3), 
           Bin.type = rep(bin.type, 3), 
           Model.type = rep(type, 3))

write.csv(DIC.res, paste("Results/Outhwaite/DIC/", bin.type, ".", res, ".", type, 
                         ".csv", sep = ""))

write.csv(results, paste("Results/Outhwaite/", bin.type, "/naive.results.", 
                         res, ".", type, ".csv", sep = ""))
write.csv(all.results[[1]], paste("Results/Outhwaite/", bin.type,
                             "/results.", res, ".", type, ".csv", sep = ""))
saveRDS(all.results[[2]], file = paste("Results/Outhwaite/Posterior_checks/",
                                       bin.type, "/", res, ".", type, ".models.rds", 
                                       sep = ""))

####################################
##### ASSESSING ALL DIC VALUES #####
####################################

(all.DIC <- do.call(rbind, lapply(list.files(path = "Results/Outhwaite/DIC/", 
                                             full.names = T), read.csv)))
