################################################################################
############ OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS #############
################################################################################

# Christopher D. Dean, Alfio Alessandro Chiarenza, Jeffrey W. Doser, Alexander
# Farnsworth, Lewis A. Jones, Sinéad Lyster, Charlotte L. Outhwaite, Paul J. 
# Valdes, Richard J. Butler, Philip D. Mannion.
# 2025
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
res <- 0.5

# Set extent
e <- extent(-155, -72, 22.5, 73)

# Set bin type
bin.type <- "scotese"

# Set model type; 1 for sparta, 2 for half cauchy random walk.
type <- 1

# Set inclusion of small bodied organisms
nomam <- F

# Load main dataset
if(nomam == F){
  master.occs <- read.csv(paste("Prepped_data/Occurrence_Data/", bin.type, "/", 
                                bin.type, "_occurrence_dataset.csv", sep = ""))
}else{
  master.occs <- read.csv(paste("Prepped_data/Occurrence_Data/", bin.type, "/", 
                                bin.type, "_no_mammal_occurrence_dataset.csv", sep = ""))
}

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
target = c("Ankylosauridae", "Hadrosauridae", "Ceratopsidae", "Tyrannosauridae")

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

# Setup phylopic
a.uuid <- "5b062105-b6a2-4405-bd75-3d0399102b9a"
c.uuid <- "4b12a77a-0a16-4b93-baf5-7d4d01d0d9bd"
t.uuid <- "f3808e65-a95f-4df5-95a0-5f5b46a221f2"
h.uuid <- "72be89b9-3f2b-4dc3-b485-e74a5f8b1fbc"
silhouette_df <- data.frame(x = c(70.5, 70, 70.5, 70.5), y = c(0.92, 0.92, 0.92, 0.92), 
                            Target = c("Ankylosauridae", "Ceratopsidae", "Hadrosauridae", 
                                      "Tyrannosauridae"),
                            name = c("Ankylosauridae", "Ceratopsidae", "Hadrosauridae", 
                                     "Tyrannosauridae"))

# Plot naive occupancy
a <- ggplot(data = subset(all.results[[1]], Data == "Naive occupancy"), aes(x = new_bins, y = value)) +
  ylab("Proportion of total sites") + 
  xlab("Time (Ma)") +
  scale_x_reverse() +
  geom_phylopic(data = silhouette_df, aes(x = x, y = y), 
                uuid = c(a.uuid, c.uuid, h.uuid, t.uuid), 
                size = c(0.15, 0.15, 0.15, 0.16), 
                alpha = 1, 
                color = "grey") +
  #deeptime::coord_geo(dat = list("stages"), 
  #                    xlim = c((max(all.results[[1]]$new_bins)+1), (min(all.results[[1]]$new_bins-1))), 
  #                    ylim = c(0, 1)) +
  coord_cartesian(ylim = c(0.01,1)) +
  geom_line(aes(x = new_bins, y = value, color = Data)) +
  scale_color_manual(breaks = c("Naïve Occupancy", "PAO", "Occupancy Probability", "Detection Probability"),
                     values=c("#252424", "#DE2D26", "#DE2D26", "#3182BD")) +
  scale_fill_manual(breaks = c("Naïve Occupancy", "PAO", "Occupancy Probability", "Detection Probability"), 
                    values=c("#FFFFFF", "white","#DE2D26", "#3182BD")) +
  # scale_color_manual(values=c("#252424")) +
  theme_few() +
  geom_smooth(method=lm, fullrange = T) +
  facet_wrap(~Target, nrow = 1) +
  theme(strip.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        plot.margin = unit(c(0.2,0.2,-0.3,1), 'lines'), 
        panel.spacing.x = unit(c(0.7, 0.7, 0.7), "lines"),
        axis.title.y = element_text(size = 14))

# Plot modelled occupancy and detection
b <- ggplot(data = subset(all.results[[1]], Data == "Mean occupancy"), aes(x = new_bins, 
                                                                   y = value)) +
  geom_blank(aes(color = Data), data = all.results[[1]]) +
  geom_ribbon(data = all.results[[1]], aes(x = new_bins, ymin = lower95CI, 
                                   ymax = upper95CI, fill = Data), alpha = 0.2) +
  ylab("Probability") + 
  xlab("Time (Ma)") +
  scale_x_reverse() +
  deeptime::coord_geo(dat = list("stages"), 
                      xlim = c((max(all.results[[1]]$new_bins)+1), (min(all.results[[1]]$new_bins-1))), 
                      ylim = c(0, 1)) +
  geom_line(data = subset(all.results[[1]], Data != "Naive occupancy"), 
            aes(x = new_bins, y = value, color = Data)) +
  scale_color_manual(values=c("#BDD7E7", "#6BAED6", "#3182BD", 
                              "#DE2D26", "#252424")) +
  scale_fill_manual(values=c("#BDD7E7", "#6BAED6", "#3182BD", 
                             "#DE2D26", "#FFFFFF")) +
  new_scale_color() +
  geom_point(data = subset(all.results[[1]], Data == "Mean occupancy"),
             aes(x = new_bins, y = value, color = rhat_threshold), size = 2) +
  scale_color_manual(name = 'Rhat', values = c('Bad (>1.1)' = 'white',
                                               'Good (<1.1)' = '#DE2D26')) +
  geom_smooth(method=lm, fullrange = T) +
  theme_few() +
  facet_wrap(~Target, nrow = 1) +
  theme(legend.position = "bottom", 
        strip.text.x = element_blank(), 
        plot.margin = unit(c(0,0.2,0,1), 'lines'), 
        panel.spacing.x = unit(c(0.7, 0.7, 0.7), "lines"), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14))

# Combine
(p <- ggarrange(a, b, nrow = 2, heights = c(0.72, 1), align = 'v'))

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
if(nomam == F){
  ggsave(paste("Figures/NEW/2.Sparta.", bin.type, ".", res, ".", type, ".png", sep = ""), plot = p, 
         device = "png")
  ggsave(paste("Figures/NEW/2.Sparta.", bin.type, ".", res, ".", type, ".pdf", sep = ""), plot = p, 
         device = "pdf")
}else{
  ggsave(paste("Figures/NEW/2.Sparta.", bin.type, ".", res, ".", type, "_no_mammals.png", sep = ""), plot = p, 
         device = "png")
  ggsave(paste("Figures/NEW/2.Sparta.", bin.type, ".", res, ".", type, "_no_mammals.pdf", sep = ""), plot = p, 
         device = "pdf")
}


########################
##### SAVE RESULTS #####
########################

DIC.res <- data.frame(Taxon = target[1:4], 
           DIC = c(all.results[[2]]$Ankylosauridae$BUGSoutput$DIC, 
                   all.results[[2]]$Hadrosauridae$BUGSoutput$DIC,
                   all.results[[2]]$Ceratopsidae$BUGSoutput$DIC, 
                   all.results[[2]]$Tyrannosauridae$BUGSoutput$DIC), 
           Res = rep(res, 4), 
           Bin.type = rep(bin.type, 4), 
           Model.type = rep(type, 4))

if(nomam == F){
  write.csv(DIC.res, paste("Results/Outhwaite/NEW/DIC/", bin.type, ".", res, ".", type, 
                           ".csv", sep = ""))
  
  write.csv(results, paste("Results/Outhwaite/NEW/", bin.type, "/naive.results.", 
                           res, ".", type, ".csv", sep = ""))
  write.csv(all.results[[1]], paste("Results/Outhwaite/NEW/", bin.type,
                                    "/results.", res, ".", type, ".csv", sep = ""))
  saveRDS(all.results[[2]], file = paste("Results/Outhwaite/NEW/Posterior_checks/",
                                         bin.type, "/", res, ".", type, ".models.rds", 
                                         sep = ""))
}else{
  write.csv(DIC.res, paste("Results/Outhwaite/NEW/DIC/", bin.type, ".", res, ".", type, 
                           ".no_mammals.csv", sep = ""))
  
  write.csv(results, paste("Results/Outhwaite/NEW/", bin.type, "/naive.results.", 
                           res, ".", type, ".no_mammals.csv", sep = ""))
  write.csv(all.results[[1]], paste("Results/Outhwaite/NEW/", bin.type,
                                    "/results.", res, ".", type, ".no_mammals.csv", sep = ""))
  saveRDS(all.results[[2]], file = paste("Results/Outhwaite/NEW/Posterior_checks/",
                                         bin.type, "/", res, ".", type, ".no_mammals.models.rds", 
                                         sep = ""))
}


####################################
##### ASSESSING ALL DIC VALUES #####
####################################

(all.DIC <- do.call(rbind, lapply(list.files(path = "Results/Outhwaite/NEW/DIC/", 
                                             full.names = T), read.csv)))
