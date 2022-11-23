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
library(devtools)
library(sparta)
library(palaeoverse)
library(ggplot2)
library(snowfall)
#library(deeptime)
library(wesanderson)
library(divDyn)

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

#==== Making formation bins====

# Get information on ICS based intervals, then trim to fit relevant time frame. 
data(stages)
stages <- stages[80:81,] # Set stages to range from Campanian to Maastrichtian

# Run formation binning
bin_limits <- c(2.5, max(formations$max_age), 66) # Set user defined bin size - change the first number to vary resolution in graphs.
Scoring_Grid_1(formations = formations, res = 0.01) # run score grid to work out appropriate bin points
#Scoring_Grid_2(formations = formations, res = 0.01)
newBins(score_grid = score_grid, formations = formations, bin_limits = bin_limits, 
        allbins = allbins, stages = stages, smallamalg = TRUE) # Run new bins to generate bins
bins <- binlist %>% # Remove non-useful bins outside of range, add range for bins
  filter(bottom < 84) %>%
  mutate(range = top-bottom)
colnames(bins) <- c("bin", "min_ma", "max_ma", "mid_ma", "range") # Rename bin names to match across
bins[1,2] <- 66 # Cap youngest bin at 66 Ma. 
master.occs.binned <- bin_time(master.occs, bins, method = 'majority') # Bin occurrences
bin.type <- "formation"

# Adding midpoint
master.occs.binned$mid_ma <- (master.occs.binned$max_ma + master.occs.binned$min_ma)/2

#============================ RUNNING MODEL ====================================

#==== Occupancy model ====
# run the model with these data for one species
formattedOccData <- formatOccData(taxa = master.occs.binned$family,
                                site = master.occs.binned$siteID,
                                survey = master.occs.binned$mid_ma,
                                replicate = master.occs.binned$collection_no,
                                closure_period = master.occs.binned$bin_assignment)

# Initiate the cluster
sfInit(parallel = TRUE, cpus = 4)

# Export data to the cluster
sfExport('formattedOccData')

# Create a function that takes a species name and runs the model
occ_mod_function <- function(taxa_name){
  occ_out <- occDetFunc(taxa_name = taxa_name,
                        n_iterations = 40000,
                        nyr = 1,
                        burnin = 20000, 
                        occDetdata = formattedOccData$occDetdata,
                        spp_vis = formattedOccData$spp_vis,
                       # modeltype = c('indran', 'halfcauchy', 'catlistlength'),
                        write_results = TRUE)  
} 

# Run the model in parallel
system.time({
  para_out <- sfClusterApplyLB(c('Ceratopsidae', "Hadrosauridae", "Tyrannosauridae"), occ_mod_function)
})
# Name each element of this output by the species
for(i in  1:length(para_out)) names(para_out)[i] <- para_out[[i]]$SPP_NAM

# Stop the cluster
sfStop()

# Plot
plot(para_out$Ceratopsidae)
plot(para_out$Hadrosauridae)
plot(para_out$Tyrannosauridae)

#==== Get raw occupancy for comparison ====

# Make a list of target taxa
target = c("Hadrosauridae", "Ceratopsidae", "Tyrannosauridae")

# Make empty results table
results <- data.frame(matrix(ncol = 6, nrow = 0))

# For each target taxon, calculate raw occupancy through time
for(i in target){
  temp <- master.occs.binned
  
  # Set everything that's not the target to 0
  temp$Target[which(temp$Target != i | is.na(temp$Target))] <- 0
  
  # Set the target to 1 and make column numeric
  temp$Target[which(temp$Target == i)] <- 1
  temp$Target <- as.numeric(temp$Target)
  
  # Summarise data to establish occupied vs non-occupied grid cells per bin
  temp <- temp %>%
    select(bin_midpoint, siteID, Target) %>%
    distinct() %>%
    group_by(bin_midpoint, siteID) %>%
    summarise(sites = ceiling(mean(Target))) 
  
  # Make into dataframe and rearrange data
  temp <- as.data.frame(table(temp$bin_midpoint, temp$sites))
  temp <- dcast(temp, Var1 ~ Var2)
  temp$Total <- temp$`0`+temp$`1`
  temp$perc <- temp$`1`/temp$Total
  temp$Target <- i
  
  # Save to results
  results <- rbind(results, temp)
}

# Make results numeric for plotting
results$Var1 <- as.numeric(as.character(results$Var1))

# Plot results
a <- ggplot(results, aes(x=Var1, y=perc)) + 
  geom_line(aes(colour = Target)) +
  #scale_x_reverse() +
  theme_bw() +
  ylab("Frequency") +
  guides(fill=guide_legend(title="Preservation score")) +
  xlab("Time (Ma)") +
  scale_fill_manual(values=(wes_palette("Zissou1")))
a
#gggeo_scale(a)