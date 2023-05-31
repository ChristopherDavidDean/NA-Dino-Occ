################################################################################
############ OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS #############
################################################################################

# Christopher D. Dean, Lewis A. Jones, Alfio A. Chiarenza, Sinéad Lyster, 
# Alex Farnsworth, Philip D. Mannion, Richard J. Butler.
# 2023
# Script written by Christopher D. Dean and Lewis A. Jones

################################################################################
#              FILE 5: VISUALISING AND ANALYSING UNMARKED RESULTS              #
################################################################################

################################################################################
# 1. INITIAL SETUP
################################################################################

# Load packages
library(ggplot2)
library(wesanderson)
library(deeptime)
library(dplyr)

# Setup
res <- 0.5
bin.type <- "scotese"

# Load Bins
bins <- read.csv("Data/Occurrences/scotesebins.csv")
bins$Bin <- bins$code 
bins <- select(bins, Bin, max_ma, min_ma, mid_ma)

# Load Results
get.results(Ceratopsidae)
get.results(Hadrosauridae)
get.results(Tyrannosauridae)

################################################################################
# 2. PLOTTING
################################################################################

# Plot modelled results
a <- plot.occ.unmarked(Hadrosauridae)
b <- plot.naive.unmarked(Hadrosauridae)
c <- plot.occ.unmarked(Tyrannosauridae)
d <- plot.naive.unmarked(Tyrannosauridae)
e <- plot.occ.unmarked(Ceratopsidae)
f <- plot.naive.unmarked(Ceratopsidae)

# Arrange
p <- ggarrange(a, c, e, b, d, f,
               nrow = 2, ncol = 3,
               align='h', labels=c('A', 'B', 'C',
                                   'D', 'E', 'F'),
               legend = "bottom",
               common.legend = T)
p

# Add phylopics
cowplot::ggdraw() +  
  cowplot::draw_plot(p) +
  cowplot::draw_image("https://images.phylopic.org/images/aeeb30a8-afdc-4e7e-9bcc-574cb290a1f6/raster/1536x575.png?build=140", 
                      x = 0.435, y = -0.25, scale = 0.1) +
  cowplot::draw_image("https://images.phylopic.org/images/f3808e65-a95f-4df5-95a0-5f5b46a221f2/raster/1536x505.png?build=140", 
                      x = 0.09, y = -0.25, scale = 0.12) +
  cowplot::draw_image("https://images.phylopic.org/images/72be89b9-3f2b-4dc3-b485-e74a5f8b1fbc/raster/1536x512.png?build=140", 
                      x = -0.245, y = -0.25, scale = 0.1) 
