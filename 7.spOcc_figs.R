################################################################################
############ OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS #############
################################################################################

# Christopher D. Dean, Lewis A. Jones, Alfio A. Chiarenza, Sinéad Lyster, 
# Alex Farnsworth, Philip D. Mannion, Richard J. Butler.
# 2023
# Script written by Christopher D. Dean

################################################################################
#                    FILE 7: MULTI-SEASON OCCUPANCY FIGURES                    #
################################################################################

################################################################################
# 1. INITIAL SETUP
################################################################################

devtools::install_github("NightingaleHealth/ggforestplot")
library(ggforestplot)
library(ggforce)

# Set resolution
res <- 0.5

# Load in files
cera <- readRDS(paste("Results/spOccupancy/", res, 
                      "/Ceratopsidae.cov.table.rds", sep = ""))
hadr <- readRDS(paste("Results/spOccupancy/", res, 
                      "/Hadrosauridae.cov.table.rds", sep = ""))
tyra <- readRDS(paste("Results/spOccupancy/", res, 
                      "/Tyrannosauridae.cov.table.rds", sep = ""))

# Combined
all.sp <- rbind(cera, hadr)

# Forest plot
ggforestplot::forestplot(
  df = all.sp,
  name = Covariate,
  estimate = Mean,
  se = SD, 
  colour = Group,
  psignif = 0.05
) +
  ggforce::facet_col(
    facets = ~Submodel,
    scales = "free_y",
    space = "free"
  )

