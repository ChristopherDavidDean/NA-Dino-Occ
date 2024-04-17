################################################################################
############ OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS #############
################################################################################

# Christopher D. Dean, Alfio Alessandro Chiarenza, Jeffrey W. Doser, Alexander
# Farnsworth, Lewis A. Jones, Sinéad Lyster, Charlotte L. Outhwaite, Richard J. 
# Butler, Philip D. Mannion.
# 2024
# Script written by Nick Isaac, adapted by Christopher D. Dean.

################################################################################
#               FILE 3: POSTERIOR CHECKS FOR SPARTA.                         #
################################################################################

################################################################################
# 1. INITIAL SETUP
################################################################################

# Load packages
library(R2jags)
library(reshape2)
library(sparta)
library(LearnBayes)
library(lattice)
library(dplyr)

# Setup
iter_update <- 100
q <- c(0.025, 0.975)
res <- 0.5
dir.create("Results/Outhwaite/Posterior_checks")
bin.type <- "scotese"

# Read in data
para_out <- readRDS(file = paste("Results/Outhwaite/Posterior_checks/", bin.type, 
                                 "/", res, ".models.rds", sep = ""))

################################################################################
# 2. RUNNING TESTS
################################################################################

# Make list for p values
bayes.p.list <- list()
# Make list for analysis
check.data <- list()

# For each taxonomic group
for(group in names(para_out)){

  # locate and load the species' model
  out <- para_out[[group]]
    
  # check whether there is enough data to proceed
  if(out$species_observations >= 50){
  
  # recompile
  out$model$recompile()
    
  # Take samples
  samp <- rjags:::coda.samples(model=out$model, 
                               parallel=TRUE, n.cores=3,
                               variable.names="Py", 
                               n.iter=iter_update, thin=3)
  
  Y_rep <- lapply(samp, function(x) apply(x, 1:2, rbinom, n=1, size=1))
  
  Y_rep <- melt(Y_rep)
  Y_rep$VisitID <- as.character(gsub(Y_rep$Var2, pa="Py\\[", repl=""))
  Y_rep$VisitID <- as.numeric(gsub(Y_rep$VisitID, pa="\\]", repl=""))
  Y_rep$iter <- with(Y_rep, Var1 + max(Var1)*(L1 - 1))
  
  # load the raw data
  temp_data <- out$model$data()
  temp_data <- as.data.frame(with(temp_data, cbind(Site, Year, y)))
  temp_data$VisitID <- 1:nrow(temp_data)
  
  Y_rep <- merge(Y_rep[,-c(1,2,4)], temp_data)
  
  # next calculate the reporting rate in differing ways
  # first, for each site:year combination report the number of sites for which a positive record was made
  RR1 <- Y_rep %>% group_by(Site, Year, iter) %>%
    summarise(obs = max(y),
              sim = max(value)) %>%
    ungroup()
  # the data in RR1 refer to whether the species was recorded in each site in each year
  
  # next, the mean of these across years 
  RR2 <- RR1 %>% group_by(Year, iter) %>%
    summarise(obs_meanSitesPerYear = mean(obs),
              sim_meanSitesPerYear = mean(sim)
    ) %>%
    ungroup()
  # the data in RR2 refer to the proportion of Sites with positive records per Year
  RR3 <- RR2 %>% group_by(iter) %>%
    summarise(obs_varSitesPerYear = var(obs_meanSitesPerYear),
              sim_varSitesPerYear = var(sim_meanSitesPerYear)
    ) %>%
    ungroup()
  
  # now the Bayesian p-value
  Bayes_pVal <- c(meanSitesPerYear = with(RR2, 
                                          mean(obs_meanSitesPerYear < sim_meanSitesPerYear)),
                  varSitesPerYear = with(RR3, 
                                         mean(obs_varSitesPerYear < sim_varSitesPerYear)))
  # save this somewhere
  bayes.p.list[[group]] <- Bayes_pVal 
  
  # save the data 
  check.data[[group]] <- RR2
  }
}

# coerce to a single dataframe
sppMeanSites <- melt(check.data, id=1:4)
names(sppMeanSites)[names(sppMeanSites) == "L1"] <- "Taxonomic group"

# Now we're ready to calculate three other summary statistics:
# sppVarSites: the variance across years in the proportion of sites with a record
# taxonMeanSites: the annual mean (across species) in the proportion of sites with a record
# taxonVarSites: the variance across years in the above

# now the variance across years 
sppVarSites <- sppMeanSites %>% group_by(`Taxonomic group`, iter) %>%
  summarise(obs_varSitesPerYear = var(obs_meanSitesPerYear),
            sim_varSitesPerYear = var(sim_meanSitesPerYear)
  ) %>%
  ungroup()

# Summarise the data, including the Bayesian p value
sppMeanSites_Summ <- sppMeanSites %>% 
  group_by(`Taxonomic group`) %>%
  summarise(Bp = mean(obs_meanSitesPerYear > sim_meanSitesPerYear),
            obs = mean(obs_meanSitesPerYear),
            sim_mean = mean(sim_meanSitesPerYear),
            sim_lower = quantile(sim_meanSitesPerYear, q[1]),
            sim_upper = quantile(sim_meanSitesPerYear, q[2]),
  ) %>%
  ungroup()

sppVarSites_Summ <- sppVarSites %>% 
  group_by(`Taxonomic group`) %>%
  summarise(Bp = mean(obs_varSitesPerYear > sim_varSitesPerYear),
            obs = mean(obs_varSitesPerYear),
            sim_mean = mean(sim_varSitesPerYear),
            sim_lower = quantile(sim_varSitesPerYear, q[1]),
            sim_upper = quantile(sim_varSitesPerYear, q[2]),
  ) %>%
  ungroup()

################################################################################
# 3. PLOTTING
################################################################################

# Now plot the results
a <- ggplot(data = sppMeanSites_Summ) +
    geom_point(aes(x=obs, y=sim_mean, col=`Taxonomic group`)) +
    geom_errorbar(aes(x=obs, ymin=sim_lower, ymax=sim_upper, col=`Taxonomic group`)) +
    geom_abline(aes(intercept=0, slope=1)) + 
    scale_x_log10() + 
    scale_y_log10() + 
    ggtitle("Mean (across years) proportion of sites with a record") +
    xlab("Observed from data") + ylab("Model prediction") + 
    theme_bw()

b <- ggplot(data = sppVarSites_Summ) +
  geom_point(aes(x=obs, y=sim_mean, col=`Taxonomic group`)) +
  geom_errorbar(aes(x=obs, ymin=sim_lower, ymax=sim_upper, col=`Taxonomic group`)) +
  geom_abline(aes(intercept=0, slope=1)) + 
  scale_x_log10() + 
  scale_y_log10() + 
  ggtitle("Variance (across years) in proportion of sites per year") +
  xlab("Observed from data") + ylab("Model prediction") +
  theme_bw()

(p1 <- ggpubr::ggarrange(a, b,
          nrow = 1, ncol = 2,
          align='hv', labels=c('A', 'B'),
          legend = "bottom",
          common.legend = T))

#ggsave("Figures/X.Sparta.Pred.Checks.png", plot = p1, 
#       device = "png", type = "cairo")

# Check p values
print(bayes.p.list)
