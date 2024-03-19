################################################################################
############ OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS #############
################################################################################

# Christopher D. Dean, Lewis A. Jones, Alfio A. Chiarenza, Sinéad Lyster, 
# Alex Farnsworth, Philip D. Mannion, Richard J. Butler.
# 2023
# Script written by Christopher D. Dean

################################################################################
#                    FILE 10: MULTI-SEASON OCCUPANCY FIGURES                   #
################################################################################

################################################################################
# 1. INITIAL SETUP
################################################################################

library(ggforestplot)
library(ggforce)
library(ggh4x)

##### Load in Functions #####
source("0.Functions.R") # Import functions from other R file (must be in same working directory)

# Setup phylopic
c.uuid <- get_uuid(name = "Ceratopsidae", n = 4)[[4]]
t.uuid <- get_uuid(name = "Tyrannosauridae", n = 4)[[4]]
h.uuid <- get_uuid(name = "Edmontosaurus", n = 3)[[3]]

# Function for cleaning figures
clean.for.fig <- function(cov.table){
  cov.table[cov.table == "scale(ann)"] <- "Ann"
  cov.table[cov.table == "scale(dry)"] <- "Dry"
  cov.table[cov.table == "scale(col)"] <- "Cold"
  cov.table[cov.table == "scale(wet)"] <- "Wet"
  cov.table[cov.table == "scale(hot)"] <- "Hot"
  cov.table[cov.table == "scale(Distance)"] <- "Distance"
  cov.table[cov.table == "scale(sedflux)"] <- "Sediment Flux" 
  cov.table[cov.table == "scale(rain)"] <- "Rainfall" 
  cov.table[cov.table == "scale(MGVF)"] <- "MGVF"
  cov.table[cov.table == "scale(occur)"] <- "Occurrences"
  cov.table[cov.table == "scale(outcrop)"] <- "Outcrop area"
  cov.table[cov.table == "factor(Year)2"] <- "Bin 2 (75 Ma)"
  cov.table[cov.table == "factor(Year)3"] <- "Bin 3 (69 Ma)"
  cov.table[cov.table == "factor(Year)4"] <- "Bin 4 (66.7 Ma)"
  cov.table[cov.table == "factor(land)2"] <- "Land cover (2)"
  cov.table[cov.table == "factor(land)3"] <- "Land cover (3)"
  cov.table[cov.table == "factor(land)4"] <- "Land cover (4)"
  cov.table[cov.table == "(Intercept)"] <- "Intercept"
  cov.table[cov.table == "scale(coll)"] <- "Collections"
  cov.table[cov.table == "Random Effect Variance; Site"] <- "REV: Site"
  cov.table[cov.table == "teyen"] <- "Bin 4 (66.7 Ma)"
  cov.table[cov.table == "teyeo"] <- "Bin 3 (69 Ma)"
  cov.table[cov.table == "teyep"] <- "Bin 2 (75 Ma)"
  cov.table[cov.table == "teyeq"] <- "Bin 1 (80.8 Ma)"
  return(cov.table)
}

################################################################################
# FIGURE 2. OCCUPANCY AND DETECTION THROUGH TIME, SPARTA
################################################################################

#################
##### PLOTS #####
#################

# Set resolution and time bin scheme
res <- 0.5
bin.type <- "scotese"

# Set target
target <- c("Ceratopsidae", 
            "Hadrosauridae",
            "Tyrannosauridae")

# Load data
all.results <- read.csv(paste("Results/Outhwaite/", bin.type,
                              "/results.", res, ".csv", sep = ""))

# Plotting occupancy (naive and modelled)
cera <- all.results %>%
  filter(Target == "Ceratopsidae")
tyran <- all.results %>%
  filter(Target == "Tyrannosauridae")
hadro <- all.results %>%
  filter(Target == "Hadrosauridae")

c.uuid <- get_uuid(name = "Ceratopsidae", n = 4)[[4]]
t.uuid <- get_uuid(name = "Tyrannosauridae", n = 4)[[4]]
h.uuid <- get_uuid(name = "Edmontosaurus", n = 3)[[3]]

# Plot modelled results
a <- plot.occ(hadro)
b <- plot.naive(hadro, h.uuid)
c <- plot.occ(tyran)
d <- plot.naive(tyran, t.uuid)
e <- plot.occ(cera)
f <- plot.naive(cera, c.uuid)

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
type <- "sp"

# Save figure
ggsave(paste("Figures/2.Sparta.", bin.type, ".", res, ".", type, ".png", sep = ""), plot = p, 
       device = "png")

ggsave(paste("Figures/2.Sparta.", bin.type, ".", res, ".", type, ".pdf", sep = ""), plot = p, 
       device = "pdf")

################################################################################
# FIGURE 3: DETECTION THROUGH TIME, SPOCCUPANCY
################################################################################

# Set extent
e <- extent(-155, -72, 22.5, 73)

# Set target
target <- c("Ceratopsidae", 
            "Hadrosauridae",
            "Tyrannosauridae")

all.p.vals <- dplyr::bind_rows(lapply(target, function(target){
  res <- c(0.5, 1)
  p.vals <- lapply(res, function(res){
    best.model <- readRDS(paste("Results/spOccupancy/", res, 
                                "/", target, ".best.model.rds", sep = ""))
    siteCoords <- siteCoordsFun(res = res, e = e, as.numeric(rownames(best.model$y))) 
    fit.mod <- fitted(best.model)
    p.samples <- fit.mod$p.samples
    p.samples <- colMeans(p.samples)
    p.vals <- as.data.frame(apply(p.samples,c(1,2),mean, na.rm = T))
    mean.p <- as.data.frame(colMeans(p.vals, na.rm = T))
    CIs <- t(apply(p.vals, 2, quantile, c(0.025, 0.975), na.rm = T))
    colnames(mean.p) <- "p"
    mean.p$Resolution <- res
    mean.p$Group <- target
    mean.p$Bin <- c("teyeq", "teyep", "teyeo", "teyen")
    mean.p <- cbind(mean.p, CIs)
    return(mean.p)
  })
  return(p.vals)
}))

all.p.vals[all.p.vals == "teyeq"] <- 80.8
all.p.vals[all.p.vals == "teyep"] <- 75
all.p.vals[all.p.vals == "teyeo"] <- 69
all.p.vals[all.p.vals == "teyen"] <- 66.7

c.uuid <- get_uuid(name = "Ceratopsidae", n = 4)[[4]]
t.uuid <- get_uuid(name = "Tyrannosauridae", n = 4)[[4]]
h.uuid <- get_uuid(name = "Edmontosaurus", n = 3)[[3]]

# Setup phylopic
silhouette_df <- data.frame(x = c(71, 71, 71), y = c(0.92, 0.92, 0.92), 
                            Group = c("Ceratopsidae", "Hadrosauridae", 
                                      "Tyrannosauridae"),
                            Resolution = c("0.5", "0.5", "0.5"),
                            name = c("Ceratopsidae", "Hadrosauridae", 
                                     "Tyrannosauridae"))

all.p.vals$Bin <- as.numeric(all.p.vals$Bin)
all.p.vals$Group <- as.factor(all.p.vals$Group)
all.p.vals$Resolution <- as.factor(all.p.vals$Resolution)

(p2 <- ggplot(data = all.p.vals, aes(x = Bin, y = p)) +
  ylab("Detection probability") + 
  xlab("Time (Ma)") +
  scale_x_reverse() +
  deeptime::coord_geo(dat = list("stages"), 
                      expand = T, 
                      ylim = c(0, 1)) +
  geom_line(aes(color = Group)) +
  theme_few() +
  geom_ribbon(aes(x = Bin, ymin = `2.5%`, 
                  ymax = `97.5%`, fill = Group), alpha = 0.2) +
  scale_color_manual(values=c("#3182BD", "#3182BD", "#3182BD")) +
  scale_fill_manual(values=c("#3182BD", "#3182BD", "#3182BD"))+
  #viridis::scale_color_viridis(discrete=TRUE) + 
  #viridis::scale_fill_viridis(discrete = TRUE) +
  theme(legend.position="none", 
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  geom_phylopic(data = silhouette_df, aes(x = x, y = y), 
                uuid = c(c.uuid, h.uuid, t.uuid), 
                size = c(0.15, 0.15, 0.15), 
                alpha = 1, 
                color = "grey") +
  facet_grid(c("Resolution", "Group"), 
             scales = "free"))

# Save figure
ggsave(paste("Figures/3.spOcc.detection.png", sep = ""), plot = p2, 
       device = "png")

ggsave(paste("Figures/3.spOcc.detection.pdf", sep = ""), plot = p2, 
       device = "pdf")

################################################################################
# 4. DETECTION PROBABILITY FIGURES
################################################################################

# Set resolution
res <- 0.5
res <- 1
# Set extent
e <- extent(-155, -72, 22.5, 73) # 0.5 degree
e <- extent(-155, -72, 23, 73) # 1 degree

# Get chosen model
target <- c("Ceratopsidae", "Hadrosauridae", "Tyrannosauridae")

test <-  lapply(target, function(target){
  best.model <- readRDS(paste("Results/spOccupancy/", res, 
                              "/", target, ".best.model.rds", sep = ""))
  siteCoords <- siteCoordsFun(res = res, e = e, as.numeric(rownames(best.model$y))) 
  fit.mod <- fitted(best.model)
  p.samples <- fit.mod$p.samples
  p.samples <- colMeans(p.samples)
  p.vals <- as.data.frame(apply(p.samples,c(1,2),mean, na.rm = T))
  
  teyen <- gen_raster(siteCoords$siteID, p.vals$V4, res = res, ext = e)
  teyeo <- gen_raster(siteCoords$siteID, p.vals$V3, res = res, ext = e)
  teyep <- gen_raster(siteCoords$siteID, p.vals$V2, res = res, ext = e)
  teyeq <- gen_raster(siteCoords$siteID, p.vals$V1, res = res, ext = e)
  
  raster.names <- c("Bin 1 (80.8 Ma)", "Bin 2 (75 Ma)", "Bin 3 (69 Ma)", "Bin 4 (66.7 Ma)")
  
  bin.list <- stack(teyeq, teyep, teyeo, teyen)
  names(bin.list) <- raster.names
  return(bin.list)
})

p.vals <- lapply(target, function(target){
  best.model <- readRDS(paste("Results/spOccupancy/", res, 
                              "/", target, ".best.model.rds", sep = ""))
  siteCoords <- siteCoordsFun(res = res, e = e, as.numeric(rownames(best.model$y))) 
  fit.mod <- fitted(best.model)
  p.samples <- fit.mod$p.samples
  p.samples <- colMeans(p.samples)
  p.vals <- as.data.frame(apply(p.samples,c(1,2),mean, na.rm = T))
  return(p.vals)
})

p.vals <- bind_rows(p.vals)

test <- stack(test[[1]], test[[2]], test[[3]])
raster.names <- c("Bin 1 (80.8 Ma)", "Bin 2 (75 Ma)", "Bin 3 (69 Ma)", "Bin 4 (66.7 Ma)", 
                  "", "", "", "", "", "", "", "")

# find map to use as backdrop
countries <- maps::map("world", plot=FALSE, fill = TRUE) 
# Turn map into spatialpolygons
countries <- maptools::map2SpatialPolygons(countries, 
                                           IDs = countries$names, 
                                           proj4string = CRS("+proj=longlat")) 
test <- crop(test, e)

teyen <- sf::st_read("Data/Covariate_Data/Outcrop/New_shapefiles",
            layer = "teyen")
teyen <- sf::as_Spatial(teyen)
teyeo <- sf::st_read("Data/Covariate_Data/Outcrop/New_shapefiles",
                     layer = "teyeo")
teyeo <- sf::as_Spatial(teyeo)
teyep <- sf::st_read("Data/Covariate_Data/Outcrop/New_shapefiles",
                     layer = "teyep")
teyep <- sf::as_Spatial(teyep)
teyeq <- sf::st_read("Data/Covariate_Data/Outcrop/New_shapefiles",
                     layer = "teyeq")
teyeq <- sf::as_Spatial(teyeq)

# Customize the colorkey
my.theme <- BuRdTheme()
my.at <- seq(min(p.vals, na.rm = T), max(p.vals, na.rm = T), length.out=length(my.theme$regions$col)-1)
my.ckey <- list(at=my.at, col=my.theme$regions$col)
p3 <- rasterVis::levelplot(test, layout = c(4, 3), margin=T, par.settings=my.theme,
                           colorkey=my.ckey, names.attr = raster.names) + 
  # Plots state lines
  latticeExtra::layer(sp.polygons(teyen, col = 0, fill = "dark grey"), packets = c(4,8,12), under = T) +
  latticeExtra::layer(sp.polygons(teyeo, col = 0, fill = "dark grey"), packets = c(3,7,11), under = T) +
  latticeExtra::layer(sp.polygons(teyep, col = 0, fill = "dark grey"), packets = c(2,6,10), under = T) +
  latticeExtra::layer(sp.polygons(teyeq, col = 0, fill = "dark grey"), packets = c(1,5,9), under = T) +
  latticeExtra::layer(sp.polygons(states, col = "white", lwd = 0.5, fill = NA), under = T)  + 
  # Plots background colour
  latticeExtra::layer(sp.polygons(countries, col = 0, fill = "#dcdcdc"), under = T)


# Make ggplot
p3 <- ggplotify::as.ggplot(p3)

# Set phylopic requirements
silhouette_df <- data.frame(x = c(0.144,0.148,0.151), y = c(0.715,0.42,0.125), 
                            Group = c("Ceratopsidae", "Hadrosauridae", 
                                      "Tyrannosauridae"),
                            Submodel = c("Detection", "Detection", "Detection"),
                            name = c("Ceratopsidae", "Hadrosauridae", 
                                     "Tyrannosauridae"))
# Combine with phylopic
p3.5 <- p3 + geom_phylopic(data = silhouette_df, aes(x = x, y = y), 
                   uuid = c(c.uuid, h.uuid, t.uuid), 
                   size = c(0.034, 0.034, 0.038), 
                   alpha = 1, 
                   color = "dark grey")

# Save figure
ggsave(paste("Figures/4.5.spOcc.detection.space.", res, ".png", sep = ""), plot = p3.5, 
       device = "png")

ggsave(paste("Figures/4.5.spOcc.detection.space.", res, ".pdf", sep = ""), plot = p3.5, 
       device = "pdf")

################################################################################
# 5. MULTI-SEASON FIGURES
################################################################################

###################
##### SPATIAL #####
###################

# Get table
cov.table <- dplyr::bind_rows(lapply(target, function(target){
  res <- c(0.5, 1)
    all.covs <- lapply(res, function(res){
      a <- readRDS(paste("Results/spOccupancy/", res, "/", target, 
                         ".best.model.rds", sep = ""))
      b <- make.table(out.sp = a, res = res, target = target)
    })
  return(all.covs)
}))

cov.table <- clean.for.fig(cov.table)
cov.table$Resolution <- as.factor(cov.table$Res)
cov.table <- cov.table[order(cov.table$Group),]
cov.table <- cov.table[order(cov.table$Resolution),]

reference <- c("Intercept", "Ann", "Cold", "Hot", "Dry", "Wet", "Collections",
               "Distance", "Occurrences", "Land cover (2)", 
               "Land cover (3)", "Land cover (4)", "MGVF", "Rainfall", "Bin 2 (75 Ma)", 
               "Bin 3 (69 Ma)", "Bin 4 (66.7 Ma)", "Outcrop area", "Sediment Flux", 
               "REV: Site")

cov.table <- cov.table %>%
  dplyr::rowwise() %>%
  dplyr::mutate(Sig = as.numeric(!between(0,`2.5%`,`97.5%`))) 
cov.table <- as.data.frame(cov.table)

cov.table <- cov.table[order(factor(cov.table$Covariate, levels = reference)),]


# Set up phylopic
silhouette_df <- data.frame(x = c(4, 4, 4), y = c(10.7, 9, 6.3), 
                            Group = c("Ceratopsidae", "Hadrosauridae", 
                                      "Tyrannosauridae"),
                            Submodel = c("Detection", "Detection", "Detection"),
                            name = c("Ceratopsidae", "Hadrosauridae", 
                                     "Tyrannosauridae"))
# Make figure
(p4 <- forestplot2(
  df = cov.table,
  name = Covariate,
  estimate = Mean,
  pvalue = Sig, 
  colour = Resolution,
  CI.max = `97.5%`, 
  CI.min = `2.5%`,
  xlab = "Beta estimates",
  title = "",
  psignif = 0.05
) +
  geom_phylopic(data = silhouette_df, aes(x = x, y = y), 
                uuid = c(c.uuid, h.uuid, t.uuid), 
                size = c(2.6, 2, 1.6), 
                alpha = 1, 
                color = "dark grey") +
  ggplot2::facet_wrap(
    facets = ~Group + Submodel,
    nrow = 3, 
    ncol = 2,
    scales = "free_y"
  ) + facetted_pos_scales(
    y = list(Submodel == "Occupancy" ~ scale_y_discrete(position = "right"))
  )  +
theme(strip.text.x = element_blank(), 
      legend.position = "bottom", 
      panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5), 
      axis.text.y.right = element_text(hjust = 0)) + 
  scale_colour_manual(values = wesanderson::wes_palette("Darjeeling1", type = 'discrete')))

# Save plot
ggsave(paste("Figures/5.Cov.plot.png", sep = ""), plot = p4, 
       device = "png")
ggsave(paste("Figures/5.Cov.plot.pdf", sep = ""), plot = p4, 
       device = "pdf")

#######################
##### NON-SPATIAL #####
#######################

# Get table
cov.table <- dplyr::bind_rows(lapply(target, function(target){
  res <- c(0.5, 1)
  all.covs <- lapply(res, function(res){
    a <- readRDS(paste("Results/spOccupancy/Non.spatial/", res, "/", target, 
                       ".best.model.rds", sep = ""))
    b <- make.table(out.sp = a, res = res, target = target)
  })
  return(all.covs)
}))

cov.table <- clean.for.fig(cov.table)
cov.table$Resolution <- as.factor(cov.table$Res)
cov.table <- cov.table[order(cov.table$Group),]
cov.table <- cov.table[order(cov.table$Resolution),]
cov.table <- cov.table %>%
  dplyr::filter(Covariate != "Intercept") %>%
  dplyr::arrange(Covariate) %>%
  rbind(dplyr::filter(cov.table, Covariate == "Intercept"), .)

# Set up phylopic
silhouette_df <- data.frame(x = c(4.2, 4.5, 4.5), y = c(10.8, 10, 4), 
                            Group = c("Ceratopsidae", "Hadrosauridae", 
                                      "Tyrannosauridae"),
                            Submodel = c("Detection", "Detection", "Detection"),
                            name = c("Ceratopsidae", "Hadrosauridae", 
                                     "Tyrannosauridae"))
# Make figure
(p4.5 <- forestplot2(
  df = cov.table,
  name = Covariate,
  estimate = Mean,
  colour = Resolution,
  CI.max = `97.5%`, 
  CI.min = `2.5%`,
  xlab = "Beta estimates",
  title = "",
  psignif = 0.05
) +
     geom_phylopic(data = silhouette_df, aes(x = x, y = y), 
                  uuid = c(c.uuid, h.uuid, t.uuid), 
                   size = c(2.2, 2, 0.8), 
                   alpha = 1, 
                   color = "dark grey") +
    ggplot2::facet_wrap(
      facets = ~Group + Submodel,
      nrow = 3, 
      ncol = 2,
      scales = "free_y"
    ) + facetted_pos_scales(
      y = list(Submodel == "Occupancy" ~ scale_y_discrete(position = "right"))
    )  +
    theme(strip.text.x = element_blank(), 
          legend.position = "bottom", 
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5), 
          axis.text.y.right = element_text(hjust = 0)) + 
    scale_colour_manual(values = wesanderson::wes_palette("Darjeeling1", type = 'discrete')))

# Save plot
ggsave(paste("Figures/5.5.Cov.plot.non.spatial.png", sep = ""), plot = p4.5, 
       device = "png")
ggsave(paste("Figures/5.5.Cov.plot.non.spatial.pdf", sep = ""), plot = p4.5, 
       device = "pdf")

#######################
##### COMPARISONS #####
#######################

# Compare WAIC for resolution
spOccupancy::waicOcc(readRDS("Results/spOccupancy/0.5/Tyrannosauridae.best.model.rds"))
spOccupancy::waicOcc(readRDS("Results/spOccupancy/1/Tyrannosauridae.best.model.rds"))

spOccupancy::waicOcc(readRDS("Results/spOccupancy/0.5/Hadrosauridae.best.model.rds"))
spOccupancy::waicOcc(readRDS("Results/spOccupancy/1/Hadrosauridae.best.model.rds"))

spOccupancy::waicOcc(readRDS("Results/spOccupancy/0.5/Ceratopsidae.best.model.rds"))
spOccupancy::waicOcc(readRDS("Results/spOccupancy/1/Ceratopsidae.best.model.rds"))

# Compare WAIC for resolution, non-spatial
spOccupancy::waicOcc(readRDS("Results/spOccupancy/Non.spatial/0.5/Tyrannosauridae.best.model.rds"))
spOccupancy::waicOcc(readRDS("Results/spOccupancy/Non.spatial/1/Tyrannosauridae.best.model.rds"))

spOccupancy::waicOcc(readRDS("Results/spOccupancy/Non.spatial/0.5/Hadrosauridae.best.model.rds"))
spOccupancy::waicOcc(readRDS("Results/spOccupancy/Non.spatial/1/Hadrosauridae.best.model.rds"))

spOccupancy::waicOcc(readRDS("Results/spOccupancy/Non.spatial/0.5/Ceratopsidae.best.model.rds"))
spOccupancy::waicOcc(readRDS("Results/spOccupancy/Non.spatial/1/Ceratopsidae.best.model.rds"))

# Compare WAIC, spatial vs. non-spatial, 0.5
spOccupancy::waicOcc(readRDS("Results/spOccupancy/Non.spatial/0.5/Tyrannosauridae.best.model.rds"))
spOccupancy::waicOcc(readRDS("Results/spOccupancy/0.5/Tyrannosauridae.best.model.rds"))

spOccupancy::waicOcc(readRDS("Results/spOccupancy/Non.spatial/0.5/Hadrosauridae.best.model.rds"))
spOccupancy::waicOcc(readRDS("Results/spOccupancy/0.5/Hadrosauridae.best.model.rds"))

spOccupancy::waicOcc(readRDS("Results/spOccupancy/Non.spatial/0.5/Ceratopsidae.best.model.rds"))
spOccupancy::waicOcc(readRDS("Results/spOccupancy/0.5/Ceratopsidae.best.model.rds"))

# Compare WAIC, spatial vs. non-spatial, 1
spOccupancy::waicOcc(readRDS("Results/spOccupancy/Non.spatial/1/Tyrannosauridae.best.model.rds"))
spOccupancy::waicOcc(readRDS("Results/spOccupancy/1/Tyrannosauridae.best.model.rds"))

spOccupancy::waicOcc(readRDS("Results/spOccupancy/Non.spatial/1/Hadrosauridae.best.model.rds"))
spOccupancy::waicOcc(readRDS("Results/spOccupancy/1/Hadrosauridae.best.model.rds"))

spOccupancy::waicOcc(readRDS("Results/spOccupancy/Non.spatial/1/Ceratopsidae.best.model.rds"))
spOccupancy::waicOcc(readRDS("Results/spOccupancy/1/Ceratopsidae.best.model.rds"))

################################################################################
# 6. SINGLE-SEASON FIGURES
################################################################################

###################
##### SPATIAL #####
###################

# Get table
cov.table <- dplyr::bind_rows(lapply(target, function(target){
  res <- c(0.5, 1)
  all.covs <- lapply(res, function(res){
    a <- readRDS(paste("Results/spOccupancy/", res, "/", target, 
                       ".best.model.rds", sep = ""))
    b <- make.table(out.sp = a, res = res, target = target)
  })
  return(all.covs)
}))

forest.covs <- function(target){
  test <-data.frame(res = c(rep(0.5, 4), rep(1, 4)), 
                    bin = rep(c("teyen", "teyeo", "teyep", "teyeq"), 2))
  all.covs <- dplyr::bind_rows(apply(test, 1, function(test){
    a <- as.numeric(test[1])
    b <- test[2]
    c <- readRDS(paste("Results/spOccupancy/single_season/", a, "/", b, ".", target, 
                         ".best.model.rds", sep = ""))
    return(make.table(out.sp = c, res = a, target = target, ss = T, bin = b))
  }))
  return(all.covs)
}

target <- c("Hadrosauridae", 
            "Ceratopsidae", 
            "Tyrannosauridae")

cov.table <- dplyr::bind_rows(lapply(target, forest.covs))

cov.table <- clean.for.fig(cov.table)

cov.table <- cov.table %>%
  dplyr::rowwise() %>%
  dplyr::mutate(Sig = as.numeric(!between(0,`2.5%`,`97.5%`))) 
cov.table <- as.data.frame(cov.table)

cov.table$Resolution <- as.factor(cov.table$Res)
cov.table <- cov.table[order(cov.table$Group),]
cov.table <- cov.table[order(cov.table$Resolution),]
cov.table <- cov.table %>%
  dplyr::filter(Covariate != "Intercept") %>%
  dplyr::arrange(Covariate) %>%
  rbind(dplyr::filter(cov.table, Covariate == "Intercept"), .)

# Setup phylopic
silhouette_df <- data.frame(x = c(5.3, 5.5, 5.5), y = c(9.8, 10.6, 10.4), 
                            Group = c("Ceratopsidae", "Hadrosauridae", 
                                      "Tyrannosauridae"),
                            Submodel = c("Detection", "Detection", "Detection"),
                            name = c("Ceratopsidae", "Hadrosauridae", 
                                     "Tyrannosauridae"))

(p5 <- forestplot2(
  df = cov.table,
  name = Covariate,
  estimate = Mean,
  shape = Resolution,
  CI.max = `97.5%`, 
  CI.min = `2.5%`,
  xlab = "Beta estimates",
  pvalue = Sig, 
  title = "",
  colour = Bin,
  psignif = 0.05
) +
  geom_phylopic(data = silhouette_df, aes(x = x, y = y), 
                uuid = c(c.uuid, h.uuid, t.uuid), 
                size = c(2.4, 2.4, 2.8), 
                alpha = 1, 
                color = "dark grey") +
    ggplot2::facet_wrap(
      facets = ~Group + Submodel,
      nrow = 3, 
      ncol = 2,
      scales = "free_y"
    ) + facetted_pos_scales(
      y = list(Submodel == "Occupancy" ~ scale_y_discrete(position = "right"))
    ) + theme(strip.text.x = element_blank(), 
              legend.position = "bottom", 
              panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5), 
              axis.text.y.right = element_text(hjust = 0)) + 
  scale_colour_manual(values = wesanderson::wes_palette("Zissou1", type = 'discrete')[c(1:2, 4:5)]))

# Save plot
ggsave(paste("Figures/6.Cov.plot.single.png", sep = ""), plot = p5, 
       device = "png")
ggsave(paste("Figures/6.Cov.plot.single.pdf", sep = ""), plot = p5, 
       device = "pdf")

###########################
##### FORMATION CELLS #####
###########################

# Get table
forest.covs <- function(target){
  test <-data.frame(res = c(rep(0.5, 4), rep(1, 4)), 
                    bin = rep(c("teyen", "teyeo", "teyep", "teyeq"), 2))
  all.covs <- dplyr::bind_rows(apply(test, 1, function(test){
    a <- as.numeric(test[1])
    b <- test[2]
    c <- readRDS(paste("Results/spOccupancy/single_season/", a, "/", b, ".", target, 
                       ".formcells.best.model.rds", sep = ""))
    return(make.table(out.sp = c, res = a, target = target, ss = T, bin = b))
  }))
  return(all.covs)
}

target <- c("Hadrosauridae", 
            "Ceratopsidae", 
            "Tyrannosauridae")

cov.table <- dplyr::bind_rows(lapply(target, forest.covs))

cov.table <- clean.for.fig(cov.table)
cov.table$Resolution <- as.factor(cov.table$Res)
cov.table <- cov.table[order(cov.table$Group),]
cov.table <- cov.table[order(cov.table$Resolution),]
cov.table <- cov.table %>%
  dplyr::filter(Covariate != "Intercept") %>%
  dplyr::arrange(Covariate) %>%
  rbind(dplyr::filter(cov.table, Covariate == "Intercept"), .)

# Setup phylopic
silhouette_df <- data.frame(x = c(7.3, 7.5, 7.5), y = c(10, 10, 10), 
                            Group = c("Ceratopsidae", "Hadrosauridae", 
                                      "Tyrannosauridae"),
                            Submodel = c("Detection", "Detection", "Detection"),
                            name = c("Ceratopsidae", "Hadrosauridae", 
                                     "Tyrannosauridae"))

(p5 <- forestplot2(
  df = cov.table,
  name = Covariate,
  estimate = Mean,
  shape = Resolution,
  CI.max = `97.5%`, 
  CI.min = `2.5%`,
  xlab = "Beta estimates",
  title = "",
  colour = Bin,
  psignif = 0.05
) +
    geom_phylopic(data = silhouette_df, aes(x = x, y = y), 
                  uuid = c(c.uuid, h.uuid, t.uuid), 
                  size = c(2.4, 2.4, 2.5), 
                  alpha = 1, 
                  color = "dark grey") +
    ggplot2::facet_wrap(
      facets = ~Group + Submodel,
      nrow = 3, 
      ncol = 2,
      scales = "free_y"
    ) + facetted_pos_scales(
      y = list(Submodel == "Occupancy" ~ scale_y_discrete(position = "right"))
    ) + theme(strip.text.x = element_blank(), 
              legend.position = "bottom", 
              panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5), 
              axis.text.y.right = element_text(hjust = 0)) + 
    scale_colour_manual(values = wesanderson::wes_palette("Zissou1", type = 'discrete')[c(1:2, 4:5)]))


################################################################################
# A1. VARIATIONS IN SPARTA COMBINED FIGURE
################################################################################

type <- "sp"
type <- "hc.rw"

# Load data
scot.5 <- read.csv(paste("Results/Outhwaite/scotese/results.0.5.", type, ".csv", sep = ""))
scot.1 <- read.csv(paste("Results/Outhwaite/scotese/results.1.", type, ".csv", sep = ""))
form.5 <- read.csv(paste("Results/Outhwaite/formation/results.0.5.", type, ".csv", sep = ""))
form.1 <- read.csv(paste("Results/Outhwaite/formation/results.1.", type, ".csv", sep = ""))

plot.combine <- function(all.results){
  # Plotting occupancy (naive and modelled)
  cera <- all.results %>%
    filter(Target == "Ceratopsidae")
  tyran <- all.results %>%
    filter(Target == "Tyrannosauridae")
  hadro <- all.results %>%
    filter(Target == "Hadrosauridae")
  # Plot modelled results
  a <- plot.occ(cera)
  b <- plot.occ(hadro)
  c <- plot.occ(tyran)
  test <- list(a, b, c)
  return(test)
}
a <- plot.combine(scot.5)
b <- plot.combine(scot.1)
c <- plot.combine(form.5)
d <- plot.combine(form.1)

# Arrange
(p <- ggarrange(a[[1]], a[[2]], a[[3]], 
                b[[1]], b[[2]], b[[3]], 
                c[[1]], c[[2]], c[[3]], 
                d[[1]], d[[2]], d[[3]], 
                
                nrow = 4, ncol = 3,
                align='h', labels=c('A', 'B', 'C',
                                    'D', 'E', 'F', 
                                    'G', 'H', 'I', 
                                    'J', 'K', 'L'),
                legend = "bottom",
                common.legend = T))
# Save figure
ggsave(paste("Figures/A1.Sparta.compare.", type, ".png", sep = ""), plot = p, 
       device = "png", width = 210, height = 297, unit = 'mm')

ggsave(paste("Figures/A1.Sparta.compare.", type, ".pdf", sep = ""), plot = p, 
       device = "pdf", width = 210, height = 297, unit = 'mm')

################################################################################
# A2. LIST LENGTH DIFFERENCES
################################################################################

test <- read.csv(paste("Results/Outhwaite/scotese/All/results", res, ".csv", sep = ""))

test <- test %>%
  dplyr::filter(Data != "Mean occupancy") %>%
  dplyr::filter(Data != "Naive occupancy")

test2 <- test %>%
  dplyr::select(-c(lower95CI, upper95CI, X, year, group, rhat_threshold)) %>%
  tidyr::pivot_wider(names_from = Data, values_from = value) %>%
  mutate(diff12 = `Mean detection (LL2)` - `Mean detection (LL1)`, 
         diff24 = `Mean detection (LL4)` - `Mean detection (LL2)`, 
         diff14 = `Mean detection (LL4)` - `Mean detection (LL1)`) %>%
  tidyr::pivot_longer(!c(Target, new_bins), names_to = "Test", 
                      values_to = "Value") %>%
  dplyr::filter(Test == c("diff12", "diff24", "diff14"))

test2[test2 == "diff12"] <- "LL2 - LL1"
test2[test2 == "diff24"] <- "LL4 - LL2"
test2[test2 == "diff14"] <- "LL4 - LL1"

test2$Difference <- factor(test2$Test, levels=c("LL2 - LL1", 
                                          "LL4 - LL2", 
                                          "LL4 - LL1"))

silhouette_df <- data.frame(x = c(69.2, 69, 69), y = c(0.86, 0.86, 0.86), 
                            Target = c("Ceratopsidae", "Hadrosauridae", 
                                      "Tyrannosauridae"),
                            name = c("Ceratopsidae", "Hadrosauridae", 
                                     "Tyrannosauridae"))

(p6 <- ggplot2::ggplot(data = test2, aes(x = new_bins, y = Value)) +
  ylab("Detection probability difference") + 
  xlab("Time (Ma)") +
  scale_x_reverse() +
  deeptime::coord_geo(dat = list("stages"), 
                      ylim = c(0, 1), 
                      expand = T) +
  geom_line(aes(color = Difference)) +
  ggthemes::theme_few() +
  rphylopic::geom_phylopic(data = silhouette_df, aes(x = x, y = y), 
                uuid = c(c.uuid, h.uuid, t.uuid), 
                size = c(0.18, 0.18, 0.18), 
                alpha = 1, 
                color = "grey") +
  facet_grid(c("Target"), 
             scales = "free") +
  scale_color_manual(values=c("#BDD7E7", "#6BAED6", "#3182BD")))

ggsave(paste("Figures/A1.LL.plot.png", sep = ""), plot = p6, 
       device = "png")
ggsave(paste("Figures/A1.LL.plot.pdf", sep = ""), plot = p6, 
       device = "pdf")

################################################################################
# A3. UNMODELLED HETEROGENEITY OF DETECTION PROBABILTY IN SPACE
################################################################################

target <- c("Ceratopsidae", "Hadrosauridae", "Tyrannosauridae")
res <- 1
# Set extent
e <- extent(-155, -72, 22.5, 73)
new.e <- extent(-122, -95, 32, 58)

raster_for_values <- stack(lapply(target, function(x){
  # Load model 
  out.sp <- readRDS(paste("Results/spOccupancy/", res, "/", x, 
                          ".best.model.rds", sep = ""))
  sp.data <- readRDS(file = paste("Prepped_data/spOccupancy/Multi_season/", res, "/", 
                                  x, "_multi_", res, ".rds", sep = ""))
  
  # Find random effects sizes
  alpha.star.means <- apply(out.sp$alpha.star.samples, 2, mean)
  
  # Create siteIDs (individual cell numbers used)
  siteIDs <- as.numeric(rownames(sp.data$y[,1,]))
  
  # Find coordinates of each site
  siteCoords <- siteCoordsFun(res = res, e = e, siteIDs)
  
  # Generate a raster using random effects
  raster_for_values <- gen_raster(siteCoords$siteID, alpha.star.means, res = res, ext = e)
  raster_for_values <- crop(raster_for_values, new.e)
  return(raster_for_values)
}))

# Find map to use as backdrop
countries <- maps::map("world", plot=FALSE, fill = TRUE) 
# Turn map into spatialpolygons
countries <<- maptools::map2SpatialPolygons(countries, 
                                            IDs = countries$names, 
                                            proj4string = CRS("+proj=longlat")) 
mapTheme <- rasterVis::rasterTheme(region=brewer.pal(8,"Reds"))

(p2 <- rasterVis::levelplot(raster_for_values,
                            margin=list(draw = T, 
                                        scales = list(y=c(0,0))), 
                            par.settings=mapTheme, 
                            at = seq(-1.5, 2, length.out=35)) + 
    # Plots state lines
    latticeExtra::layer(sp.polygons(states, col = "white", fill = NA), under = T)  + 
    # Plots background colour
    latticeExtra::layer(sp.polygons(countries, col = 0, fill = "light grey"), under = T))


raster_for_values[[1]]
raster_for_values[[2]]
raster_for_values[[3]]

################################################################################
# TABLE S3-4. TOP MODEL RESULTS
################################################################################





################################################################################
# A4. OCCUPANCY PREDICTION
################################################################################

# Load palaeoclimatic rasters
wc <- list.files(paste("Prepped_data/Covariate_Data/All_data/", 
                       res, "deg/Palaeo/", sep = ""), 
                 pattern=paste0("^", "teyen", ".*", sep = ""))
stacked <- raster::stack(paste("Prepped_data/Covariate_Data/All_data/", 
                               res, "deg/Palaeo/", wc, 
                               sep =""))
stacked <- dropLayer(stacked, c(1,2, 3,4, 5, 7, 9, 10))
names(stacked) <- gsub(".[[:digit:]]", "", names(stacked))
stacked <- crop(stacked, e)

covariates <- as.data.frame(stacked)
full.coords <- as.data.frame(xyFromCell(stacked, 1:16600))
covs <- cbind(full.coords, covariates)

# Number of prediction sites.
J.pred <- nrow(covs)

# Number of prediction years.
n.years.pred <- 1
# Number of predictors (including intercept)
p.occ <- ncol(out.ar1$beta.samples)
# Get covariates and standardize them using values used to fit the model
hot <- (covs$hot_mean - mean(revi.data$occ.covs$hot$teyen)) / sd(revi.data$occ.covs$hot$teyen)
wet <- (covs$wet_mean - mean(revi.data$occ.covs$wet$teyen)) / sd(revi.data$occ.covs$wet$teyen)
dry <- (covs$dry_mean - mean(revi.data$occ.covs$dry$teyen)) / sd(revi.data$occ.covs$dry$teyen)

# Create three-dimensional array
X.0 <- array(1, dim = c(J.pred, n.years.pred, p.occ))
# Fill in the array
# Years
X.0[, , 2] <- hot
# Elevation
X.0[, , 3] <- wet
# Elevation^2
X.0[, , 4] <- dry
# Check out the structure
str(X.0)

# Indicate which primary time periods (years) we are predicting for
t.cols <- c(4)
# Approx. run time: < 30 sec
out.pred <- predict(out.ar1, X.0, coords.0 = full.coords, t.cols = t.cols, ignore.RE = TRUE, type = 'occupancy')

# Plotting
plot.dat <- data.frame(x = full.coords$x, 
                       y = full.coords$y, 
                       mean.psi = apply(out.pred$psi.0.samples[, , 1], 2, mean), 
                       sd.psi = apply(out.pred$psi.0.samples[, , 1], 2, sd), 
                       stringsAsFactors = FALSE)
# Make a species distribution map showing the point estimates,
# or predictions (posterior means)
dat.stars <- st_as_stars(plot.dat, dims = c('x', 'y'))
# 2009
ggplot() + 
  geom_stars(data = dat.stars, aes(x = x, y = y, fill = mean.psi)) +
  scale_fill_viridis_c(na.value = 'transparent') +
  labs(x = 'Easting', y = 'Northing', fill = '', 
       title = 'Mean REVI occurrence probability 2009') +
  theme_bw()
