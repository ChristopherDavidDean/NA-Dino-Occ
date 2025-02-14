################################################################################
############ OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS #############
################################################################################

# Christopher D. Dean, Alfio Alessandro Chiarenza, Jeffrey W. Doser, Alexander
# Farnsworth, Lewis A. Jones, Sinéad Lyster, Charlotte L. Outhwaite, Paul J. 
# Valdes, Richard J. Butler, Philip D. Mannion.
# 2025
# Script written by Christopher D. Dean

################################################################################
#                    FILE 7: MULTI-SEASON OCCUPANCY FIGURES                   #
################################################################################

################################################################################
# 1. INITIAL SETUP
################################################################################

library(ggforestplot)
library(ggforce)
library(ggh4x)
library(rnaturalearthdata)
library(rnaturalearth)
library(sf)

##### Load in Functions #####
source("0.Functions.R") # Import functions from other R file (must be in same working directory)

# Setup phylopic
a.uuid <- get_uuid(name = "Ankylosauridae", n = 6)[[6]]
c.uuid <- get_uuid(name = "Ceratopsidae", n = 5)[[5]]
t.uuid <- get_uuid(name = "Tyrannosauridae", n = 5)[[5]]
h.uuid <- get_uuid(name = "Edmontosaurus", n = 3)[[3]]

a.uuid <- "5b062105-b6a2-4405-bd75-3d0399102b9a"
c.uuid <- "4b12a77a-0a16-4b93-baf5-7d4d01d0d9bd"
t.uuid <- "f3808e65-a95f-4df5-95a0-5f5b46a221f2"
h.uuid <- "72be89b9-3f2b-4dc3-b485-e74a5f8b1fbc"
  
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
target <- c("Ankylosauridae", 
            "Ceratopsidae", 
            "Hadrosauridae",
            "Tyrannosauridae")

# Choose type of model
type <- "hc.rw"

# Load data
all.results <- read.csv(paste("Results/Outhwaite/", bin.type,
                              "/results.", res, ".", type,".csv", sep = ""))

# Plotting occupancy (naive and modelled)
cera <- all.results %>%
  filter(Target == "Ceratopsidae")
tyran <- all.results %>%
  filter(Target == "Tyrannosauridae")
hadro <- all.results %>%
  filter(Target == "Hadrosauridae")
ankyl <- all.results %>%
  filter(Target == "Ankylosauridae")

# Plot modelled results
a <- plot_occ(hadro)
b <- plot_naive(hadro, h.uuid)
c <- plot_occ(tyran)
d <- plot_naive(tyran, t.uuid)
e <- plot_occ(cera)
f <- plot_naive(cera, c.uuid)
g <- plot_occ(ankyl)
h <- plot_naive(ankyl, a.uuid)

# Arrange
(p <- ggarrange(h, f, b, d, g, e, a, c, 
                nrow = 2, ncol = 4,
                align='h', labels=c('A', 'B', 'C',
                                    'D', 'E', 'F', 
                                    'G', 'H'),
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


# Save figure
ggsave(paste("Figures/2.Sparta.", bin.type, ".", res, ".", type, ".png", sep = ""), plot = p, 
       device = "png")

ggsave(paste("Figures/2.Sparta.", bin.type, ".", res, ".", type, ".pdf", sep = ""), plot = p, 
       device = "pdf")


#########################
##### UPDATED PLOTS #####
#########################

# Set variables
bin.type <- "scotese"
type <- "hc.rw"

# Get results
res <- c(0.5, 1)
occ.table <- dplyr::bind_rows(lapply(res, function(res){
    a <- read.csv(paste("Results/Outhwaite/NEW/", bin.type ,"/results.", res, ".", type ,".csv", sep = ""))
    a$Res <- res
    return(a)
}))

# Setup phylopic
a.uuid <- "5b062105-b6a2-4405-bd75-3d0399102b9a"
c.uuid <- "4b12a77a-0a16-4b93-baf5-7d4d01d0d9bd"
t.uuid <- "f3808e65-a95f-4df5-95a0-5f5b46a221f2"
h.uuid <- "72be89b9-3f2b-4dc3-b485-e74a5f8b1fbc"


silhouette_df <- data.frame(x = c(72, 71.5, 72, 72), y = c(0.94, 0.94, 0.94, 0.94), 
                            Target = c("Ankylosauridae", "Ceratopsidae", "Hadrosauridae", 
                                      "Tyrannosauridae"),
                            Res = c("1", "1", "1", "1"),
                            name = c("Ankylosauridae", "Ceratopsidae", "Hadrosauridae", 
                                     "Tyrannosauridae"))

occ.table$Res <- as.factor(occ.table$Res)

# Plot naive occupancy
a <- ggplot(data = subset(occ.table, Data == "Naive occupancy"), aes(x = new_bins, 
                                                                y = value)) +
  ylab("Proportion of total sites") + 
  xlab("Time (Ma)") +
  scale_x_reverse() +
  geom_phylopic(data = silhouette_df, aes(x = x, y = y), 
                uuid = c(a.uuid, c.uuid, h.uuid, t.uuid), 
                size = c(0.14, 0.15, 0.14, 0.15), 
                alpha = 1, 
                color = "grey") +
  #deeptime::coord_geo(dat = list("stages"), 
  #                    xlim = c((max(all.results[[1]]$new_bins)+1), (min(all.results[[1]]$new_bins-1))), 
  #                    ylim = c(0, 1)) +
  coord_cartesian(ylim = c(0.01,1)) +
  geom_line(aes(x = new_bins, y = value)) +
 scale_color_manual(breaks = c("Naïve Occupancy", "PAO", "Occupancy Probability", "Detection Probability"),
                    values=c("#252424", "#DE2D26", "#DE2D26", "#3182BD")) +
 scale_fill_manual(breaks = c("Naïve Occupancy", "PAO", "Occupancy Probability", "Detection Probability"), 
                   values=c("#FFFFFF", "white","#DE2D26", "#3182BD")) +
   scale_color_manual(values=c("#252424")) +
  theme_bw() +
  geom_smooth(method=lm, fullrange = T, alpha=0.3) +
  facet_nested(. ~ Target + Res, 
               scales = "free_y",
               space = "free_y",
               labeller = label_wrap_gen(40),
               strip = strip_nested(size = "variable")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        plot.margin = unit(c(0.2,0.2,-0.1,1), 'lines'), 
        panel.spacing.x = unit(c(-0.1, 0.5, -0.1, 0.5, -0.1, 0.5, -0.1), "lines"), 
        axis.title.y = element_text(size = 12))


# Plot modelled occupancy and detection
b <- ggplot(data = subset(occ.table, Data == "Mean occupancy"), aes(x = new_bins, 
                                                               y = value)) +
  geom_blank(aes(color = Data), data = occ.table) +
  geom_ribbon(data = occ.table, aes(x = new_bins, ymin = lower95CI, 
                                           ymax = upper95CI, fill = Data), alpha = 0.2) +
  ylab("Probability") + 
  xlab("Time (Ma)") +
  scale_x_reverse() +
  deeptime::coord_geo(dat = list("stages"), 
                      xlim = c((max(occ.table$new_bins)+1), (min(occ.table$new_bins-1))), 
                      ylim = c(0, 1)) +
  geom_line(data = subset(occ.table, Data != "Naive occupancy"), 
            aes(x = new_bins, y = value, color = Data)) +
 scale_color_manual(values=c("#BDD7E7", "#6BAED6", "#3182BD", 
                             "#DE2D26", "#252424", "#DE2D26", "#252424")) +
 scale_fill_manual(values=c("#BDD7E7", "#6BAED6", "#3182BD", 
                            "#DE2D26", "#FFFFFF", "#DE2D26", "#252424")) +
  new_scale_color() +
  geom_point(data = subset(occ.table, Data == "Mean occupancy"),
             aes(x = new_bins, y = value, color = rhat_threshold), size = 2) +
  #scale_color_manual(name = 'Rhat', values = c('Bad (>1.1)' = 'white',
  #                                             'Good (<1.1)' = '#DE2D26')) +
  geom_smooth(method=lm, fullrange = T) +
  theme_bw() +
  facet_nested(. ~ Target + Res, 
               scales = "free_y",
               space = "free_y",
               labeller = label_wrap_gen(40),
               strip = strip_nested(size = "variable")) +
  theme(legend.position = "bottom", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.background = element_blank(),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        strip.text.x = element_blank(),
        panel.spacing.x = unit(c(-0.1, 0.5, -0.1, 0.5, -0.1, 0.5, -0.1), "lines"), 
        plot.margin = unit(c(-0.1,0.2, 1, 0.2), 'lines'), 
        axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12),)

# Combine
(p <- ggarrange(a, b, nrow = 2, heights = c(0.72, 1), align = 'v'))

# Save figure
ggsave(paste("Figures/NEW/2.Sparta.", bin.type, ".", type, ".png", sep = ""), plot = p, 
       device = "png")
ggsave(paste("Figures/NEW/2.Sparta.", bin.type, ".", type, ".pdf", sep = ""), plot = p, 
       device = "pdf")

################################################################################
# FIGURE 3: DETECTION THROUGH TIME, SPOCCUPANCY
################################################################################

# Set extent
e <- extent(-155, -72, 22.5, 73)
max_val <- 40

# Set target
target <- c("Ankylosauridae", 
            "Ceratopsidae", 
            "Hadrosauridae",
            "Tyrannosauridae")

all.p.vals <- dplyr::bind_rows(lapply(target, function(target){
  res <- c(0.5, 1)
  p.vals <- lapply(res, function(res){
    best.model <- readRDS(paste("Results/spOccupancy/NEW/", res, 
                                "/", target, "_", max_val, ".best.model.rds", sep = ""))
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

# Setup phylopic
silhouette_df <- data.frame(x = c(71, 71, 71, 71), y = c(0.92, 0.92, 0.92, 0.92), 
                            Group = c("Ankylosauridae", "Ceratopsidae", "Hadrosauridae", 
                                      "Tyrannosauridae"),
                            Resolution = c("0.5", "0.5", "0.5", "0.5"),
                            name = c("Ankylosauridae", "Ceratopsidae", "Hadrosauridae", 
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
  scale_color_manual(values=c("#3182BD", "#3182BD", "#3182BD", "#3182BD")) +
  scale_fill_manual(values=c("#3182BD", "#3182BD", "#3182BD", "#3182BD"))+
  #viridis::scale_color_viridis(discrete=TRUE) + 
  #viridis::scale_fill_viridis(discrete = TRUE) +
  theme(legend.position="none", 
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  geom_phylopic(data = silhouette_df, aes(x = x, y = y), 
                uuid = c(a.uuid, c.uuid, h.uuid, t.uuid), 
                size = c(0.15, 0.15, 0.15, 0.17), 
                alpha = 1, 
                color = "grey") +
  facet_grid(c("Resolution", "Group"), 
             scales = "free"))

# Save figure
ggsave(paste("Figures/3.spOcc.detection", "_", max_val,".png", sep = ""), plot = p2, 
       device = "png")

ggsave(paste("Figures/3.spOcc.detection","_", max_val, ".pdf", sep = ""), plot = p2, 
       device = "pdf")

##### INCLUDING ALL MAX VALS #####
# Set extent
e <- extent(-155, -72, 22.5, 73)

# Set target
target <- c("Ankylosauridae", 
            "Ceratopsidae", 
            "Hadrosauridae",
            "Tyrannosauridae")

all.p.vals <- data.frame()
for(t in c(10, 40)){
  temp.all.p.vals <- dplyr::bind_rows(lapply(target, function(target){
    res <- c(0.5, 1)
    p.vals <- lapply(res, function(res){
      best.model <- readRDS(paste("Results/spOccupancy/NEW/", res, 
                                  "/", target, "_", t, ".best.model.rds", sep = ""))
      siteCoords <- siteCoordsFun(res = res, e = e, as.numeric(rownames(best.model$y))) 
      fit.mod <- fitted(best.model)
      p.samples <- fit.mod$p.samples
      p.samples <- colMeans(p.samples)
      p.vals <- as.data.frame(apply(p.samples,c(1,2),mean, na.rm = T))
      mean.p <- as.data.frame(colMeans(p.vals, na.rm = T))
      CIs <- t(apply(p.vals, 2, quantile, c(0.025, 0.975), na.rm = T))
      colnames(mean.p) <- "p"
      mean.p$Resolution <- res
      mean.p$Sample_cap <- t
      mean.p$Group <- target
      mean.p$Bin <- c("teyeq", "teyep", "teyeo", "teyen")
      mean.p <- cbind(mean.p, CIs)
      return(mean.p)
    })
    return(p.vals)
  }))
  all.p.vals <- rbind(temp.all.p.vals, all.p.vals)
}

all.p.vals[all.p.vals == "teyeq"] <- 80.8
all.p.vals[all.p.vals == "teyep"] <- 75
all.p.vals[all.p.vals == "teyeo"] <- 69
all.p.vals[all.p.vals == "teyen"] <- 66.7

# Setup phylopic
silhouette_df <- data.frame(x = c(73, 72, 73, 73), y = c(0.92, 0.92, 0.92, 0.92), 
                            Group = c("Ankylosauridae", "Ceratopsidae", "Hadrosauridae", 
                                      "Tyrannosauridae"),
                            Resolution = c("1", "1", "1", "1"),
                            Sample_cap = c("10", "10", "10", "10"),
                            name = c("Ankylosauridae", "Ceratopsidae", "Hadrosauridae", 
                                     "Tyrannosauridae"))

all.p.vals$Bin <- as.numeric(all.p.vals$Bin)
all.p.vals$Group <- as.factor(all.p.vals$Group)
all.p.vals$Resolution <- as.factor(all.p.vals$Resolution)
all.p.vals$Sample_cap <- as.factor(all.p.vals$Sample_cap)

(p2 <- ggplot(data = all.p.vals, aes(x = Bin, y = p)) +
    ylab("Detection probability") + 
    xlab("Time (Ma)") +
    scale_x_reverse() +
    deeptime::coord_geo(dat = list("stages"), 
                        expand = T, 
                        ylim = c(0, 1)) +
    geom_line(aes(color = Group)) +
    theme_bw() +
    geom_ribbon(aes(x = Bin, ymin = `2.5%`, 
                    ymax = `97.5%`, fill = Group), alpha = 0.2) +
    scale_color_manual(values=c("#3182BD", "#3182BD", "#3182BD", "#3182BD")) +
    scale_fill_manual(values=c("#3182BD", "#3182BD", "#3182BD", "#3182BD"))+
    #viridis::scale_color_viridis(discrete=TRUE) + 
    #viridis::scale_fill_viridis(discrete = TRUE) +
    #theme(legend.position="none", 
    #      strip.background = element_blank(),
    #      strip.text.x = element_blank()) +
    geom_phylopic(data = silhouette_df, aes(x = x, y = y), 
                  uuid = c(a.uuid, c.uuid, h.uuid, t.uuid), 
                  size = c(0.12, 0.12, 0.12, 0.12), 
                  alpha = 1, 
                  color = "grey") +
    facet_nested(Sample_cap ~ Group + Resolution, 
                 scales = "free_y",
                 space = "free_y",
                 labeller = label_wrap_gen(40),
                 strip = strip_nested(size = "variable")) +
  theme(legend.position="none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.spacing.x = unit(c(-0.1, 0.5, -0.1, 0.5, -0.1, 0.5, -0.1), "lines"), 
        panel.spacing.y = unit(c(-0.1), "lines"), 
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10), 
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12)))

# Save figure
ggsave(paste("Figures/NEW/3.spOcc.detection.png", sep = ""), plot = p2, 
       device = "png")

ggsave(paste("Figures/NEW/3.spOcc.detection.pdf", sep = ""), plot = p2, 
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
# Set max_val
max_val <- 10

# Get chosen model
target <- c("Ankylosauridae", "Ceratopsidae", "Hadrosauridae", "Tyrannosauridae")

test <-  lapply(target, function(target){
  best.model <- readRDS(paste("Results/spOccupancy/NEW/", res, 
                              "/", target, "_", max_val, ".best.model.rds", sep = ""))
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
  best.model <- readRDS(paste("Results/spOccupancy/NEW/", res, 
                              "/", target, "_", max_val, ".best.model.rds", sep = ""))
  siteCoords <- siteCoordsFun(res = res, e = e, as.numeric(rownames(best.model$y))) 
  fit.mod <- fitted(best.model)
  p.samples <- fit.mod$p.samples
  p.samples <- colMeans(p.samples)
  p.vals <- as.data.frame(apply(p.samples,c(1,2),mean, na.rm = T))
  return(p.vals)
})

p.vals <- bind_rows(p.vals)

test <- stack(test[[1]], test[[2]], test[[3]], test[[4]])
raster.names <- c("Bin 1 (80.8 Ma)", "Bin 2 (75 Ma)", "Bin 3 (69 Ma)", "Bin 4 (66.7 Ma)", 
                  "", "", "", "", "", "", "", "", "", "", "", "")

# Turn map into spatialpolygons
us <- as_Spatial(ne_states(country = "united states of america", returnclass = "sf"))
ca <- as_Spatial(ne_states(country = "canada", returnclass = "sf"))
mx <- as_Spatial(ne_states(country = "mexico", returnclass = "sf"))

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
p3 <- rasterVis::levelplot(test, layout = c(4, 4), margin=T, par.settings=my.theme,
                           colorkey=my.ckey, names.attr = raster.names) + 
  # Plots state lines
  latticeExtra::layer(sp.polygons(teyen, col = 0, fill = "dark grey"), packets = c(4,8,12,16), under = T) +
  latticeExtra::layer(sp.polygons(teyeo, col = 0, fill = "dark grey"), packets = c(3,7,11,15), under = T) +
  latticeExtra::layer(sp.polygons(teyep, col = 0, fill = "dark grey"), packets = c(2,6,10,14), under = T) +
  latticeExtra::layer(sp.polygons(teyeq, col = 0, fill = "dark grey"), packets = c(1,5,9,13), under = T) +
  latticeExtra::layer(sp.polygons(us, col = "white", fill = "#dcdcdc", lwd = 0.5), under = T) +
  latticeExtra::layer(sp.polygons(mx, col = "white", fill = "#dcdcdc", lwd = 0.5), under = T) +
  latticeExtra::layer(sp.polygons(ca, col = "white", fill = "#dcdcdc", lwd = 0.5), under = T)

# Make ggplot
p3 <- ggplotify::as.ggplot(p3)

# Set phylopic requirements
silhouette_df <- data.frame(x = c(0.17,0.165,0.17,0.17), y = c(0.78, 0.555,0.33,0.105), 
                            Group = c("Ankylosauridae", "Ceratopsidae", "Hadrosauridae", 
                                      "Tyrannosauridae"),
                            Submodel = c("Detection", "Detection", "Detection", "Detection"),
                            name = c("Ankylosauridae", "Ceratopsidae", "Hadrosauridae", 
                                     "Tyrannosauridae"))
# Combine with phylopic
p3.5 <- p3 + geom_phylopic(data = silhouette_df, aes(x = x, y = y), 
                   uuid = c(a.uuid, c.uuid, h.uuid, t.uuid), 
                   size = c(0.033, 0.033, 0.034, 0.037), 
                   alpha = 1, 
                   color = "dark grey")

# Save figure
ggsave(paste("Figures/NEW/4.5.spOcc.detection.space.", res, "_", max_val,".png", sep = ""), plot = p3.5, 
       device = "png")

ggsave(paste("Figures/NEW/4.5.spOcc.detection.space.", res, "_", max_val,".pdf", sep = ""), plot = p3.5, 
       device = "pdf")

################################################################################
# 5. MULTI-SEASON FIGURES
################################################################################

###################
##### SPATIAL #####
###################

# Set max_val
max_val <- 40

target <- c("Ankylosauridae", 
            "Ceratopsidae", 
            "Hadrosauridae",
            "Tyrannosauridae")

# Get table
cov.table <- dplyr::bind_rows(lapply(target, function(target){
  res <- c(0.5, 1)
    all.covs <- lapply(res, function(res){
      a <- readRDS(paste("Results/spOccupancy/NEW/", res, "/", target, "_", max_val,
                         ".best.model.rds", sep = ""))
      b <- make_table(out.sp = a, res = res, target = target)
    })
  return(all.covs)
}))

# Clean and prep table
cov.table <- clean_for_fig(cov.table)
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

##### 40 cap #####
# Set up phylopic
silhouette_df <- data.frame(x = c(3.8, 3.8, 3.8, 3.8), y = c(4.5, 11, 5.5, 3.7), 
                            Group = c("Ankylosauridae", "Ceratopsidae", "Hadrosauridae", 
                                      "Tyrannosauridae"),
                            Submodel = c("Detection","Detection", "Detection", "Detection"),
                            name = c("Ankylosauridae", "Ceratopsidae", "Hadrosauridae", 
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
                  uuid = c(a.uuid, c.uuid, h.uuid, t.uuid), 
                  size = c(1, 3, 1.2, 0.95), 
                  alpha = 1, 
                  color = "dark grey") +
    ggplot2::facet_wrap(
      facets = ~Group + Submodel,
      nrow = 4, 
      ncol = 2,
      scales = "free_y"
    ) + facetted_pos_scales(
      y = list(Submodel == "Occupancy" ~ scale_y_discrete(position = "right"))
    )  +
    theme(strip.text.x = element_blank(), 
          legend.position = "bottom", 
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5), 
          axis.text.y.right = element_text(hjust = 0, size = 13), 
          axis.text.y.left = element_text(size = 13),
          axis.text.x = element_text(size = 13),
          axis.title.x = element_text(size = 14),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 14)) + 
    scale_colour_manual(values = wesanderson::wes_palette("Darjeeling1", type = 'discrete')))

##### 10 cap #####
# Set up phylopic
silhouette_df <- data.frame(x = c(3.8, 3.8, 3.8, 3.8), y = c(10, 8, 8, 3), 
                            Group = c("Ankylosauridae", "Ceratopsidae", "Hadrosauridae", 
                                      "Tyrannosauridae"),
                            Submodel = c("Detection","Detection", "Detection", "Detection"),
                            name = c("Ankylosauridae", "Ceratopsidae", "Hadrosauridae", 
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
                uuid = c(a.uuid, c.uuid, h.uuid, t.uuid), 
                size = c(2, 2, 1.8, 0.7), 
                alpha = 1, 
                color = "dark grey") +
  ggplot2::facet_wrap(
    facets = ~Group + Submodel,
    nrow = 4, 
    ncol = 2,
    scales = "free_y"
  ) + facetted_pos_scales(
    y = list(Submodel == "Occupancy" ~ scale_y_discrete(position = "right"))
  )  +
  theme(strip.text.x = element_blank(), 
        legend.position = "bottom", 
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5), 
        axis.text.y.right = element_text(hjust = 0, size = 13), 
        axis.text.y.left = element_text(size = 13),
        axis.text.x = element_text(size = 13),
        axis.title.x = element_text(size = 14),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14)) + 
  scale_colour_manual(values = wesanderson::wes_palette("Darjeeling1", type = 'discrete')))


# Save plot
ggsave(paste("Figures/NEW/5.Cov.plot", "_", max_val,".png", sep = ""), plot = p4, 
       device = "png")
ggsave(paste("Figures/NEW/5.Cov.plot", "_", max_val,".pdf", sep = ""), plot = p4, 
       device = "pdf")

#####################################################################
##### CHECKING PROBABILITY OF MARGINALLY SIGNIFICANT COVARIATES #####
#####################################################################

# Set variables
marg.table <- data.frame()
for(t in c(10,40)){
  marg.table.temp <- dplyr::bind_rows(lapply(target, function(target){
    res <- c(0.5, 1)
    all.covs <- lapply(res, function(res){
      a <- readRDS(paste("Results/spOccupancy/NEW/", res, "/", target, "_", t,
                         ".best.model.rds", sep = ""))
      b <- apply(a$alpha.samples, 2, function(a) mean(a > 0))
      c <- data.frame(Target = rep(target, length(b)), Res = res, Sample_cap = t, 
                      Covariate = clean_for_fig(names(b)), 
                      Value = b)
      rownames(c) <- NULL
      return(c)
    })
    return(all.covs)
  }))
  marg.table <- rbind(marg.table.temp, marg.table)
}
marg.table <- subset(marg.table, Covariate != "Intercept")

# Get full covariate table
cov.table.2 <- data.frame()
for(t in c(10,40)){
  cov.table.temp <- dplyr::bind_rows(lapply(target, function(target){
    res <- c(0.5, 1)
    all.covs <- lapply(res, function(res){
      a <- readRDS(paste("Results/spOccupancy/NEW/", res, "/", target, "_", t,
                         ".best.model.rds", sep = ""))
      b <- make_table(out.sp = a, res = res, target = target)
      b$Sample_cap <- rep(t, nrow(b))
      return(b)
    })
    return(all.covs)
  }))
  cov.table.2 <- rbind(cov.table.temp, cov.table.2)
}

# Clean table
cov.table.2 <- clean_for_fig(cov.table.2)
cov.table.2$Resolution <- as.factor(cov.table.2$Res)
cov.table.2 <- cov.table.2[order(cov.table.2$Group),]
cov.table.2 <- cov.table.2[order(cov.table.2$Resolution),]
cov.table.2 <- cov.table.2 %>%
  dplyr::rowwise() %>%
  dplyr::mutate(Sig = as.numeric(!between(0,`2.5%`,`97.5%`))) 
cov.table.2 <- as.data.frame(cov.table.2)

# Keep detection covariates, remove ones not necessary for final table
cov.table.2 <- cov.table.2 %>%
  dplyr::filter(Submodel == "Detection") %>%
  dplyr::select(Res, Group, Covariate, Sample_cap, Sig) %>%
  filter(Covariate != "Intercept" & Covariate != "REV: Site") %>%
  rename(Target = Group)
marg.sig <- merge(cov.table.2, marg.table)

# Select only non-significant covariates
marg.sig <- marg.sig %>%
  dplyr::filter(Sig == 0)

write.csv(x = marg.sig, file = "Results/spOccupancy/NEW/Marginal_signiciance.csv")

#######################
##### NON-SPATIAL #####
#######################

# Get table
cov.table <- dplyr::bind_rows(lapply(target, function(target){
  res <- c(0.5, 1)
  all.covs <- lapply(res, function(res){
    a <- readRDS(paste("Results/spOccupancy/Non.spatial/", res, "/", target, 
                       ".best.model.rds", sep = ""))
    b <- make_table(out.sp = a, res = res, target = target)
  })
  return(all.covs)
}))

cov.table <- clean_for_fig(cov.table)
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
# A1. VARIATIONS IN SPARTA COMBINED FIGURE
################################################################################

type <- "sp"
type <- "hc.rw"

# Load data
scot.5 <- read.csv(paste("Results/Outhwaite/NEW/scotese/results.0.5.", type, ".csv", sep = ""))
scot.1 <- read.csv(paste("Results/Outhwaite/NEW/scotese/results.1.", type, ".csv", sep = ""))
form.5 <- read.csv(paste("Results/Outhwaite/NEW/formation/results.0.5.", type, ".csv", sep = ""))
form.1 <- read.csv(paste("Results/Outhwaite/NEW/formation/results.1.", type, ".csv", sep = ""))

plot.combine <- function(all.results){
  # Plotting occupancy (naive and modelled)
  ankyl <- all.results %>%
    filter(Target == "Ankylosauridae")
  cera <- all.results %>%
    filter(Target == "Ceratopsidae")
  tyran <- all.results %>%
    filter(Target == "Tyrannosauridae")
  hadro <- all.results %>%
    filter(Target == "Hadrosauridae")
  # Plot modelled results
  a <- plot_occ(ankyl, x = F)
  b <- plot_occ(cera, x = F)
  c <- plot_occ(hadro, x = F)
  d <- plot_occ(tyran, x = F)
  test <- list(a, b, c, d)
  return(test)
}
a <- plot.combine(scot.5)
b <- plot.combine(scot.1)
c <- plot.combine(form.5)
d <- plot.combine(form.1)

# Arrange
(p <- ggarrange(a[[1]], a[[2]], a[[3]], a[[4]], 
                b[[1]], b[[2]], b[[3]], b[[4]], 
                c[[1]], c[[2]], c[[3]], c[[4]], 
                d[[1]], d[[2]], d[[3]], d[[4]], 
                
                nrow = 4, ncol = 4,
                align='h', labels=c('A', 'B', 'C', 'D', 
                                    'E', 'F', 'G', 'H', 
                                    'I', 'J', 'K', 'L', 
                                    'M', 'N', 'O', 'P'),
                legend = "bottom",
                common.legend = T))
# Save figure
ggsave(paste("Figures/NEW/A1.Sparta.compare.", type, ".png", sep = ""), plot = p, 
       device = "png", width = 250, height = 297, unit = 'mm')

ggsave(paste("Figures/NEW/A1.Sparta.compare.", type, ".pdf", sep = ""), plot = p, 
       device = "pdf", width = 250, height = 297, unit = 'mm')

################################################################################
# A2. LIST LENGTH DIFFERENCES
################################################################################

res <- 0.5
type <- "sp"

test <- read.csv(paste("Results/Outhwaite/NEW/scotese/results.", res, ".", type, ".csv", sep = ""))

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

silhouette_df <- data.frame(x = c(68, 68, 68, 68 ), y = c(0.86, 0.86, 0.86, 0.86), 
                            Target = c("Ankylosauridae", "Ceratopsidae", "Hadrosauridae", 
                                      "Tyrannosauridae"),
                            name = c("Ankylosauridae" , "Ceratopsidae", "Hadrosauridae", 
                                     "Tyrannosauridae"))

(p7 <- ggplot2::ggplot(data = test2, aes(x = new_bins, y = Value)) +
  ylab("Detection probability difference") + 
  xlab("Time (Ma)") +
  scale_x_reverse() +
  deeptime::coord_geo(dat = list("stages"), 
                      ylim = c(0, 1), 
                      expand = T) +
  geom_line(aes(color = Difference)) +
  ggthemes::theme_few() +
  rphylopic::geom_phylopic(data = silhouette_df, aes(x = x, y = y), 
                uuid = c(a.uuid, c.uuid, h.uuid, t.uuid), 
                size = c(0.18, 0.18, 0.18, 0.18), 
                alpha = 1, 
                color = "grey") +
  facet_grid(c("Target"), 
             scales = "free") +
  scale_color_manual(values=c("#BDD7E7", "#6BAED6", "#3182BD")))

ggsave(paste("Figures/NEW/A1.LL.", res, ".", type, ".plot.png", sep = ""), plot = p7, 
       device = "png")
ggsave(paste("Figures/NEW/A1.LL.", res, ".", type, ".plot.pdf", sep = ""), plot = p7, 
       device = "pdf")

################################################################################
# A3. UNMODELLED HETEROGENEITY OF DETECTION PROBABILTY IN SPACE
################################################################################

target <- c("Ankylosauridae", "Ceratopsidae", "Hadrosauridae", "Tyrannosauridae")
res <- 1
# Set extent
e <- extent(-155, -72, 22.5, 73)
new.e <- extent(-122, -95, 32, 58)

raster_for_values <- stack(lapply(target, function(target){
  max_val <- c(10, 40)
  stacked <- stack(lapply(max_val, function(max_val){
    # Load model 
    out.sp <- readRDS(paste("Results/spOccupancy/NEW/", res, "/", target, "_", max_val,
                            ".best.model.rds", sep = ""))
    sp.data <- readRDS(file = paste("Prepped_data/spOccupancy/Updated/Multi_season/", res, "/", 
                                    target, "_multi_", res, "_", max_val, ".rds", sep = ""))
    
    # Find random effects sizes
    alpha.star.means <- apply(out.sp$alpha.star.samples, 2, mean)
    
    # Create siteIDs (individual cell numbers used)
    siteIDs <- as.numeric(rownames(sp.data$y[,1,]))
    # Find coordinates of each site
    siteCoords <- siteCoordsFun(res = res, e = e, siteIDs)
    
    # Generate a raster using random effects
    raster_for_values <- gen_raster(siteCoords$siteID, alpha.star.means, res = res, ext = e)
    raster_for_values <- crop(raster_for_values, new.e)
    names(raster_for_values) <- paste(target, ":", max_val, sep = "")
    return(raster_for_values)
  }))
  return(stacked)
}))

mapTheme <- rasterVis::rasterTheme(region=brewer.pal(8,"Reds"))
us <- as_Spatial(ne_states(country = "united states of america", returnclass = "sf"))
ca <- as_Spatial(ne_states(country = "canada", returnclass = "sf"))
mx <- as_Spatial(ne_states(country = "mexico", returnclass = "sf"))

raster_for_values <- subset(raster_for_values, (c(1,5,2,4,3,7,6,8)))

(p2 <- rasterVis::levelplot(raster_for_values,
                            layout = c(4, 2),
                            margin=list(draw = T, 
                                        scales = list(y=c(0,0))), 
                            par.settings=mapTheme, 
                            at = seq(-1.5, 2, length.out=35)) + 
    latticeExtra::layer(sp.polygons(us, col = "white", fill = "#dcdcdc", lwd = 0.5), under = T) +
    latticeExtra::layer(sp.polygons(mx, col = "white", fill = "#dcdcdc", lwd = 0.5), under = T) +
    latticeExtra::layer(sp.polygons(ca, col = "white", fill = "#dcdcdc", lwd = 0.5), under = T))

# Make ggplot
p2 <- ggplotify::as.ggplot(p2)

# Set phylopic requirements
silhouette_df <- data.frame(x = c(0.28, 0.47, 0.64, 0.81), y = c(0.92, 0.92, 0.92, 0.92), 
                            Group = c("Ankylosauridae", "Ceratopsidae", "Hadrosauridae", 
                                      "Tyrannosauridae"),
                            name = c("Ankylosauridae", "Ceratopsidae", "Hadrosauridae", 
                                     "Tyrannosauridae"))
# Combine with phylopic
(p2.5 <- p2 + geom_phylopic(data = silhouette_df, aes(x = x, y = y), 
                           uuid = c(a.uuid, c.uuid, h.uuid, t.uuid), 
                           size = c(0.04, 0.04, 0.04, 0.042), 
                           alpha = 1, 
                           color = "dark grey"))

ggsave(paste("Figures/NEW/A4.unmodelled.png", sep = ""), plot = p2.5, 
       device = "png")
ggsave(paste("Figures/NEW/A4.unmodelled.pdf", sep = ""), plot = p2.5, 
       device = "pdf")

levelplot(raster_for_values[[3]], margin= T)

################################################################################
# A4. OCCUPANCY PREDICTION
################################################################################
#
## Load palaeoclimatic rasters
#wc <- list.files(paste("Prepped_data/Covariate_Data/All_data/", 
#                       res, "deg/Palaeo/", sep = ""), 
#                 pattern=paste0("^", "teyen", ".*", sep = ""))
#stacked <- raster::stack(paste("Prepped_data/Covariate_Data/All_data/", 
#                               res, "deg/Palaeo/", wc, 
#                               sep =""))
#stacked <- dropLayer(stacked, c(1,2, 3,4, 5, 7, 9, 10))
#names(stacked) <- gsub(".[[:digit:]]", "", names(stacked))
#stacked <- crop(stacked, e)
#
#covariates <- as.data.frame(stacked)
#full.coords <- as.data.frame(xyFromCell(stacked, 1:16600))
#covs <- cbind(full.coords, covariates)
#
## Number of prediction sites.
#J.pred <- nrow(covs)
#
## Number of prediction years.
#n.years.pred <- 1
## Number of predictors (including intercept)
#p.occ <- ncol(out.ar1$beta.samples)
## Get covariates and standardize them using values used to fit the model
#hot <- (covs$hot_mean - mean(revi.data$occ.covs$hot$teyen)) / sd(revi.data$occ.covs$hot$teyen)
#wet <- (covs$wet_mean - mean(revi.data$occ.covs$wet$teyen)) / sd(revi.data$occ.covs$wet$teyen)
#dry <- (covs$dry_mean - mean(revi.data$occ.covs$dry$teyen)) / sd(revi.data$occ.covs$dry$teyen)
#
## Create three-dimensional array
#X.0 <- array(1, dim = c(J.pred, n.years.pred, p.occ))
## Fill in the array
## Years
#X.0[, , 2] <- hot
## Elevation
#X.0[, , 3] <- wet
## Elevation^2
#X.0[, , 4] <- dry
## Check out the structure
#str(X.0)
#
## Indicate which primary time periods (years) we are predicting for
#t.cols <- c(4)
## Approx. run time: < 30 sec
#out.pred <- predict(out.ar1, X.0, coords.0 = full.coords, t.cols = t.cols, ignore.RE = TRUE, type = 'occupancy')
#
## Plotting
#plot.dat <- data.frame(x = full.coords$x, 
#                       y = full.coords$y, 
#                       mean.psi = apply(out.pred$psi.0.samples[, , 1], 2, mean), 
#                       sd.psi = apply(out.pred$psi.0.samples[, , 1], 2, sd), 
#                       stringsAsFactors = FALSE)
## Make a species distribution map showing the point estimates,
## or predictions (posterior means)
#dat.stars <- st_as_stars(plot.dat, dims = c('x', 'y'))
## 2009
#ggplot() + 
#  geom_stars(data = dat.stars, aes(x = x, y = y, fill = mean.psi)) +
#  scale_fill_viridis_c(na.value = 'transparent') +
#  labs(x = 'Easting', y = 'Northing', fill = '', 
#       title = 'Mean REVI occurrence probability 2009') +
#  theme_bw()
#