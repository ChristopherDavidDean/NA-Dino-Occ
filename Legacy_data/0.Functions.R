#======================== FUNCTIONS FOR OCCUPANCY AND RANGE SIZE =========================#
#                                                                                         #
#      AUTHOR: CHRISTOPHER D. DEAN                                                        #
#                                                                                         #
#=========================================================================================#

#=========================== iPAK AND REQUIRED PACKAGES ==================================

# function that automatically installs necessary packages that the user is lacking.

#===== iPAK =====
ipak <- function(pkg){ # Function to install packages. Read in character vector of any packages required. 
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
} 

#===== REQUIRED PACKAGES =====

library(tidyr)
library(tibble)
library(raster)
library(plyr)
library(dplyr)
library(latticeExtra)
library(rasterVis)
library(sp)
library(maps)
library(maptools)
library(stringr)
library(ggplot2)
library(patchwork)
library(sf)
library(rgeos)
library(maps)
library(mapview)
library(mapdata)
library(ggplot2)

#==================================== GET_GRID ===========================================

# Creates a raster of chosen resolution, and attaches associated grid cell IDs to occurrences/collections
get_grid <- function(data, res, e){ # data is first output from combine_data (fossil.colls). Res is chosen resolution in degrees
  r <- raster(res = res, val = 1, ext = e) # Value must be added because extract uses values
  r <<- r
  crs <- r@crs
  xy <- SpatialPointsDataFrame(cbind.data.frame(data$lng, data$lat), data, proj4string = crs)
  Final <- raster::extract(r, xy, sp = TRUE, cellnumbers = TRUE)
  Final <<- as.data.frame(Final)
}

#=================================== GET_EXTENT ==========================================

# Setup raster for resolution and extent. Note: these values should be the same ones used for the file 1.Setup_occupancy_DR.
get_extent <- function(data){
  maxLat <- round_any((max(data$lat) + 1), 0.5)  #get value for defining extent, and increase by x for visualisation purposes
  minLat <- round_any((min(data$lat) - 1), 0.5)  #get value for defining extent, and increase by x for visualisation purposes
  maxLng <- round_any((max(data$lng) + 1), 0.5)  #get value for defining extent, and increase by x for visualisation purposes
  minLng <- round_any((min(data$lng) - 1), 0.5) #get value for defining extent, and increase by x for visualisation purposes
  e <<- extent(minLng, maxLng, minLat, maxLat) # build extent object
}

#=================================== FORMBIN_M1 ==========================================

FormBin_M1 <- function(occs, formations, bins) {

  M1_List <- list()# make an empty list of the occurrences in each bin
  Form_list <- split(occs, occs$formation)
  # Code for assigning formations and associated occurrences to bins
  for(b in 1:nrow(bins)){ # for each new bin
    temp_recs <- data.frame()
    for (f in 1:nrow(formations)){ 
      if (formations$max_age[f] >= bins[b,3] && formations$min_age[f] <= bins[b,4]){ # If Formation max. age is greater than Bin min. age AND if Formation min. age is less then Bin max. age (i.e. falls within bin at some point)
        temp_recs <- rbind(temp_recs, Form_list[[f]]) # Add occurrences from that formation to occurrence list
      }
    }
    if (nrow(temp_recs) > 0){
      temp_recs$bin_no <- b
    }
    M1_List[[b]] <- temp_recs
  }
  df <- do.call("rbind", M1_List)
  binned_occs <<- df
}

#============================= PRESENCE.ABSENCE.RASTER ===================================

# Function from https://amywhiteheadresearch.wordpress.com/2013/05/27/creating-a-presence-absence-raster-from-point-data/
# Credit to Amy Whitehead

presence.absence.raster <- function (mask.raster,species.data,raster.label="") {
  require(raster)
  
  # set the background cells in the raster to 0
  mask.raster[!is.na(mask.raster)] <- 0
  
  #set the cells that contain points to 1
  speciesRaster <- rasterize(species.data,mask.raster,field=1)
  speciesRaster <- merge(speciesRaster,mask.raster)
  
  #label the raster
  names(speciesRaster) <- raster.label
  return(speciesRaster)
}

#================================== OCCUPANCY_EST ========================================

occupancy_est <- function(fossils, res, bins, formations, level = "genus"){
  e <- get_extent(fossils) # Get extent
  binned_occs <- FormBin_M1(fossils, formations, bins) # Bin occurrences
  if(level == "family"){
    total_family_occupancy <- data.frame(stringsAsFactors = FALSE)
  }
  if(level == "genus"){
    total_genus_occupancy <- data.frame(stringsAsFactors = FALSE)
  }
  
  all_total_cells <- c()
  
  for (t in rev(unique(binned_occs$bin_no))){ # For each bin
    temp_data <- binned_occs %>%
      dplyr::filter(bin_no == t) 
    if(level == "family"){
      names <- unique(temp_data$family_name)
    }
    if (level == "genus"){
      names <- unique(temp_data$occurrence.genus_name)
    }
    
    grid_data <- get_grid(temp_data, res, e) # Make raster and grid up cells
    total_cells <- length(unique(grid_data$cells))
    all_total_cells <- c(all_total_cells, total_cells)
    
    temp_stack <- stack()

    if (level == "family"){
      family_occupancy <- data.frame(stringsAsFactors = FALSE)
      for (f in names){
        temp_occ_data <- grid_data %>%
          dplyr::filter(family_name == f)
        temp_vec <- c(t, f, length(unique(temp_occ_data$cells)), length(unique(temp_occ_data$cells))/total_cells*100)
        family_occupancy <- rbind(family_occupancy, temp_vec, stringsAsFactors = FALSE)
        
        r <- raster(res = res, val = 1, ext = e)
        xy <- cbind(as.double(temp_occ_data$lng), as.double(temp_occ_data$lat))  
        occ_raster <- presence.absence.raster(mask.raster = r, species.data = xy, raster.label = f)
        temp_stack <- stack(temp_stack, occ_raster)
      }
      colnames(family_occupancy) <- c("Bin", "Family", "Occupied_cells", "Perc_Occupied_cells")
      total_family_occupancy <- rbind(total_family_occupancy, family_occupancy, stringsAsFactors = FALSE)
    }
    if (level == "genus"){
      genus_occupancy <- data.frame(stringsAsFactors = FALSE)
      for (f in names){
        temp_occ_data <- grid_data %>%
          dplyr::filter(occurrence.genus_name == f)
        temp_vec <- c(t, as.character(temp_occ_data$family_name[[1]]), f, length(unique(temp_occ_data$cells)), length(unique(temp_occ_data$cells))/total_cells*100)
        genus_occupancy <- rbind(genus_occupancy, temp_vec, stringsAsFactors = FALSE)
        
        r <- raster(res = res, val = 1, ext = e)
        xy <- cbind(as.double(temp_occ_data$lng), as.double(temp_occ_data$lat))  
        occ_raster <- presence.absence.raster(mask.raster = r, species.data = xy, raster.label = f)
        temp_stack <- stack(temp_stack, occ_raster)
      }
      colnames(genus_occupancy) <- c("Bin", "Family", "Genus", "Occupied_cells", "Perc_Occupied_cells")
      total_genus_occupancy <- rbind(total_genus_occupancy, genus_occupancy, stringsAsFactors = FALSE)
    }
    stack_name <- paste("Bin", t, "Stack", sep = "_")
    assign(stack_name, temp_stack, envir = .GlobalEnv)
  }
  if (level == "family"){
    reshaped <- total_family_occupancy[,-3]
    perc_fam_occ <- reshape(reshaped, idvar = "Bin", timevar = "Family", direction = "wide")

    reshaped <- total_family_occupancy[,-4]
    total_fam_occ <- reshape(reshaped, idvar = "Bin", timevar = "Family", direction = "wide")
    
    f.names <- colnames(perc_fam_occ)[2:length(perc_fam_occ)]
    f.names <- c("Bin", gsub(".*\\.","",f.names))
    colnames(perc_fam_occ) <- f.names
    colnames(total_fam_occ) <- f.names
    perc_fam_occ <<- perc_fam_occ
    total_fam_occ <<- total_fam_occ

    
    total_family_occupancy$Bin <- as.numeric(total_family_occupancy$Bin)
    total_family_occupancy$Perc_Occupied_cells <- as.numeric(total_family_occupancy$Perc_Occupied_cells)
    total_family_occupancy <<- total_family_occupancy

  }
  if (level == "genus"){
    names <- unique(total_genus_occupancy$Family)
    family_grouped_occupancy <- list()
    for (f in names){
      temp_fam_split <- total_genus_occupancy %>%
        dplyr::filter(Family == f) %>%
        dplyr::select(Bin, Genus, Perc_Occupied_cells) 
      temp_fam_split <- reshape(temp_fam_split, idvar = "Bin", timevar = "Genus", direction = "wide")
      temp_fam_split$Bin <- as.numeric(temp_fam_split$Bin)
      if (nrow(temp_fam_split) < nrow(bins)){
        temp_fam_split <- temp_fam_split %>%
          tidyr::complete(Bin = seq(1, nrow(bins), 1)) %>%
          purrr::map_df(rev)
      }
      g.names <- colnames(temp_fam_split)[2:length(temp_fam_split)]
      g.names <- c("Bin", gsub(".*\\.","",g.names))
      colnames(temp_fam_split) <- g.names
      family_grouped_occupancy[[f]] <- temp_fam_split
    }
    names(family_grouped_occupancy) <- names
    family_grouped_occupancy <<- family_grouped_occupancy
    
    total_genus_occupancy$Bin <- as.numeric(total_genus_occupancy$Bin)
    total_genus_occupancy$Perc_Occupied_cells <- as.numeric(total_genus_occupancy$Perc_Occupied_cells)
    total_genus_occupancy <<- total_genus_occupancy
  }
  all_total_cells <<- all_total_cells
}

#=================================== OCC_RASTER ==========================================

# Set up background info

occ_raster <- function(tax_name){
  regexp <- "[[:digit:]]+"
  stack_names <- rev(sort(grep("Stack",names(.GlobalEnv),value=TRUE)))
  tax_stack <- stack()
  for (t in stack_names){
    temp_stack <- get(t)
    if (tax_name %in% names(temp_stack) == TRUE){
      named_layer <- raster::subset(temp_stack, grep(tax_name, names(temp_stack), value = T))
      values(named_layer)[values(named_layer) == 0] = NA
      bin <- str_extract(t, regexp)
      names(named_layer) <- paste(tax_name, "_Bin_", bin, sep = "")
      tax_stack <- stack(tax_stack, named_layer)
    }
  }
  tax_stack <<- tax_stack
}


occ_raster_plotter <- function(tax_stack){
  countries <- maps::map("world", plot=FALSE, fill = TRUE) # find map to use as backdrop
  states <- maps::map("state", plot = FALSE, fill = TRUE)
  countries <<- maptools::map2SpatialPolygons(countries, IDs = countries$names, proj4string = CRS("+proj=longlat")) # Turn map into spatialpolygons
  states <<- maptools::map2SpatialPolygons(states, IDs = states$names, proj4string = CRS("+proj=longlat")) # Turn map into spatialpolygons
  mapTheme <- rasterVis::rasterTheme(region=brewer.pal(8,"Greens")) # MGVF
  
  title_name <- gsub("_.*","",names(tax_stack)[1])
  subbing <- paste(".*", title_name, "_", sep = "")
  names(tax_stack) <- sub(subbing, '', names(tax_stack))
  rasterVis::levelplot(tax_stack, margin=F, par.settings=mapTheme,  main = title_name) + #create levelplot for raster
    latticeExtra::layer(sp.polygons(states, col = "white", fill = NA), under = T)  +                   
    latticeExtra::layer(sp.polygons(countries, col = 0, fill = "light grey"), under = T)
}



#================================== OCC_PLOTTER =========================================

occplot <- function (data, title, bins){
  time <- read.csv("pbdb_time_int.csv")
  stages <- time[34:35,]
  stages$mid_ma <- (stages$max_ma + stages$min_ma)/2
  periods <- time[46,]
  s <- seq(1, nrow(bins), 2) #generate evens sequence for subsetting stages for plotting grey shading
  throwing_shade <- bins[s,]
  g.plot <- ggplot(data = data, aes()) +
    #geom_segment(data = periods, mapping=aes(x = min_ma, xend = min_ma, y =0, yend = 2), linetype = 2, size = 1, color = "grey90") + #draw period boundaries
    geom_rect(data = throwing_shade, mapping=aes(xmin=bottom, xmax=top, ymin=0, ymax= 50), linetype = 0, color="grey90", alpha=0.1)  + #draw stage shading
    #geom_rect(data = periods, mapping=aes(xmin=113, xmax=66, ymin=-0.1, ymax= 0), linetype = 1, colour = "black", fill="black", alpha=1)  + #draw black box underneath stages, personal preference, not needed
    geom_rect(data = periods, mapping=aes(xmin=83.6, xmax=66, ymin= -5, ymax= -2.5), linetype = 1, colour = "black", fill=periods$color, alpha=1)  + #draw period boxes below 0, ymin should be 5% of yend in line 14
    geom_rect(data = stages, mapping=aes(xmin=min_ma, xmax=max_ma, ymin= -2.5, ymax= 0.0), linetype = 1, colour = "black", fill=stages$color, alpha=1)  +
    geom_text(data = periods, mapping=aes(x=(66+83.6)/2, y= -3.75, label = abbrev), colour = "black", alpha=1)  + #draw period labels below 0, ymin should be 2.5% of yend in line 14
    geom_text(data = stages, mapping=aes(x=mid_ma, y= -1.25, label = abbrev), colour = "black", alpha=1)  +
    geom_line(aes(x = x, y = y), size = 0.8, alpha = 0.9) + # Data
    geom_point(aes(x=x, y=y), size = 1.2) +
    scale_x_reverse(expand=c(0,0), limits = c(83.6, 66)) +
    scale_y_continuous(expand=c(0,0), limits = c(-5, 50)) +
    labs(x = "Time (Ma)", y = "Occupancy", title = title) +
    theme(panel.background = element_blank(),
          plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),
          axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
          axis.text.y = element_text(size = 12, angle = 0, hjust = 0),
          axis.title.y = element_text(size = 12, face = "bold", vjust = 4),
          axis.title.y.right = element_text(size = 12, face = "bold", vjust = 4),
          axis.title.x = element_text(size = 12, face = "bold", vjust = -1),
          plot.title = element_text(size = 12, face = "bold", hjust = 0, vjust = 6),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          aspect.ratio = 0.6)
}

occ_plotter <- function(name, level){
  if (level == "family"){ # plots  occupancy for whole family individually
    family_data <- family_grouped_occupancy[[name]]
    plot.list <- list()
    for (f in 2:ncol(family_data)){
        temp_data <- data.frame(rev(bins$mid), as.numeric(family_data[,f]))
        colnames(temp_data) <- c("x", "y")
        g.plot <- occplot(data = temp_data, title = colnames(family_data)[f], bins = bins)
    plot.list[[f-1]] <- g.plot
    }
    print(wrap_plots(plot.list))
  }
  if (level == "genus"){ #  plots occupancy for chosen genera
    genus_data <- total_genus_occupancy %>%
      filter(Genus == name) %>% # filter data to just chosen genus
      tidyr::complete(Bin = seq(1, nrow(bins), 1)) %>% # Add any bins that the taxon isn't found in
      purrr::map_df(rev)  %>% # reverse order for plotting
      add_column(x = rev(bins$mid)) %>% # add column of bin midpoints for plotting
      select(x, Perc_Occupied_cells) %>% # only select columns for plotting
      rename(y = Perc_Occupied_cells) # rename to y
    g.plot <- occplot(data = genus_data, title = name, bins = bins)
      print(g.plot)
  }
  if (level == "family_total"){ # plots occupancy of the whole family
    family_data <- total_family_occupancy %>%
      filter(Family == name) %>% # filter data to just chosen genus
      tidyr::complete(Bin = seq(1, nrow(bins), 1)) %>% # Add any bins that the taxon isn't found in
      purrr::map_df(rev)  %>% # reverse order for plotting
      add_column(x = rev(bins$mid)) %>% # add column of bin midpoints for plotting
      select(x, 'Perc_Occupied_cells') %>% # only select columns for plotting
      rename(y = 'Perc_Occupied_cells') # rename to y
    f.plot <- occplot(data = family_data, title = name, bins = bins)
    print(f.plot)
  }
}

get_hulls <- function(bin_stack){
  df.all <- data.frame()
  for(s in 1:length(bin_stack@layers)){
    r <- bin_stack[[s]]
    r[r==0] <- NA
    test <- rasterToPoints(r, spatial = TRUE)
    df1 <- as.data.frame(test)
    df1[,1] <- names(bin_stack[[s]])
    colnames(df1)[1] <- "Species"
    df.all <- rbind(df.all, df1)
  }
  
  df.all <- df.all %>%
    st_as_sf( coords = c( "x", "y" ), crs = 4326 )
  
  hulls <- df.all %>%
    group_by(Species) %>%
    summarise( geometry = st_combine( geometry ) ) %>%
    st_convex_hull()
  
  hulls$area <- st_area(hulls) #
  
  plot(hulls[1])
}

#====================================== RUNNING ======================================


# Read in occs
master.occs <- read.csv("NA-Dino-Occ/Data/Occurrences/Master_spreadsheet_v1_051219.csv")
camp.occs <- read.csv("NA-Dino-Occ/Data/Occurrences/Camp_data_V1_Species_removed.csv", stringsAsFactors = FALSE) # Load in occurrences
maas.occs <- read.csv("NA-Dino-Occ/Data/Occurrences/Maas_data_V1_Species_removed.csv", stringsAsFactors = FALSE) # Load in occurrences

master.occs$mid_ma = (master.occs$max_ma + master.occs$min_ma)/2

occs <- read.csv("NA-Dino-Occ/Data/Occurrences_Final.csv")

# Read in bins
bins <- read.csv("NA-Dino-Occ/Data/Bins_SG2_res2.csv")
bins <- read.csv("NA-Dino-Occ/Data/Bins_SG2_res3.csv")

# Read in formations
formations <- read.csv("NA-Dino-Occ/Data/Formations_Final.csv")
formations$max_age <- as.numeric(as.character(formations$max_age)) # Make Numeric
formations$min_age <- as.numeric(as.character(formations$min_age)) # Make Numeric 
formations <- formations[order(formations$Formation),]
myformations <- sort(as.vector(formations$Formation))

# Select appropriate occurrences
testoccs <- occs[occs$formation %in% myformations,] # Only include occurrences from formation list
testoccs <- droplevels.data.frame(testoccs) # Remove old levels
formations <- formations[-(35),]


# Read in bins
res <- 0.5
get_extent(master.occs)

FormBin_M1(occs = testoccs, formations, bins)

#occupancy_est(fossils = binned_occs, res = res, bins = bins, formations = formations, level = "genus")
occupancy_est(fossils = binned_occs, res = res, bins = bins, formations = formations, level = "family")

occ_raster("Hadrosauridae")
occ_raster_plotter(tax_stack)

occ_plotter("Hadrosauridae", level = "family_total")
occ_plotter("Tyrannosauridae", level = "family_total")
occ_plotter("Ceratopsidae", level = "family_total")
occ_plotter("Hadrosauridae", level = "family")
occ_plotter("Edmontosaurus", level = "genus")

get_hulls(Bin_6_Stack)
get_hulls(Bin_5_Stack)
get_hulls(Bin_4_Stack)
get_hulls(Bin_3_Stack)
get_hulls(Bin_2_Stack)
get_hulls(Bin_1_Stack)

