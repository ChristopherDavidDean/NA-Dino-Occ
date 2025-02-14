################################################################################
############ OCCUPANCY OF LATE CRETACEOUS NORTH AMERICAN DINOSAURS #############
################################################################################

# Christopher D. Dean, Alfio Alessandro Chiarenza, Jeffrey W. Doser, Alexander
# Farnsworth, Lewis A. Jones, Sinéad Lyster, Charlotte L. Outhwaite, Paul J. 
# Valdes, Richard J. Butler, Philip D. Mannion.
# 2025
# Script written by Christopher D. Dean

################################################################################
#              FILE 0: FINDING GEOLOGICAL OUTCROP FROM MACROSTRAT              #
################################################################################

################################################################################
# 1. INITIAL SETUP
################################################################################

##### Load packages #####

library(geojsonsf)
library(sf)
library(tictoc)

##### Load in Functions #####
source("0.Functions.R") # Import functions from other R file (must be in same working directory)

# Set resolution
res <- 0.5
res <- 1
bin.res <- NA

# Set extent
e <- extent(-155, -72, 22.5, 73)
e <- extent(-155, -72, 23, 73)

# Set bin type
bin.type <- "scotese"

# Set formation based cells or not
form_cell <- "N"

# Load main dataset
master.occs <- read.csv(paste("Prepped_data/Occurrence_Data/", bin.type, "/", 
                              res, "_", bin.type, "_occurrence_dataset.csv", sep = ""))
# Load formations
formations <- read.csv("Data/Occurrences/Formations.csv") # Load formations

################################################################################
# 2. GET LIST OF FORMATION STRAT NAME IDS
################################################################################

# Make formations good for searching
test <- gsub(' ', '%20', gsub('\\.', '', unique(master.occs$formation)))

#find strat IDs for formations
test.forms <- lapply(test, function(x){
  URL <- paste("https://macrostrat.org/api/units?strat_name=", x, "&response=long&format=csv", sep = "")
  print(x)
  if(getURL(URL, .opts = curlOptions(ssl.verifypeer=FALSE, verbose=TRUE)) == "\r\n" |
     grepl("error", getURL(URL)) == T){
    strat_dframe = NULL
  }else{
    strat_dframe <- read.csv(URL, stringsAsFactors = FALSE)
  }
  print(x)
  strat_dframe$units_above <- as.character(strat_dframe$units_above)
  strat_dframe$units_below <- as.character(strat_dframe$units_below)
  return(strat_dframe)
})

# Discard null rows, bind together
test.forms <- test.forms %>% purrr::discard(is.null)
test.forms <- bind_rows(test.forms)

######### FILTERING ###########
test.forms2 <- test.forms %>%
  dplyr::filter(b_age > 66.6) %>%
  dplyr::filter(t_age < 83.7)

######### BINNING ###########
bins <- read.csv("Data/Occurrences/scotesebins.csv")
bins[1,5] <- 45
bins[4,4] <- 110
bins$bin <- bins$code
test.forms2$Formation <- test.forms2$Fm

combined.test <- merge(test.forms2, formations, by = "Formation")

names(combined.test)[names(combined.test) == 'min_age'] <- "min_ma"
names(combined.test)[names(combined.test) == 'max_age'] <- "max_ma"

combined.test <- bin_time(combined.test, bins, method = "all")

teyen.matched <- combined.test %>%
  dplyr::filter(bin_assignment == "teyen")
teyeo.matched <- combined.test %>%
  dplyr::filter(bin_assignment == "teyeo")
teyep.matched <- combined.test %>%
  dplyr::filter(bin_assignment == "teyep")
teyeq.matched <- combined.test %>%
  dplyr::filter(bin_assignment == "teyeq")

################################################################################
# 3. GET MAP IDS FOR NON-MATCHING UNITS
################################################################################

# Get list of strat names that couldn't be matched
miss.list <- c()
miss.list <- lapply(test, function(x){
  URL <- paste("https://macrostrat.org/api/units?strat_name=", x, "&response=long&format=csv", sep = "")
  if(getURL(URL, .opts = curlOptions(ssl.verifypeer=FALSE, verbose=TRUE)) == "\r\n" |
     grepl("error", getURL(URL)) == T){
    miss.list <- x
  }
  return(miss.list)
})

# Discard nulls, return to normal
miss.list <- miss.list %>% purrr::discard(is.null)
miss.list <- gsub('%20', ' ', unlist(miss.list))

# Find associated row that match missing list
miss.data <- master.occs[master.occs$formation %in% miss.list,]
miss.data <- miss.data %>%
  dplyr::select(formation, lat, lng) %>%
  dplyr::distinct() %>%
  dplyr::mutate(ID = 1:nrow(.))

# Get map IDS
map.ids <- lapply(miss.data$ID, function(x){
  test <- miss.data %>%
    filter(ID == x)
  lat <- test$lat
  lng <- test$lng
  URL <- paste("https://macrostrat.org/api/v2/geologic_units/map?lat=", lat, 
               "&lng=", lng, "&format=csv", sep = "")
  map_ids <- read.csv(URL, stringsAsFactors = FALSE)
  print(test)
  return(map_ids)
})

# Remerge together
test <- rapply(map.ids, as.character, how = "replace")
test <- bind_rows(test)
test$t_age <- as.numeric(test$t_age)
test$b_age <- as.numeric(test$b_age)

###### FILTERING WRONG VALUES ######
test2 <- test %>%
  dplyr::distinct() %>%
  dplyr::filter(strat_name != "") %>%
  dplyr::filter(b_age > 66.6) %>%
  dplyr::filter(t_age < 83.7)

# Clean
test2$Formation <- gsub(" Formation.*","",test2$strat_name)
test2$Formation <- gsub("Bearpaw","Bearpaw Shale",test2$Formation)
test2$Formation <- gsub("Lennep Sandstone Member","Sedan",test2$Formation)

# conjoin with formations dataset
test2 <- left_join(test2, formations, "Formation")

test3 <- test2 %>%
  #dplyr::filter(Depositional.environment != "Marine" | is.na(Depositional.environment) == T) %>%
  dplyr::filter(!grepl("volcanic",name)) %>%
  dplyr::filter(!grepl("Group",strat_name))

test3$max_age[is.na(test3$max_age) == T] <- test3$b_age[is.na(test3$max_age) == T]
test3$min_age[is.na(test3$min_age) == T] <- test3$t_age[is.na(test3$min_age) == T]

names(test3)[names(test3) == 'min_age'] <- "min_ma"
names(test3)[names(test3) == 'max_age'] <- "max_ma"

combined.test.2 <- bin_time(test3, bins, method = "all")
combined.test.2$map_id <- as.numeric(combined.test.2$map_id)

teyen.unmatched <- combined.test.2 %>%
  dplyr::filter(bin_assignment == "teyen")
teyeo.unmatched <- combined.test.2 %>%
  dplyr::filter(bin_assignment == "teyeo")
teyep.unmatched <- combined.test.2 %>%
  dplyr::filter(bin_assignment == "teyep")
teyeq.unmatched <- combined.test.2 %>%
  dplyr::filter(bin_assignment == "teyeq")

################################################################################
# 4. GETTING SHAPEFILES
################################################################################

type <- "strat_name_id"
type <- "map_id"

get.shapefiles <- function(ID.list, type){
  test.shape <- lapply(ID.list, function(a){
    # Get map units using strat_name
    print(a)
    y <- getURL(paste("https://macrostrat.org/api/v2/geologic_units/map?", type,
                        "=", a, "&format=geojson_bare", sep = ""))
    if(y == "{\"type\":\"FeatureCollection\",\"features\":[]}"){
      sf_map <- NULL
    }else{
      sf_temp <- geojsonsf::geojson_wkt(y)
      sf_map <- sf::st_as_sf(sf_temp, wkt = 'geometry', crs = 4326) 
    }
    return(sf_map)
  })
  test.shape <- test.shape %>% purrr::discard(is.null)
  return(test.shape)
}

combine.shapefiles <- function(matched, unmatched){
  matched.shp <- get.shapefiles(sort(unique(matched$strat_name_id)), type = "strat_name_id")
  matched.shp <- do.call(rbind, matched.shp)
  unmatched.shp <- get.shapefiles(unmatched$map_id, type = "map_id")
  unmatched.shp <- do.call(rbind, unmatched.shp)
  total.shp <- rbind(matched.shp, unmatched.shp)
  sf::st_crs(total.shp) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  
  return(total.shp)
}

# Get shapefiles for each bin
teyen.out <- combine.shapefiles(matched = teyen.matched, unmatched = teyen.unmatched)
teyeo.out <- combine.shapefiles(matched = teyeo.matched, unmatched = teyeo.unmatched)
teyep.out <- combine.shapefiles(matched = teyep.matched, unmatched = teyep.unmatched)
teyeq.out <- combine.shapefiles(matched = teyeq.matched, unmatched = teyeq.unmatched)

# Teyen
r1 <- raster(ext = e, res = res)
newData.n <- rasterize(teyen.out, r1, getCover = TRUE, progress = "window")
plot(newData.n)
teyen.occs <- master.occs %>%
  dplyr::filter(bin_assignment == "teyen")
points(teyen.occs$lng, teyen.occs$lat, pch = 3)
writeRaster(newData.n, paste("Prepped_data/Covariate_Data/All_data/", res, 
                           "deg/Out_teyen_", res, ".asc", sep = ""), 
            pattern = "ascii", overwrite = TRUE)
st_write(teyen.out, "Data/Covariate_Data/Outcrop/New_Shapefiles/teyen.shp")

# Teyeo
r1 <- raster(ext = e, res = res)
newData.o <- rasterize(teyeo.out, r1, getCover = TRUE, progress = "window")
plot(newData.o)
teyeo.occs <- master.occs %>%
  dplyr::filter(bin_assignment == "teyeo")
points(teyeo.occs$lng, teyeo.occs$lat, pch = 3)
writeRaster(newData.o, paste("Prepped_data/Covariate_Data/All_data/", res, 
                             "deg/Out_teyeo_", res, ".asc", sep = ""), 
            pattern = "ascii", overwrite = TRUE)
st_write(teyeo.out, "Data/Covariate_Data/Outcrop/New_Shapefiles/teyeo.shp")

# Teyep
r1 <- raster(ext = e, res = res)
newData.p <- rasterize(teyep.out, r1, getCover = TRUE, progress = "window")
plot(newData.p)
teyep.occs <- master.occs %>%
  dplyr::filter(bin_assignment == "teyep")
points(teyep.occs$lng, teyep.occs$lat, pch = 3)
writeRaster(newData.p, paste("Prepped_data/Covariate_Data/All_data/", res, 
                             "deg/Out_teyep_", res, ".asc", sep = ""), 
            pattern = "ascii", overwrite = TRUE)
st_write(teyep.out, "Data/Covariate_Data/Outcrop/New_Shapefiles/teyep.shp")

# Teyeq
r1 <- raster(ext = e, res = res)
newData.q <- rasterize(teyeq.out, r1, getCover = TRUE, progress = "window")
plot(newData.q)
teyeq.occs <- master.occs %>%
  dplyr::filter(bin_assignment == "teyeq")
points(teyeq.occs$lng, teyeq.occs$lat, pch = 3)
writeRaster(newData.q, paste("Prepped_data/Covariate_Data/All_data/", res, 
                             "deg/Out_teyeq_", res, ".asc", sep = ""), 
            pattern = "ascii", overwrite = TRUE)
st_write(teyeq.out, "Data/Covariate_Data/Outcrop/New_Shapefiles/teyeq.shp")

################################################################################
# 5. ADDING TO RELEVANT SP DATA FOR MODELLING
################################################################################

res <-  0.5
# Set extent
e <- extent(-155, -72, 22.5, 73)
# Set max limit value
max_val <- 40
max_val_on <- TRUE
bin.type <- "scotese"
target <- "Tyrannosauridae"
bin <- "teyeq" # Change this if you want a single season model
form_cells <- "Y" # Change this to "Y" for single season models with cells grouped 

# Multi-season
y.data <- readRDS(file = paste("Prepped_data/spOccupancy/Multi_season/", res, "/", target, 
                     "_multi_", res, ".rds", sep = ""))

# Single-season
y.data <- readRDS(file = paste("Prepped_data/spOccupancy/Single_season/", res, "/", bin, "/", 
                               target, "_single_", res, ".rds", sep = ""))

# Single-season Form Cells
y.data <- readRDS(file = paste("Prepped_data/spOccupancy/Single_season/", res, "/", bin, "/", 
                               target, "_single_", res, ".formcells.rds", sep = ""))

if(form_cells == "Y"){
  site_IDs <- as.numeric(gsub("([0-9]+).*$", "\\1", rownames(y.data$y)))
}else{
  site_IDs <- as.numeric(rownames(y.data$y))
}

siteCoords <- siteCoordsFun(res, e, site_IDs)

# Modern det
wc <- list.files(paste("Prepped_data/Covariate_Data/All_data/", 
                       res, "deg/", sep = ""), 
                 pattern=".asc")
stacked <- raster::stack(paste("Prepped_data/Covariate_Data/All_data/", 
                               res, "deg/", wc, 
                               sep =""))
covs <- get_cov(siteCoords, stacked, colls = FALSE)

names(covs) <- gsub(".[[:digit:]]", "", names(covs))

# Outcrop
outcrop <- data.frame(teyeq = covs$Out_teyeq, 
                      teyep = covs$Out_teyep,
                      teyeo = covs$Out_teyeo,
                      teyen = covs$Out_teyen
)
if(is.na(bin) == F){
  outcrop <- as.vector(outcrop[bin])
  outcrop <- outcrop[[1]]
}

y.data$det.covs$outcrop <- outcrop

if(is.na(bin) == T){
  saveRDS(y.data,
          file = paste("Prepped_data/spOccupancy/Multi_season/", res, "/", target, 
                       "_multi_", res, ".rds", sep = ""))
}else{
  if(form_cells == "Y"){
    saveRDS(y.data, 
            file = paste("Prepped_data/spOccupancy/Single_season/", res, "/", 
                         bin, "/", target, "_single_", res, ".formcells.rds", sep = ""))
  }else{
    saveRDS(y.data, 
            file = paste("Prepped_data/spOccupancy/Single_season/", res, "/", 
                         bin, "/", target, "_single_", res, ".rds", sep = ""))
  }
}