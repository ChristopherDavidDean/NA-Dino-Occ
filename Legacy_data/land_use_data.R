install.packages("tmap")
library(tmap)
library(sp)
library(raster)
data(land) #get land cover data
lc <- land$cover #get land cover data
plot(lc) #plot land cover data
lc[lc == 20] <- NA #set ocean and water bodies as NA
#https://globalmaps.github.io/glcnmo.html
#Code	Class Name	                      Code  Class Name
#1	  Broadleaf Evergreen Forest	      11	  Cropland
#2	  Broadleaf Deciduous Forest	      12	  Paddy field
#3	  Needleleaf Evergreen Forest	      13	  Cropland / Other Vegetation Mosaic
#4	  Needleleaf Deciduous Forest	      14	  Mangrove
#5	  Mixed Forest	                    15	  Wetland
#6	  Tree Open	                        16	  Bare area,consolidated(gravel,rock)
#7	  Shrub	                            17	  Bare area,unconsolidated (sand)
#8	  Herbaceous	                      18	  Urban
#9	  Herbaceous with Sparse Tree/Shrub	19	  Snow / Ice
#10	  Sparse vegetation	                20    Water bodies
lc[lc == 10] <- 100#reclassify 
lc[lc == 16] <- 100#reclassify 
lc[lc == 17] <- 100#reclassify 
lc[lc <= 20] <- 0 #remove all values less than or equal to 20
plot(lc)

setwd("C://Users/deancd/Documents/RESEARCH/PROJECTS/OCC1/Occ_Datasets/")
hires <- raster("Land_Cover/NLCD_2001_Land_Cover_L48_20190424.ige")
test10 <- raster::select(hires)
plot(testd)
testd <- raster("Elevation_GRID/Elevation_GRID/NA_Elevation/data/NA_Elevation/na_elevation/w001001.adf")
testlc <- raster("Land_Cover/GlobalLandCover_tif/LCType.tif")
test <- raster("Precip/precip.mon.ltm.v501.nc")
test2 <- raster::select(testlc)
test3 <- select(test)
test4 <- select(testd)
test2[test2 == 0] <- NA
test2[test2 == 16] <- 100
test2[test2 <= 15] <- 0
plot(test2)
plot(test3)
plot(test4)
test5 <- terrain(test4, opt='slope', unit='degrees', neighbors=8)
plot(test5)
test6 <- aggregate(test2, fact=3)
plot(test6)
plot(testlc)
testlc$landcvi020l
plot(testlc[testlc == 19])
