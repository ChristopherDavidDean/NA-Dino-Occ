#========================================= CLEANING DINOSAUR OCCURRENCE DATA  ====================================================#
#                                                                                                                                 #
#      AUTHOR: CHRISTOPHER D. DEAN                                                                                                #
#                                                                                                                                 #    
#=================================================================================================================================#

#========================================== REQUIRED PACKAGES & SET-UP =================================================

library(dplyr)

# Set your working directory
setwd("C:/Users/deancd/Documents/RESEARCH/PROJECTS/DINO_RANGE/NA-Dino-Occ/") 

# Load in Occurrences and Formations
CampanianOccs <- read.csv("Data/Occurrences/Broad_Colls_Data_Campanian_edited.csv") 
MaastrichtianOccs <- read.csv("Data/Occurrences/Broad_Colls_Data_Maastrichtian.csv")
Formations <- read.csv("Data/Occurrences/Formations.csv")

# Set Bins
bins <- c(84.3, 81, 76.3, 72.1, 69.9, 66)
bin_names <- c("Early Campanian", "Middle Campanian", "Late Campanian", "Early Maastrichtian", "Late Maastrichtian")
bins <- data.frame(bin = bin_names, # Combines bin data to make dataframe of minimum, maximum and mid point of each new bin
                      bottom = as.numeric(bins[1:(length(bins)-1)]), 
                      top = as.numeric(bins[2:(length(bins))]))
bins$mid <- (bins$bottom + bins$top) / 2
# Campanian
test1 <- CampanianOccs %>%  
  mutate(early_interval_PBDB = case_when(grepl("Early Campanian", early_interval) ~ bins[1,2],
                           grepl("Middle Campanian", early_interval) ~ bins[2,2],
                           grepl("Late Campanian", early_interval) ~ bins[3,2],
                           grepl("Late Santonian", early_interval) ~ 84.5,
                           grepl("Late Coniacian", early_interval) ~ 88))

test2 <- test1 %>%  
  mutate(late_interval_PBDB = case_when(grepl("Early Campanian", late_interval) ~ bins[1,3],
                           grepl("Middle Campanian", late_interval) ~ bins[2,3],
                           grepl("Late Campanian", late_interval) ~ bins[3,3],
                           grepl("Early Maastrichtian", late_interval) ~ bins[4,3],
                           grepl("Late Maastrichtian", late_interval) ~ bins[5,3]))

for (r in 1:nrow(test2)){
  if (!is.na(test2$early_interval_PBDB[[r]]) &&
      is.na(test2$late_interval_PBDB[[r]])){
    test2$late_interval_PBDB[[r]] <- bins[(grep(test2$early_interval_PBDB[[r]], bins$bottom)),3]
  }
}
                          
sum(is.na(test2$early_interval_PBDB))
sum(is.na(test2$late_interval_PBDB))

for (n in 1:nrow(test2)){
  for(f in 1:nrow(Formations)){
    if(test2$formation[[n]] == Formations$Formation[[f]]){
      test2$early_interval_updated[[n]] <- Formations$max_age[[f]]
      test2$late_interval_updated[[n]] <- Formations$min_age[[f]]
    }
  }
}

sum(is.na(test2$early_interval_updated))
sum(is.na(test2$late_interval_updated))


Form_list <- split(test2, test2$formation)

Occ_List <- list()# make an empty list of the occurrencess in each bin

for(n in 1:nrow(test2)){
  for(b in 1:nrow(bins)){
    if(is.na(test2$late_interval_PBDB[n]) | is.na(test2$early_interval_PBDB[n])){
      next
    }
    if (test2$late_interval_PBDB[n] > bins[b,2] |
        bins[b,3] > test2$early_interval_PBDB[n]) { # If the occurrence DOES NOT sit in this bin 
      next # Just print, don't do anything else.
    }
    else { # Otherwise (i.e. if a formation DOES sit in/cross this bin in any way)
      if (test2$early_interval_PBDB[n] <= bins[b,2] && 
          test2$late_interval_PBDB[n] >= bins[b,3]){# If occurrence sits within boundaries
        test2$PBDB_bin[n] <- as.character(bins[b,1])
        next
      }
      if (test2$early_interval_PBDB[n] >= bins[b,3] && 
          test2$late_interval_PBDB[n] <= bins[b,3] && 
          test2$early_interval_PBDB[n] <= bins[b,2]){ #If formation crosses upper age limit (minimum age)
        x <- as.numeric(test2$early_interval_PBDB[n] - as.numeric(bins[b,3]))
        y <- as.numeric(as.numeric(bins[b,3]) - test2$late_interval_PBDB[n])
        if (x > y){
          test2$PBDB_bin[n] <- as.character(bins[b,1]) # Add occurrences from that formation to occurrence list
          next
        }
      }
      if (test2$early_interval_PBDB[n] >= bins[b,2] && 
          test2$late_interval_PBDB[n] <= bins[b,2] &&
          test2$late_interval_PBDB[n] >= bins[b,3]){ #If formation crosses lower age limit (maximum age)
        x <- as.numeric(test2$early_interval_PBDB[n] - as.numeric(bins[b,2]))
        y <- as.numeric(as.numeric(bins[b,2]) - test2$late_interval_PBDB[n])
        if (y > x){
          test2$PBDB_bin[n] <- as.character(bins[b,1]) # Add occurrences from that formation to occurrence list
          next
        }
      }
      if (test2$early_interval_PBDB[n] > bins[b,2] &&
          test2$late_interval_PBDB[n] < bins[b,3]){
        next
      }
    }
  }
}


for(n in 1:nrow(test2)){
  for(b in 1:nrow(bins)){
    if(is.na(test2$late_interval_updated[n]) | is.na(test2$early_interval_updated[n])){
      next
    }
    if (test2$late_interval_updated[n] > bins[b,2] |
        bins[b,3] > test2$early_interval_updated[n]) { # If the occurrence DOES NOT sit in this bin 
      next # Just print, don't do anything else.
    }
    else { # Otherwise (i.e. if a formation DOES sit in/cross this bin in any way)
      if (test2$early_interval_updated[n] <= bins[b,2] && 
          test2$late_interval_updated[n] >= bins[b,3]){# If occurrence sits within boundaries
        test2$updated_bin[n] <- as.character(bins[b,1])
        test2$updated_fit[n] <- "Single Stage"
        next
      }
      if (test2$early_interval_updated[n] >= bins[b,3] && 
          test2$late_interval_updated[n] <= bins[b,3] && 
          test2$early_interval_updated[n] <= bins[b,2]){ #If formation crosses upper age limit (minimum age)
        x <- as.numeric(test2$early_interval_updated[n] - as.numeric(bins[b,3]))
        y <- as.numeric(as.numeric(bins[b,3]) - test2$late_interval_updated[n])
        if (x > y){
          test2$updated_bin[n] <- as.character(bins[b,1]) # Add occurrences from that formation to occurrence list
          test2$updated_fit[n] <- "Top Boundary Crosser"
          next
        }
      }
      if (test2$early_interval_updated[n] >= bins[b,2] && 
          test2$late_interval_updated[n] <= bins[b,2] &&
          test2$late_interval_updated[n] >= bins[b,3]){ #If formation crosses lower age limit (maximum age)
        x <- as.numeric(test2$early_interval_updated[n] - as.numeric(bins[b,2]))
        y <- as.numeric(as.numeric(bins[b,2]) - test2$late_interval_updated[n])
        if (y > x){
          test2$updated_bin[n] <- as.character(bins[b,1]) # Add occurrences from that formation to occurrence list
          test2$updated_fit[n] <- "Bottom Boundary Crosser"
          next
        }
      }
      if (test2$early_interval_updated[n] > bins[b,2] &&
          test2$late_interval_updated[n] < bins[b,3]){
        test2$updated_fit[n] <- "Top and Bottom Boundary Crosser"
        next
      }
    }
  }
}

PBDB_data <- test2 %>%
  dplyr::select(accepted_name, collection_no, formation, early_interval, late_interval, PBDB_bin)
colnames(PBDB_data)[colnames(PBDB_data)=="PBDB_bin"] <- "updated_bin"


Data_for_revision <- test2 %>%
  dplyr::select(accepted_name, collection_no, formation, early_interval_PBDB, late_interval_PBDB, early_interval_updated, late_interval_updated, early_interval, late_interval, PBDB_bin, updated_bin)

write.csv(test2, file = "Camp_data_for_revision_Boundaries_added.csv")

colnames(data)[colnames(data)=="old_name"] <- "new_name"
Updated_data <- test2 %>%
  dplyr::select(accepted_name, collection_no, formation, early_interval, late_interval, updated_bin)

test3 <- anti_join(PBDB_data, Updated_data, by=NULL)
test3<-test3[!(test3$updated_bin==81),]

test4 <- test2[!(test2$PBDB_bin %in% test2$updated_bin),]

# STEPS
# 1. Set bins
# 2. Apply binning technique using midpoints of PBDB data
# 3. Apply binning technique using midpoints from Formations data
# 4. Find occurrences that don't match up, investigate further.


# Maastrichtian
test1 <- MaastrichtianOccs %>%  
  mutate(early_interval_PBDB = case_when(grepl("Early Campanian", early_interval) ~ bins[1,2],
                                         grepl("Middle Campanian", early_interval) ~ bins[2,2],
                                         grepl("Late Campanian", early_interval) ~ bins[3,2],
                                         grepl("Early Maastrichtian", early_interval) ~ bins[4,2],
                                         grepl("Late Maastrichtian", early_interval) ~ bins[5,2]))

test2 <- test1 %>%  
  mutate(late_interval_PBDB = case_when(grepl("Early Maastrichtian", late_interval) ~ bins[4,3],
                                        grepl("Late Maastrichtian", late_interval) ~ bins[5,3]))


for (r in 1:nrow(test2)){
  if (!is.na(test2$early_interval_PBDB[[r]]) &&
      is.na(test2$late_interval_PBDB[[r]])){
    test2$late_interval_PBDB[[r]] <- bins[(grep(test2$early_interval_PBDB[[r]], bins$bottom)),3]
  }
}

sum(is.na(test2$early_interval_PBDB))
sum(is.na(test2$late_interval_PBDB))

for (n in 1:nrow(test2)){
  for(f in 1:nrow(Formations)){
    if(test2$formation[[n]] == Formations$Formation[[f]]){
      test2$early_interval_updated[[n]] <- Formations$max_age[[f]]
      test2$late_interval_updated[[n]] <- Formations$min_age[[f]]
    }
  }
}

sum(is.na(test2$early_interval_updated))
sum(is.na(test2$late_interval_updated))


Form_list <- split(test2, test2$formation)

Occ_List <- list()# make an empty list of the occurrencess in each bin

test2$PBDB_bin <- NA

for(n in 1:nrow(test2)){
  for(b in 1:nrow(bins)){
    if(is.na(test2$late_interval_PBDB[n]) | is.na(test2$early_interval_PBDB[n])){
      next
    }
    if (test2$late_interval_PBDB[n] > bins[b,2] |
        bins[b,3] > test2$early_interval_PBDB[n]) { # If the occurrence DOES NOT sit in this bin 
      next # Just print, don't do anything else.
    }
    else { # Otherwise (i.e. if a formation DOES sit in/cross this bin in any way)
      if (test2$early_interval_PBDB[n] <= bins[b,2] && 
          test2$late_interval_PBDB[n] >= bins[b,3]){# If occurrence sits within boundaries
        test2$PBDB_bin[n] <- as.character(bins[b,1])
        next
      }
      if (test2$early_interval_PBDB[n] >= bins[b,3] && 
          test2$late_interval_PBDB[n] <= bins[b,3] && 
          test2$early_interval_PBDB[n] <= bins[b,2]){ #If formation crosses upper age limit (minimum age)
        x <- as.numeric(test2$early_interval_PBDB[n] - as.numeric(bins[b,3]))
        y <- as.numeric(as.numeric(bins[b,3]) - test2$late_interval_PBDB[n])
        if (x > y){
          test2$PBDB_bin[n] <- as.character(bins[b,1]) # Add occurrences from that formation to occurrence list
          next
        }
      }
      if (test2$early_interval_PBDB[n] >= bins[b,2] && 
          test2$late_interval_PBDB[n] <= bins[b,2] &&
          test2$late_interval_PBDB[n] >= bins[b,3]){ #If formation crosses lower age limit (maximum age)
        x <- as.numeric(test2$early_interval_PBDB[n] - as.numeric(bins[b,2]))
        y <- as.numeric(as.numeric(bins[b,2]) - test2$late_interval_PBDB[n])
        if (y > x){
          test2$PBDB_bin[n] <- as.character(bins[b,1]) # Add occurrences from that formation to occurrence list
          next
        }
      }
      if (test2$early_interval_PBDB[n] > bins[b,2] &&
          test2$late_interval_PBDB[n] < bins[b,3]){
        next
      }
    }
  }
}
test2$updated_bin <- NA
test2$Updated_Fit <- NA

for(n in 1:nrow(test2)){
  for(b in 1:nrow(bins)){
    if(is.na(test2$late_interval_updated[n]) | is.na(test2$early_interval_updated[n])){
      next
    }
    if (test2$late_interval_updated[n] > bins[b,2] |
        bins[b,3] > test2$early_interval_updated[n]) { # If the occurrence DOES NOT sit in this bin 
      next # Just print, don't do anything else.
    }
    else { # Otherwise (i.e. if a formation DOES sit in/cross this bin in any way)
      if (test2$early_interval_updated[n] <= bins[b,2] && 
          test2$late_interval_updated[n] >= bins[b,3]){# If occurrence sits within boundaries
        test2$updated_bin[n] <- as.character(bins[b,1])
        test2$Updated_Fit[n] <- "Single Stage"
        next
      }
      if (test2$early_interval_updated[n] >= bins[b,3] && 
          test2$late_interval_updated[n] <= bins[b,3] && 
          test2$early_interval_updated[n] <= bins[b,2]){ #If formation crosses upper age limit (minimum age)
        x <- as.numeric(test2$early_interval_updated[n] - as.numeric(bins[b,3]))
        y <- as.numeric(as.numeric(bins[b,3]) - test2$late_interval_updated[n])
        if (x > y){
          test2$updated_bin[n] <- as.character(bins[b,1])
          test2$Updated_Fit[n] <- "Top Boundary Crosser" # Add occurrences from that formation to occurrence list
          next
        }
      }
      if (test2$early_interval_updated[n] >= bins[b,2] && 
          test2$late_interval_updated[n] <= bins[b,2] &&
          test2$late_interval_updated[n] >= bins[b,3]){ #If formation crosses lower age limit (maximum age)
        x <- as.numeric(test2$early_interval_updated[n] - as.numeric(bins[b,2]))
        y <- as.numeric(as.numeric(bins[b,2]) - test2$late_interval_updated[n])
        if (y > x){
          test2$updated_bin[n] <- as.character(bins[b,1])
          test2$Updated_Fit[n] <- "Bottom Boundary Crosser" # Add occurrences from that formation to occurrence list
          next
        }
      }
      if (test2$early_interval_updated[n] > bins[b,2] &&
          test2$late_interval_updated[n] < bins[b,3]){
        test2$Updated_Fit[n] <- "Top and Bottom Boundary Crosser"
        next
      }
    }
  }
}

PBDB_data <- test2 %>%
  dplyr::select(accepted_name, collection_no, formation, early_interval, late_interval, PBDB_bin)
colnames(PBDB_data)[colnames(PBDB_data)=="PBDB_bin"] <- "updated_bin"


Data_for_revision <- test2 %>%
  dplyr::select(accepted_name, collection_no, formation, early_interval_PBDB, late_interval_PBDB, early_interval_updated, late_interval_updated, early_interval, late_interval, PBDB_bin, updated_bin)

write.csv(test2, file = "Maas_data_for_revision.csv")

