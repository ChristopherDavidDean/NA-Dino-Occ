library(ggplot2)
library(dplyr)

# Setup
fam <- "Hadrosauridae"
res <- 0.5
sub <- 5

S.1 <- read.csv(paste("Results/stage/S.1/", res, "/Model_results/", fam, ".combined.results.SS", sub, ".csv", sep = ""))
S.2 <- read.csv(paste("Results/stage/S.2/", res, "/Model_results/", fam, ".combined.results.SS", sub, ".csv", sep = ""))
other.res.S.1 <- read.csv(paste("Results/stage/S.1/Targeted_res_stats.csv"))
other.res.S.2 <- read.csv(paste("Results/stage/S.2/Targeted_res_stats.csv"))

transform.results <- function(data, fam){
  temp <- data.frame(Bin = data[1,6],
                     null.occ = data[1,2],
                     null.occ_CI_5 = data[1,3],
                     null.occ_CI_95 = data[1,4],
                     null.det = data[2,2],
                     null.det_CI_5 = data[2,3],
                     null.det_CI_95 = data[2,4],
                     det = data[4,2],
                     det_SD = data[4,5],
                     PAO = data[3,2],
                     PAO_CI_5 = data[3,3],
                     PAO_CI_95 = data[3,4]
  )
  temp_name <- paste(deparse(substitute(data)),".", deparse(fam), sep = "")
  assign(temp_name, temp, envir = .GlobalEnv)
}

results <- rbind(transform.results(data = S.1, fam), transform.results(data = S.2, fam))

test <- filter(other.res.S.1, X == paste(fam,".",res, sep = ""))
results$naive[1] <- (test[1,3]+test[1,11])/(test[1,2]+test[1,10])
test <- filter(other.res.S.2, X == paste(fam,".",res, sep = ""))
results$naive[2] <- (test[1,3]+test[1,11])/(test[1,2]+test[1,10])

data(stages)
stages <- stages[80:81,]
results$bin[2]<- stages[1,7]
results$bin[1]<- stages[2,7]

data <- results

#==========
detplus <- data$det + data$det_SD
detminus <- data$det - data$det_SD

par(mar = c(4.1, 4.1, 4, 2.1))

data(stages)

# Null model
tsplot(stages[80:81,], boxes=c("short","system"), ylab = "Occupancy and Detection Probability", # Creates plot using data from DivDyn package. 
       ylim=c(0,1), 
       shading="stage", boxes.col=c("col","systemCol"), labels.args=list(cex=0.75))  

polygon(x = c(data$bin, rev(data$bin)), y = c(data$null.occ_CI_5, rev(data$null.occ_CI_95)), col = adjustcolor("blue", alpha.f = 0.40), border = NA)
lines(data$bin, data$null.occ, type = "o", pch = 21, col = adjustcolor("blue", alpha.f = 0.4), bg = "blue", lwd = 1)

polygon(x = c(data$bin, rev(data$bin)), y = c(data$null.det_CI_5, rev(data$null.det_CI_95)), col = adjustcolor("red", alpha.f = 0.40), border = NA)
lines(data$bin, data$null.det, type = "o", pch = 21, col = adjustcolor("red", alpha.f = 0.4), bg = "red", lwd = 1)

lines(data$bin, data$naive, type = "o", pch = 21, col = "black", bg = "grey", lwd = 1, lty = 1)

# Model results
tsplot(stages[80:81,], boxes=c("short","system"), ylab = "Occupancy and Detection Probability", # Creates plot using data from DivDyn package. 
       ylim=c(0,1), 
       shading="stage", boxes.col=c("col","systemCol"), labels.args=list(cex=0.75))  

polygon(x = c(data$bin, rev(data$bin)), y = c(data$PAO_CI_5, rev(data$PAO_CI_95)), col = adjustcolor("blue", alpha.f = 0.40), border = NA)
lines(data$bin, data$PAO, type = "o", pch = 21, col = adjustcolor("blue", alpha.f = 0.4), bg = "blue", lwd = 1)

polygon(x = c(data$bin, rev(data$bin)), y = c(detminus, rev(detplus)), col = adjustcolor("red", alpha.f = 0.40), border = NA)
lines(data$bin, data$det, type = "o", pch = 21, col = adjustcolor("red", alpha.f = 0.4), bg = "red", lwd = 1)

lines(data$bin, data$naive, type = "o", pch = 21, col = "black", bg = "grey", lwd = 1, lty = 1)



#======================
res <- 0.5
sub <- 20

# DETECTION
Hadro <- read.csv(paste("Results/stage/S.1/", res, "/Model_results/Hadrosauridae.model.results.SS", sub, ".csv", sep = ""))
Tyran <- read.csv(paste("Results/stage/S.1/", res, "/Model_results/Tyrannosauridae.model.results.SS", sub, ".csv", sep = ""))
Cerat <- read.csv(paste("Results/stage/S.1/", res, "/Model_results/Ceratopsidae.model.results.SS", sub, ".csv", sep = ""))

Hadro <- Hadro[order(Hadro$name),]
Tyran <- Tyran[order(Tyran$name),]
Cerat <- Cerat[order(Cerat$name),]

Hadro <- Hadro[1:7,]
Tyran <- Tyran[1:7,]
Cerat <- Cerat[1:7,]

barplot(Hadro$value, names.arg = Hadro$name, ylim = c(-1,1))
barplot(Cerat$value, names.arg = Cerat$name, ylim = c(-1,1))
barplot(Tyran$value, names.arg = Tyran$name, ylim = c(-1,1))

Hadro <- read.csv(paste("Results/stage/S.2/", res, "/Model_results/Hadrosauridae.model.results.SS", sub, ".csv", sep = ""))
Tyran <- read.csv(paste("Results/stage/S.2/", res, "/Model_results/Tyrannosauridae.model.results.SS", sub, ".csv", sep = ""))
Cerat <- read.csv(paste("Results/stage/S.2/", res, "/Model_results/Ceratopsidae.model.results.SS", sub, ".csv", sep = ""))

Hadro <- Hadro[order(Hadro$name),]
Tyran <- Tyran[order(Tyran$name),]
Cerat <- Cerat[order(Cerat$name),]

Hadro <- Hadro[1:7,]
Tyran <- Tyran[1:7,]
Cerat <- Cerat[1:7,]

barplot(Hadro$value, names.arg = Hadro$name, ylim = c(-1,1))
barplot(Cerat$value, names.arg = Cerat$name, ylim = c(-1,1))
barplot(Tyran$value, names.arg = Tyran$name, ylim = c(-1,1))


# OCCUPANCY
Hadro <- read.csv(paste("Results/stage/S.1/", res, "/Model_results/Hadrosauridae.model.results.SS", sub, ".csv", sep = ""))
Tyran <- read.csv(paste("Results/stage/S.1/", res, "/Model_results/Tyrannosauridae.model.results.SS", sub, ".csv", sep = ""))
Cerat <- read.csv(paste("Results/stage/S.1/", res, "/Model_results/Ceratopsidae.model.results.SS", sub, ".csv", sep = ""))

Hadro <- Hadro[order(Hadro$name),]
Tyran <- Tyran[order(Tyran$name),]
Cerat <- Cerat[order(Cerat$name),]

Hadro <- Hadro[8:12,]
Tyran <- Tyran[8:12,]
Cerat <- Cerat[8:12,]

barplot(Hadro$value, names.arg = Hadro$name)
barplot(Cerat$value, names.arg = Cerat$name)
barplot(Tyran$value, names.arg = Tyran$name)

Hadro <- read.csv(paste("Results/stage/S.2/", res, "/Model_results/Hadrosauridae.model.results.SS", sub, ".csv", sep = ""))
Tyran <- read.csv(paste("Results/stage/S.2/", res, "/Model_results/Tyrannosauridae.model.results.SS", sub, ".csv", sep = ""))
Cerat <- read.csv(paste("Results/stage/S.2/", res, "/Model_results/Ceratopsidae.model.results.SS", sub, ".csv", sep = ""))

Hadro <- Hadro[order(Hadro$name),]
Tyran <- Tyran[order(Tyran$name),]
Cerat <- Cerat[order(Cerat$name),]

Hadro <- Hadro[8:12,]
Tyran <- Tyran[8:12,]
Cerat <- Cerat[8:12,]

barplot(Hadro$value, names.arg = Hadro$name)
barplot(Cerat$value, names.arg = Cerat$name)
barplot(Tyran$value, names.arg = Tyran$name)


#===============================================================================

# Setup
fam <- "Tyrannosauridae"
res <- 0.1
sub <- 5

S.1 <- read.csv(paste("Results/scotese/SC.1/", res, "/Model_results/", fam, ".combined.results.SS", sub, ".csv", sep = ""))
S.2 <- read.csv(paste("Results/scotese/SC.2/", res, "/Model_results/", fam, ".combined.results.SS", sub, ".csv", sep = ""))
S.3 <- read.csv(paste("Results/scotese/SC.3/", res, "/Model_results/", fam, ".combined.results.SS", sub, ".csv", sep = ""))
S.4 <- read.csv(paste("Results/scotese/SC.4/", res, "/Model_results/", fam, ".combined.results.SS", sub, ".csv", sep = ""))

other.res.S.1 <- read.csv(paste("Results/scotese/SC.1/Targeted_res_stats.csv"))
other.res.S.2 <- read.csv(paste("Results/scotese/SC.2/Targeted_res_stats.csv"))
other.res.S.3 <- read.csv(paste("Results/scotese/SC.3/Targeted_res_stats.csv"))
other.res.S.4 <- read.csv(paste("Results/scotese/SC.4/Targeted_res_stats.csv"))


transform.results <- function(data, fam){
  temp <- data.frame(Bin = data[1,6],
                     null.occ = data[1,2],
                     null.occ_CI_5 = data[1,3],
                     null.occ_CI_95 = data[1,4],
                     null.det = data[2,2],
                     null.det_CI_5 = data[2,3],
                     null.det_CI_95 = data[2,4],
                     det = data[4,2],
                     det_SD = data[4,5],
                     PAO = data[3,2],
                     PAO_CI_5 = data[3,3],
                     PAO_CI_95 = data[3,4]
  )
  temp_name <- paste(deparse(substitute(data)),".", deparse(fam), sep = "")
  assign(temp_name, temp, envir = .GlobalEnv)
}

results <- rbind(transform.results(data = S.1, fam), transform.results(data = S.2, fam), 
                 transform.results(data = S.3, fam), transform.results(data = S.4, fam))

test <- filter(other.res.S.1, X == paste(fam,".",res, sep = ""))
results$naive[1] <- (test[1,3]+test[1,11])/(test[1,2]+test[1,10])
test <- filter(other.res.S.2, X == paste(fam,".",res, sep = ""))
results$naive[2] <- (test[1,3]+test[1,11])/(test[1,2]+test[1,10])
test <- filter(other.res.S.3, X == paste(fam,".",res, sep = ""))
results$naive[3] <- (test[1,3]+test[1,11])/(test[1,2]+test[1,10])
test <- filter(other.res.S.4, X == paste(fam,".",res, sep = ""))
results$naive[4] <- (test[1,3]+test[1,11])/(test[1,2]+test[1,10])

bins <- read.csv("Data/Occurrences/scotesebins.csv")
bins$bin <- c("SC.1","SC.2","SC.3","SC.4")
results$bin[1]<- bins[1,6]
results$bin[2]<- bins[2,6]
results$bin[3]<- bins[3,6]
results$bin[4]<- bins[4,6]

data <- results

#==========
detplus <- data$det + data$det_SD
detminus <- data$det - data$det_SD

par(mar = c(4.1, 4.1, 4, 2.1))

data(stages)

# Null model
tsplot(stages[80:81,], boxes=c("short","system"), ylab = "Occupancy and Detection Probability", # Creates plot using data from DivDyn package. 
       ylim=c(0,1), 
       shading="stage", boxes.col=c("col","systemCol"), labels.args=list(cex=0.75))  

polygon(x = c(data$bin, rev(data$bin)), y = c(data$null.occ_CI_5, rev(data$null.occ_CI_95)), col = adjustcolor("blue", alpha.f = 0.40), border = NA)
lines(data$bin, data$null.occ, type = "o", pch = 21, col = adjustcolor("blue", alpha.f = 0.4), bg = "blue", lwd = 1)

polygon(x = c(data$bin, rev(data$bin)), y = c(data$null.det_CI_5, rev(data$null.det_CI_95)), col = adjustcolor("red", alpha.f = 0.40), border = NA)
lines(data$bin, data$null.det, type = "o", pch = 21, col = adjustcolor("red", alpha.f = 0.4), bg = "red", lwd = 1)

lines(data$bin, data$naive, type = "o", pch = 21, col = "black", bg = "grey", lwd = 1, lty = 1)

# Model results
tsplot(stages[80:81,], boxes=c("short","system"), ylab = "Occupancy and Detection Probability", # Creates plot using data from DivDyn package. 
       ylim=c(0,1), 
       shading="stage", boxes.col=c("col","systemCol"), labels.args=list(cex=0.75))  

polygon(x = c(data$bin, rev(data$bin)), y = c(data$PAO_CI_5, rev(data$PAO_CI_95)), col = adjustcolor("blue", alpha.f = 0.40), border = NA)
lines(data$bin, data$PAO, type = "o", pch = 21, col = adjustcolor("blue", alpha.f = 0.4), bg = "blue", lwd = 1)

polygon(x = c(data$bin, rev(data$bin)), y = c(detminus, rev(detplus)), col = adjustcolor("red", alpha.f = 0.40), border = NA)
lines(data$bin, data$det, type = "o", pch = 21, col = adjustcolor("red", alpha.f = 0.4), bg = "red", lwd = 1)

lines(data$bin, data$naive, type = "o", pch = 21, col = "black", bg = "grey", lwd = 1, lty = 1)

