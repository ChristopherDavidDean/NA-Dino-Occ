library(ggplot2)
library(dplyr)

# Setup
data <- read.csv("Results/First_results.csv")
fam <- "Ceratopsidae"
fam <- "Hadrosauridae"
fam <- "Tyrannosauridae"

data <- data %>%
  filter(family == fam)

detplus <- data$det + data$det_SD
detminus <- data$det - data$det_SD

par(mar = c(4.1, 4.1, 4, 2.1))

# Null model
tsplot(stages, boxes=c("short","system"), ylab = "Occupancy and Detection Probability", # Creates plot using data from DivDyn package. 
       ylim=c(0,1), 
       shading="stage", boxes.col=c("col","systemCol"), labels.args=list(cex=0.75))  

polygon(x = c(data$bin, rev(data$bin)), y = c(data$null.occ_CI_5, rev(data$null.occ_CI_95)), col = adjustcolor("blue", alpha.f = 0.40), border = NA)
lines(data$bin, data$null.occ, type = "o", pch = 21, col = "black", bg = "grey", lwd = 1)

polygon(x = c(data$bin, rev(data$bin)), y = c(data$null.det_CI_5, rev(data$null.det_CI_95)), col = adjustcolor("red", alpha.f = 0.40), border = NA)
lines(data$bin, data$null.det, type = "o", pch = 21, col = "black", bg = "grey", lwd = 1)

lines(data$bin, data$naive, type = "o", pch = 21, col = "black", bg = "grey", lwd = 1, lty = 1)

# Model results
tsplot(stages, boxes=c("short","system"), ylab = "Occupancy and Detection Probability", # Creates plot using data from DivDyn package. 
       ylim=c(0,1), 
       shading="stage", boxes.col=c("col","systemCol"), labels.args=list(cex=0.75))  

polygon(x = c(data$bin, rev(data$bin)), y = c(data$PAO_CI_5, rev(data$PAO_CI_95)), col = adjustcolor("blue", alpha.f = 0.40), border = NA)
lines(data$bin, data$PAO, type = "o", pch = 21, col = adjustcolor("blue", alpha.f = 0.4), bg = "blue", lwd = 1)

polygon(x = c(data$bin, rev(data$bin)), y = c(detminus, rev(detplus)), col = adjustcolor("red", alpha.f = 0.40), border = NA)
lines(data$bin, data$det, type = "o", pch = 21, col = adjustcolor("red", alpha.f = 0.4), bg = "red", lwd = 1)

lines(data$bin, data$naive, type = "o", pch = 21, col = "black", bg = "grey", lwd = 1, lty = 1)
