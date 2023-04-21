library(ecodist)
library(vegan)

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv(file = "experimental_design.csv",header=T)
DESIGN.matrix <- cbind(DESIGN$DFB, DESIGN$DFE)
DESIGN.dist <- dist(DESIGN.matrix)

setwd("~/R/Analysis/2_UNE/Autocorrelation")
distance.matrix <- read.csv("dd.distances.csv",header=T, row.names=1)
dist.dist <- as.dist(distance.matrix)

# For ITS
setwd("~/R/Analysis/2_UNE/ITS")
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T)
ASV <- ASV.table [,1:(ncol(ASV.table)-7)]
ASV.t <- t(ASV)
ASV.t <- as.matrix(ASV.t)

ASV.dist <- vegdist(ASV.t,method="bray")

# Run MRM analysis
MRM <- lm(as.numeric(ASV.dist)~ as.numeric(DESIGN.dist)+as.numeric(dist.dist))
anova <- anova(MRM)

# Save
write.csv(anova, "MRM_ITS_spatial.vs.experiment.csv")

# For 16S
setwd("~/R/Analysis/2_UNE/16S")
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T)
ASV <- ASV.table [,1:(ncol(ASV.table)-6)]
ASV.t <- t(ASV)
ASV.t <- as.matrix(ASV.t)

ASV.dist<-vegdist(ASV.t,method="bray")

# Run MRM analysis
MRM <- lm(as.numeric(ASV.dist)~ as.numeric(DESIGN.dist)+as.numeric(dist.dist))
anova <- anova(MRM)

# Save
write.csv(anova, "MRM_16S_spatial.vs.experiment.csv")

