library(maptools)
library(ggplot2)
library(ggrepel)

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv(file = "experimental_design.csv",header=T)
setwd("~/R/Analysis/2_UNE/16s/Tax4Fun2")
stress.16S <- read.csv (file = "pathway_prediction.csv",row.names=1,header=T,sep="\t",quote = "") 
data <- t(stress.16S[,1:(ncol(stress.16S)-3)])

# cor.test
Results <- NULL
for (i in 1:ncol(data)){
cor.DFB <- cor.test(-DESIGN$DFB, data[,i])
cor.DFE <- cor.test(-DESIGN$DFE, data[,i])
bind <- cbind(cor.DFB$estimate,cor.DFE$estimate)
Results <- rbind (Results,bind)}
rownames (Results) <- colnames(data)[1:ncol(data)]
colnames(Results) <- c("Urbanization","Fragmentation")
Results <- cbind(Results,stress.16S[,(ncol(stress.16S)-2):ncol(stress.16S)])

# Save
write.csv(Results,"correlation.withpathway.csv")
