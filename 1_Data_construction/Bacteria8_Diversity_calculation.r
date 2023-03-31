library(vegan)
library(ggplot2)
library(lme4)
library(lmerTest)

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv("experimental_design.csv",header=T)
setwd("~/R/Analysis/2_UNE/16S")
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T)

ASV <- ASV.table [,1:(ncol(ASV.table)-6)] 
ASV.t <- t(ASV)
shannon <- diversity(ASV.t, index="shannon",base=2)
simpson <- diversity(ASV.t, index="simpson")
invsimpson <- diversity(ASV.t, index="invsimpson")
fisher <- fisher.alpha(ASV.t)
data <- cbind(shannon, simpson, invsimpson,fisher)

# Save
write.csv (data, "diversity.csv")
