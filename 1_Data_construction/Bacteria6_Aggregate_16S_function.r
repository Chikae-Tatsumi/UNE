library(tidyverse)
library (dplyr)
library(ggplot2)

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv("experimental_design.csv",header=T)
DESIGN <- na.omit(DESIGN)
setwd("~/R/Analysis/2_UNE/16S/function")

function.table <- read.csv(file="rarefied_16S_function.csv",row.names = 1,header=T)
ASV <- function.table [,1:(ncol(function.table)-19)] 
fgs <- function.table [,(ncol(function.table)-12):ncol(function.table)] 
percent <- ASV / mean(colSums(ASV)) *100

name <- colnames(fgs) 
colnames(fgs) <- c("Cellulolytic bacteria", "Assimiratory nitrite reducing bacteria", "Dissimiratory nitrite reducing bacteria", 
"Assimiratory nitrate reducing bacteria", "N fixing bacteria", "Dissimiratory nitrate reducing bacteria",
"Nitrifying bacteria", "Denitrifying bacteria", "Chitinolytic bacteria", "Lignolytic bacteria",              
"Methanotroph bacteria", "Copiotroph bacteria", "Oligotroph bacteria")

Result <- NULL
for (i in 1:ncol(fgs)){
aggregated <- aggregate(percent, by=list(fgs[,i]),FUN = sum,na.rm=F) 
rownames(aggregated) <- aggregated [,1]
aggregated <- aggregated [,-1]
aggregated.t <- t(aggregated)
Result <- cbind (Result, aggregated.t)
Result <- Result [,-ncol(Result)]}

# Save
colnames(Result) <- name
write.csv (Result, "aggregated.function.table.csv")