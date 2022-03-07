library(tidyverse)
library (dplyr)
library(lme4)
library(lmerTest)
library(vegan)

setwd("~/R/Analysis/2_UNE/ITS/FungalTrait")

fungaltrait.table <- read.csv(file="fungaltrait.table.csv",row.names = 1)
ASV <- fungaltrait.table [,1:(ncol(fungaltrait.table)-31)] 
guild <- fungaltrait.table [,(ncol(fungaltrait.table)-24):ncol(fungaltrait.table)] 
percent <- ASV / mean(colSums(ASV)) *100

percent.table <- cbind (percent, guild)
ECM.table<- percent.table[grep("ectomycorrhizal", percent.table$primary_lifestyle),]
write.csv(ECM.table, "ASV.table.ECMonly.FungalTrait.csv")