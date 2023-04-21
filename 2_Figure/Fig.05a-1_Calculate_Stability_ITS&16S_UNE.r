library(vegan)
library(reshape2)
library(ggplot2)
library(robCompositions)
library(coda.base)
library(lme4)
library(lmerTest)

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <-  read.csv("experimental_design.csv",header=T)
metadata <- read.csv("metadata.csv",header=T)

# For ITS
setwd("~/R/Analysis/2_UNE/ITS")
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T)
ASV <- ASV.table [,1:(ncol(ASV.table)-7)]
ASV.t <- t(ASV)
ASV.t <- as.matrix(ASV.t)
setwd("effect_of_year")

# Get Aitchison distance 
ASV.dist<-dist(ASV.t+1, method="aitchison")
ASV.dist.matrix <- as.matrix(ASV.dist)

# Remove not corresponding data
ASV.dist.matrix.add <- rbind(ASV.dist.matrix[1:197,],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,ASV.dist.matrix[198:nrow(ASV.dist.matrix),]) # HF05 in 2021
ASV.dist.matrix.add <- rbind(ASV.dist.matrix.add[1:177,],NA,NA,ASV.dist.matrix.add[178:nrow(ASV.dist.matrix.add),]) # BH02 90m in 2021
ASV.dist.matrix.add <- rbind(ASV.dist.matrix.add[1:68,],NA,ASV.dist.matrix.add[69:nrow(ASV.dist.matrix.add),]) # HW07 90m A in 2018

ASV.dist.matrix.add <- cbind(ASV.dist.matrix.add[,1:197],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,ASV.dist.matrix.add[,198:ncol(ASV.dist.matrix.add)])
ASV.dist.matrix.add <- cbind(ASV.dist.matrix.add[,1:177],NA,NA,ASV.dist.matrix.add[,178:ncol(ASV.dist.matrix.add)])
ASV.dist.matrix.add <- cbind(ASV.dist.matrix.add[,1:68],NA,ASV.dist.matrix.add[,69:ncol(ASV.dist.matrix.add)])

# Get the Aitchison distance between years in each location
Result <- NULL
for (i in 1:80){
    result <- NULL
dist1819<- ASV.dist.matrix.add[i+80,i]
dist1821<- ASV.dist.matrix.add[i+160,i]
result <- cbind(dist1819,dist1821)
Result <- rbind(Result,result)}

Result <-data.frame(Result)

ag <- aggregate (metadata, by=list(DESIGN$UID),FUN = mean,na.rm=T)
DESIGN.ag <- DESIGN[80:159,] # Because 2019 has a full data (n=80)

data1 <- cbind(Result$dist1819, DESIGN.ag[,-3:-4], ag)
data2 <- cbind(Result$dist1821, DESIGN.ag[,-3:-4], ag)
colnames(data1)[1] <- "Dissimilarity"
colnames(data2)[1] <- "Dissimilarity"
data <- rbind(data1, data2)
data$Year <- "2019"
data$Year[81:160] <- "2021"
write.csv(data, "ITS_2018baseline_withAveragedMetadata.csv")

# For 16S
setwd("~/R/Analysis/2_UNE/16S")
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T)
ASV <- ASV.table [,1:(ncol(ASV.table)-6)]
ASV.t <- t(ASV)
ASV.t <- as.matrix(ASV.t)
setwd("effect_of_year")

# Get the Aitchison distance 
ASV.dist<-dist(ASV.t+1, method="aitchison")
ASV.dist.matrix <- as.matrix(ASV.dist)

# Remove not corresponding data
ASV.dist.matrix.add <- rbind(ASV.dist.matrix[1:197,],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,ASV.dist.matrix[198:nrow(ASV.dist.matrix),])
ASV.dist.matrix.add <- rbind(ASV.dist.matrix.add[1:177,],NA,NA,ASV.dist.matrix.add[178:nrow(ASV.dist.matrix.add),])
ASV.dist.matrix.add <- rbind(ASV.dist.matrix.add[1:68,],NA,ASV.dist.matrix.add[69:nrow(ASV.dist.matrix.add),])

ASV.dist.matrix.add <- cbind(ASV.dist.matrix.add[,1:197],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,ASV.dist.matrix.add[,198:ncol(ASV.dist.matrix.add)])
ASV.dist.matrix.add <- cbind(ASV.dist.matrix.add[,1:177],NA,NA,ASV.dist.matrix.add[,178:ncol(ASV.dist.matrix.add)])
ASV.dist.matrix.add <- cbind(ASV.dist.matrix.add[,1:68],NA,ASV.dist.matrix.add[,69:ncol(ASV.dist.matrix.add)])

# Get the Aitchison distance between years in each location
Result <- NULL
for (i in 1:80){
    result <- NULL
dist1819<- ASV.dist.matrix.add[i+80,i]
dist1821<- ASV.dist.matrix.add[i+160,i]
result <- cbind(dist1819,dist1821)
Result <- rbind(Result,result)}

Result <-data.frame(Result)

ag <- aggregate (metadata, by=list(DESIGN$UID),FUN = mean,na.rm=T)
DESIGN.ag <- DESIGN[80:159,]

data1 <- cbind(Result$dist1819, DESIGN.ag[,-3:-4], ag)
data2 <- cbind(Result$dist1821, DESIGN.ag[,-3:-4], ag)
colnames(data1)[1] <- "Dissimilarity"
colnames(data2)[1] <- "Dissimilarity"
data <- rbind(data1, data2)
data$Year <- "2019"
data$Year[81:160] <- "2021"
write.csv(data, "16S_2018baseline_withAveragedMetadata.csv")
