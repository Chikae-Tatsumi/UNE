library(vegan)
library(ggplot2)
library(tidyr)
library(tidyr)
library(tidyverse)
library(reshape2)

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv("experimental_design.csv",header=T)

# For 16S
setwd("~/R/Analysis/2_UNE/16S")
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T)
ASV <- ASV.table [,1:(ncol(ASV.table)-6)]
ASV.t <- t(ASV)
ASV.t <- as.matrix(ASV.t)

# Categolize plots to four location groups
DESIGN$Group4 <- NA
for (i in 1: nrow(DESIGN)){
if (DESIGN$Urban.DFE[i]=="Urban0"|DESIGN$Urban.DFE[i]=="Urban15") {DESIGN$Group4[i]<-"Urban.edge"} else 
if (DESIGN$Urban.DFE[i]=="Urban30"|DESIGN$Urban.DFE[i]=="Urban60"|DESIGN$Urban.DFE[i]=="Urban90") {DESIGN$Group4[i]<-"Urban.interior"} else 
if (DESIGN$Urban.DFE[i]=="Rural0"|DESIGN$Urban.DFE[i]=="Rural15") {DESIGN$Group4[i]<-"Rural.edge"} else 
if (DESIGN$Urban.DFE[i]=="Rural30"|DESIGN$Urban.DFE[i]=="Rural60"|DESIGN$Urban.DFE[i]=="Rural90") {DESIGN$Group4[i]<-"Rural.interior"}  
}

subset <- subset(ASV.t, DESIGN$Year=="Y2019") # Because 2019 has a full soil data
DESIGN.subset <- subset(DESIGN, DESIGN$Year=="Y2019")
permdisp2019.16S <- betadisper(vegdist(subset, method="bray"), DESIGN.subset$Group4)
dist.16S <- permdisp2019.16S$distances

# For ITS
setwd("~/R/Analysis/2_UNE/ITS")
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T)
ASV <- ASV.table [,1:(ncol(ASV.table)-7)]
ASV.t <- t(ASV)
ASV.t <- as.matrix(ASV.t)

# Categolize plots to four location groups
DESIGN$Group4 <- NA
for (i in 1: nrow(DESIGN)){
if (DESIGN$Urban.DFE[i]=="Urban0"|DESIGN$Urban.DFE[i]=="Urban15") {DESIGN$Group4[i]<-"Urban.edge"} else 
if (DESIGN$Urban.DFE[i]=="Urban30"|DESIGN$Urban.DFE[i]=="Urban60"|DESIGN$Urban.DFE[i]=="Urban90") {DESIGN$Group4[i]<-"Urban.interior"} else 
if (DESIGN$Urban.DFE[i]=="Rural0"|DESIGN$Urban.DFE[i]=="Rural15") {DESIGN$Group4[i]<-"Rural.edge"} else 
if (DESIGN$Urban.DFE[i]=="Rural30"|DESIGN$Urban.DFE[i]=="Rural60"|DESIGN$Urban.DFE[i]=="Rural90") {DESIGN$Group4[i]<-"Rural.interior"}  
}

subset <- subset(ASV.t, DESIGN$Year=="Y2019")
DESIGN.subset <- subset(DESIGN, DESIGN$Year=="Y2019")
permdisp2019.ITS <- betadisper(vegdist(subset, method="bray"), DESIGN.subset$Group4)
dist.ITS <- permdisp2019.ITS$distances

# For soil properties
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv("experimental_design.csv",header=T)
METADATA <- read.csv("metadata.csv", header=T)
Soil <- cbind(METADATA$Moist,METADATA$SOM,METADATA$Temp,METADATA$pH,METADATA$NH4, METADATA$NO3)

# Categolize plots to four location groups
DESIGN$Group4 <- NA
for (i in 1: nrow(DESIGN)){
if (DESIGN$Urban.DFE[i]=="Urban0"|DESIGN$Urban.DFE[i]=="Urban15") {DESIGN$Group4[i]<-"Urban Edge"} else 
if (DESIGN$Urban.DFE[i]=="Urban30"|DESIGN$Urban.DFE[i]=="Urban60"|DESIGN$Urban.DFE[i]=="Urban90") {DESIGN$Group4[i]<-"Urban Interior"} else 
if (DESIGN$Urban.DFE[i]=="Rural0"|DESIGN$Urban.DFE[i]=="Rural15") {DESIGN$Group4[i]<-"Rural Edge"} else 
if (DESIGN$Urban.DFE[i]=="Rural30"|DESIGN$Urban.DFE[i]=="Rural60"|DESIGN$Urban.DFE[i]=="Rural90") {DESIGN$Group4[i]<-"Rural Interior"}  
}

subset <- subset(Soil, DESIGN$Year=="Y2019")
DESIGN.subset <- subset(DESIGN, DESIGN$Year=="Y2019")
bind <- cbind(subset,DESIGN.subset)
bind <- na.omit(bind)
subset <- bind[,1:ncol(subset)]
DESIGN.subset <- bind[,(ncol(subset)+1):ncol(bind)]
permdisp2019.soil <- betadisper(vegdist(subset, method="bray"), DESIGN.subset$Group4)

dist.soil <- permdisp2019.soil$distances
dist.soil<- c(dist.soil[1:45],NA,dist.soil[46:79]) # Becasue one temperature data is missing
Group4 <- c(DESIGN.subset$Group4[1:45],NA,DESIGN.subset$Group4[46:79])

# cor.test
data <- cbind(dist.soil, dist.ITS, dist.16S, Group4)
data <- data.frame(data)
data$dist.soil <- as.numeric(data$dist.soil)
data$dist.16S <- as.numeric(data$dist.16S)
data$dist.ITS <- as.numeric(data$dist.ITS)

setwd("~/R/Analysis/2_UNE/Others")
cor <- cor.test (data$dist.soil, data$dist.ITS,method="p")
sink("Homogenization_ITS_cor.txt")
cor
sink()

cor <- cor.test (data$dist.soil, data$dist.16S,method="p")
sink("Homogenization_16S_cor.txt")
cor
sink()

# Format data
data.16S <- cbind(data$dist.soil,data$dist.16S, data$Group4,"Bacteria")
colnames(data.16S) <- c("soil","microbe", "Group4","Type")
data.ITS <- cbind(data$dist.soil,data$dist.ITS,data$Group4,"Fungi")
colnames(data.ITS) <- c("soil","microbe", "Group4","Type")
vdata <- rbind(data.ITS, data.16S)

vdata  <- data.frame(vdata )
vdata  <- na.omit(vdata)

vdata$soil <- as.numeric(vdata$soil)
vdata$microbe <- as.numeric(vdata$microbe)

# Visualize
vdata$Group4 <- factor (vdata$Group4, levels=c("Urban Edge","Urban Interior","Rural Edge","Rural Interior"))
vdata$Type <- factor (vdata$Type, levels=c("Fungi","Bacteria"))

ggplot(vdata)+
geom_point(aes(x=soil, y=microbe,color=Type, shape=Group4))+
geom_smooth(method="lm", aes(x=soil, y=microbe, group=Type, color=Type))+
scale_shape_manual(values = c(15, 16, 17, 18))+
scale_color_manual(values = c("black", "darkgray"))+  
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
labs (x="",y="microbial community",color="", shape="") + 
scale_x_continuous("Distance to centroids \nof soil properties")+
scale_y_continuous("Distance to centroids \nof microbial communities")

ggsave("ScatterPlot_Homogenization.png", height=4,width=5)
