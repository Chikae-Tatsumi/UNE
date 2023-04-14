library(vegan)
library(ggplot2)
library(lme4)
library(lmerTest)

# Import files
setwd("~/R/Analysis/2_UNE")
design <- read.csv("experimental_design.csv",header=T)
setwd("~/R/Analysis/2_UNE/ITS/FungalTrait")
fungaltrait <- read.csv("aggregated.fungaltrait.table.csv",header=T)

# ITS
ASV.table <- read.table(file="~/R/Analysis/2_UNE/ITS/rarefied_ASV_table.txt",header=T)
ASV <- ASV.table [,1:(ncol(ASV.table)-7)] 
ASV.t <- t(ASV)
shannon <- diversity(ASV.t, index="shannon",base=2)
data.ITS <- cbind(design,fungaltrait,shannon)
data.ITS$Type <- "Fungi"

cor <- cor.test (data.ITS$ECM, data.ITS$shannon,method="p")
sink("Diversity.ITS.vs.ECM_cor.txt")
cor
sink()

# 16S
ASV.table <- read.table(file="~/R/Analysis/2_UNE/16S/rarefied_ASV_table.txt",header=T)
ASV <- ASV.table [,1:(ncol(ASV.table)-6)] 
ASV.t <- t(ASV)
shannon <- diversity(ASV.t, index="shannon",base=2)
data.16S <- cbind(design,fungaltrait,shannon)
data.16S$Type <- "Bacteria"

cor <- cor.test (data.16S$ECM, data.16S$shannon,method="p")
sink("Diversity.16S.vs.ECM_cor.txt")
cor
sink()

# Format data
data <- rbind(data.ITS,data.16S)
data$Type_Urban <- paste(data$Type, data$Urban, sep="-")

# Visualize
vdata <- data

vdata$Type_Urban <- factor (vdata$Type_Urban, levels=c("Fungi-Urban","Fungi-Rural","Bacteria-Urban","Bacteria-Rural"))
vdata$Urban <- factor (vdata$Urban, levels=c("Urban","Rural"))
vdata$Type <- factor (vdata$Type, levels=c("Fungi","Bacteria"))

ggplot(vdata)+
geom_point(aes(x=ECM, y=shannon,shape=Type), color="black")+
geom_smooth(method="lm", aes(x=ECM, y=shannon, group=Type, linetype=Type),color="black")+
scale_shape_manual(values = c(16,2))+ 
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
labs (x="ECM fungi (%)",y="Shannon diversity index",shape=NULL, linetype=NULL) +
ylim(2,10)

ggsave("Plot_ECM.vs.Diversity.png",width = 4, height = 3)
