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
METADATA<- read.csv("metadata.csv",header=T)

setwd("~/R/Analysis/2_UNE/ITS/effect_of_year")
sim.ITS <- read.csv("ITS_Stability_withAveragedMetadata.csv", header=T, row.names=1)
setwd("~/R/Analysis/2_UNE/16S/effect_of_year")
sim.16S <- read.csv("16S_Stability_withAveragedMetadata.csv", header=T, row.names=1)
setwd("~/R/Analysis/2_UNE/Others")

# Format data
data.ITS.1921 <- data.frame(cbind(sim.ITS$Similarity, sim.ITS$DFB,sim.ITS$DFE, sim.ITS$Urban, sim.ITS$Site, sim.ITS$Year, "Fungi"))
data.16S.1921 <- data.frame(cbind(sim.16S$Similarity, sim.16S$DFB,sim.16S$DFE, sim.16S$Urban, sim.16S$Site, sim.16S$Year, "Bacteria"))
colnames(data.ITS.1921) <- c("Similarity","DFB","DFE","Urban","Site","Year","Type")
colnames(data.16S.1921) <- colnames(data.ITS.1921) 

data.ITS.1921$Similarity <- as.numeric(data.ITS.1921$Similarity)
data.ITS.1921$DFB <- as.numeric(data.ITS.1921$DFB)
data.ITS.1921$DFE <- as.numeric(data.ITS.1921$DFE)
data.ITS.1921$Year <- as.numeric(data.ITS.1921$Year)

data.16S.1921$Similarity <- as.numeric(data.16S.1921$Similarity)
data.16S.1921$DFB <- as.numeric(data.16S.1921$DFB)
data.16S.1921$DFE <- as.numeric(data.16S.1921$DFE)
data.16S.1921$Year <- as.numeric(data.16S.1921$Year)

data.ITS.18 <- cbind(1,data.ITS.1921[,2:5], 2018, "Fungi")
data.16S.18 <- cbind(1,data.16S.1921[,2:5], 2018, "Bacteria")
colnames(data.ITS.18) <- colnames(data.ITS.1921) 
colnames(data.16S.18) <- colnames(data.ITS.1921) 

data.ITS <- data.frame(rbind(data.ITS.18 , data.ITS.1921))
data.16S <- data.frame(rbind(data.16S.18, data.16S.1921))

data <- rbind(data.ITS, data.16S)

# Categorize locations to edge or interior
data$Edge <- NA
for (i in 1: nrow(data)){
if (data$DFE[i]==0|data$DFE[i]==15) {data$Edge[i]<-"Edge"} else 
if (data$DFE[i]==30|data$DFE[i]==60|data$DFE[i]==90) {data$Edge[i]<-"Interior"} }

# ANOVA for ITS
anova <- anova(lmer(scale(Similarity)~ scale(DFB)*scale(DFE)+(1|Site),data=data.ITS.1921))
write.csv(anova, "Similarity_ITS.csv")

# ANOVA for 16S
anova <- anova(lmer(scale(Similarity)~ scale(DFB)*scale(DFE)+(1|Site),data=data.16S.1921))
write.csv(anova, "Similarity_16S.csv")

# Visualize
ag <- aggregate(data, by=list(data$Urban,data$Edge, data$Year, data$Type),FUN= "mean", na.rm=TRUE)

vdata <- cbind(ag[,1:4], ag$Similarity)
colnames(vdata) <- c("Urban","Edge","Year","Type","Similarity")
vdata$Urban.Edge <- paste(vdata$Urban,vdata$Edge, sep=".")

vdata$Year <- as.character(vdata$Year)
vdata$Urban <- factor (vdata$Urban, levels=c("Urban","Rural"))
vdata$Type <- factor (vdata$Type, levels=c("Fungi","Bacteria"))

ggplot(vdata)+
    geom_line(aes(x=Year, y=Similarity, color=Urban, linetype=Edge, group=Urban.Edge))+
    geom_point(aes(x=Year, y=Similarity, color=Urban, shape=Edge),size=2)+
    scale_color_manual(values = c("red","black"))+  
    theme_classic()+
    theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
    labs (y="",x="Year",color ="", shape="", linetype="") +
    scale_y_continuous("Community similarity \ncompared to 2018")+
    facet_wrap(~Type)

ggsave ("Line_Similarity.tif",height=4, width=5,device='tiff', dpi=1200)
