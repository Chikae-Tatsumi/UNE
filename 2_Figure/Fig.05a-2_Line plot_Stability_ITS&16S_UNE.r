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
dissim.ITS <- read.csv("ITS_2018baseline_withAveragedMetadata.csv", header=T, row.names=1)
setwd("~/R/Analysis/2_UNE/16S/effect_of_year")
dissim.16S <- read.csv("16S_2018baseline_withAveragedMetadata.csv", header=T, row.names=1)
setwd("~/R/Analysis/2_UNE/Others")

# Format data
data.ITS.1921 <- data.frame(cbind(dissim.ITS$Dissimilarity, dissim.ITS$DFB,dissim.ITS$DFE, dissim.ITS$Urban, dissim.ITS$Site, dissim.ITS$Year, "Fungi"))
data.16S.1921 <- data.frame(cbind(dissim.16S$Dissimilarity, dissim.16S$DFB,dissim.16S$DFE, dissim.16S$Urban, dissim.16S$Site, dissim.16S$Year, "Bacteria"))
colnames(data.ITS.1921) <- c("Dissimilarity","DFB","DFE","Urban","Site","Year","Type")
colnames(data.16S.1921) <- colnames(data.ITS.1921) 

data.ITS.1921$Dissimilarity <- as.numeric(data.ITS.1921$Dissimilarity)
data.ITS.1921$DFB <- as.numeric(data.ITS.1921$DFB)
data.ITS.1921$DFE <- as.numeric(data.ITS.1921$DFE)
data.ITS.1921$Year <- as.numeric(data.ITS.1921$Year)

data.16S.1921$Dissimilarity <- as.numeric(data.16S.1921$Dissimilarity)
data.16S.1921$DFB <- as.numeric(data.16S.1921$DFB)
data.16S.1921$DFE <- as.numeric(data.16S.1921$DFE)
data.16S.1921$Year <- as.numeric(data.16S.1921$Year)

data.ITS.18 <- cbind(0,data.ITS.1921[,2:5], 2018, "Fungi")
data.16S.18 <- cbind(0,data.16S.1921[,2:5], 2018, "Bacteria")
colnames(data.ITS.18) <- colnames(data.ITS.1921) 
colnames(data.16S.18) <- colnames(data.ITS.1921) 

data.ITS <- data.frame(rbind(data.ITS.18 , data.ITS.1921))
data.16S <- data.frame(rbind(data.16S.18, data.16S.1921))

# Categorize locations to edge or interior
data.ITS$Edge <- NA
for (i in 1: nrow(data.ITS)){
if (data.ITS$DFE[i]==0|data.ITS$DFE[i]==15) {data.ITS$Edge[i]<-"Edge"} else 
if (data.ITS$DFE[i]==30|data.ITS$DFE[i]==60|data.ITS$DFE[i]==90) {data.ITS$Edge[i]<-"Interior"} }

data.16S$Edge <- NA
for (i in 1: nrow(data.16S)){
if (data.16S$DFE[i]==0|data.16S$DFE[i]==15) {data.16S$Edge[i]<-"Edge"} else 
if (data.16S$DFE[i]==30|data.16S$DFE[i]==60|data.16S$DFE[i]==90) {data.16S$Edge[i]<-"Interior"} }

# Convert dissimilarity to similarity
data <- rbind(data.ITS, data.16S)
data.16S.1921$Similarity <- (100-data.16S.1921$Dissimilarity)/100
data.ITS.1921$Similarity <- (100-data.ITS.1921$Dissimilarity)/100
data$Similarity <- (100-data$Dissimilarity)/100

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

ggsave ("Line_Similarity.png",height=4, width=5)
