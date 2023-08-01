library(vegan)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(reshape2)
library(multcomp)

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv("experimental_design.csv",header=T)

setwd("~/R/Analysis/2_UNE/ITS")
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T)
ASV <- ASV.table [,1:(ncol(ASV.table)-7)]
ASV.t <- t(ASV)
ASV.t.ITS <- as.matrix(ASV.t)

setwd("~/R/Analysis/2_UNE/16S")
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T)
ASV <- ASV.table [,1:(ncol(ASV.table)-6)]
ASV.t <- t(ASV)
ASV.t.16S <- as.matrix(ASV.t)

ASV.t <-cbind(ASV.t.ITS, ASV.t.16S)
setwd("~/R/Analysis/2_UNE/Others")

# Categolize plots to four location groups
DESIGN$Group4 <- NA
for (i in 1: nrow(DESIGN)){
if (DESIGN$Urban.DFE[i]=="Urban0"|DESIGN$Urban.DFE[i]=="Urban15") {DESIGN$Group4[i]<-"Urban.edge"} else 
if (DESIGN$Urban.DFE[i]=="Urban30"|DESIGN$Urban.DFE[i]=="Urban60"|DESIGN$Urban.DFE[i]=="Urban90") {DESIGN$Group4[i]<-"Urban.interior"} else 
if (DESIGN$Urban.DFE[i]=="Rural0"|DESIGN$Urban.DFE[i]=="Rural15") {DESIGN$Group4[i]<-"Rural.edge"} else 
if (DESIGN$Urban.DFE[i]=="Rural30"|DESIGN$Urban.DFE[i]=="Rural60"|DESIGN$Urban.DFE[i]=="Rural90") {DESIGN$Group4[i]<-"Rural.interior"}  
}

# Permdisp analysis
subset <- subset(ASV.t, DESIGN$Year=="Y2018")
DESIGN.subset <- subset(DESIGN, DESIGN$Year=="Y2018")
permdisp2018 <- betadisper(vegdist(subset, method="bray"), DESIGN.subset$Group4)

subset <- subset(ASV.t, DESIGN$Year=="Y2019")
DESIGN.subset <- subset(DESIGN, DESIGN$Year=="Y2019")
permdisp2019 <- betadisper(vegdist(subset, method="bray"), DESIGN.subset$Group4)

subset <- subset(ASV.t, DESIGN$Year=="Y2021")
DESIGN.subset <- subset(DESIGN, DESIGN$Year=="Y2021")
permdisp2021 <- betadisper(vegdist(subset, method="bray"), DESIGN.subset$Group4)

# Multiple comparison between locations
dist <- NULL
dist <- c(permdisp2018$distance,permdisp2019$distance,permdisp2021$distance)
data <- cbind (dist, DESIGN$Group4)
colnames(data) <- c("dist", "Group4")
data <- data.frame(data)
data$dist <- as.numeric(data$dist)
data$Group4 <- as.factor(data$Group4)

amod<-aov(dist~Group4,data=data) 
tuk<-glht(amod,linfct=mcp(Group4= "Tukey"))
tuk.cld<-cld(tuk,level=0.05,decreasing=TRUE)
letter <- tuk.cld$mcletters$Letters
letter.t <- t(letter)
letter  <- data.frame(letter.t)

# Format data
data$Urban <- NA
data$Edge <- NA
for (i in 1: nrow(data)){
if (data$Group4[i]=="Urban.edge") {data$Urban[i]<-"Urban" 
data$Edge[i]<-"Edge"} else 
if (data$Group4[i]=="Urban.interior") {data$Urban[i]<-"Urban" 
data$Edge[i]<-"Interior"} else 
if (data$Group4[i]=="Rural.edge") {data$Urban[i]<-"Rural" 
data$Edge[i]<-"Edge"} else 
if (data$Group4[i]=="Rural.interior") {data$Urban[i]<-"Rural" 
data$Edge[i]<-"Interior"} 
}

data$resU <- NA
data$resR <- NA

for (k in 1:nrow(data)){
if (data$Edge[k] == "Edge") {data$resU[k] <- letter$Urban.edge} else{data$resU[k] <- letter$Urban.interior}
if (data$Edge[k] == "Edge") {data$resR[k] <- letter$Rural.edge} else{data$resR[k] <- letter$Rural.interior}}

# Visualize
data$Urban <- factor (data$Urban, levels=c("Urban","Rural"))

ggplot(data)+
geom_boxplot(aes(x=Urban, y=dist, fill =Urban, group=Urban))+
scale_fill_manual(values = c("red","black"))+  
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=14,color="black"), legend.position="none")+
labs (y="Distance to centroids",x=" ",fill="")+
facet_wrap(~Edge)+
geom_text(data,mapping=aes(label=resU),x=-Inf, y=Inf, hjust=-4.5,  vjust =9, size=5)+
geom_text(data,mapping=aes(label=resR),x=-Inf, y=Inf, hjust=-11.5, vjust =9, size=5)

ggsave("Boxplot_Homogenization_ITS&16S.tif", height=4,width=4,device='tiff', dpi=1200)

# ANOVA
anova <- anova(lm(dist~Urban*Edge, data=data))
write.csv(anova,"anova.homogenization.ITS&16S.csv")
