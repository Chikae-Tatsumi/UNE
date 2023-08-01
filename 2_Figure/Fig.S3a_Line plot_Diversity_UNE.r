library(vegan)
library(ggplot2)
library(lme4)
library(lmerTest)

setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv("experimental_design.csv",header=T)

# ITS
setwd("~/R/Analysis/2_UNE/ITS")
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T)
ASV <- ASV.table [,1:(ncol(ASV.table)-7)] 
ASV.t <- t(ASV)
shannon.ITS <- diversity(ASV.t, index="shannon",base=2)

data.ITS <- cbind(DESIGN, shannon.ITS)
colnames(data.ITS) [ncol(data.ITS)] <- "shannon"
data.ITS$Type <- "Fungi"

anova <- anova(lmer(scale(shannon)~ scale(DFB)*scale(DFE)+(1|Site),data=data.ITS))
write.csv(anova, "Diversity_ITS_anova.csv")

# 16S
setwd("~/R/Analysis/2_UNE/16S")
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T)
ASV <- ASV.table [,1:(ncol(ASV.table)-6)] 
ASV.t <- t(ASV)
shannon.16S <- diversity(ASV.t, index="shannon",base=2)

data.16S <- cbind(DESIGN, shannon.16S)
colnames(data.16S) [ncol(data.16S)] <- "shannon"
data.16S$Type <- "Bacteria"

anova <- anova(lmer(scale(shannon)~ scale(DFB)*scale(DFE)+(1|Site),data=data.16S))
write.csv(anova, "Diversity_16S_anova.csv")

# Format data
setwd("~/R/Analysis/2_UNE/Others")
bind <- rbind(data.ITS,data.16S)
diversity <- as.numeric(bind$shannon)

standard_error <- function(x) sd(x) / sqrt(length(x)) 
mean <- aggregate(diversity, by=list(bind$Urban,bind$DFE,bind$Type),FUN= "mean")
se <- aggregate(diversity, by=list(bind$Urban,bind$DFE,bind$Type),FUN= standard_error)

vdata <- cbind(mean[,1:3], mean$x, se$x)
colnames(vdata) <- c("Urban","DFE","Type","shannon","se")

vdata$Type_Urban <- paste(vdata$Type, vdata$Urban, sep="-")

# Visualize
vdata$Type_Urban <- factor (vdata$Type_Urban, levels=c("Fungi-Urban","Fungi-Rural","Bacteria-Urban","Bacteria-Rural"))
vdata$Urban <- factor (vdata$Urban, levels=c("Urban","Rural"))
vdata$Type <- factor (vdata$Type, levels=c("Fungi","Bacteria"))

ggplot(vdata)+
geom_line(aes(x=DFE, y=shannon, color=Urban, linetype=Type, group=Type_Urban))+
geom_point(aes(x=DFE, y=shannon, color=Urban, shape=Type, group=Type_Urban), 
size=2,position=position_dodge(0.2))+
geom_errorbar(aes(x=DFE, ymin=shannon-se, ymax=shannon+se, color=Urban, group=Type_Urban),position=position_dodge(0.2)) +
scale_color_manual(values = c("red","black"))+  
scale_shape_manual(values = c(16,2))+  
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
labs (y="Shannon diversity index", x="Distance from edge (m)",color=NULL, shape =NULL, linetype=NULL)  

ggsave("Plot_Diversity.png",width = 4, height = 3)
