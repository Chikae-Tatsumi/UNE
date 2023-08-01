library(ggplot2)
library(dplyr)
library(lme4)
library(lmerTest)

# Katie's data
setwd("~/R/Analysis/2_UNE/Others")
METADATA <- read.csv("UNE root data.csv", header=T)
data <- filter(METADATA , Horizon == "O")

data$DFB <- NA
for (i in 1: nrow(data)){
if (data$Site[i]=="Arnold Arboretum") {data$DFB[i]<-8.01} else 
if (data$Site[i]=="Blue Hills") {data$DFB[i]<-16.34} else 
if (data$Site[i]=="Broad Meadow") {data$DFB[i]<-60.41} else 
if (data$Site[i]=="Harvard Forest 4") {data$DFB[i]<-95.16} else 
if (data$Site[i]=="Harvard Forest 6") {data$DFB[i]<-92.63} else 
if (data$Site[i]=="Hammond Woods") {data$DFB[i]<-9.63} else 
if (data$Site[i]=="Sutherland Woods") {data$DFB[i]<-13.2}}

data$Urban<- NA
for (i in 1: nrow(data)){
if (data$DFB[i]<50) {data$Urban[i]<-"Urban"} else  {data$Urban[i]<-"Rural"} }

data$Edge <- NA
for (i in 1: nrow(data)){
if (data$Site[i]=="Blue Hills") {
if (data$Pit_Type[i]=="Urban Edge"|data$Pit_Type[i]=="Rural Edge") {data$DFE[i]<-0} else {data$DFE[i]<-60}} else{
if (data$Pit_Type[i]=="Urban Edge"|data$Pit_Type[i]=="Rural Edge") {data$DFE[i]<-0} else {data$DFE[i]<-90}} }

for (i in 1: nrow(data)){
if (data$Site[i]=="Arnold Arboretum") {data$Site[i]<-"AA01"} else 
if (data$Site[i]=="Blue Hills") {data$Site[i]<-"BH02"} else 
if (data$Site[i]=="Broad Meadow") {data$Site[i]<-"BM03"} else 
if (data$Site[i]=="Harvard Forest 4") {data$Site[i]<-"HF04"} else 
if (data$Site[i]=="Harvard Forest 6") {data$Site[i]<-"HF06"} else 
if (data$Site[i]=="Hammond Woods") {data$Site[i]<-"HW07"} else 
if (data$Site[i]=="Sutherland Woods") {data$Site[i]<-"SW08"}}

rootka <- data
  
# Chikae's data
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv("experimental_design.csv",header=T)
METADATA <- read.csv("metadata.csv", header=T)
setwd("~/R/Analysis/2_UNE/Others")

rootct <- cbind(METADATA, DESIGN)

dataka <- cbind(rootka$RootDensity, rootka$DFB,data$DFE, rootka$Urban, rootka$Site)
datact <- cbind(rootct$Root_density, rootct$DFB, rootct$DFE, rootct$Urban, rootct$Site)

data <- rbind(dataka, datact)
colnames(data) <- c("RootDensity","DFB","DFE","Urban","Site")
data <- data.frame(data)
data$DFB <- as.numeric(data$DFB)
data$DFE <- as.numeric(data$DFE)
data$RootDensity <- as.numeric(data$RootDensity)
data <- na.omit(data)

# ANOVA
anova <- anova(lmer(scale(RootDensity)~ scale(DFB)*scale(DFE)+(1|Site),data=data))
write.csv(anova, "Root_Density_anova.csv")

# Visualize
vdata <- data
vdata$Urban <- factor (vdata$Urban, levels=c("Urban","Rural"))
ylab=expression(paste("Root density (g ",
                      {cm^-3}, ")", sep=""))

ggplot(vdata)+
geom_point(aes(x=DFE, y=RootDensity, color=Urban),position=position_jitter( width=2, height=0))+ 
geom_smooth(method="lm", aes(x=DFE, y=RootDensity, color=Urban, group=Urban))+  
scale_color_manual(values = c("red","black"))+ 
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
labs (y=ylab, x="Distance from edge (m)", color ="") 

ggsave("Plot_Root_Density.png",width = 5, height = 4)
