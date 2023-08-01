library(vegan)
library(reshape2)
library(ggplot2)

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
data <- cbind(sim.ITS$Similarity, sim.16S$Similarity, sim.16S[,2:ncol(sim.16S)])
colnames(data)[1:2] <-c("fungi","bacteria")

# cor.test
cor <- cor.test (data$fungi, data$Moist,method="p")
sink("Similarity.vs.Moist_ITS_cor.txt")
cor
sink()

cor <- cor.test (data$bacteria, data$Moist,method="p")
sink("Similarity.vs.Moist_16S_cor.txt")
cor
sink()

# Format data
sim.ITS$Type <- "Fungi"
sim.16S$Type <- "Bacteria"
vdata <- rbind(sim.ITS,sim.16S)
vdata$Similarity <- c(data$fungi, data$bacteria)

# Visualize
vdata$Urban <- factor (vdata$Urban, levels=c("Urban","Rural"))
vdata$Type <- factor (vdata$Type, levels=c("Fungi","Bacteria"))

ggplot(vdata)+
geom_point(aes(x=Moist, y=Similarity,color=Type))+
geom_smooth(method="lm", aes(x=Moist, y=Similarity, color=Type))+
scale_color_manual(values = c("black", "darkgray"))+  
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
labs (y="",x="Soil moisture (%)",color ="", shape="") +
scale_y_continuous("Community similarity \ncompared to 2018")

ggsave("ScatterPlot_Similarity.vs.Moist.tif", height=4,width=5,device='tiff', dpi=1200)

