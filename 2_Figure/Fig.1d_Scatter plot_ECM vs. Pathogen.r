library(ggplot2)
library(lme4)
library(lmerTest)

# Import files
setwd("~/R/Analysis/2_UNE")
design <- read.csv("experimental_design.csv",header=T)
setwd("~/R/Analysis/2_UNE/ITS/FungalTrait")
fungaltrait <- read.csv("aggregated.fungaltrait.table.csv",header=T)

data <- cbind(design,fungaltrait)

# Cor.test
cor <- cor.test (data$ECM, data$Plant_pathogenic_capacity,method="p")
sink("Plant_pathogen_cor.txt")
cor
sink()

cor <- cor.test (data$ECM, data$Animal_parasite,method="p")
sink("Animal_parasite_cor.txt")
cor
sink()

# Format data
Plant <- data.frame(cbind(data$DFB, data$DFE, data$ECM, data$Plant_pathogenic_capacity,"Plant"))
Animal <- data.frame(cbind(data$DFB, data$DFE, data$ECM, data$Animal_parasite, "Animal"))

vdata <- rbind(Plant, Animal)
colnames(vdata) <- c("DFB","DFE","ECM","Pathogenic_capacity","Type")

vdata <- data.frame(vdata)
vdata$Pathogenic_capacity <- as.numeric(vdata$Pathogenic_capacity)
vdata$DFE <- as.numeric(vdata$DFE)
vdata$ECM <- as.numeric(vdata$ECM)

# Visualize
vdata$Type <- factor (vdata$Type, levels=c("Plant","Animal"))

ggplot(vdata)+
geom_point(aes(x=ECM, y=Pathogenic_capacity, shape= Type),color="black", position=position_jitter(width=2, height=0))+
geom_smooth(method="lm", aes(x=ECM, y=Pathogenic_capacity, group=Type, linetype=Type),color="black")+
scale_shape_manual(values = c(16,2))+ 
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=14,color="black"))+
labs (x="Soil ECM fungi (%)",y="Soil pathogenic fungi (%)",shape=NULL, linetype=NULL)  

ggsave("ECM.vs.Pathogen_UNE.png",width = 4, height = 3)
