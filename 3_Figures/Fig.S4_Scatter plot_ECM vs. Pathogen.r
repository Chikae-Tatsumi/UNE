library(vegan)
library(ggplot2)
library(lme4)
library(lmerTest)

setwd("~/R/Analysis/2_UNE")
design <- read.csv("experimental_design.csv",header=T)

setwd("~/R/Analysis/2_UNE/ITS/FungalTrait")
fungaltrait <- read.csv("aggregated.fungaltrait.table.csv",header=T)

data <- cbind(design,fungaltrait)
data$Urban <- factor (data$Urban, levels=c("Urban","Rural"))

########### Plant Pathogen ###########

cor <- cor.test (data$ECM, data$Plant_pathogenic_capacity,method="p")
Pval <- cor$p.value
estimate <- round (cor$estimate, digits =2)
if (Pval[1] > 0.05) {result <- paste ("R= ",estimate, sep="")
} else if (Pval[1] > 0.01) {result <- paste ("R= ",estimate," *", sep="")
} else if (Pval[1] > 0.001) {result <- paste("R= ",estimate," **", sep="")
} else {result <- paste("R= ",estimate," ***", sep="")}

ggplot(data)+
geom_point(aes(x=ECM, y=Plant_pathogenic_capacity),position=position_jitter(width=2, height=0))+
geom_smooth(method="lm", aes(x=ECM, y=Plant_pathogenic_capacity), color ="black")+
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
labs (x="ECM fungi (%)",y="Fungal plant-pathogenic capacity (%)") + 
annotate("text", x=Inf, y=Inf, hjust=2, vjust =2, label=result, size=5) 

ggsave("Plot_ECMvs.PlantPatho.png",width = 4, height = 4)

########### Animal parasite ###########

cor <- cor.test (data$ECM, data$Animal_parasite,method="p")
Pval <- cor$p.value
estimate <- round (cor$estimate, digits =2)
if (Pval[1] > 0.05) {result <- paste ("R= ",estimate, sep="")
} else if (Pval[1] > 0.01) {result <- paste ("R= ",estimate," *", sep="")
} else if (Pval[1] > 0.001) {result <- paste("R= ",estimate," **", sep="")
} else {result <- paste("R= ",estimate," ***", sep="")}

ggplot(data)+
geom_point(aes(x=ECM, y=Animal_parasite),position=position_jitter(width=2, height=0))+
geom_smooth(method="lm", aes(x=ECM, y=Animal_parasite), color ="black")+
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
labs (x="ECM fungi (%)",y="Fungal animal parasite (%)") + 
annotate("text", x=Inf, y=Inf, hjust=2, vjust =2, label=result, size=5) 

ggsave("Plot_ECMvs.AnimalParasite.png",width = 4, height = 4)
