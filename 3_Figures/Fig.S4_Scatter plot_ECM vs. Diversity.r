library(vegan)
library(ggplot2)
library(lme4)
library(lmerTest)

setwd("~/R/Analysis/2_UNE")
design <- read.csv("experimental_design.csv",header=T)

setwd("~/R/Analysis/2_UNE/ITS/FungalTrait")
fungaltrait <- read.csv("aggregated.fungaltrait.table.csv",header=T)

########### Fungal Diversity ###########

ASV.table <- read.table(file="~/R/Analysis/2_UNE/ITS/rarefied_ASV_table.txt",header=T)

ASV <- ASV.table [,1:(ncol(ASV.table)-7)] 
ASV.t <- t(ASV)
shannon <- diversity(ASV.t, index="shannon",base=2)
data <- cbind(design,fungaltrait,shannon)

cor <- cor.test (data$ECM, data$shannon,method="p")
Pval <- cor$p.value
estimate <- round (cor$estimate, digits =2)
if (Pval[1] > 0.05) {result <- paste ("R= ",estimate, sep="")
} else if (Pval[1] > 0.01) {result <- paste ("R= ",estimate," *", sep="")
} else if (Pval[1] > 0.001) {result <- paste("R= ",estimate," **", sep="")
} else {result <- paste("R= ",estimate," ***", sep="")}

ggplot(data)+
geom_point(aes(x=ECM, y=shannon),position=position_jitter(width=2, height=0))+
geom_smooth(method="lm", aes(x=ECM, y=shannon), color ="black")+
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
labs (x="ECM fungi (%)",y="Fungal total diversity") + 
annotate("text", x=Inf, y=Inf, hjust=2, vjust =2, label=result, size=5) 

ggsave("Plot_ECMvs.FungalDiversity.png",width = 4, height = 4)

########### Bacterial Diversity ###########

ASV.table <- read.table(file="~/R/Analysis/2_UNE/16S/rarefied_ASV_table.txt",header=T)

ASV <- ASV.table [,1:(ncol(ASV.table)-6)] 
ASV.t <- t(ASV)
shannon <- diversity(ASV.t, index="shannon",base=2)
data <- cbind(design,fungaltrait)
data <- na.omit(data)
data <- cbind(data,shannon)

data$Urban <- factor (data$Urban, levels=c("Urban","Rural"))

cor <- cor.test (data$ECM, data$shannon,method="p")
Pval <- cor$p.value
estimate <- round (cor$estimate, digits =2)
if (Pval[1] > 0.05) {result <- paste ("R= ",estimate, sep="")
} else if (Pval[1] > 0.01) {result <- paste ("R= ",estimate," *", sep="")
} else if (Pval[1] > 0.001) {result <- paste("R= ",estimate," **", sep="")
} else {result <- paste("R= ",estimate," ***", sep="")}

ggplot(data)+
geom_point(aes(x=ECM, y=shannon),position=position_jitter(width=2, height=0))+
geom_smooth(method="lm", aes(x=ECM, y=shannon), color ="black")+
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
labs (x="ECM fungi (%)",y="Bacterial total diversity") + 
annotate("text", x=Inf, y=Inf, hjust=2, vjust =2, label=result, size=5) 

ggsave("Plot_ECMvs.BacterialDiversity.png",width = 4, height = 4)
