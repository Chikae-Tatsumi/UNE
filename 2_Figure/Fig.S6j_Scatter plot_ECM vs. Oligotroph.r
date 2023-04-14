library(ggplot2)
library(lme4)
library(lmerTest)

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv("experimental_design.csv",header=T)
setwd("~/R/Analysis/2_UNE/ITS/FungalTrait")
fungaltrait <- read.csv("aggregated.fungaltrait.table.csv",header=T, row.names=1)
setwd("~/R/Analysis/2_UNE/16S/function")
func <- read.csv(file="aggregated.function.table.csv",row.names = 1,header=T)

# cor.test
data <- cbind (func, fungaltrait, DESIGN)

cor <- cor.test (data$ECM, data$Oligotroph,method="p")
sink("ECM.vs.Oligotroph_cor.txt")
cor
sink()

# Visualize
vdata <- data

ggplot(vdata)+
geom_point(aes(x=ECM, y=Oligotroph),color="black")+
geom_smooth(method="lm", aes(x=ECM, y=Oligotroph),color="black")+
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
theme(legend.position="none")+
labs (x="ECM fungi (%)",y="Oligotroph bacteria (%)") + 
xlim(0,110)

ggsave("Plot_ECM.vs.Oligotroph.png",width = 3, height = 3)
