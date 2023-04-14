library(ggplot2)
library(lme4)
library(lmerTest)

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv("experimental_design.csv",header=T)
setwd("~/R/Analysis/2_UNE/ITS/FungalTrait")
fungaltrait <- read.csv("aggregated.fungaltrait.table.csv",header=T, row.names=1)
setwd("~/R/Analysis/2_UNE/qPCR")
qPCR <- read.csv (file = "qPCR.calc.result.UNE.csv",row.names=1,header=T) 

# cor.test
data <- cbind (log10(qPCR$nosZ.copies.gsoil+1), fungaltrait, DESIGN)
colnames(data)[1] <- "Denitrifier"

cor <- cor.test (data$ECM, data$Denitrifier,method="p")
sink("ECM.vs.Denitrifier_cor.txt")
cor
sink()

# Visualize
vdata <- data

ggplot(vdata)+
geom_point(aes(x=ECM, y=Denitrifier),color="black")+
geom_smooth(method="lm", aes(x=ECM, y=Denitrifier),color="black")+
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
theme(legend.position="none")+
labs (x="ECM fungi (%)",y="")+
scale_y_continuous("Log number of \ndenitrifying bacteria \n(copies g-1)")+
xlim(0,110)

ggsave("Plot_ECM.vs.Denitrifier.png",width = 3, height = 3)
