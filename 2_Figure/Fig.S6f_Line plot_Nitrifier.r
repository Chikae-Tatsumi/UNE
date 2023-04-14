library(lme4)
library(lmerTest)
library(ggplot2)

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv("experimental_design.csv",header=T)
setwd("~/R/Analysis/2_UNE/qPCR")

qPCR <- read.csv (file = "qPCR.calc.result.UNE.csv",row.names=1,header=T) 
data <- cbind (log10(qPCR$amoA.copies.gsoil+1), DESIGN)
colnames(data)[1] <- "Nitrifier"

# ANOVA
anova <- anova(lmer(scale(Nitrifier)~ scale(DFB)*scale(DFE)+(1|Site),data=data))
write.csv(anova, "Nitrifier_rel.abund_anova.csv")

# Format data
nitrifier <- data$Nitrifier

standard_error <- function(x) sd(x) / sqrt(length(x)) 
mean <- aggregate(nitrifier, by=list(DESIGN$Urban,DESIGN$DFE),FUN= "mean")
se <- aggregate(nitrifier, by=list(DESIGN$Urban,DESIGN$DFE),FUN= standard_error)

vdata <- cbind(mean[,1:2], mean$x, se$x)
colnames(vdata) <- c("Urban","DFE","Nitrifier","se")

# Visualize
vdata$Urban <- factor (vdata$Urban, levels=c("Urban","Rural"))

ggplot(vdata)+
geom_line(aes(x=DFE, y=Nitrifier, color=Urban, group=Urban))+
geom_point(aes(x=DFE, y=Nitrifier, color=Urban, group=Urban), 
size=2,position=position_dodge(0.2))+
geom_errorbar(aes(x=DFE, ymin=Nitrifier-se, ymax=Nitrifier+se, color=Urban, group=Urban),position=position_dodge(0.2)) +
scale_color_manual(values = c("red","black"))+  
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
theme(legend.position="none")+
labs (y="", x="Distance from edge (m)")+ 
scale_y_continuous("Log number of \nnitrifying bacteria \n(copies g-1)")

# Save
ggsave("Nitrifier_abund.png",width = 3, height = 3)
