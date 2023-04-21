library(ggplot2)
library(lme4)
library(lmerTest)

setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv("experimental_design.csv",header=T)
METADATA <- read.csv(file = "metadata.csv",header=T)
setwd("~/R/Analysis/2_UNE/Others")

# ANOVA
data <- cbind (METADATA, DESIGN)
anova <- anova(lmer(scale(ECM.tree)~ scale(DFB)*scale(DFE)+(1|Site),data=data))
write.csv(anova, "ECMtree_anova.csv")

# Visualize
vdata <- data
vdata$Urban <- factor (vdata$Urban, levels=c("Urban","Rural"))

ggplot(vdata)+
geom_point(aes(x=DFE, y=ECM.tree, color=Urban),position=position_jitter(width=2, height=0))+ 
geom_smooth(method="lm", aes(x=DFE, y=ECM.tree, group=Urban, color=Urban))+  
scale_color_manual(values = c("red","black"))+  # if you want to change the colors
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
labs (y= "ECM tree (%)", x="Distance from edge (m)", color="")  # if you want to change the axis titles

ggsave("Plot_ECMtree.png",width=5, height=4)
