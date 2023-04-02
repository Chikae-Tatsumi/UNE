library(ggplot2)
library(lme4)
library(lmerTest)

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv("experimental_design.csv",header=T)
setwd("~/R/Analysis/2_UNE/ITS/FungalTrait")
fungaltrait <- read.csv(file="aggregated.fungaltrait.table.csv",header=T,row.names = 1)

data <- cbind (fungaltrait, DESIGN)

# ANOVA
anova <- anova(lmer(scale(ECM)~ scale(DFB)*scale(DFE)+(1|Site),data=data))
write.csv(anova, "ECM_rel.abund_anova.csv")

# Visualize
data$Urban <- factor (data$Urban, levels=c("Urban","Rural"))

ggplot(data)+
geom_point(aes(x=DFE, y=ECM, color=Urban),position=position_jitter( width=2, height=0))+ 
geom_smooth(method="lm", aes(x=DFE, y=ECM, group=Urban, color=Urban))+  
scale_color_manual(values = c("red","black"))+ 
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
labs (y="Soil ECM fungi (%)", x="Distance from edge (m)",color="") 

ggsave("ECM_rel.abund_UNE.png",width = 5, height = 4)
