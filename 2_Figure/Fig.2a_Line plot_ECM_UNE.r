library(ggplot2)
library(lme4)
library(lmerTest)

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv("experimental_design.csv",header=T)
setwd("~/R/Analysis/2_UNE/ITS/FungalTrait")
fungaltrait <- read.csv(file="aggregated.fungaltrait.table.csv",header=T,row.names = 1)

# ANOVA
data <- cbind (fungaltrait, DESIGN)
anova <- anova(lmer(scale(ECM)~ scale(DFB)*scale(DFE)+(1|Site),data=data))
write.csv(anova, "ECM_rel.abund_anova.csv")

# Format data
standard_error <- function(x) sd(x) / sqrt(length(x)) 
mean <- aggregate(fungaltrait, by=list(DESIGN$Urban,DESIGN$DFE),FUN= "mean")
se <- aggregate(fungaltrait, by=list(DESIGN$Urban,DESIGN$DFE),FUN= standard_error)

vdata <- cbind(mean[,1:2], mean$ECM, se$ECM)
colnames(vdata) <- c("Urban","DFE","ECM","se")

# Visualize
vdata$Urban <- factor (vdata$Urban, levels=c("Urban","Rural"))

ggplot(vdata)+
geom_line(aes(x=DFE, y=ECM, color=Urban, group=Urban))+
geom_point(aes(x=DFE, y=ECM, color=Urban, group=Urban), 
size=2,position=position_dodge(0.2))+
geom_errorbar(aes(x=DFE, ymin=ECM-se, ymax=ECM+se, color=Urban, group=Urban),position=position_dodge(0.2)) +
scale_color_manual(values = c("red","black"))+ 
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
labs (y="Soil ECM fungi (%)", x="Distance from edge (m)") 

ggsave("Plot_ECM_rel.abund.tif",width = 4, height = 3, device='tiff', dpi=1200)
