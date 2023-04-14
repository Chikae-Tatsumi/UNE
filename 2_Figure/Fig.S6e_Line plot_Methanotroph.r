library(lme4)
library(lmerTest)
library(ggplot2)

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv("experimental_design.csv",header=T)
setwd("~/R/Analysis/2_UNE/16S/function")

func <- read.csv(file="aggregated.function.table.csv",row.names = 1,header=T)
data <- cbind (func$Methanotroph, DESIGN)
colnames(data)[1] <- "Methanotroph"

# ANOVA
anova <- anova(lmer(scale(Methanotroph)~ scale(DFB)*scale(DFE)+(1|Site),data=data))
write.csv(anova, "Methanotroph_rel.abund_anova.csv")

# Format data
standard_error <- function(x) sd(x) / sqrt(length(x)) 
mean <- aggregate(func, by=list(DESIGN$Urban,DESIGN$DFE),FUN= "mean")
se <- aggregate(func, by=list(DESIGN$Urban,DESIGN$DFE),FUN= standard_error)

vdata <- cbind(mean[,1:2], mean$Methanotroph, se$Methanotroph)
colnames(vdata) <- c("Urban","DFE","Methanotroph","se")

# Visualize
vdata$Urban <- factor (vdata$Urban, levels=c("Urban","Rural"))

ggplot(vdata)+
geom_line(aes(x=DFE, y=Methanotroph, color=Urban, group=Urban))+
geom_point(aes(x=DFE, y=Methanotroph, color=Urban, group=Urban), 
size=2,position=position_dodge(0.2))+
geom_errorbar(aes(x=DFE, ymin=Methanotroph-se, ymax=Methanotroph+se, color=Urban, group=Urban),position=position_dodge(0.2)) +
scale_color_manual(values = c("red","black"))+  
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
theme(legend.position="none")+
labs (y="Methanotroph bacteria (%)", x="Distance from edge (m)") 

ggsave("Plot_Methanotroph_rel.abund.png",width = 3, height = 3)
