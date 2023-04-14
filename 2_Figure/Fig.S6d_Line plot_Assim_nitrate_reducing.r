library(lme4)
library(lmerTest)
library(ggplot2)

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv("experimental_design.csv",header=T)
setwd("~/R/Analysis/2_UNE/16S/function")

func <- read.csv(file="aggregated.function.table.csv",row.names = 1,header=T)
data <- cbind (func$Assim_nitrate_reduction, DESIGN)
colnames(data)[1] <- "Assim_nitrate_reduction"

# ANOVA
anova <- anova(lmer(scale(Assim_nitrate_reduction)~ scale(DFB)*scale(DFE)+(1|Site),data=data))
write.csv(anova, "Assim_nitrate_reduction_rel.abund_anova.csv")

# Format data
standard_error <- function(x) sd(x) / sqrt(length(x)) 
mean <- aggregate(func, by=list(DESIGN$Urban,DESIGN$DFE),FUN= "mean")
se <- aggregate(func, by=list(DESIGN$Urban,DESIGN$DFE),FUN= standard_error)

vdata <- cbind(mean[,1:2], mean$Assim_nitrate_reduction, se$Assim_nitrate_reduction)
colnames(vdata) <- c("Urban","DFE","Assim_nitrate_reduction","se")

# Visualize
vdata$Urban <- factor (vdata$Urban, levels=c("Urban","Rural"))

ggplot(vdata)+
geom_line(aes(x=DFE, y=Assim_nitrate_reduction, color=Urban, group=Urban))+
geom_point(aes(x=DFE, y=Assim_nitrate_reduction, color=Urban, group=Urban), 
size=2,position=position_dodge(0.2))+
geom_errorbar(aes(x=DFE, ymin=Assim_nitrate_reduction-se, ymax=Assim_nitrate_reduction+se, color=Urban, group=Urban),position=position_dodge(0.2)) +
scale_color_manual(values = c("red","black"))+  
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
labs (y="Assim_nitrate_reduction bacteria (%)", x="Distance from edge (m)") +
theme(legend.position="none")+
scale_y_continuous("Assimiratory nitrate-reducing \nbacteria (%)")

ggsave("Plot_Assim_nitrate_reduction_rel.abund.png",width = 3, height = 3)
