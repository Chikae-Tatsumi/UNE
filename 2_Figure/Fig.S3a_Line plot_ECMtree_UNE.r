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
standard_error <- function(x, na.rm = TRUE) {
  sd(x, na.rm = TRUE) / sqrt(length(na.omit(x)))}
mean <- aggregate(METADATA, by=list(DESIGN$Urban,DESIGN$DFE),FUN= "mean", na.rm=TRUE)
se <- aggregate(METADATA, by=list(DESIGN$Urban,DESIGN$DFE),FUN= standard_error)

data <- cbind(mean[,1:2], mean$ECM.tree, se$ECM.tree)
colnames(data) <- c("Urban","DFE","ECM.tree","se")

data$Urban <- factor (data$Urban, levels=c("Urban","Rural"))

ggplot(data)+
geom_line(aes(x=DFE, y=ECM.tree, color=Urban, group=Urban))+
geom_point(aes(x=DFE, y=ECM.tree, color=Urban, group=Urban), 
size=2,position=position_dodge(0.2))+
geom_errorbar(aes(x=DFE, ymin=ECM.tree-se, ymax=ECM.tree+se, color=Urban, group=Urban),position=position_dodge(0.2)) +
scale_color_manual(values = c("red","black"))+  
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
labs (y= "ECM tree (%)", x="Distance from edge (m)", color="")   

ggsave("Plot_ECMtree.png",width=5, height=4)
