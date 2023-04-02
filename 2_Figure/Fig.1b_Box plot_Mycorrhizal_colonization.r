library(ggplot2)
library(lme4)
library(lmerTest)

# Import files
setwd("~/R/Analysis/2_UNE/Others")
data <- read.csv("Mycorrhizal_colonization.csv",header=T)

# ANOVA
anova <- anova(lm(scale(Percent_Mycorrhizal_Colonization)~ scale(DFB),data=data))
write.csv(anova, "Mycorrhizal_Colonization_anova.csv")

# Visualize
data$Urban <- factor (data$Urban, levels=c("Urban","Rural"))

ggplot(data)+
geom_boxplot(aes(x=Urban, y=Percent_Mycorrhizal_Colonization, fill = Urban))+
scale_fill_manual(values = c("red","black"))+ 
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
labs (y="ECM colonization (%)", x="", fill="") 

ggsave("ECM_Colonization_UNE.png",width = 5, height = 4)

