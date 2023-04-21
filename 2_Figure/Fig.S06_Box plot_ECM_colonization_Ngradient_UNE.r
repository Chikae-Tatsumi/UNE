library(ggplot2)
library(lme4)
library(lmerTest)
library(multcomp)

# Import files
setwd("~/R/Analysis/2_UNE/Others")
data <- read.csv("Mycorrhizal_colonization_Ngradient.csv",header=T)
data$Fertilization <- factor (data$Fertilization, levels=c("High","Low","Control"))

# ANOVA
anova <- anova(lm(scale(Percent_Mycorrhizal_Colonization)~ Fertilization,data=data))
write.csv(anova, "ECM_Colonization_Ngradient_anova.csv")

# Multiple comparison
amod<-aov(Percent_Mycorrhizal_Colonization~Fertilization,data=data) 
tuk<-glht(amod,linfct=mcp(Fertilization= "Tukey"))
tuk.cld<-cld(tuk,level=0.05,decreasing=TRUE)
letter <- tuk.cld$mcletters$Letters
letter.t <- t(letter)
letter  <- data.frame(letter.t)

# Visualize
ggplot(data)+
geom_boxplot(aes(x=Fertilization, y=Percent_Mycorrhizal_Colonization, fill = Fertilization))+
scale_fill_manual(values = c("black","gray","white"))+ 
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
labs (y="ECM colonization (%)", x="N fertilization", fill="") +
geom_text(vdata,mapping=aes(label=letter$High),x=-Inf, y=Inf, hjust=-3.5,  vjust =19, size=5)+
geom_text(vdata,mapping=aes(label=letter$Low),x=-Inf, y=Inf, hjust=-12.5, vjust =5, size=5)+
geom_text(vdata,mapping=aes(label=letter$Control),x=-Inf, y=Inf, hjust=-21.5,  vjust =2, size=5)

ggsave("Plot_ECM_Colonization_Ngradeint.png",width = 5, height = 4)

