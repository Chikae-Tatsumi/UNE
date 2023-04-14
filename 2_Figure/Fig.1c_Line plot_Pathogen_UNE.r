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

anova <- anova(lmer(scale(Plant_pathogenic_capacity)~ scale(DFB)*scale(DFE)+(1|Site),data=data))
write.csv(anova, "Plant_pathogen_anova.csv")

anova <- anova(lmer(scale(Animal_parasite)~ scale(DFB)*scale(DFE)+(1|Site),data=data))
write.csv(anova, "Animal_parasite_anova.csv")

# Format data
Plant <- data.frame(cbind(DESIGN$Urban, DESIGN$DFE, "Plant", fungaltrait$Plant_pathogenic_capacity))
Animal <- data.frame(cbind(DESIGN$Urban, DESIGN$DFE, "Animal",fungaltrait$Animal_parasite))

pathogen <- data.frame(rbind(Plant, Animal))
colnames(pathogen) <- c("Urban","DFE","Type","Pathogenic_capacity")
pathogen.numeric <- as.numeric(pathogen$Pathogenic_capacity)

standard_error <- function(x) sd(x) / sqrt(length(x)) 
mean <- aggregate(pathogen.numeric, by=list(pathogen$Urban,pathogen$DFE,pathogen$Type),FUN= "mean")
se <- aggregate(pathogen.numeric, by=list(pathogen$Urban,pathogen$DFE,pathogen$Type),FUN= standard_error)

vdata <- cbind(mean[,1:3], mean$x, se$x)
colnames(vdata) <- c("Urban","DFE","Type","pathogen","se")

vdata$DFE <- as.numeric(vdata$DFE)
vdata$Type_Urban <- paste(vdata$Type, vdata$Urban, sep="-")

# Visualize
vdata$Urban <- factor (vdata$Urban, levels=c("Urban","Rural"))
vdata$Type <- factor (vdata$Type, levels=c("Plant","Animal"))
vdata$Type_Urban <- factor (vdata$Type_Urban, levels=c("Plant-Urban","Plant-Rural","Animal-Urban","Animal-Rural"))

ggplot(vdata)+
geom_line(aes(x=DFE, y=pathogen, color=Urban, linetype=Type, group=Type_Urban))+
geom_point(aes(x=DFE, y=pathogen, color=Urban, shape=Type, group=Type_Urban), 
size=2,position=position_dodge(0.2))+
geom_errorbar(aes(x=DFE, ymin=pathogen-se, ymax=pathogen+se, color=Urban, group=Type_Urban),position=position_dodge(0.2)) +
scale_color_manual(values = c("red","black"))+  
scale_shape_manual(values = c(16,2))+  
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=14,color="black"))+
labs (y="Soil pathogenic fungi (%)", x="Distance from edge (m)", color=NULL, shape=NULL, linetype=NULL)  

ggsave("Plot_Pathogen_rel.abund.png",width = 4, height = 3)
