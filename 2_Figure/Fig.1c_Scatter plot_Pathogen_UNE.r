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
anova <- anova(lmer(scale(Plant_pathogenic_capacity)~ scale(DFB)*scale(DFE)+(1|Site),data=data))
write.csv(anova, "Plant_pathogen_anova.csv")

anova <- anova(lmer(scale(Animal_parasite)~ scale(DFB)*scale(DFE)+(1|Site),data=data))
write.csv(anova, "Animal_parasite_anova.csv")

# Format data
Plant <- data.frame(cbind(data$DFB, data$DFE, data$Plant_pathogenic_capacity,"Plant"))
Animal <- data.frame(cbind(data$DFB, data$DFE, data$Animal_parasite, "Animal"))
Plant <- cbind(Plant, data$Urban)
Animal <- cbind(Animal, data$Urban)

Pathogen <- rbind(Plant, Animal)
colnames(Pathogen) <- c("DFB","DFE","Pathogenic_capacity","Type","Urban")

Pathogen <- data.frame(Pathogen)
Pathogen$Pathogenic_capacity <- as.numeric(Pathogen$Pathogenic_capacity)
Pathogen$DFE <- as.numeric(Pathogen$DFE)

Pathogen$Type_Urban <- paste(Pathogen$Type, Pathogen$Urban, sep="-")

# Visualize
Pathogen$Urban <- factor (Pathogen$Urban, levels=c("Urban","Rural"))
Pathogen$Type <- factor (Pathogen$Type, levels=c("Plant","Animal"))
Pathogen$Type_Urban <- factor (Pathogen$Type_Urban, levels=c("Plant-Urban","Plant-Rural","Animal-Urban","Animal-Rural"))

ggplot(Pathogen)+
geom_point(aes(x=DFE, y=Pathogenic_capacity, color=Urban, shape=Type),position=position_jitter(width=2, height=0))+ 
geom_smooth(method="lm", aes(x=DFE, y=Pathogenic_capacity, group=Type_Urban, color=Urban))+  
scale_color_manual(values = c("red","black"))+  
scale_shape_manual(values = c(16,2))+  
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=14,color="black"))+
labs (y="Soil pathogenic fungi (%)", x="Distance from edge (m)",color="", shape ="")  

ggsave("Pathogen_rel.abund_UNE.png",width = 5, height = 4)
