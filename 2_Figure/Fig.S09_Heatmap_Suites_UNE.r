library(vegan)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(reshape2)
library(wesanderson)

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv("experimental_design.csv",header=T)
func.16S <- read.csv (file = "~/R/Analysis/2_UNE/16S/function/aggregated.function.table.csv",row.names=1,header=T) 
func.ITS <- read.csv (file = "~/R/Analysis/2_UNE/ITS/FungalTrait/aggregated.fungaltrait.table.csv",row.names=1,header=T) 
qPCR <- read.csv (file = "~/R/Analysis/2_UNE/qPCR/qPCR.calc.result.UNE.csv",row.names=1,header=T) 
picrust <- read.csv (file = "~/R/Analysis/2_UNE/16S/PICRUSt2/aggregated.picrust2.function.table.csv",row.names=1,header=T) 
setwd("~/R/Analysis/2_UNE/Others")

# Categolize plots to four location groups
DESIGN$Group4 <- NA
for (i in 1: nrow(DESIGN)){
if (DESIGN$Urban.DFE[i]=="Urban0"|DESIGN$Urban.DFE[i]=="Urban15") {DESIGN$Group4[i]<-"Urban.edge"} else 
if (DESIGN$Urban.DFE[i]=="Urban30"|DESIGN$Urban.DFE[i]=="Urban60"|DESIGN$Urban.DFE[i]=="Urban90") {DESIGN$Group4[i]<-"Urban.interior"} else 
if (DESIGN$Urban.DFE[i]=="Rural0"|DESIGN$Urban.DFE[i]=="Rural15") {DESIGN$Group4[i]<-"Rural.edge"} else 
if (DESIGN$Urban.DFE[i]=="Rural30"|DESIGN$Urban.DFE[i]=="Rural60"|DESIGN$Urban.DFE[i]=="Rural90") {DESIGN$Group4[i]<-"Rural.interior"}  
}

# Make dataset
func <- cbind (func.ITS$ECM, log10(qPCR$AM.copies.gsoil+1), 
func.ITS$soil_saprotroph, func.ITS$wood_saprotroph, func.ITS$litter_saprotroph, func.ITS$dung_saprotroph, 
func.ITS$Plant_pathogenic_capacity, func.ITS$Animal_parasite,
func.16S$Copiotroph, func.16S$Oligotroph, 
func.16S$Cellulolytic, func.16S$Lignolytic,
func.16S$Methanotroph, func.16S$Chitinolytic, 
func.16S$Assim_nitrate_reduction, func.16S$Dissim_nitrate_reduction,
log10(qPCR$nifH.copies.gsoil+1), log10(qPCR$amoA.copies.gsoil+1),log10(qPCR$nosZ.copies.gsoil+1),
picrust$Xenobiotics_degradation.count, picrust$Plant_pathogen_interaction.count, picrust$Infectious_disease_bacterial.count
)
colnames(func) <- c("ECM fungi","AM fungi", 
"Soil saprotroph fungi","Wood saprotroph fungi","Litter saprotroph fungi","Dung saprotroph fungi",
"Plant-pathogenic fungi", "Animal-parasitic fungi",
"Copiotroph bacteria", "Oligotroph bacteria",
"Cellulolytic bacteria","Lignolytic bacteria",
"Methanotroph bacteria", "Chitinolytic bacteria",    
"Assimiratory nitrate reducing bacteria", "Dissimiratory nitrate reducing bacteria", 
"N fixing bacteria","Nitrifying bacteria", "Denitrifying bacteria",
"Xenobiotics-degrading bacteria", "Plant-pathogenic bacteria", "Human-pathogenic bacteria"
)

ag <- aggregate(func, by=list(DESIGN$Group4),FUN = mean,na.rm=F) 
rownames(ag) <- ag[,1]
ag <- ag[,-1]
ag <- as.matrix(ag)
scale <- scale (ag)

# Visualize
vdata <- melt(scale)
vdata$Var1 <- factor (vdata$Var1, levels=c("Urban.edge","Urban.interior","Rural.edge","Rural.interior"))
pal <- wes_palette("Zissou1", 100, type = "continuous")

ggplot (vdata)+
geom_tile(aes(x=Var1, y=Var2, fill =value))+
scale_fill_gradientn(colours = pal, name = "Normalized \nabundance")+
labs(x="",y="")+
theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.3)) 

# Save
ggsave(file = "Heatmap_Suites.png", height=5, width=6)
