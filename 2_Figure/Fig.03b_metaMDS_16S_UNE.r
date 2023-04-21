library(ggplot2)
library(vegan)
library(ggrepel)

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv(file = "experimental_design.csv",header=T)
setwd("~/R/Analysis/2_UNE/16S")
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T)

setwd("~/R/Analysis/2_UNE/16S/function")
func.16S <- read.csv(file="aggregated.function.table.csv",row.names = 1,header=T)
qPCR <- read.csv (file = "~/R/Analysis/2_UNE/qPCR/qPCR.calc.result.UNE.csv",row.names=1,header=T) 
picrust <- read.csv (file = "~/R/Analysis/2_UNE/16S/PICRUSt2/aggregated.picrust2.function.table.csv",row.names=1,header=T) 

func <- cbind (func.16S$Copiotroph, func.16S$Oligotroph, func.16S$Cellulolytic, func.16S$Lignolytic,
func.16S$Methanotroph, 
log10(qPCR$nifH.copies.gsoil+1), log10(qPCR$amoA.copies.gsoil+1),log10(qPCR$nosZ.copies.gsoil+1),
picrust$Xenobiotics_degradation.count, picrust$Plant_pathogen_interaction.count, picrust$Infectious_disease_bacterial.count)
colnames(func) <- c("Copiotroph", "Oligotroph",
"Cellulolytic","Lignolytic","Methanotroph",
"N-fixing","Nitrifying", "Denitrifying",
"Xenobiotics-degrading", "Plant-pathogen", "Human-pathogen")

env <- cbind(-DESIGN$DFB,-DESIGN$DFE,func)
colnames(env)[1:2]<-c("Urbanization","Fragmentation")

# % table
ASV <- ASV.table [,1:(ncol(ASV.table)-6)] 
taxonomy <- ASV.table [,(ncol(ASV.table)-5):ncol(ASV.table)] 
percent <- ASV / mean(colSums(ASV)) *100

# nmds & envfit
percent.t <- t (percent)
nmds <-metaMDS(percent.t, trace=F, distance="bray", perm=100000)
env.fit <- envfit(nmds, env, perm=100000, na.rm=TRUE)

# Format data
vdata <- cbind (nmds$points, DESIGN)
env.values <- (env.fit$vectors[1]$arrows)
pval <- (env.fit$vectors[4]$pvals)
env.arrows <- cbind (env.values, pval)
env.arrows <- data.frame(subset (env.arrows, pval < 0.05))

# Vizualize
rate = 0.4
col = c(rep("gray40",2),rep("black",(nrow(env.arrows)-2) ))

ggplot(vdata)+
geom_point(data=vdata, aes(x=MDS1,y=MDS2, color = DFB , size = DFE))+
scale_colour_gradient(low="red",high="black","Distance \nfrom \nBoston (km)")+         
scale_size_continuous("Distance \nfrom \nedge (m)")+         
geom_segment(data=env.arrows, aes(x = 0, y = 0, xend = (NMDS1*rate), yend = (NMDS2*rate)), arrow = arrow(length = unit(0.2,"cm")),color=col)+
geom_text_repel(data=env.arrows, aes(x=(NMDS1*rate), y=(NMDS2*rate), label=rownames(env.arrows)),  size=3.5, color=col) +
theme_classic()+
theme(text=element_text(size=10,color="black"),axis.text=element_text(size=10,color="black"))

# Save
ggsave(file = "NMDS_16S.png",width =5, height =4)
