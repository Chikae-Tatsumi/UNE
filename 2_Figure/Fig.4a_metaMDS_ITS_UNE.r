library(ggplot2)
library(vegan)
library(ggrepel)

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv(file = "experimental_design.csv",header=T)
setwd("~/R/Analysis/2_UNE/ITS")
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T)
setwd("~/R/Analysis/2_UNE/ITS/FungalTrait")
fungaltrait <- read.csv(file="aggregated.fungaltrait.table.csv",header=T,row.names=1)

env <- cbind(fungaltrait$ECM, fungaltrait$Plant_pathogenic_capacity, fungaltrait$Animal_parasite,
fungaltrait$soil_saprotroph, fungaltrait$wood_saprotroph,fungaltrait$litter_saprotroph)
colnames(env) <- c("ECM","Plant-pathogen","Animal-pathogen","Soil-saportroph","Wood-saprotroph","Litter-saprotroph")
env <- cbind(-DESIGN$DFB,-DESIGN$DFE,env)
colnames(env)[1:2]<-c("Urbanization","Fragmentation")

# % table
ASV <- ASV.table [,1:(ncol(ASV.table)-7)] 
taxonomy <- ASV.table [,(ncol(ASV.table)-6):ncol(ASV.table)]
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

# Visualize
rate = 0.3
col = c(rep("gray40",2),rep("black",(nrow(env.arrows)-2) ))

ggplot()+
geom_point(data=vdata, aes(x=MDS1,y=MDS2, color = DFB , size = DFE))+
scale_colour_gradient(low="red",high="black","Distance \nfrom \nBoston (km)")+         
scale_size_continuous("Distance \nfrom \nedge (m)")+         
geom_segment(data=env.arrows, aes(x = 0, y = 0, xend = (NMDS1*rate), yend = (NMDS2*rate)), arrow = arrow(length = unit(0.2,"cm")),color=col)+
geom_text_repel(data=env.arrows, aes(x=(NMDS1*rate), y=(NMDS2*rate), label=rownames(env.arrows)),  size=4, color=col) +
theme_classic()+
theme(text=element_text(size=10,color="black"),axis.text=element_text(size=10,color="black"))

# Save
ggsave(file = "NMDS_ITS.tif",width =5, height =4, device='tiff', dpi=1200)
