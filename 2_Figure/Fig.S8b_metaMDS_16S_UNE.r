library(ggplot2)
library(vegan)
library(ggrepel)

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv(file = "experimental_design.csv",header=T)
setwd("~/R/Analysis/2_UNE/16S")
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T)

# % table
ASV <- ASV.table [,1:(ncol(ASV.table)-6)] 
taxonomy <- ASV.table [,(ncol(ASV.table)-5):ncol(ASV.table)] 
percent <- ASV / mean(colSums(ASV)) *100

# nmds
percent.t <- t (percent)
nmds <-metaMDS(percent.t, trace=F, distance="bray", perm=100000)

# Visualize
vdata <- cbind (nmds$points, DESIGN)
rate = 0.7
col = c(rep("darkgray",2),rep("black",(nrow(env.arrows)-2) ))

ggplot(vdata)+
geom_point(data=vdata, aes(x=MDS1,y=MDS2, color = ISA , size = DFE))+
scale_colour_gradient(low="black",high="red","Impervious \nsurface \narea (%)")+         
scale_size_continuous("Distance \nfrom \nedge (m)")+         
theme_classic()+
theme(text=element_text(size=10,color="black"),axis.text=element_text(size=10,color="black"))

# Save
ggsave(file = "NMDS_16S_ISA.png",width =5, height =4)

