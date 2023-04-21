library(vegan)
library(geosphere)

# Import data
setwd("~/R/Analysis/2_UNE")
distance.matrix <- read.csv("Autocorrelation/dd.distances.csv",header=T, row.names=1)
distance.matrix <- as.matrix(distance.matrix)

setwd("~/R/Analysis/2_UNE/ITS")
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T)
ASV <- ASV.table [,1:(ncol(ASV.table)-7)]
ASV.t <- t(ASV)
ASV.t <- as.matrix(ASV.t)

# Run Mantel test
ASV.dist <- vegdist(ASV.t,method="bray")
dist.dist <- as.dist(distance.matrix)

mantel <- mantel(ASV.dist, dist.dist, method = "spearman", permutations = 99999, na.rm = TRUE)
sink("Mantel_ITS.txt")
mantel
sink()

# Visualize
ggplot()+
geom_point(aes(x=dist.dist, y=ASV.dist),shape=1)+ 
geom_smooth(method="lm", aes(x=dist.dist, y=ASV.dist),colour="black")+  
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
labs (y="Fungal community difference", x="Geological distance (km)")  

ggsave("Mantel.ITS.png",width=5,height=4)
