library(vegan)
library(ggplot2)

# Import files
setwd("~/R/Analysis/2_UNE")
distance.matrix <- read.csv("autocorrelation/dd.distances.csv",header=T, row.names=1)
distance.matrix <- as.matrix(distance.matrix)

setwd("~/R/Analysis/2_UNE/ITS")
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T)
ASV <- ASV.table [,1:(ncol(ASV.table)-7)]
ASV.t <- t(ASV)
ASV.t <- as.matrix(ASV.t)

# Calculate
ASV.dist<-vegdist(ASV.t,method="bray")
dist.dist <- as.dist(distance.matrix)
mantel <- mantel(ASV.dist, dist.dist, method = "spearman", permutations = 99999, na.rm = TRUE)

# ggplot
ggplot()+
geom_point(aes(x=dist.dist, y=ASV.dist),shape=1)+ 
geom_smooth(method="lm", aes(x=dist.dist, y=ASV.dist),colour="black")+  
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
labs (y="Fungal community difference", x="Geological distance (km)") + 
annotate("text", x=-Inf, y=-Inf, hjust=-.5, vjust = -3, label=paste("R=",round(mantel$statistic,digit=4)), size=5) +
annotate("text", x=-Inf, y=-Inf, hjust=-.5, vjust = -1, label=paste("P=",round(mantel$signif,digit=8)), size=5) 

# Save
ggsave("Mantel.ITS.png",width=5,height=4)