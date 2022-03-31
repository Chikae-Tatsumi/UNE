library(vegan)
library(geosphere)

setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv(file = "experimental_design.csv",header=T)
distance.matrix <- read.csv("autocorrelation/dd.distances.csv",header=T, row.names=1)

setwd("~/R/Analysis/2_UNE/16S")
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T)
ASV <- ASV.table [,1:(ncol(ASV.table)-6)]
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
labs (y="Bacterial community difference", x="Geological distance (km)") +  # if you want to change the axis titles
annotate("text", x=-Inf, y=-Inf, hjust=-.5, vjust = -3, label=paste("R=",round(mantel$statistic,digit=4)), size=5) +
annotate("text", x=-Inf, y=-Inf, hjust=-.5, vjust = -1, label=paste("P=",round(mantel$signif,digit=8)), size=5) 

ggsave("Mantel.16S.png",width=5,height=4)