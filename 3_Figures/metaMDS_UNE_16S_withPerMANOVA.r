library(ggplot2)
library(vegan)
library(ggrepel)

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv(file = "experimental_design.csv",header=T)
DESIGN <- na.omit(DESIGN)
setwd("~/R/Analysis/2_UNE/16S")
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T)

# % table
ASV <- ASV.table [,1:(ncol(ASV.table)-6)] # I have changed 7 into 6 because the table did not contain the Species column.
taxonomy <- ASV.table [,(ncol(ASV.table)-5):ncol(ASV.table)] # Also, I have changed 6 into 5.
percent <- ASV / mean(colSums(ASV)) *100
# Remove "k__","p__", "c__"  before phylum name
taxonomy <- data.frame(lapply(taxonomy, function(x){gsub(pattern="k__", replacement = "", x)}),stringsAsFactors = FALSE)
taxonomy <- data.frame(lapply(taxonomy, function(x){gsub(pattern="p__", replacement = "", x)}),stringsAsFactors = FALSE)
taxonomy <- data.frame(lapply(taxonomy, function(x){gsub(pattern="c__", replacement = "", x)}),stringsAsFactors = FALSE)

# Make Kingdom table
kingdom <- aggregate(percent, by=list(taxonomy$Kingdom),FUN = sum,na.rm=F) 
row.names(kingdom)<-kingdom [,1]
kingdom <- kingdom[,-1]
kingdom <- data.frame(kingdom)
kingdom.t <- t (kingdom)

# Make >1% Phylum table
phylum <- aggregate(percent, by=list(taxonomy$Phylum),FUN = sum,na.rm=F) 
row.names(phylum)<-phylum[,1]
phylum <- phylum[,-1]
phylum <- data.frame(phylum)
rowMeans <- rowMeans(phylum) 
phylum <- cbind(phylum,rowMeans)
major.phylum <- phylum[phylum[,"rowMeans"] > 5,] # Change
major.phylum <- major.phylum[,-ncol(major.phylum)] 
major.phylum.t <- t (major.phylum)

# Make >1% Class (in >10% phylum) table
percent.table <- cbind (percent, taxonomy)
more.major.phylum <- phylum[phylum[,"rowMeans"] > 10,]
more.major.phylum.t <- t(more.major.phylum)
class.table <- data.frame()
for (i in 1: ncol (more.major.phylum.t)){
class.subset <- subset (percent.table, Phylum = colnames(more.major.phylum.t)[i])
class.table <- rbind (class.table, class.subset)}
class.ASV <- class.table[,1:(ncol(class.table)-6)] # I have changed 7 into 6 because the table did not contain the Species column.
class.taxonomy <- class.table[,(ncol(class.table)-5):ncol(class.table)] # Also, I have changed 6 into 5.
class <- aggregate(class.ASV, by=list(class.taxonomy$Class),FUN = sum,na.rm=F) 
row.names(class)<-class[,1]
class <- class[,-1]
class <- data.frame(class)
rowMeans <- rowMeans(class) 
class <- cbind(class,rowMeans)
major.class <- class[class[,"rowMeans"] > 5,] # Change
major.class <- major.class[,-ncol(major.class)] 
major.class.t <- t (major.class)

# nmds
percent.t <- t (percent)
nmds <-metaMDS(percent.t, trace=F, distance="bray", perm=100000)
env <- cbind(-DESIGN$DFB,-DESIGN$DFE)
colnames(env)<-c("Urbanization","Fragmentation")
env.fit <- envfit(nmds, env, perm=100000, na.rm=TRUE)

# Make dataset
data <- cbind (nmds$points, DESIGN)
env.values <- (env.fit$vectors[1]$arrows)
pval <- (env.fit$vectors[4]$pvals)
env.arrows <- cbind (env.values, pval)
env.arrows <- data.frame(subset (env.arrows, pval < 0.05))

# PerMANOVA
ASV.t <- t(ASV)
adonis <- adonis(ASV.t ~ DFB*DFE,  data=DESIGN, permutations=10000)
Pval <- adonis[[1]][,6]
if (adonis[[1]][1,6] > 0.05) {DFB.result <- ""
} else if (adonis[[1]][1,6] > 0.01) {DFB.result <- "        Urbanization *      "
} else if (adonis[[1]][1,6] > 0.001) {DFB.result <- "        Urbanization **     "
} else {DFB.result <- "        Urbanization ***    "}
if (adonis[[1]][2,6] > 0.05) {DFE.result <- ""
}else if (adonis[[1]][2,6] > 0.01) {DFE.result <- "        Fragmentation *       "
}else if (adonis[[1]][2,6] > 0.001) {DFE.result <- "        Fragmentation **     "
}else {DFE.result <- "        Fragmentation ***    "}
if (adonis[[1]][3,6] > 0.05) {DFBDFE.result <- ""
}else if (adonis[[1]][3,6] > 0.01) {DFBDFE.result <- "        U×F *      "
}else if (adonis[[1]][3,6] > 0.001) {DFBDFE.result <- "        U×F**     "
}else {DFBDFE.result <- "        U×F ***    "}

# ggplot
rate = 0.5 # how many times shorter than acctual arrows
ggplot()+
geom_point(data=data, aes(x=MDS1,y=MDS2, color = DFB , size = DFE))+
scale_colour_gradient(low="#f6766d",high="navy")+            
geom_segment(data=env.arrows, aes(x = 0, y = 0, xend = (NMDS1*rate), yend = (NMDS2*rate)), arrow = arrow(length = unit(0.3,"cm")),color="black")+
geom_text_repel(data=env.arrows, aes(x=(NMDS1*rate), y=(NMDS2*rate), label=rownames(env.arrows)),  size=6, color="black") +
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
annotate("text", x=Inf, y=-Inf, hjust=0.9, vjust =-6, label=DFB.result, size=5) +
annotate("text", x=Inf, y=-Inf, hjust=0.9, vjust =-4, label=DFE.result, size=5) +
annotate("text", x=Inf, y=-Inf, hjust=0.9, vjust =-2, label=DFBDFE.result, size=5) 

# Save
ggsave(file = "NMDS_16S.png",width =5, height =4)