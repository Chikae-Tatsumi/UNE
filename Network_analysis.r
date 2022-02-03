library(igraph)  
library(Hmisc)  
library(Matrix)  

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv("experimental_design.csv")
setwd("~/R/Analysis/2_UNE/ITS")
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T)
ASV <- ASV.table [,1:(ncol(ASV.table)-7)]
taxonomy <- ASV.table [,(ncol(ASV.table)-6):ncol(ASV.table)]
percent <- ASV / mean(colSums(ASV)) *100
percent.t <- t(percent)

# To make minor phylum "Others"
taxonomy.minusp <- data.frame(lapply(taxonomy, function(x){gsub(pattern="p__", replacement = "", x)}),stringsAsFactors = FALSE)
rownames(taxonomy.minusp) <- rownames(taxonomy)
taxonomy <- taxonomy.minusp 
phylum <- aggregate(percent, by=list(taxonomy$Phylum),FUN = sum,na.rm=F) 
row.names(phylum)<-phylum[,1]
phylum <- phylum[,-1]
phylum <- data.frame(phylum)
rowMeans <- rowMeans(phylum) 
phylum <- cbind(phylum,rowMeans)
minor.phylum <- phylum[phylum[,"rowMeans"] < 1,] # Change
minor.phylum.list <- rownames(minor.phylum)

for (i in 1:length (minor.phylum.list)){
taxonomy$Phylum <- gsub(minor.phylum.list[i],"Others",taxonomy$Phylum)}

# Subset
bind <- cbind (percent.t,DESIGN)
subset <- subset(bind, bind$Urban=="Urban") # Change
subset <- subset(subset, subset$Edge=="Edge") # Change

subset <- subset[,1:(ncol(subset)-ncol(DESIGN))]

# Filter to pick up parts of ASVs
subset.t.filter <- subset[ ,colMeans(percent.t) >= 0.05] # To pick up >0.05% ASVs
subset.t.filter <- subset.t.filter[,colSums(subset.t.filter)>0]
print(c(ncol(subset),"versus",ncol(subset.t.filter)))

# Calculate network
percent.cor <- rcorr(as.matrix(subset.t.filter), type="spearman")
percent.pval <- forceSymmetric(percent.cor$P) # Self-correlation as NA
#Select only the taxa for the filtered ASVs by using rownames of percent.pval
sel.tax <- taxonomy[rownames(percent.pval),,drop=FALSE]
#Sanity check --> should be "[1] TRUE"
all.equal(rownames(sel.tax), rownames(percent.pval))

p.yes <- percent.cor$P<0.05
r.yes <- percent.cor$r>0
r.high <- percent.cor$r>0.7
r.val <- percent.cor$r # select all the correlation values 
p.r.yes = p.yes*r.yes*r.val*r.high
adjm<-p.r.yes 

net.grph=graph.adjacency(adjm,mode="undirected",weighted=TRUE,diag=FALSE)

# hs <- hub_score(net.grph, weights=NA)$vector　
# as <- authority_score(net.grph, weights=NA)$vector
# pr <- page.rank(net.grph,directed=F)$vector　
deg <- degree(net.grph, mode="all")　

# Align taxonomy names
sel.tax$Phylum <- factor(sel.tax$Phylum)
others.n <- which(levels(sel.tax$Phylum)=="Others")
levels <- NULL
if (others.n == 1){
    for (i in 1:(length(levels(sel.tax$Phylum)))){
    levels <- c(levels, levels(sel.tax$Phylum)[i])}
    levels <- c(levels, "Others")}else{
for (i in 1:(others.n-1)){
    levels <- c(levels, levels(sel.tax$Phylum)[i])}
for (i in (others.n+1):(length(levels(sel.tax$Phylum)))){
    levels <- c(levels, levels(sel.tax$Phylum)[i])}
levels <- c(levels, "Others")
levels(sel.tax$Phylum) <- levels}

# Illustrate
col=rainbow(length(levels(sel.tax$Phylum)))
plot.igraph(net.grph, vertex.size=deg*0.15,vertex.label=NA, vertex.color=col[unclass(sel.tax$Phylum)],layout=layout.kamada.kawai)
# plot(net.grph, vertex.size=deg*0.15,vertex.label=NA,vertex.color=col[unclass(sel.tax$Phylum)],layout=layout.random)
# plot(net.grph, vertex.size=deg*0.15,vertex.label=NA,vertex.color=col[unclass(sel.tax$Phylum)],layout=layout.fruchterman.reingold)
legend(x = "bottomleft", legend = levels(sel.tax$Phylum), pch = 19, col = col, bty = "n", pt.cex=2, title = "Color legend")

gsize <- gsize(net.grph)
edge_density <- round(edge_density(net.grph),digit=5)
text(x=1,y=-1,paste("The number of edge = ", gsize))
text(x=1,y=-1.1,paste("edge density = ", edge_density))
title("XXXXXX") #Change

# Save 
dev.copy(pdf, file="~/R/Analysis/2_UNE/ITS/Network_Urban_Edge.pdf", height=5, width=10)
dev.off()
