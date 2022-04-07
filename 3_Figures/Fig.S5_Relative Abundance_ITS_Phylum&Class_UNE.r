library(ggplot2)
library(stringr)
library(stringi)

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv("experimental_design.csv",header=T)
setwd("~/R/Analysis/2_UNE/ITS")
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T)

# % table
ASV <- ASV.table [,1:(ncol(ASV.table)-7)]
taxonomy <- ASV.table [,(ncol(ASV.table)-6):ncol(ASV.table)]
percent <- ASV / mean(colSums(ASV)) *100
# Remove "p__" and "c__" before phylum name
taxonomy <- data.frame(lapply(taxonomy, function(x){gsub(pattern="p__", replacement = "", x)}),stringsAsFactors = FALSE)
taxonomy <- data.frame(lapply(taxonomy, function(x){gsub(pattern="c__", replacement = "", x)}),stringsAsFactors = FALSE)

# >1% phylum table
phylum.ag <- aggregate(percent, by=list(taxonomy$Phylum),FUN = sum,na.rm=F)
row.names(phylum.ag)<- phylum.ag[,1]
phylum.ag <- phylum.ag[,-1]
phylum.ag <- data.frame(phylum.ag)
rowMeans <- rowMeans(phylum.ag)
phylum.ag <- cbind(phylum.ag,rowMeans)
major.phylum <- phylum.ag[phylum.ag[,"rowMeans"] > 1,]
major.phylum <- major.phylum[order(major.phylum$rowMeans,decreasing = T),]
minor.phylum <- phylum.ag[phylum.ag[,"rowMeans"] < 1,]
others <- colSums(percent)-colSums (major.phylum[,-ncol(major.phylum)])
selected.phylum <- rbind (major.phylum, others) 
rownames (selected.phylum)[nrow(selected.phylum)] <- "Others"
selected.phylum <- selected.phylum[,-ncol(selected.phylum)] 

# >1% class in >10% class tables
more.major.phylum <- phylum.ag[phylum.ag[,"rowMeans"] > 10,]
more.major.phylum <- more.major.phylum[order(more.major.phylum$rowMeans,decreasing = T),]
percent.table <- cbind (percent, taxonomy)
class.list <- list()
for (i in 1:nrow(more.major.phylum)){
subset <-subset(percent.table,percent.table$Phylum==rownames(more.major.phylum)[i])
class.list <- c(class.list, list(subset))}

selected.class.list <- list()
for (i in 1:length(class.list)){
class.ASV <- class.list[[i]][,1:(ncol(ASV.table)-7)]
class.taxonomy <- class.list[[i]][,(ncol(ASV.table)-6):ncol(ASV.table)]
class.ag <- aggregate(class.ASV, by=list(class.taxonomy$Class),FUN = sum,na.rm=F)
row.names(class.ag)<- class.ag[,1]
class.ag <- class.ag[,-1]
class.ag <- data.frame(class.ag)
rowMeans <- rowMeans(class.ag)
class.ag <- cbind(class.ag,rowMeans)
major.class <- class.ag[class.ag[,"rowMeans"] > 1,]
major.class <- major.class[order(major.class$rowMeans,decreasing = T),]
other.class <- colSums(class.ASV)-colSums (major.class[,-ncol(major.class)])
selected.class <- rbind (major.class, other.class) 
rownames (selected.class)[nrow(selected.class)] <- paste ("Other", rownames(more.major.phylum)[i],sep="_")
selected.class <- selected.class[,-ncol(selected.class)]
selected.class.list <- c(selected.class.list, list(selected.class))}

# Make dataset
phylum.class.table <- data.frame()
for (i in 1:nrow(more.major.phylum)){
phylum.class.table <- rbind(phylum.class.table,selected.class.list[[i]])}
phylum.class.table <- rbind (phylum.class.table, selected.phylum [((nrow(more.major.phylum)+1):nrow(selected.phylum)),])
selected.t <- t (phylum.class.table)
write.csv(selected.t, "aggregated.phylum.class.table.csv")

bind <- cbind (selected.t, DESIGN)
category <- paste (bind$Urban, bind$DFE,sep="&") 
bind <- cbind (bind, category)
table <- aggregate(selected.t, by=list(bind$category),FUN = mean) 
split <- str_split (table[,1], "&",simplify = TRUE)
colnames(split) <- c("Urban","DFE") 
table <- cbind (split,table[2:(ncol(table))])
Urban <- as.character(table[,1])
DFE <- as.character(table[,2])
table <- data.frame(table)
data <- data.frame()
for (i in 3:(ncol(table))){
Abundance <-table[,i]
Example <- colnames (table)[i]
bind <- cbind(Urban, DFE, Abundance, Example) 
data<-rbind (data,bind)}
data$Abundance <- as.numeric(as.character(data$Abundance))

rownames<-rownames(phylum.class.table)[nrow(phylum.class.table):1]
data$Example <- factor(data$Example, levels = rownames)
data$Urban <- factor (data$Urban, levels=c("Urban","Rural"))

# ggplot
ggplot (data,  mapping=aes(x=DFE, y=Abundance, fill=Example))+ 
geom_bar(aes(), stat="identity", position="stack",color="black")+
scale_fill_manual(values = c("gray","#8900ae","#7CAE00","#e9f7f8","#c3eaeb","#81d5d7","#00BFC4","#1d7b7e","#135355","#ffbcb3","#F8766D","#753c37"))+  
theme_classic()+
theme(text=element_text(size=12,color="black"),axis.text=element_text(size=12,color="black"))+
labs (x="Distance from edge (m)",y="Abundance (%)",fill="")+
facet_wrap(~Urban,ncol=5) 

# Save
ggsave(file = "Rel.Abund.ITS.Phylum.png", width=8, height=5 )
