library(indicspecies)
library(tidyr)
library(stringr)
library(seqinr)

# Import files
setwd("~/R/Analysis/2_UNE/16S")
seqs <- read.table(file="rarefied_seq.txt", header=T, row.names=1)
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T,row.names=1)
tax <- ASV.table [,(ncol(ASV.table)-5):ncol(ASV.table)]

database <- cbind(rownames(ASV.table),seqs$x, tax)
colnames(database) <- c("ASV_ID","seqs",colnames(tax))

setwd("~/R/Analysis/2_UNE/16S/indispecies")

# Rural.Edge
data <- read.csv(file="Rural.Edge.csv",header=T)
merged <- merge (data, database, by.x='ASV_ID', by.y='ASV_ID')
rownames(merged) <- merged[,1]
merged$Group <- "Rural.edge"
Rural.edge <- merged[,-1:-3]

# Urban.Edge
data <- read.csv(file="Urban.Edge.csv",header=T)
merged <- merge (data, database, by.x='ASV_ID', by.y='ASV_ID')
rownames(merged) <- merged[,1]
merged$Group <- "Urban.edge"
Urban.edge <- merged[,-1:-3]

# Rural.Interior
data <- read.csv(file="Rural.Interior.csv",header=T)
merged <- merge (data, database, by.x='ASV_ID', by.y='ASV_ID')
rownames(merged) <- merged[,1]
merged$Group <- "Rural.interior"
Rural.interior <- merged[,-1:-3]

# Urban.Interior
data <- read.csv(file="Urban.Interior.csv",header=T)
merged <- merge (data, database, by.x='ASV_ID', by.y='ASV_ID')
rownames(merged) <- merged[,1]
merged$Group <- "Urban.interior"
Urban.interior <- merged[,-1:-3]

# Save sequences
data <- rbind(Urban.edge, Rural.edge, Urban.interior, Rural.interior)
write.fasta (sequences = as.list(data[,1]),names = rownames(data),file.out="seqs.fasta")

write.csv(data,"phylogenetic_info.csv")

# To visualize in iTOL
# Make a binary file to show their present or absent
data$Urban.edge <- 0
data$Rural.edge <- 0
data$Urban.interior <- 0
data$Rural.interior <- 0

for(i in 1:nrow(data)){
if(data$Group[i]=="Urban.edge"){data$Urban.edge[i]<-1}else{data$Urban.edge[i]<-0}
if(data$Group[i]=="Rural.edge"){data$Rural.edge[i]<-1}else{data$Rural.edge[i]<-0}
if(data$Group[i]=="Urban.interior"){data$Urban.interior[i]<-1}else{data$Urban.interior[i]<-0}
if(data$Group[i]=="Rural.interior"){data$Rural.interior[i]<-1}else{data$Rural.interior[i]<-0}}

heatmap <- cbind(rownames(data),data$Urban.edge,data$Rural.edge, data$Urban.interior,data$Rural.interior)
write.csv(heatmap,"heatmap_iTOL.csv")
