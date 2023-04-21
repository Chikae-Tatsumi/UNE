library(indicspecies)
library(tidyr)
library(stringr)
library(seqinr)

# Import files
setwd("~/R/Analysis/2_UNE/ITS")
seqs <- read.table(file="rarefied_seq.txt", header=T, row.names=1)
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T,row.names=1)
tax <- ASV.table [,(ncol(ASV.table)-6):ncol(ASV.table)]

database <- cbind(rownames(ASV.table),colnames(seqs), tax)
colnames(database) <- c("ASV_ID","seqs",colnames(tax))

setwd("~/R/Analysis/2_UNE/ITS/indispecies")

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
