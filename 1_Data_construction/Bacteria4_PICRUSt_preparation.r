library(seqinr)
library(ape)

setwd("~/R/Analysis/2_UNE/16S")
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T,row.names=1)
ASV <- ASV.table [,1:(ncol(ASV.table)-6)]

seq <- read.table(file="rarefied_seq.txt",header=T, row.names=1)
write.fasta (sequences = as.list(seq[,1]),names = rownames(seq),file.out="rarefied_seqs.fasta")

ASV <- cbind (rownames(ASV),ASV)
ASV <- rbind (colnames(ASV),ASV)
ASV[1,1] <- "ID"
rownames(ASV) <- NULL 
colnames(ASV) <- NULL
write.table(ASV, "rarefied_ASV.txt",sep="\t",col.names = F, row.names = F,quote=F)
