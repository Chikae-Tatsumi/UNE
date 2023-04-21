# Download 28S sequence file from https://www.arb-silva.de/no_cache/download/archive/release_138_1/Exports/

library(indicspecies)
library(tidyr)
library(stringr)
library(seqinr)
library(tidyverse)

# Import files 
setwd("~/R/Analysis/2_UNE/ITS/indispecies")
sample <- read.csv ("phylogenetic_info.csv", header=T)
sample <- data.frame(lapply(sample, function(x){gsub(pattern="o__", replacement = "", x)}),stringsAsFactors = FALSE) 
sample <- data.frame(lapply(sample, function(x){gsub(pattern="f__", replacement = "", x)}),stringsAsFactors = FALSE) 
sample <- data.frame(lapply(sample, function(x){gsub(pattern="g__", replacement = "", x)}),stringsAsFactors = FALSE) 
colnames(sample)[1] <- "ASV_ID"

DATABASE <- read.fasta("~/R/Database/tax/SILVA_138.1_LSURef_NR99_tax_silva.fasta", as.string=TRUE)
at <- data.frame(attributes(DATABASE))

# Format database
tax.list <- NULL
id.list <- NULL
seq.list <- NULL
for (i in 1:nrow(at)){
    description <- attributes(DATABASE[[at[i,]]])[["Annot"]]
    annot <- strsplit(description, " +")[[1]][2]
    if (grepl("Fungi", annot, fixed = TRUE)){
        id <- at[i,]
        id.list <- c(id.list, id)
        seq <- DATABASE[[at[i,]]]
        seq.list <- rbind(seq.list, seq)
        tax.name <- strsplit(annot, ";")[[1]]
        tax.name <- c(tax.name, rep(NA,(14-length(tax.name))))
        tax.list <- rbind(tax.list, tax.name)}
}
tax.data <- data.frame (tax.list)
seq.data <- data.frame(seq.list)

database <- tax.data[,6:14]
database$ID <- id.list
database$seq <- seq.data$seq.list
database <- data.frame(database)

colnames(database) <- c("Kingdom","Subkingdom","Phylum","Subphylum","Class","Order","Family","Genus","Species","ID","Sequence")


# Assign 28S sequences to ASVs referring to Genus
database.genus <- na.omit(cbind(database$Genus, database$Sequence))
colnames(database.genus) <- c("Genus","Sequence")
database.genus <- data.frame(database.genus)
db.genus <- database.genus[!duplicated(database.genus$Genus), ]

merge.genus <- merge (sample, db.genus, by.x='Genus', by.y='Genus', all.x=TRUE)
merge.genus <- merge.genus[order(merge.genus$ASV_ID,decreasing = FALSE),]

# Assign 28S sequences to ASVs referring to Family
database.family <- na.omit(cbind(database$Family, database$Sequence))
colnames(database.family) <- c("Family","Sequence")
database.family <- data.frame(database.family)
db.family <- database.family[!duplicated(database.family$Family), ]

merge.family <- merge (sample, db.family, by.x='Family', by.y='Family', all.x=TRUE)
merge.family <- merge.family[order(merge.family$ASV_ID,decreasing = FALSE),]

# Assign 28S sequences to ASVs referring to Order
database.order <- na.omit(cbind(database$Order, database$Sequence))
colnames(database.order) <- c("Order","Sequence")
database.order <- data.frame(database.order)
db.order <- database.order[!duplicated(database.order$Order), ]

merge.order <- merge (sample, db.order, by.x='Order', by.y='Order', all.x=TRUE)
merge.order <- merge.order[order(merge.order$ASV_ID,decreasing = FALSE),]


# Combine the taxonomy information
data <- sample[order(sample$ASV_ID,decreasing = FALSE),]
data$seqs <- merge.genus$Sequence

for (i in 1: nrow(data)){
if (is.na(data$seqs[i])){data$seqs[i] <- merge.family$Sequence[i]} 
if (is.na(data$seqs[i])){data$seqs[i] <- merge.order$Sequence[i]}}

data <- data[!is.na(data$seqs),]
rownames(data) <- data$ASV_ID

write.csv(data[,-1],"phylogenetic_info_28S.csv")
write.fasta (sequences = as.list(data$seqs),names = data$ASV_ID, file.out="seqs_28S.fasta")


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

heatmap <- cbind(data$ASV_ID,data$Urban.edge,data$Rural.edge, data$Urban.interior,data$Rural.interior)
write.csv(heatmap,"heatmap_iTOL_28S.csv")
