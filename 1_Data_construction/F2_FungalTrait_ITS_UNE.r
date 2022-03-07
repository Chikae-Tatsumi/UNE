# Import files
setwd("~/R/Analysis/2_UNE/ITS")
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T)
DATABASE <- read.csv(file="~/R/Database/FungalTrait.csv",header=T)

ASV <- ASV.table [,1:(ncol(ASV.table)-7)] 
taxonomy <- ASV.table [,(ncol(ASV.table)-6):ncol(ASV.table)]
taxonomy <- data.frame(lapply(taxonomy, function(x){gsub(pattern="g__", replacement = "", x)}),stringsAsFactors = FALSE)
DATABASE$Genus <- DATABASE$GENUS
taxonomy$id  <- 1:nrow(taxonomy)

merge <- merge(taxonomy, DATABASE, by="Genus", all.x=TRUE)
merge.tidy <- merge[order(merge$id), ]
merge.tidy <- merge.tidy[,-2:-(ncol(taxonomy))]
DATABASE <- DATABASE[,-ncol(DATABASE)]
colnames (merge.tidy) [2:ncol(merge.tidy)] <- colnames (DATABASE)
rownames(merge.tidy) <- rownames (ASV.table)

write.csv(merge.tidy, "FungalTrait/fungaltrait.csv")

fungaltrait.table <- cbind(ASV.table, merge.tidy[,-1])
rownames(fungaltrait.table) <- rownames (ASV.table)
write.csv(fungaltrait.table, "FungalTrait/fungaltrait.table.csv")
