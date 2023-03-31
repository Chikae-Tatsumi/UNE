setwd("~/R/Analysis/2_UNE/ITS/FungalTrait")

# For rarefied table
fungaltrait.table <- read.csv(file="rarefied_fungaltrait.table.csv",row.names = 1)
ECM.table<- fungaltrait.table[grep("ectomycorrhizal", fungaltrait.table$primary_lifestyle),]
write.csv(ECM.table, "rarefied_ECM_table_FungalTrait.csv")
