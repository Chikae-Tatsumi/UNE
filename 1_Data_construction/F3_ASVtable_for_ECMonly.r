setwd("~/R/Analysis/2_UNE/ITS/FungalTrait")

# For rarefied table
fungaltrait.table <- read.csv(file="rarefied_fungaltrait.table.csv",row.names = 1)
ASV <- fungaltrait.table [,1:(ncol(fungaltrait.table)-31)] 
guild <- fungaltrait.table [,(ncol(fungaltrait.table)-24):ncol(fungaltrait.table)] 
percent <- ASV / mean(colSums(ASV)) *100

percent.table <- cbind (percent, guild)
ECM.table<- percent.table[grep("ectomycorrhizal", percent.table$primary_lifestyle),]
write.csv(ECM.table, "rarefied_ECM_table_FungalTrait.csv")


# For full table
fungaltrait.table <- read.csv(file="fungaltrait.table.csv",row.names = 1)
ASV <- fungaltrait.table [,1:(ncol(fungaltrait.table)-31)] 
guild <- fungaltrait.table [,(ncol(fungaltrait.table)-24):ncol(fungaltrait.table)] 
percent <- ASV / mean(colSums(ASV)) *100

percent.table <- cbind (percent, guild)
ECM.table<- percent.table[grep("ectomycorrhizal", percent.table$primary_lifestyle),]
write.csv(ECM.table, "ECM_table_FungalTrait.csv")
