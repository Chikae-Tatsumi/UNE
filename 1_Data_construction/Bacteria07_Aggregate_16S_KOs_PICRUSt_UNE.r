library(tidyverse)
library (dplyr)

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv("experimental_design.csv",header=T)
setwd("~/R/Analysis/2_UNE/16S/PICRUSt2")
fun <- read.table("pred_metagenome_unstrat.tsv",header=T,row.names=1)
DATABASE.xeno <- read.csv("~/R/Database/Xenobiotic_degradation_KOs.csv", header=T)
DATABASE.plant <- read.csv("~/R/Database/Plant_pathogen_interaction_KOs.csv", header=T)
DATABASE.inf <- read.csv("~/R/Database/Infectious_disease_bacterial_KOs.csv", header=T)

fun$ID <- rownames(fun)

# Xenobiotic degradation
DATABASE.xeno.unique <- DATABASE.xeno[!duplicated(DATABASE.xeno$ID),]
erged <- merge(DATABASE.xeno.unique, fun, by="ID", all.x=TRUE)
result <- na.omit(merged[,-1:-3])
Xenobiotics_degradation.abundance <- colSums(result)
result.t <- t(result)
Xenobiotics_degradation.diversity <- diversity(result.t, index="shannon", base=2)
result[(result>0)]<- 1
Xenobiotics_degradation.count <- colSums(result)

# Platn pathogen
DATABASE.plant.unique <- DATABASE.plant[!duplicated(DATABASE.plant$ID),]
merged <- merge(DATABASE.plant.unique, fun, by="ID", all.x=TRUE)
result <- na.omit(merged[,-1:-2])
Plant_pathogen_interaction.abundance <- colSums(result)
result.t <- t(result)
Plant_pathogen_interaction.diversity <- diversity(result.t, index="shannon", base=2)
result[(result>0)]<- 1
Plant_pathogen_interaction.count <- colSums(result)

# Human desease infection
DATABASE.inf.unique <- DATABASE.inf[!duplicated(DATABASE.inf$ID),]
merged <- merge(DATABASE.inf.unique, fun, by="ID", all.x=TRUE)
result <- na.omit(merged[,-1:-3])
Infectious_disease_bacterial.abundance <- colSums(result)
result.t <- t(result)
Infectious_disease_bacterial.diversity <- diversity(result.t, index="shannon", base=2)
result[(result>0)]<- 1
Infectious_disease_bacterial.count <- colSums(result)

# Save
Result <- cbind (Xenobiotics_degradation.abundance,Plant_pathogen_interaction.abundance, Infectious_disease_bacterial.abundance,
Xenobiotics_degradation.diversity,Plant_pathogen_interaction.diversity, Infectious_disease_bacterial.diversity,
Xenobiotics_degradation.count,Plant_pathogen_interaction.count, Infectious_disease_bacterial.count
)
write.csv (Result, "aggregated.picrust2.function.table.csv")

