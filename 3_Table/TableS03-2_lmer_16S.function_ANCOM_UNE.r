library (lme4)
library (lmerTest)
library(data.table)
library(tidyverse)
library (dplyr)

### Assign_function_otu-tax ###
# Copied from https://github.com/zoey-rw/NEFI_microbe/blob/master/16S/data_construction/1._create_tax_to_function_reference.r

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv("experimental_design.csv",header=T)
setwd("~/R/Analysis/2_UNE/16S")
ASV.table <- read.table(file="ancom_ASV_table.txt",header=T) 
setwd("~/R/Analysis/2_UNE/16S/function") 

# % table 
ASV <- ASV.table [,1:(ncol(ASV.table)-6)]  
taxonomy <- ASV.table [,(ncol(ASV.table)-5):ncol(ASV.table)]  

# Load database
tax_fun <- readRDS("bacteria_tax_to_function.rds")
otu <- ASV
tax <- taxonomy

# Assign function to taxonomy
pathway_names <- colnames(tax_fun)[3:15]
tax[, pathway_names] <- "other"

for (i in 1:length(pathway_names)) {
  p <- pathway_names[i]
  
# Classifications from literature search (multiple taxon levels)
  has_pathway <- tax_fun[tax_fun[,p] == 1,]
  levels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  for (j in 1:length(levels)) {
    taxon_level <- levels[j]
    has_pathway_taxon_level <- has_pathway[has_pathway$Taxonomic.level==taxon_level,]
    if (taxon_level == "Species") {
      if(nrow(tax[which(paste(tax$Genus, tax$Species) %in% has_pathway_taxon_level$Taxon),]) > 0) {
      tax[which(paste(tax$Genus, tax$Species) %in% has_pathway_taxon_level$Taxon),][,p] <- p
      }
    } else {
    if (nrow(tax[tax[[taxon_level]] %in% has_pathway_taxon_level$Taxon,]) > 0){
      tax[tax[[taxon_level]] %in% has_pathway_taxon_level$Taxon,][,p] <- p
    }
    }
  }
}

# Save output
saveRDS(tax, "tax_with_function_ancom.rds") 

bind <-cbind (ASV,tax)
write.csv(bind, "16S_function_ancom.csv")


# Aggregate the function table
function.table <- bind
ASV <- function.table [,1:(ncol(function.table)-19)] 
fgs <- function.table [,(ncol(function.table)-12):ncol(function.table)] 

name <- colnames(fgs) 

Result <- NULL
for (i in 1:ncol(fgs)){
aggregated <- aggregate(ASV, by=list(fgs[,i]),FUN = sum,na.rm=F) 
rownames(aggregated) <- aggregated [,1]
aggregated <- aggregated [,-1]
aggregated.t <- t(aggregated)
data <- cbind (aggregated.t, DESIGN)

Result <- cbind (Result, aggregated.t)
Result <- Result [,-ncol(Result)]}
colnames(Result) <- name
write.csv (Result, "aggregated.function.table_ancom.csv")

### Run lmer ###
objective <- Result

NAME1 <- "Urbanization"
NAME2 <- "Fragmentation"

# ANOVA
Results <- NULL
for (i in 1:ncol(objective)){
data <- cbind(DESIGN, objective[,i])
colnames (data)[ncol(data)] <- colnames(objective)[i]
data <- na.omit (data)
model <- lmer(scale(data[,ncol(data)])~ scale(-DFB)*scale(-DFE)+(1|Site),data=data)
anova <- anova(model)
summary <- summary(model)

estimate <- summary$coefficients[,1]
Fval <- anova[,5]
Pval <-anova[,6]
Bind <- c(estimate[2], Fval[1], Pval[1], estimate[3], Fval[2], Pval[2], estimate[4],Fval[3], Pval[3])
names(Bind) <- c(paste(NAME1,".R",sep=""),paste(NAME1,".F",sep=""),paste(NAME1,".P",sep=""),paste(NAME2,".R",sep=""),paste(NAME2,".F",sep=""),paste(NAME2,".P",sep=""),paste(NAME1,"*",NAME2,".R",sep=""),paste(NAME1,"*",NAME2,".F",sep=""),paste(NAME1,"*",NAME2,".P",sep=""))
Results <- rbind(Results, Bind)}

rownames(Results) <- colnames(objective)

# P.val --> Asterisk
Results.asterisk <- Results
for (k in 1:ncol(Results)){
if(k%%3 == 0){
for (i in 1:nrow(Results)){
if (Results[i,k] < 0.001) {Results.asterisk[i,k] <- "***"
} else if (Results[i,k] < 0.01) {Results.asterisk[i,k] <- "**"
} else if (Results[i,k] < 0.05) {Results.asterisk[i,k] <- "*"
} else {Results.asterisk[i,k]<- "n.s"}}
}}

# Save
write.csv(Results.asterisk, "lmer.16S.function_ancom.csv") 
