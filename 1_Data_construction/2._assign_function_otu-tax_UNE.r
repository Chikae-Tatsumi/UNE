# bin phylo groups for Delgado 

setwd("~/R/Analysis/2_UNE/16S") # Added
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T) # Added
# % table # Added
ASV <- ASV.table [,1:(ncol(ASV.table)-6)]  # I have changed 7 into 6 because the table did not contain the Species column. # Added
taxonomy <- ASV.table [,(ncol(ASV.table)-5):ncol(ASV.table)]  # Also, I have changed 6 into 5. # Added
setwd("~/R/Analysis/2_UNE/16S/function") # Added


#clear environment, source paths, packages and functions.
#rm(list=ls())
library(data.table)

#load data.----
#otu <- readRDS(delgado_dada2_SV_table.path)
#tax <- readRDS(delgado_dada2_tax_table.path)
tax_fun <- readRDS("bacteria_tax_to_function.rds")
otu <- ASV
tax <- taxonomy


#### prep taxonomic data. ####

# remove leading "k__" in taxonomy.
#for (i in 1:ncol(tax)) {
#  tax[, i] <- substring(tax[, i], 4)
#}
#tax <- as.data.frame(tax)

# remove taxa that do not assign to a kingdom from tax and otu table.
#tax <- tax[tax$Kingdom == 'Bacteria',] 
#otu <- otu[, colnames(otu) %in% rownames(tax)]
#tax <- tax[rownames(tax) %in% colnames(otu),]

# rarefy otu table
#set.seed(5) # so that rarefaction is repeatable.
#otu <- otu[rowSums(otu) >= 10000,]
#otu <- vegan::rrarefy(otu, 10000)

# assign function to taxonomy
pathway_names <- colnames(tax_fun)[3:15]
tax[, pathway_names] <- "other"

# taxon assignments
for (i in 1:length(pathway_names)) {
  p <- pathway_names[i]
  
  # Classifications from literature search (multiple taxon levels)
  has_pathway <- tax_fun[tax_fun[,p] == 1,]
  levels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
 # levels <- c("phylum", "class", "order", "family", "genus", "species")
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

#save output.----
saveRDS(tax, "tax_with_function.rds") 

bind <-cbind (ASV,tax)
write.csv(bind, "rarefied_16S_function.csv")