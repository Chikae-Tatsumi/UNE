library(ape)
library(Biostrings)
library(ggplot2)
library(ggtree)
library(ggnewscale)

# Import files
setwd("~/R/Analysis/2_UNE/ITS/indispecies")
tree <- read.tree("Phylogenetic_tree_28S_Newick")
tree

data <- read.csv("phylogenetic_info_28S.csv", header=T)
data <- data.frame(lapply(data, function(x){gsub(pattern="p__", replacement = "", x)}),stringsAsFactors = FALSE) 

# To change color by Phylum
groupInfo <- split(data$X,data$Phylum)
tree <- groupOTU(tree, groupInfo, group_name = "Phylum")

# Save
ggtree(tree, aes(color=Phylum), layout='circular')  
ggsave("tree_phylum.png", width=20,height=20)
