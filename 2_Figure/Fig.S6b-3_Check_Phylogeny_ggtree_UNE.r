library(ape)
library(Biostrings)
library(ggplot2)
library(ggtree)
library(ggnewscale)

# Import files
setwd("~/R/Analysis/2_UNE/16S/indispecies")
tree <- read.tree("Phylogenetic_tree_Newick_16S_NJ")
tree

data <- read.csv("phylogenetic_info.csv", header=T)
data <- data.frame(lapply(data, function(x){gsub(pattern="c__", replacement = "", x)}),stringsAsFactors = FALSE) # Change p__ --> the first alphabet of the level you analyze

# To change color by Phylum
groupInfo <- split(data$X,data$Phylum)

# To pick up frequent groups
list <- NULL
for (i in 1:length(groupInfo)){
if (length(groupInfo[[i]]) > 20){list[[length(list)+1]] <- groupInfo[[i]]}}
names <- NULL
for (i in 1:length(groupInfo)){
if (length(groupInfo[[i]]) > 20){names <- c(names,names(groupInfo[i]))}}
names(list) <- names

tree <- groupOTU(tree, list,group_name = "Phylum")

# Save
ggtree(tree, aes(color=Phylum), layout='circular')  
ggsave("tree_Phylum.png", width=20,height=20)
