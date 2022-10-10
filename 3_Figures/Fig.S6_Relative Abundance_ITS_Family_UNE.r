library(ggplot2)
library(stringr)
library(stringi)

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv("experimental_design.csv",header=T)
setwd("~/R/Analysis/2_UNE/ITS")
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T)

# % table
ASV <- ASV.table [,1:(ncol(ASV.table)-7)] 
taxonomy <- ASV.table [,(ncol(ASV.table)-6):ncol(ASV.table)] 
percent <- ASV / mean(colSums(ASV)) *100
# Remove "p__" before phylum name
taxonomy <- data.frame(lapply(taxonomy, function(x){gsub(pattern="f__", replacement = "", x)}),stringsAsFactors = FALSE)

# Aggregate
agrregated <- aggregate(percent, by=list(taxonomy$Family),FUN = sum,na.rm=F)
row.names(agrregated)<-agrregated[,1]
agrregated <- agrregated[,-1]
agrregated <- data.frame(agrregated)
rowMeans <- rowMeans(agrregated)
agrregated <- cbind(agrregated,rowMeans)

# >2% abund 
majors <- agrregated[agrregated[,"rowMeans"] > 2,]
majors <- majors[order(majors$rowMeans,decreasing = T),]
selected <- majors
selected <- selected[,-ncol(selected)] 

# Make dataset
selected.t <- t (selected)
write.csv(selected.t, "aggregated.family.table.csv")
bind <- cbind (selected.t, DESIGN)
category <- paste (bind$Urban, bind$DFE,sep="&") 
bind <- cbind (bind, category)
table <- aggregate(selected.t, by=list(bind$category),FUN = mean) 
split <- str_split (table[,1], "&",simplify = TRUE)
colnames(split) <- c("Urban","DFE")
table <- cbind (split,table[2:(ncol(table))])
Urban <- as.character(table[,1]) 
DFE <- as.character(table[,2]) 
table <- data.frame(table)
data <- data.frame()
for (i in 3:(ncol(table))){
Abundance <-table[,i]
Example <- colnames (table)[i]
bind <- cbind(Urban, DFE, Abundance, Example) 
data<-rbind (data,bind)}
data$Abundance <- as.numeric(as.character(data$Abundance))

rownames<-rownames(selected)[nrow(selected):1]
data$Example <- factor(data$Example, levels = rownames)
data$Urban <- factor (data$Urban, levels=c("Urban","Rural"))
data <- na.omit(data)

# ggplot
ggplot (data,  mapping=aes(x=DFE, y=Abundance, fill=Example))+ 
geom_bar(aes(), stat="identity", position="stack",color="black")+
theme_classic()+
theme(text=element_text(size=12,color="black"),axis.text=element_text(size=12,color="black"))+
labs (x="Distance from edge (m)",y="Abundance (%)",fill = "")+
facet_wrap(~Urban,ncol=5)

# Save
ggsave(file = "Rel.Abund.Family.ITS.png", width = 8,height=5)
