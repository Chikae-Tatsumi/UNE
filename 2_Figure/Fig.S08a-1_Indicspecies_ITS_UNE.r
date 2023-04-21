library(indicspecies)
library(tidyr)
library(stringr)

# Import files & Format data
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv("experimental_design.csv",header=T)
for (i in 1: nrow(DESIGN)){
if (DESIGN$DFE[i]==0|DESIGN$DFE[i]==15) {DESIGN$Edge[i]<-"Edge"} else 
if (DESIGN$DFE[i]==30|DESIGN$DFE[i]==60|DESIGN$DFE[i]==90) {DESIGN$Edge[i]<-"Interior"} }

DESIGN$Group4 <- paste (DESIGN$Urban,DESIGN$Edge,sep=".")

setwd("~/R/Analysis/2_UNE/ITS")
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T, row.names=1)
ASV <- ASV.table [,1:(ncol(ASV.table)-7)] 
taxonomy <- ASV.table [,(ncol(ASV.table)-6):ncol(ASV.table)] 
ASV.t <- t(ASV)
percent <- ASV / mean(colSums(ASV)) *100
percent.t <- t(percent)

setwd("indispecies")

# Run indicspecies
inv = multipatt(ASV.t, cluster=factor(DESIGN$Group4)
,func="r.g",duleg=TRUE)

# Save information
options(max.print=100000)
sink("summary.inv.Group4.csv")  
summary(inv)
sink()     

str <- inv$str
write.csv(str,"correlation.inv.Group4.csv")

# Get indicator ASVs of each location group
summary <- read.csv("summary.inv.Group4.csv")
rows <- grep (" Group", summary$Multilevel.pattern.analysis) 

Rural.Edge <- data.frame(summary[(rows[1]+1):(rows[2]-1),])
Rural.Interior <- data.frame(summary[(rows[2]+1):(rows[3]-1),])
Urban.Edge <- data.frame(summary[(rows[3]+1):(rows[4]-1),])
Urban.Interior <- data.frame(summary[(rows[4]+1):(nrow(summary)-2),])

NAME <- c("Rural.Edge","Rural.Interior","Urban.Edge", "Urban.Interior")

sample.list <- list()
sample.list [[1]] <- Rural.Edge
sample.list [[2]] <- Rural.Interior
sample.list [[3]] <- Urban.Edge
sample.list [[4]] <- Urban.Interior

# Save
for (i in 1:4){
sample.list [[i]] <- str_split_fixed(sample.list [[i]][2:nrow(sample.list [[i]]),], ' 0.',3)
rownames(sample.list [[i]]) <- sample.list [[i]][,1]
rownames(sample.list [[i]]) <- gsub(" ","",rownames(sample.list [[i]]))
sample.list [[i]][,1] <- rownames(sample.list [[i]])
colnames(sample.list [[i]]) <- c("ASV_ID","stat","p.value")
write.csv(sample.list [[i]], paste (NAME[i],".csv", sep=""),row.names=F)
}
