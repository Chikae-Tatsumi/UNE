# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv(file = "experimental_design.csv",header=T)
setwd("~/R/Analysis/2_UNE/16S")
Pros <- read.csv(file = "function/aggregated.function.table.csv",header=T, row.names=1)

# Prepare Cons
pathway <- read.table(file = "Tax4Fun2/pathway_prediction.txt",row.names=1,header=T,sep="\t",dec = ".",quote="")
ASV <- pathway[,1:(ncol(pathway)-3)]
description <- pathway[,(ncol(pathway)-2):ncol(pathway)]
ASV.t <- t(ASV)

# Aggregate
level3 <- aggregate(ASV, by=list(description$level3),FUN = sum,na.rm=F) 
rownames(level3) <- level3 [,1]
level3 <- level3 [,-1]
level3.t <- data.frame(t(level3))
level2 <- aggregate(ASV, by=list(description$level2),FUN = sum,na.rm=F) 
rownames(level2) <- level2 [,1]
level2 <- level2 [,-1]
level2.t <- data.frame(t(level2))

Cons <- cbind(level3.t$Human.Disease, level2.t$Xenobiotics.biodegradation.and.metabolism)
colnames(Cons) <- c("Human Disease", "Xenobiotics biodegradation and metabolism")

# Cor.test
Result <- NULL
for (i in 1:ncol(Pros)){
res <- NULL
for (j in 1:ncol(Cons)){
cor<- cor.test(Pros[,i],Cons[,j], method="p")
res <- c (res, cor$estimate,cor$p.value)}
Result <- rbind (Result, res)}

colnames (Result) <- c(paste(colnames(Cons)[1],".R"),paste(colnames(Cons)[1],".P"),paste(colnames(Cons)[2],".R"),paste(colnames(Cons[2]),".P"))
rownames(Result) <- colnames (Pros) 

# Save
write.csv(Result, "lm.tradeoff.16S.csv")
