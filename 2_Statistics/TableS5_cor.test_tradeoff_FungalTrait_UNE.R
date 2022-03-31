# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv(file = "experimental_design.csv",header=T)
setwd("~/R/Analysis/2_UNE/ITS/FungalTrait")
data <- read.csv(file = "aggregated.fungaltrait.table.csv",header=T, row.names=1)

Pros <- cbind(data$ECM.fungi, data$Saprotroph.fungi, 
data$Soil.saportroph.fungi, data$Wood.saprotroph.fungi, data$Litter.saprotroph.fungi, data$Dung.saprotroph.fungi)
Cons <- cbind(data$Plant.pathogenic.capacity,	data$Animal.parasite)
colnames(Pros) <- c("ECM fungi","Saprotroph fungi",	"Soil saportroph fungi",
"Wood saprotroph fungi","Litter saprotroph fungi","Dung saprotroph fungi")
colnames(Cons) <- c("Plant pathogenic capacity",	"Animal parasite")

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
write.csv(Result, "lm.tradeoff.FungalTrait.csv")
