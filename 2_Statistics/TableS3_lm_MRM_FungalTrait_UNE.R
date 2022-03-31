library(ecodist)

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv(file = "experimental_design.csv",header=T)
objective <- read.csv(file = "ITS/FungalTrait/aggregated.fungaltrait.table.csv",header=T, row.names=1)
setwd("~/R/Analysis/2_UNE/Autocorrelation")
distance.matrix <- read.csv("dd.distances.csv",header=T, row.names=1)

NAME1 <- "Urbanization"
NAME2 <- "Fragmentation"
NAME3 <- "Spatial.distance"

# Calculation
Results <- NULL
for (i in 1:ncol(objective)){
MRM <- lm(dist(objective[,i])~ dist(DESIGN$DFB)*dist(DESIGN$DFE)+as.dist(distance.matrix))
anova <- anova(MRM)
Fval <- anova[,4]
Pval <-anova[,5]
Bind <- c(Fval[1], Pval[1], Fval[2], Pval[2], Fval[4], Pval[4], Fval[3], Pval[3])
Results <- rbind(Results, Bind)}
colnames(Results) <- c(paste(NAME1,".F",sep=""),paste(NAME1,".P",sep=""),paste(NAME2,".F",sep=""),paste(NAME2,".P",sep=""),paste(NAME1,"*",NAME2,".F",sep=""),paste(NAME1,"*",NAME2,".P",sep=""),paste(NAME3,".F",sep=""),paste(NAME3,".P",sep=""))
rownames(Results) <- colnames(METADATA)

# P.val --> Asterisk
Results.asterisk <- Results
for (k in 1:ncol(Results)){
if(k%%2 == 0){
for (i in 1:nrow(Results)){
if (Results[i,k] < 0.001) {Results.asterisk[i,k] <- "***"
} else if (Results[i,k] < 0.01) {Results.asterisk[i,k] <- "**"
} else if (Results[i,k] < 0.05) {Results.asterisk[i,k] <- "*"
} else {Results.asterisk[i,k]<- "n.s"}}
}}

# Save
write.csv(Results.asterisk, "lm.MRM.FungalTrait.csv")

