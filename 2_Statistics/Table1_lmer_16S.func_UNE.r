library (lme4)
library (lmerTest)

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv("experimental_design.csv")
setwd("~/R/Analysis/2_UNE/16S/function")  
objective <- read.csv("aggregated.function.table.csv", header=T, row.names=1) 

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
write.csv(Results.asterisk, "lmer.16S.function.csv") 
