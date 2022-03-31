library(ggplot2)
library(vegan)
library(ggrepel)

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv(file = "experimental_design.csv",header=T)
setwd("~/R/Analysis/2_UNE/ITS/FungalTrait")
objective <- read.csv("aggregated.fungaltrait.table.csv", header=T, row.names=1)



###### For Urban ######
bind <- cbind (objective,DESIGN)
subset <- subset(bind, bind$Urban=="Urban") 
subset.ASV <- subset[,1:(ncol(subset)-ncol(DESIGN))]
DESIGN.subset <- subset[,(ncol(subset)-ncol(DESIGN)+1):ncol(subset)]

NAME1 <- "DFE"
NAME2 <- "YearnestedbyDFE"

# ANOVA
Results <- NULL
for (i in 1:ncol(subset.ASV )){
data <- cbind(DESIGN.subset, subset.ASV [,i])
colnames (data)[ncol(data)] <- colnames(subset.ASV )[i]
data <- na.omit (data)
anova <- anova(lmer(scale(data[,ncol(data)])~ DFE/Year+(1|Site),data=data))

Fval <- anova[,5]
Pval <-anova[,6]
Bind <- c(Fval[1], Pval[1], Fval[2], Pval[2])
names(Bind) <- c(paste(NAME1,".F",sep=""),paste(NAME1,".P",sep=""),paste(NAME2,".F",sep=""),paste(NAME2,".P",sep=""))
Results <- rbind(Results, Bind)}

rownames(Results) <- colnames(subset.ASV )

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
write.csv(Results.asterisk, "lmer.ITS.Urban.Year.csv") 



###### For Rural ######
bind <- cbind (objective,DESIGN)
subset <- subset(bind, bind$Urban=="Rural") 
subset.ASV <- subset[,1:(ncol(subset)-ncol(DESIGN))]
DESIGN.subset <- subset[,(ncol(subset)-ncol(DESIGN)+1):ncol(subset)]

NAME1 <- "DFE"
NAME2 <- "YearnestedbyDFE"

# ANOVA
Results <- NULL
for (i in 1:ncol(subset.ASV )){
data <- cbind(DESIGN.subset, subset.ASV [,i])
colnames (data)[ncol(data)] <- colnames(subset.ASV )[i]
data <- na.omit (data)
anova <- anova(lmer(scale(data[,ncol(data)])~ DFE/Year+(1|Site),data=data))

Fval <- anova[,5]
Pval <-anova[,6]
Bind <- c(Fval[1], Pval[1], Fval[2], Pval[2])
names(Bind) <- c(paste(NAME1,".F",sep=""),paste(NAME1,".P",sep=""),paste(NAME2,".F",sep=""),paste(NAME2,".P",sep=""))
Results <- rbind(Results, Bind)}

rownames(Results) <- colnames(subset.ASV )

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
write.csv(Results.asterisk, "lmer.ITS.Rural.Year.csv") 



###### For Edge ######
bind <- cbind (objective,DESIGN)
subset <- subset(bind, bind$Edge=="Edge") 
subset.ASV <- subset[,1:(ncol(subset)-ncol(DESIGN))]
DESIGN.subset <- subset[,(ncol(subset)-ncol(DESIGN)+1):ncol(subset)]

NAME1 <- "DFB"
NAME2 <- "YearnestedbyDFB"

# ANOVA
Results <- NULL
for (i in 1:ncol(subset.ASV )){
data <- cbind(DESIGN.subset, subset.ASV [,i])
colnames (data)[ncol(data)] <- colnames(subset.ASV )[i]
data <- na.omit (data)
anova <- anova(lmer(scale(data[,ncol(data)])~ DFB/Year+(1|Site),data=data))

Fval <- anova[,5]
Pval <-anova[,6]
Bind <- c(Fval[1], Pval[1], Fval[2], Pval[2])
names(Bind) <- c(paste(NAME1,".F",sep=""),paste(NAME1,".P",sep=""),paste(NAME2,".F",sep=""),paste(NAME2,".P",sep=""))
Results <- rbind(Results, Bind)}

rownames(Results) <- colnames(subset.ASV )

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
write.csv(Results.asterisk, "lmer.ITS.Edge.Year.csv") 



###### For Interior ######
bind <- cbind (objective,DESIGN)
subset <- subset(bind, bind$Edge=="Interior") 
subset.ASV <- subset[,1:(ncol(subset)-ncol(DESIGN))]
DESIGN.subset <- subset[,(ncol(subset)-ncol(DESIGN)+1):ncol(subset)]

NAME1 <- "DFB"
NAME2 <- "YearnestedbyDFB"

# ANOVA
Results <- NULL
for (i in 1:ncol(subset.ASV )){
data <- cbind(DESIGN.subset, subset.ASV [,i])
colnames (data)[ncol(data)] <- colnames(subset.ASV )[i]
data <- na.omit (data)
anova <- anova(lmer(scale(data[,ncol(data)])~ DFB/Year+(1|Site),data=data))

Fval <- anova[,5]
Pval <-anova[,6]
Bind <- c(Fval[1], Pval[1], Fval[2], Pval[2])
names(Bind) <- c(paste(NAME1,".F",sep=""),paste(NAME1,".P",sep=""),paste(NAME2,".F",sep=""),paste(NAME2,".P",sep=""))
Results <- rbind(Results, Bind)}

rownames(Results) <- colnames(subset.ASV )

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
write.csv(Results.asterisk, "lmer.ITS.Interior.Year.csv") 
