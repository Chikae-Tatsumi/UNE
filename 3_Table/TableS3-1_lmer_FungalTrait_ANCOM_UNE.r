library (lme4)
library (lmerTest)

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv("experimental_design.csv",header=T)
setwd("~/R/Analysis/2_UNE/ITS")
ASV.table <- read.table(file="ancom_ASV_table.txt",header=T)
DATABASE <- read.csv(file="~/R/Database/FungalTrait.csv",header=T)
setwd("~/R/Analysis/2_UNE/ITS/FungalTrait")

# FungalTrait
ASV <- ASV.table [,1:(ncol(ASV.table)-7)] 
taxonomy <- ASV.table [,(ncol(ASV.table)-6):ncol(ASV.table)]
taxonomy <- data.frame(lapply(taxonomy, function(x){gsub(pattern="g__", replacement = "", x)}),stringsAsFactors = FALSE)
DATABASE$Genus <- DATABASE$GENUS
taxonomy$id  <- 1:nrow(taxonomy)

merge <- merge(taxonomy, DATABASE, by="Genus", all.x=TRUE)
merge.tidy <- merge[order(merge$id), ]
merge.tidy <- merge.tidy[,-2:-(ncol(taxonomy))]
DATABASE <- DATABASE[,-ncol(DATABASE)]
colnames (merge.tidy) [2:ncol(merge.tidy)] <- colnames (DATABASE)
rownames(merge.tidy) <- rownames (ASV.table)

fungaltrait.table <- cbind(ASV, taxonomy[,-ncol(taxonomy)], merge.tidy[,-1])
rownames(fungaltrait.table) <- rownames (ASV)
write.csv(fungaltrait.table, "fungaltrait.table_ancom.csv")

# Aggregate fungaltrait table
ASV <- fungaltrait.table [,1:(ncol(fungaltrait.table)-31)] 
guild <- fungaltrait.table [,(ncol(fungaltrait.table)-24):ncol(fungaltrait.table)] 
guild$lifestyle <- paste(guild$primary_lifestyle, "_",guild$secondary_lifestyle)

# For ECM
aggregated <- aggregate(ASV, by=list(guild$lifestyle),FUN = sum,na.rm=F) 
row.names(aggregated)<-aggregated[,1]
aggregated <- aggregated[,-1]
aggregated <- data.frame(aggregated)
rows <- grep ("ectomycorrhizal", rownames (aggregated)) 
subset <- aggregated[rows,]
subset.t <- t(subset)
ECM <- data.frame(rowSums (subset.t))
colnames (ECM) [1] <- "ECM"

# For Saprotroph
rows <- grep ("saprotroph", rownames (aggregated)) 
subset <- aggregated[rows,]
subset.t <- t(subset)
Saprotroph <- data.frame(rowSums (subset.t))
colnames (Saprotroph) [1] <- "Saprotroph" 

# For Pathotroph
rows <- grep ("patho", rownames (aggregated)) 
subset <- aggregated[rows,]
subset.t <- t(subset)
Pathotroph <- data.frame(rowSums (subset.t))
colnames (Pathotroph) [1] <- "Pathotroph"

# Plant ptatogenic capacity
pat.cap <- cbind(ASV, guild$Plant_pathogenic_capacity_template)
colnames(pat.cap)[ncol(pat.cap)] <- "Plant_pathogenic_capacity_template"
pat.cap <- na.omit(pat.cap)
pat.cap[pat.cap$Plant_pathogenic_capacity_template == "",] <-NA
pat.cap <- na.omit (pat.cap)
colSums <- colSums(pat.cap[,-ncol(pat.cap)])
Plant_pathogenic_capacity <- data.frame(colSums)
colnames (Plant_pathogenic_capacity) [1] <- "Plant_pathogenic_capacity"

# Animal_biotrophic_capacity
pat.cap <- cbind(ASV, guild$Animal_biotrophic_capacity_template)
colnames(pat.cap)[ncol(pat.cap)] <- "Animal_biotrophic_capacity_template"
pat.cap <- na.omit(pat.cap)
pat.cap[pat.cap$Animal_biotrophic_capacity_template == "",] <-NA
pat.cap <- na.omit (pat.cap)
rows <- grep ("parasite", pat.cap$Animal_biotrophic_capacity_template)
subset <- pat.cap[rows,]
subset <- subset[,-ncol(subset)]
Animal_parasite <- data.frame(colSums (subset))
colnames (Animal_parasite) [1] <- "Animal_parasite"

# For soil saprotroph
rows <- grep ("soil_saprotroph", rownames (aggregated)) 
subset <- aggregated[rows,]
subset.t <- t(subset)
soil_saprotroph <- data.frame(rowSums (subset.t))
colnames (soil_saprotroph) [1] <- "soil_saprotroph"

# For wood saprotroph
rows <- grep ("wood_saprotroph", rownames (aggregated)) 
subset <- aggregated[rows,]
subset.t <- t(subset)
wood_saprotroph <- data.frame(rowSums (subset.t))
colnames (wood_saprotroph) [1] <- "wood_saprotroph"

# For litter saprotroph
rows <- grep ("litter_saprotroph", rownames (aggregated)) 
subset <- aggregated[rows,]
subset.t <- t(subset)
litter_saprotroph <- data.frame(rowSums (subset.t))
colnames (litter_saprotroph) [1] <- "litter_saprotroph"

# For dung saprotroph
rows <- grep ("dung_saprotroph", rownames (aggregated)) 
subset <- aggregated[rows,]
subset.t <- t(subset)
dung_saprotroph <- data.frame(rowSums (subset.t))
colnames (dung_saprotroph) [1] <- "dung_saprotroph"

# Summarize
Result <- cbind (ECM, Saprotroph, Pathotroph,Plant_pathogenic_capacity, Animal_parasite,soil_saprotroph, wood_saprotroph, litter_saprotroph, dung_saprotroph)
write.csv(Result, "aggregated.fungaltrait.table_ancom.csv")

# Run lmer
objective <- Result
NAME1 <- "Urbanization"
NAME2 <- "Fragmentation"

objective <- objective[,-ncol(objective)] # Becasue of no dung-saprotoph ASVs

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
write.csv(Results.asterisk, "lmer.ITS.FungalTrait_ancom.csv") 
