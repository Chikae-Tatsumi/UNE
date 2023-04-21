# Import files
setwd("~/R/Analysis/2_UNE")
METADATA <- read.csv(file = "metadata.csv",header=T)
DESIGN <- read.csv(file = "experimental_design.csv",header=T)
fungaltrait.table <- read.csv(file="~/R/Analysis/2_UNE/ITS/FungalTrait/rarefied_fungaltrait.table.csv",header=T,row.names = 1)
function.table <- read.csv(file="~/R/Analysis/2_UNE/16S/function/rarefied_16S_function.csv",row.names = 1,header=T)
func.16S <- read.csv (file = "~/R/Analysis/2_UNE/16S/function/aggregated.function.table.csv",row.names=1,header=T) 
func.ITS <- read.csv (file = "~/R/Analysis/2_UNE/ITS/FungalTrait/aggregated.fungaltrait.table.csv",row.names=1,header=T) 
qPCR <- read.csv (file = "~/R/Analysis/2_UNE/qPCR/qPCR.calc.result.UNE.csv",row.names=1,header=T) 

# Quantify "Other" Fungi
ASV <- fungaltrait.table [,1:(ncol(fungaltrait.table)-31)] 
guild <- fungaltrait.table [,(ncol(fungaltrait.table)-24):ncol(fungaltrait.table)] 
percent <- ASV / mean(colSums(ASV)) *100

aggregated <- aggregate(percent, by=list(guild$primary_lifestyle),FUN = sum,na.rm=F) 
row.names(aggregated)<-aggregated[,1]
aggregated <- aggregated[,-1]
aggregated <- data.frame(aggregated)
Other_fungi <- 100 - colSums(aggregated) + aggregated[1,]
Other_fungi <- data.frame(t(Other_fungi))

func.ITS$Other_fungi <- Other_fungi[,1]
colnames(func.ITS)[ncol(func.ITS)] <- "Other_fungi"

# Quantify "Other" Bacteria
ASV <- function.table [,1:(ncol(function.table)-19)] 
fgs <- function.table [,(ncol(function.table)-12):ncol(function.table)] 
percent <- ASV / mean(colSums(ASV)) *100

fgs$Other <- paste(fgs[,1],fgs[,2],fgs[,3],fgs[,4],fgs[,5],fgs[,6],fgs[,7],
fgs[,8],fgs[,9],fgs[,10],fgs[,11],fgs[,12],fgs[,13],sep="")

subset <- subset(percent, fgs$Other=="otherotherotherotherotherotherotherotherotherotherotherotherother")
Other_bacteria <- data.frame(colSums(subset))

func.16S$Other_bacteria <- Other_bacteria[,1]
colnames(func.16S)[ncol(func.16S)] <- "Other_bacteria"

# Get absolute abundance
qPCR.func.ITS.full <- qPCR$ITS.copies.gsoil * func.ITS/100
qPCR.func.16S.full <- qPCR$X16S.copies.gsoil * func.16S/100
qPCR.nifamonos <- cbind(qPCR$nifH.copies.gsoil, qPCR$amoA.copies.gsoil,qPCR$nosZ.copies.gsoil)
colnames(qPCR.nifamonos) <- c("N-fixation","Nitrification","Denitrification")

# Format data
qPCR.func.ITS <- cbind(qPCR.func.ITS.full$ECM, qPCR.func.ITS.full$Plant_pathogenic_capacity,qPCR.func.ITS.full$Animal_parasite,
qPCR.func.ITS.full$soil_saprotroph,qPCR.func.ITS.full$wood_saprotroph,qPCR.func.ITS.full$litter_saprotroph, qPCR.func.ITS.full$Other_fungi)
colnames(qPCR.func.ITS) <- c("ECM", "Plant_pathogen","Animal_parasite",
"soil_saprotroph","wood_saprotroph","litter_saprotroph", "Other_fungi")

qPCR.func.16S <- cbind(qPCR.func.16S.full$Copiotroph, qPCR.func.16S.full$Oligotroph,
qPCR.func.16S.full$Cellulolytic, qPCR.func.16S.full$Lignolytic, qPCR.func.16S.full$Methanotroph,
qPCR.func.16S.full$Chitinolytic, qPCR.func.16S.full$Assim_nitrate_reduction, qPCR.func.16S.full$Dissim_nitrate_reduction, qPCR.func.16S.full$Other_bacteria)
colnames(qPCR.func.16S) <- c("Copiotroph", "Oligotroph",
"Cellulolytic","Lignolytic","Methanotroph",
"Chitinolytic","Assim_nitrate_reduction","Dissim_nitrate_reduction", "Other_bacteria")

data <- cbind(qPCR.func.ITS,qPCR.func.16S,qPCR.nifamonos)

# Save
write.csv(data, "Network/function.absolute.abundance.csv")