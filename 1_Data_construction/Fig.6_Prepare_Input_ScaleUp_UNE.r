# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv("experimental_design.csv",header=T)
fungaltrait <- read.csv (file = "~/R/Analysis/2_UNE/ITS/FungalTrait/aggregated.fungaltrait.table.csv",row.names=1,header=T) 
func.16S <- read.csv (file = "~/R/Analysis/2_UNE/16S/function/aggregated.function.table.csv",row.names=1,header=T) 
qPCR <- read.csv (file = "~/R/Analysis/2_UNE/qPCR/qPCR.calc.result.UNE.csv",row.names=1,header=T) 

# Assign forest type to UNE site
# See the forest code: https://www.fia.fs.usda.gov/library/database-documentation/#FIADB
DESIGN$Forest_type <- NA
for (i in 1: nrow(DESIGN)){
if (DESIGN$Site[i]=="AA01") {DESIGN$Forest_type[i]<-"100"} else 
if (DESIGN$Site[i]=="HF06") {DESIGN$Forest_type[i]<-"400"} else 
if (DESIGN$Site[i]=="BH02"|DESIGN$Site[i]=="HW07"|DESIGN$Site[i]=="SW08") {DESIGN$Forest_type[i]<-"500"} else 
if (DESIGN$Site[i]=="BM03"|DESIGN$Site[i]=="HF04"|DESIGN$Site[i]=="HF05") {DESIGN$Forest_type[i]<-"800"} }

# Format data
data <- cbind(fungaltrait$ECM, fungaltrait$Plant_pathogenic_capacity, fungaltrait$Animal_parasite,(qPCR$nosZ.copies.gsoil+1)/(qPCR$X16S.copies.gsoil+1)*100)
colnames(data) <- c("ECM_fungi","Plant_pathogenic_fungi","Animal_parasite_fungi","Denitifying_bacteria")

aggregated <- aggregate(data, by=list(DESIGN$DFE, DESIGN$Forest_type),FUN = mean,na.rm=F) 
colnames(aggregated)[1:2] <- c("Distance_from_edge","Forest_type")

# Save
setwd("~/R/Analysis/2_UNE/FIA")
write.csv(aggregated, "Aggregated_function_by_foresttype_DFE.csv", row.names=F)