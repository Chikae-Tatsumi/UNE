library(ggplot2)
library(reshape2)
library(wesanderson)

# Import files
setwd("~/R/Analysis/2_UNE")
METADATA <- read.csv(file = "metadata.csv",header=T)
DESIGN <- read.csv(file = "experimental_design.csv",header=T)
func.16S <- read.csv (file = "~/R/Analysis/2_UNE/16S/function/aggregated.function.table.csv",row.names=1,header=T) 
func.ITS <- read.csv (file = "~/R/Analysis/2_UNE/ITS/FungalTrait/aggregated.fungaltrait.table.csv",row.names=1,header=T) 
qPCR <- read.csv (file = "~/R/Analysis/2_UNE/qPCR/qPCR.calc.result.UNE.csv",row.names=1,header=T) 
fun.div <- read.csv (file = "~/R/Analysis/2_UNE/ITS/diversity.csv",row.names=1,header=T) 
bac.div <- read.csv (file = "~/R/Analysis/2_UNE/16S/diversity.csv",row.names=1,header=T) 
network <- read.csv (file = "~/R/Analysis/2_UNE/Network/05_genus_combForest_227samples.csv",row.names=1,header=T) 
setwd("~/R/Analysis/2_UNE/Others")

# Make dataset
func <- cbind (
bac.div$shannon, fun.div$shannon, 
qPCR$ITS.copies.gsoil/qPCR$X16S.copies.gsoil,
log10(qPCR$X16S.copies.gsoil+1),log10(qPCR$ITS.copies.gsoil+1), 
log10(qPCR$nosZ.copies.gsoil+1),log10(qPCR$amoA.copies.gsoil+1),log10(qPCR$nifH.copies.gsoil+1),
func.16S$Dissim_nitrate_reduction, func.16S$Assim_nitrate_reduction, 
func.16S$Chitinolytic, func.16S$Methanotroph, 
func.16S$Lignolytic,func.16S$Cellulolytic,
func.16S$Oligotroph, func.16S$Copiotroph, 
func.ITS$Saprotroph, log10(qPCR$AM.copies.gsoil+1), func.ITS$ECM
)
colnames(func) <- c(
"Bacterial diversity", "Fungal diversity",
"Fungi: bacteria ratio",
"Total bacterial abundance", "Total fungal abundance", 
"Denitrifying bacteria","Nitrifying bacteria", "N fixing bacteria",
"Dissimiratory nitrate reducing bacteria", "Assimiratory nitrate reducing bacteria",
"Chitinolytic bacteria", "Methanotroph bacteria",
"Lignolytic bacteria", "Cellulolytic bacteria",
"Oligotroph bacteria", "Copiotroph bacteria",
"Saprotroph fungi","AM fungi", "ECM fungi")

metadata <- cbind(METADATA$BA,METADATA$MR,METADATA$AR,METADATA$NR,METADATA$Soil.C,METADATA$Soil.CN,
METADATA$Respiration,METADATA$leaf.C_copied,METADATA$leaf.N_copied,
METADATA$Moist,METADATA$pH,METADATA$SOM,METADATA$Temp,
METADATA$NH4,METADATA$NO3)
colnames(metadata) <- c("Basal area", "N mineralization","Ammonification", "Nitrification","Soil C","Soil C:N", 
"Soil respiration", "Foliar C","Foliar N",
"Soil moisture","Soil pH","Soil organic matter","Soil temperature",
"Soil NH4+","Soil NO3-")

# Make coefficient matrix
Result <- NULL
for (i in 1:ncol(func)){
estimate <- NULL
for (j in 1:ncol(metadata)){
cor<- cor.test(func[,i],metadata[,j], method="p")
estimate <- c (estimate, cor$estimate)}
Result <- rbind (Result, estimate)}
colnames (Result) <- colnames (metadata)
rownames(Result) <- colnames (func) 

# Make Pval matrix
Result.Pval <- NULL
for (i in 1:ncol(func)){
Pval <- NULL
for (j in 1:ncol(metadata)){
cor<- cor.test(func[,i],metadata[,j], method="p")
Pval <- c (Pval, cor$p.value)}
Result.Pval <- rbind (Result.Pval, Pval)}
colnames (Result.Pval) <- colnames (metadata)
rownames(Result.Pval) <- colnames (func) 

# For network properties
# Make dataset
func <- cbind (
network$G_assort, network$K_assort,
network$Degree,network$Betweenness, network$Diameter, network$Short_Path, 
network$Transitivity, network$Complexity, network$Density, network$Edges)
colnames(func) <- c(
"Guild assortativity", "Kingdom assortativity",
"Degree centrality", "Betweenness centrality", "Diameter", "Shortest path length",
"Transitivity", "Complexity", "Edge density", "The number of Edges")

metadata <- cbind(METADATA$BA,METADATA$MR,METADATA$AR,METADATA$NR,METADATA$Soil.C,METADATA$Soil.CN,
METADATA$Respiration,METADATA$leaf.C,METADATA$leaf.N,
METADATA$Moist,METADATA$pH,METADATA$SOM,METADATA$Temp,
METADATA$NH4,METADATA$NO3)
colnames(metadata) <- c("Basal area", "N mineralization","Ammonification", "Nitrification","Soil C","Soil C:N", 
"Soil respiration", "Foliar C","Foliar N",
"Soil moisture","Soil pH","Soil organic matter","Soil temperature",
"Soil NH4+","Soil NO3-")

# Make coefficient matrix
Result.net <- NULL
for (i in 1:ncol(func)){
estimate <- NULL
for (j in 1:ncol(metadata)){
cor<- cor.test(func[,i],metadata[,j], method="p")
estimate <- c (estimate, cor$estimate)}
Result.net <- rbind (Result.net, estimate)}
colnames (Result.net) <- colnames (metadata)
rownames(Result.net) <- colnames (func) 

# Make Pval matrix
Result.Pval.net <- NULL
for (i in 1:ncol(func)){
Pval <- NULL
for (j in 1:ncol(metadata)){
cor<- cor.test(func[,i],metadata[,j], method="p")
Pval <- c (Pval, cor$p.value)}
Result.Pval.net <- rbind (Result.Pval.net, Pval)}
colnames (Result.Pval.net) <- colnames (metadata)
rownames(Result.Pval.net) <- colnames (func) 

# Make Asterisk matrix
Result.Pval.combined <- rbind(Result.Pval.net, Result.Pval)
Result.asterisk <- Result.Pval.combined
for (k in 1:ncol(Result.Pval.combined)){
for (i in 1:nrow(Result.Pval.combined)){
if (Result.Pval.combined[i,k] < 0.001) {Result.asterisk[i,k] <- "***"
} else if (Result.Pval.combined[i,k] < 0.01) {Result.asterisk[i,k] <- "**"
} else if (Result.Pval.combined[i,k] < 0.05) {Result.asterisk[i,k] <- "*"
} else {Result.asterisk[i,k]<- " "}}
}

# Combine netowrk propeties + others
Result.combined <- rbind(Result.net, Result)
melt <- melt(Result.combined)
melt.asterisk <- melt(Result.asterisk)

# Visualize
vdata <- cbind(melt, melt.asterisk$value)
colnames(vdata) <- c("microbe","environment","coef","Pval")
pal <- wes_palette("Zissou1", 100, type = "continuous")

ggplot (vdata)+
geom_tile(aes(x=environment, y=microbe, fill =coef))+
scale_fill_gradientn(colours = pal, name = "Correlation \ncoefficient")+
geom_text(aes(x=environment, y=microbe, label = Pval))+
labs(x="",y="")+
theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.3)) 

ggsave(file = "Heatmap.png", height=6, width=8)
