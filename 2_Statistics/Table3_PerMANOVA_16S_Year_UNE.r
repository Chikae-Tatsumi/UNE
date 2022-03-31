library(ggplot2)
library(vegan)
library(ggrepel)

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv(file = "experimental_design.csv",header=T)
setwd("~/R/Analysis/2_UNE/16S")
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T)
dir.create("effect_of_year")

ASV <- ASV.table [,1:(ncol(ASV.table)-6)] 
taxonomy <- ASV.table [,(ncol(ASV.table)-5):ncol(ASV.table)]
ASV.t <- t (ASV)



###### For Urban ######
bind <- cbind (ASV.t,DESIGN)
subset <- subset(bind, bind$Urban=="Urban") 
subset.ASV <- subset[,1:(ncol(subset)-ncol(DESIGN))]
DESIGN.subset <- subset[,(ncol(subset)-ncol(DESIGN)+1):ncol(subset)]

sink("effect_of_year/Urban.PerMANOVA_YearnestedbyDFE.txt")
# PerMANOVA for Urban
adonis(subset.ASV ~ DFE/Year,  data=DESIGN.subset, 
strata=DESIGN.subset$Site, 
permutations=10000)
sink()



###### For Rural ######
bind <- cbind (ASV.t,DESIGN)
subset <- subset(bind, bind$Urban=="Rural") 
subset.ASV <- subset[,1:(ncol(subset)-ncol(DESIGN))]
DESIGN.subset <- subset[,(ncol(subset)-ncol(DESIGN)+1):ncol(subset)]

sink("effect_of_year/Rural.PerMANOVA_YearnestedbyDFE.txt")
# PerMANOVA for Rural
adonis(subset.ASV ~ DFE/Year,  data=DESIGN.subset, 
strata=DESIGN.subset$Site, 
permutations=10000)
sink()



###### For Edge ######
bind <- cbind (ASV.t,DESIGN)
subset <- subset(bind, bind$Edge=="Edge") 
subset.ASV <- subset[,1:(ncol(subset)-ncol(DESIGN))]
DESIGN.subset <- subset[,(ncol(subset)-ncol(DESIGN)+1):ncol(subset)]

sink("effect_of_year/Edge.PerMANOVA_YearnestedbyDFB.txt")
# PerMANOVA for Edge
adonis(subset.ASV ~ DFB/Year,  data=DESIGN.subset, 
strata=DESIGN.subset$Site, 
permutations=10000)
sink()



###### For Interior ######
bind <- cbind (percent.t,DESIGN)
subset <- subset(bind, bind$Edge=="Interior") 
subset.ASV <- subset[,1:(ncol(subset)-ncol(DESIGN))]
DESIGN.subset <- subset[,(ncol(subset)-ncol(DESIGN)+1):ncol(subset)]

sink("effect_of_year/Interior.PerMANOVA_YearnestedbyDFB.txt")
# PerMANOVA for Interior
adonis(subset.ASV ~ DFB/Year,  data=DESIGN.subset, 
strata=DESIGN.subset$Site, 
permutations=10000)
sink()
