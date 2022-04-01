library(ggplot2)
library(lme4)
library(lmerTest)

setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv("experimental_design.csv",header=T)
setwd("~/R/Analysis/2_UNE/ITS/FungalTrait")

fungaltrait.table <- read.csv(file="rarefied_fungaltrait.table.csv",header=T,row.names = 1)
ASV <- fungaltrait.table [,1:(ncol(fungaltrait.table)-31)] 
guild <- fungaltrait.table [,(ncol(fungaltrait.table)-24):ncol(fungaltrait.table)] 
percent <- ASV / mean(colSums(ASV)) *100
guild$lifestyle <- paste(guild$primary_lifestyle, "_",guild$secondary_lifestyle)

# For ECM
aggregated <- aggregate(percent, by=list(guild$lifestyle),FUN = sum,na.rm=F) 
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
pat.cap <- cbind(percent, guild$Plant_pathogenic_capacity_template)
colnames(pat.cap)[ncol(pat.cap)] <- "Plant_pathogenic_capacity_template"
pat.cap <- na.omit(pat.cap)
pat.cap[pat.cap$Plant_pathogenic_capacity_template == "",] <-NA
pat.cap <- na.omit (pat.cap)
colSums <- colSums(pat.cap[,-ncol(pat.cap)])
Plant_pathogenic_capacity <- data.frame(colSums)
colnames (Plant_pathogenic_capacity) [1] <- "Plant_pathogenic_capacity"

# Animal_biotrophic_capacity
pat.cap <- cbind(percent, guild$Animal_biotrophic_capacity_template)
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
write.csv(Result, "aggregated.fungaltrait.table.csv")
colnames(Result) <- c("ECM fungi","Saprotroph fungi", "Pathotroph fungi", "Plant pathogenic capacity", "Animal parasite", "Soil saportroph fungi", "Wood saprotroph fungi", "Litter saprotroph fungi", "Dung saprotroph fungi")

data <- cbind (Result, DESIGN)
data$Urban <- factor (data$Urban, levels=c("Urban","Rural"))
for (i in 1:ncol(Result)){

# To show the result of lmer in the graph
anova <- anova(lmer(scale(data[,i])~ scale(DFB)*scale(DFE)+(1|Site),data=data))
Pval <- anova[,6]
if (Pval[1] > 0.05) {DFB.result <- ""
} else if (Pval[1] > 0.01) {DFB.result <- "        Urbanization *      "
} else if (Pval[1] > 0.001) {DFB.result <- "        Urbanization **     "
} else {DFB.result <- "        Urbanization ***    "}
if (Pval[2] > 0.05) {DFE.result <- ""
}else if (Pval[2] > 0.01) {DFE.result <- "        Fragmentation *       "
}else if (Pval[2] > 0.001) {DFE.result <- "        Fragmentation **     "
}else {DFE.result <- "        Fragmentation ***    "}
if (Pval[3] > 0.05) {DFBDFE.result <- ""
}else if (Pval[3] > 0.01) {DFBDFE.result <- "        U×F *      "
}else if (Pval[3] > 0.001) {DFBDFE.result <- "        U×F **     "
}else {DFBDFE.result <- "        U×F ***    "}

ggplot(data)+
geom_point(aes(x=DFE, y=data[,i], color=Urban),position=position_jitter( width=2, height=0))+ 
geom_smooth(method="lm", aes(x=DFE, y=data[,i], group=Urban, color=Urban))+  
# scale_fill_manual(values = c("#C77CFF","#7CAE00","#00BFC4","#F8766D"))+  # if you want to change the colors
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
labs (y=paste(colnames(data)[i], "(%)",sep = " "), x="Distance from edge (m)") + # if you want to change the axis titles
annotate("text", x=-Inf, y=Inf, hjust=0, vjust =2, label=DFB.result, size=5) +
annotate("text", x=-Inf, y=Inf, hjust=0, vjust =4, label=DFE.result, size=5) +
annotate("text", x=-Inf, y=Inf, hjust=0, vjust =6, label=DFBDFE.result, size=5) 

ggsave(paste(colnames(data)[i],".png"),width = 5, height = 4)}
