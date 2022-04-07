library(ggplot2)
library(lme4)
library(lmerTest)

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv("experimental_design.csv",header=T)
METADATA <- read.csv("metadata.csv",header=T)
setwd("~/R/Analysis/2_UNE/ITS/FungalTrait")

fungaltrait.table <- read.csv(file="rarefied_fungaltrait.table.csv",header=T,row.names = 1)
ASV <- fungaltrait.table [,1:(ncol(fungaltrait.table)-31)] 
guild <- fungaltrait.table [,(ncol(fungaltrait.table)-24):ncol(fungaltrait.table)] 
percent <- ASV / mean(colSums(ASV)) *100
guild$lifestyle <- paste(guild$primary_lifestyle, "_",guild$secondary_lifestyle)

# Relative abundance of ECM fungi
aggregated <- aggregate(percent, by=list(guild$lifestyle),FUN = sum,na.rm=F) 
row.names(aggregated)<-aggregated[,1]
aggregated <- aggregated[,-1]
aggregated <- data.frame(aggregated)
rows <- grep ("ectomycorrhizal", rownames (aggregated)) 
subset <- aggregated[rows,]
subset.t <- t(subset)
ECM <- data.frame(rowSums (subset.t))
colnames (ECM) [1] <- "ECM"

# Normalizing ECM fungi % by ECM tree %
Result <- ECM
data <- cbind (Result, DESIGN,METADATA)
data$Urban <- factor (data$Urban, levels=c("Urban","Rural"))
data$ECM <- log10(data$ECM)-(log10(data$ECM.tree+0.1)) # To use 0% data

# To show the result of lmer in the graph
anova <- anova(lmer(scale(ECM)~ scale(DFB)*scale(DFE)+(1|Site),data=data))
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

# ggplot
ggplot(data)+
geom_point(aes(x=DFE, y=ECM, color=Urban),position=position_jitter( width=2, height=0))+ 
geom_smooth(method="lm", aes(x=DFE, y=ECM, group=Urban, color=Urban))+  
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
labs (y="Normalized abundance of ECM fungi", x="Distance from edge (m)") + 
annotate("text", x=-Inf, y=Inf, hjust=0, vjust =2, label=DFB.result, size=5) +
annotate("text", x=-Inf, y=Inf, hjust=0, vjust =4, label=DFE.result, size=5) +
annotate("text", x=-Inf, y=Inf, hjust=0, vjust =6, label=DFBDFE.result, size=5)

# Save
ggsave("normalized.ECM.abund.png",width = 5, height = 4)
