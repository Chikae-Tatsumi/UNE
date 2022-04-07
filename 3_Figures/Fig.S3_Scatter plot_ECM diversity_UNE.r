library(ggplot2)
library(lme4)
library(lmerTest)
library(vegan)

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv("experimental_design.csv",header=T)
setwd("~/R/Analysis/2_UNE/ITS/FungalTrait")
ECM.table<- read.csv(file="rarefied_ECM_table_FungalTrait.csv",header=T,row.names = 1)
ASV.ECM <- ECM.table [,1:(ncol(ECM.table)-31)] 
guild.ECM <- ECM.table [,(ncol(ECM.table)-24):ncol(ECM.table)] 
ASV.t <- t(ASV.ECM)

# Diversity index calculation
shannon <- diversity(ASV.t, index="shannon",base=2)
simpson <- diversity(ASV.t, index="simpson")
invsimpson <- diversity(ASV.t, index="invsimpson")
data <- cbind(shannon, simpson, invsimpson)
write.csv (data, "ECM.diversity.FungalTrait.csv")

data<-cbind(DESIGN, data)
data$Urban <- factor (data$Urban, levels=c("Urban","Rural"))

# To show the result of lmer in the graph
anova <- anova(lmer(scale(shannon)~ scale(DFB)*scale(DFE)+(1|Site),data=data))
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
geom_point(aes(x=DFE, y=shannon, color=Urban),position=position_jitter(width=2, height=0))+
geom_smooth(method="lm", aes(x=DFE, y=shannon, group=Urban, color=Urban))+
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
labs (y="Shannon's diversity index",x="Distance from edge (m)") + 
annotate("text", x=-Inf, y=Inf, hjust=0, vjust =2, label=DFB.result, size=5) +
annotate("text", x=-Inf, y=Inf, hjust=0, vjust =4, label=DFE.result, size=5) +
annotate("text", x=-Inf, y=Inf, hjust=0, vjust =6, label=DFBDFE.result, size=5) 

# Save
ggsave("Plot_Shannon.ECM.FungalTrait.png",width = 5, height = 4)
