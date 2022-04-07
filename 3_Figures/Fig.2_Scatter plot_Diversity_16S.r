library(vegan)
library(ggplot2)
library(lme4)
library(lmerTest)

setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv("experimental_design.csv",header=T)
DESIGN <- na.omit (DESIGN)
setwd("~/R/Analysis/2_UNE/16S")
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T)

ASV <- ASV.table [,1:(ncol(ASV.table)-6)] 
ASV.t <- t(ASV)
shannon <- diversity(ASV.t, index="shannon",base=2)
simpson <- diversity(ASV.t, index="simpson")
invsimpson <- diversity(ASV.t, index="invsimpson")
fisher <- fisher.alpha(ASV.t)
data <- cbind(shannon, simpson, invsimpson,fisher)
write.csv (data, "diversity.csv")
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

ggplot(data)+
geom_point(aes(x=DFE, y=shannon, color=Urban),position=position_jitter( width=2, height=0))+
geom_smooth(method="lm", aes(x=DFE, y=shannon, group=Urban, color=Urban))+
# scale_fill_manual(values = c("#C77CFF","#7CAE00","#00BFC4","#F8766D"))+  # if you want to change the colors
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
labs (y="Shannon's diversity index",x="Distance from edge (m)") + # if you want to change the axis titles
annotate("text", x=-Inf, y=Inf, hjust=0, vjust =2, label=DFB.result, size=5) +
annotate("text", x=-Inf, y=Inf, hjust=0, vjust =4, label=DFE.result, size=5) +
annotate("text", x=-Inf, y=Inf, hjust=0, vjust =6, label=DFBDFE.result, size=5) 

ggsave("Plot_Shannon.png",width = 5, height = 4)
