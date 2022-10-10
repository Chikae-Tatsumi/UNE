library(tidyverse)
library (dplyr)
library(ggplot2)

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv("experimental_design.csv",header=T)
DESIGN <- na.omit(DESIGN)
setwd("~/R/Analysis/2_UNE/16S/function")

function.table <- read.csv(file="rarefied_16S_function.csv",row.names = 1,header=T)
ASV <- function.table [,1:(ncol(function.table)-19)] 
fgs <- function.table [,(ncol(function.table)-12):ncol(function.table)] 
percent <- ASV / mean(colSums(ASV)) *100

name <- colnames(fgs) 
colnames(fgs) <- c("Cellulolytic bacteria", "Assimiratory nitrite reducing bacteria", "Dissimiratory nitrite reducing bacteria", 
"Assimiratory nitrate reducing bacteria", "N fixing bacteria", "Dissimiratory nitrate reducing bacteria",
"Nitrifying bacteria", "Denitrifying bacteria", "Chitinolytic bacteria", "Lignolytic bacteria",              
"Methanotroph bacteria", "Copiotroph bacteria", "Oligotroph bacteria")

Result <- NULL
for (i in 1:ncol(fgs)){
aggregated <- aggregate(percent, by=list(fgs[,i]),FUN = sum,na.rm=F) 
rownames(aggregated) <- aggregated [,1]
aggregated <- aggregated [,-1]
aggregated.t <- t(aggregated)
data <- cbind (aggregated.t, DESIGN)
data$Urban <- factor (data$Urban, levels=c("Urban","Rural"))

# To show the result of lmer in the graph
anova <- anova(lmer(scale(data[,1])~ scale(DFB)*scale(DFE)+(1|Site),data=data))
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
geom_point(aes(x=DFE, y=data[,1], color=Urban),position=position_jitter( width=2, height=0))+ 
geom_smooth(method="lm", aes(x=DFE, y=data[,1], group=Urban, color=Urban))+  
# scale_fill_manual(values = c("#C77CFF","#7CAE00","#00BFC4","#F8766D"))+  # if you want to change the colors
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
labs (y=paste(colnames(fgs)[i], "(%)",sep = " "), x="Distance from edge (m)") + # if you want to change the axis titles
annotate("text", x=-Inf, y=Inf, hjust=0, vjust =2, label=DFB.result, size=5) +
annotate("text", x=-Inf, y=Inf, hjust=0, vjust =4, label=DFE.result, size=5) +
annotate("text", x=-Inf, y=Inf, hjust=0, vjust =6, label=DFBDFE.result, size=5) 

# Save
ggsave(paste(colnames(fgs)[i],".png"),width = 5, height = 4)

Result <- cbind (Result, aggregated.t)
Result <- Result [,-ncol(Result)]}
colnames(Result) <- name
write.csv (Result, "aggregated.function.table.csv")
