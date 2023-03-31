library(ggplot2)
library(lme4)
library(lmerTest)

setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv("experimental_design.csv",header=T)
setwd("~/R/Analysis/2_UNE/ITS/FungalTrait")
fungaltrait <- read.csv(file="aggregated.fungaltrait.table.csv",header=T,row.names = 1)

data <- cbind (fungaltrait, DESIGN)
data$Urban <- factor (data$Urban, levels=c("Urban","Rural"))

# To show the result of lmer in the graph
anova <- anova(lmer(scale(ECM)~ scale(DFB)*scale(DFE)+(1|Site),data=data))
# Pval <- anova[,6]
# if (Pval[1] > 0.05) {DFB.result <- ""
# } else if (Pval[1] > 0.01) {DFB.result <- "        Urbanization *      "
# } else if (Pval[1] > 0.001) {DFB.result <- "        Urbanization **     "
# } else {DFB.result <- "        Urbanization ***    "}
# if (Pval[2] > 0.05) {DFE.result <- ""
# }else if (Pval[2] > 0.01) {DFE.result <- "        Fragmentation *       "
# }else if (Pval[2] > 0.001) {DFE.result <- "        Fragmentation **     "
# }else {DFE.result <- "        Fragmentation ***    "}
# if (Pval[3] > 0.05) {DFBDFE.result <- ""
# }else if (Pval[3] > 0.01) {DFBDFE.result <- "        U×F *      "
# }else if (Pval[3] > 0.001) {DFBDFE.result <- "        U×F **     "
# }else {DFBDFE.result <- "        U×F ***    "}

write.csv(anova, "ECM_rel.abund_anova.csv")

ggplot(data)+
geom_point(aes(x=DFE, y=ECM, color=Urban),position=position_jitter( width=2, height=0))+ 
geom_smooth(method="lm", aes(x=DFE, y=ECM, group=Urban, color=Urban))+  
scale_color_manual(values = c("red","black"))+  # if you want to change the colors
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
labs (y="Soil ECM fungi (%)", x="Distance from edge (m)",color="")  # if you want to change the axis titles
# annotate("text", x=-Inf, y=Inf, hjust=0, vjust =2, label=DFB.result, size=5) +
# annotate("text", x=-Inf, y=Inf, hjust=0, vjust =4, label=DFE.result, size=5) +
# annotate("text", x=-Inf, y=Inf, hjust=0, vjust =6, label=DFBDFE.result, size=5) 

ggsave("ECM_rel.abund_UNE_nostats.png",width = 5, height = 4)
