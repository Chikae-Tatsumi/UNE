library(ggplot2)
library(vegan)
library(ggrepel)

# Import files
setwd("~/R/Analysis/2_UNE")
DESIGN <- read.csv("experimental_design.csv")
setwd("~/R/Analysis/2_UNE/Network")
func <- read.csv("function.absolute.abundance.csv",header=T,row.names=1)

DFBDFE <- cbind(DESIGN$DFB, DESIGN$DFE)
colnames(DFBDFE) <- c("DFB","DFE")
ag <- aggregate(DFBDFE, by=list(DESIGN$Plot40),FUN = mean,na.rm=F) 
Plot40 <- ag[,1]
DFB <- ag$DFB
DFE <- ag$DFE
Urban <- c(rep("Urban",10),rep("Rural",20),rep("Urban",10))

# Calculate the data characteristics
Result<-NULL
for (i in 1:40){
data <- subset(func, DESIGN$Plot40==Plot40[i])
data <- t(data)
proportion <- colSums(data)/sum(data)*100
var <- var(proportion)
sd <- sd(proportion)
median <- median(proportion)
shannon <- diversity(rowMeans(data), index="shannon",base=2)
bind<-cbind(shannon, var,sd,median)
Result <- rbind (Result, bind)
}
rownames(Result) <- Plot40
colnames(Result) <-c("shannon","variance","sd", "median")

write.csv(Result, "proportion.csv")

# ANOVA
Result <- data.frame(Result)
Result$DFB <- DFB
Result$DFE <- DFE
Result$Urban <- Urban

stats <- NULL
for (i in 1:4){
    anova<- anova(lm(Result[,i]~DFB*DFE, data=Result))
    Fval <- anova[,4]
    Pval <- anova[,5]
    stats <- cbind(stats,Pval)
}
colnames(stats) <- colnames(Result)[1:4]
rownames(stats) <- rownames(anova)
write.csv(stats, "Pval.csv")

# Visualization
Result$Urban <- factor (Result$Urban, levels=c("Urban","Rural"))

ggplot(Result)+
geom_point(aes(x=DFE, y=shannon, color=Urban),position=position_jitter(width=2, height=0))+
geom_smooth(method="lm", aes(x=DFE, y=shannon, color=Urban, group=Urban))+
scale_color_manual(values = c("#f6766d","#4B0082"))+  # if you want to change the colors
labs (y="Shannon's diversity index", x="Distance from edge (m)",color="") + # if you want to change the axis titles
theme_classic()

ggsave("shannon.png",width = 4, height = 4)
