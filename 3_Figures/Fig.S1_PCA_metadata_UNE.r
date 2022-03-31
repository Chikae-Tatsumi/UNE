library(maptools)
library(ggplot2)
library(ggrepel)

# Import files
setwd("~/R/Analysis/2_UNE")
METADATA <- read.csv(file = "metadata.csv",header=T)
DESIGN <- read.csv(file = "experimental_design.csv",header=T)

data <- cbind (-(DESIGN$DFB),-(DESIGN$DFE), METADATA$BA, METADATA$Moist,METADATA$pH,METADATA$SOM,METADATA$Temp,METADATA$NH4,METADATA$NO3,METADATA$NOx,METADATA$O3,METADATA$freeze_to_thaw_transition_0.5C,METADATA$wet_to_dry_transition_10perc)
colnames(data) <- c("Urbanization","Fragmentation","Basal area","Moisture","pH","SOM","Temperature","NH4","NO3","atm NOx","atm O3","F/T cycle","W/D cycle")

# Make dataset
data<- na.omit(data)
pilots.pca <- prcomp(data,scale=TRUE, center = TRUE)  #standardized
biplot(pilots.pca) # Check the result
loading <- sweep(pilots.pca$rotation,MARGIN=2,pilots.pca$sdev,FUN="*")
loading <- data.frame(loading)

# ggplot
ggplot(loading) + 
geom_segment(aes(xend=PC1, yend=PC2), x=0, y=0, arrow=arrow(length = unit(0.1,"cm"))) + 
geom_text_repel(aes(x=PC1, y=PC2, label=rownames(loading)),  size=4, color='black') +
xlim(-1,1) + 
ylim(-1,1) +
theme_classic()+
theme(text=element_text(size=10,color="black"),
axis.text=element_text(size=10,color="black"))+
coord_fixed()

# Save
dir.create("PCA")
ggsave(file = "PCA/PCA.metadata.png",height=4,width=4)